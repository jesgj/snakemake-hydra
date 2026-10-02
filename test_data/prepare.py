#!/usr/bin/env python3
"""Fetch real libraries, sample synchronized pairs, and remove temporary originals."""
import argparse
import gzip
import hashlib
import json
import random
import subprocess
import tempfile
from itertools import zip_longest
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def digest(path, algorithm="sha256"):
    result = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            result.update(chunk)
    return result.hexdigest()


def fetch(url, path):
    # Catalog sources are under 100 MiB per file. Reject accidental large downloads.
    subprocess.run([
        "curl", "-L", "-f", "-sS", "--connect-timeout", "20", "--max-time", "180",
        "--max-filesize", str(100 * 1024 * 1024), "--retry", "2", "-o", str(path), url,
    ], check=True)


def records(path):
    with gzip.open(path, "rb") as handle:
        while header := handle.readline():
            seq, plus, quality = (handle.readline() for _ in range(3))
            if (not header.startswith(b"@") or not plus.startswith(b"+")
                    or not quality.endswith(b"\n") or len(seq.rstrip()) != len(quality.rstrip())):
                raise ValueError(f"Invalid FASTQ record in {path}")
            yield (header, seq, plus, quality)


def read_id(record):
    name = record[0].split()[0]
    return name[:-2] if name.endswith((b"/1", b"/2")) else name


def write_fastq(path, records_to_write):
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as handle:
            for record in records_to_write:
                handle.writelines(record)


def sample_pair(paths, count, seed, outpaths):
    rng = random.Random(seed)
    reservoir = []
    total = 0
    for index, pair in enumerate(zip_longest(*(records(path) for path in paths))):
        if None in pair or read_id(pair[0]) != read_id(pair[1]):
            raise ValueError(f"Unsynchronized mates at record {index + 1}: {paths}")
        total = index + 1
        if index < count:
            reservoir.append((index, pair))
        else:
            slot = rng.randrange(index + 1)
            if slot < count:
                reservoir[slot] = (index, pair)
    if total <= count:
        raise ValueError(f"Need more than {count} source pairs, found {total}")
    reservoir.sort(key=lambda item: item[0])
    for mate, path in enumerate(outpaths):
        write_fastq(path, (pair[mate] for _, pair in reservoir))
    return total


def make_configs(catalog):
    import copy
    import yaml  # Already provided by the committed Pixi environment.

    base = yaml.safe_load((ROOT.parent / "config/config.yaml").read_text())
    configs = ROOT / "configs"
    configs.mkdir(exist_ok=True)
    variants = [("rnaseq", "rnaseq"), ("wgbs", "wgbs"), ("wgbs_full", "wgbs"), ("chipseq", "chip_cr"),
                ("chipseq_se", "chip_cr"), ("cutrun", "chip_cr"),
                ("atacseq_macs3", "atacseq"), ("atacseq_genrich", "atacseq")]
    for variant, mode in variants:
        category = mode if mode in ("atacseq", "wgbs") else variant.removesuffix("_se")
        config = {"pipeline": mode, **{key: {} for key in ("rnaseq", "wgbs", "chip_cr", "atacseq")}}
        section = copy.deepcopy(base[mode])
        for key in section:
            if key.endswith("_dir"):
                section[key] = str(ROOT / "work" / variant / key)
        refcategory = "yeast" if mode in ("chip_cr", "atacseq") else mode
        section["ref_genome"] = str(ROOT / "references" / refcategory / "genome.fa")
        if variant == "wgbs_full":
            # Bismark indexes every FASTA in a directory: keep full and small genomes apart.
            section["ref_genome"] = str(ROOT / "references/tair10/genome.fa")
        section["raw_fastqs_dir"] = str(ROOT / "fastqs" / category)
        section["samples_info"] = {}
        for sample in catalog["samples"]:
            if sample["category"] != category:
                continue
            stem = ROOT / "fastqs" / category / sample["sample"]
            info = {"R1": str(stem) + "_R1.fastq.gz", "type": "SE" if variant.endswith("_se") else "PE"}
            if info["type"] == "PE":
                info["R2"] = str(stem) + "_R2.fastq.gz"
            section["samples_info"][sample["sample"]] = info
        section.setdefault("picard", {})["java_opts"] = "-Xmx1g"
        if mode == "rnaseq":
            section["transcriptome_fasta"] = str(ROOT / "references/rnaseq/transcriptome.fa")
            section["kallisto_index"] = str(ROOT / "work" / variant / "kallisto_index/transcriptome.idx")
            Path(section["kallisto_index"]).parent.mkdir(parents=True, exist_ok=True)
            section["kallisto"] = {"extra_args": "--bootstrap-samples=0"}
            section["gene_body_coverage"].update(enabled=True, refgene_bed=str(ROOT / "references/rnaseq/genes.bed12"))
        if mode == "wgbs":
            section["bismark"].update(threads=2, parallel=1)
        if mode == "chip_cr":
            section["gene_bed"] = str(ROOT / "references/yeast/genes.bed")
            if variant == "cutrun":
                # This genuine early CUT&RUN library has 25 bp reads.
                section["fastp"]["extra_args"] = "-q 20 -l 20"
        if mode == "atacseq":
            section["peak_caller"] = variant.split("_")[-1]
            section["macs3"]["genome_size"] = "1.2e7"
        config[mode] = section
        (configs / (variant + ".yaml")).write_text(yaml.safe_dump(config, sort_keys=False))


def bed12_from_gtf():
    import re
    transcripts = {}
    for line in (ROOT / "references/rnaseq/genes.gtf").read_text().splitlines():
        if line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 9 or fields[2] != "exon":
            continue
        match = re.search(r'transcript_id "([^"]+)"', fields[8])
        if match:
            key = (fields[0], match[1], fields[6])
            transcripts.setdefault(key, set()).add((int(fields[3]) - 1, int(fields[4])))
    lines = []
    for (chrom, name, strand), exons in sorted(transcripts.items()):
        blocks = sorted(exons)
        start, end = blocks[0][0], max(e for _, e in blocks)
        sizes = ",".join(str(e - s) for s, e in blocks) + ","
        offsets = ",".join(str(s - start) for s, _ in blocks) + ","
        lines.append(f"{chrom}\t{start}\t{end}\t{name}\t0\t{strand}\t{start}\t{end}\t0\t{len(blocks)}\t{sizes}\t{offsets}\n")
    if not lines:
        raise ValueError("No transcript exons in RNA-seq GTF")
    (ROOT / "references/rnaseq/genes.bed12").write_text("".join(lines))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--configs-only", action="store_true", help="Regenerate isolated configs without downloads")
    args = parser.parse_args()
    catalog = json.loads((ROOT / "datasets.json").read_text())
    if args.configs_only:
        make_configs(catalog)
        return
    report = {"seed": catalog["seed"], "method": "paired reservoir sampling without replacement over every source record", "samples": [], "references": []}
    for sample in catalog["samples"]:
        print(f"Preparing {sample['category']}: {sample['sample']}", flush=True)
        outdir = ROOT / "fastqs" / sample["category"]
        outdir.mkdir(parents=True, exist_ok=True)
        outpaths = [outdir / (sample["sample"] + f"_R{i}.fastq.gz") for i in (1, 2)]
        if any(path.exists() for path in outpaths):
            raise FileExistsError(f"Fixture already exists: {outpaths}; move it aside before regenerating")
        with tempfile.TemporaryDirectory(prefix="hydra-fastqs-") as temporary:
            originals = [Path(temporary) / f"R{i}.fastq.gz" for i in (1, 2)]
            for index, (url, path) in enumerate(zip(sample["urls"], originals)):
                fetch(url, path)
                if "md5" in sample and digest(path, "md5") != sample["md5"][index]:
                    raise ValueError(f"Source MD5 mismatch: {url}")
            source_hashes = [digest(path) for path in originals]
            source_bytes = [path.stat().st_size for path in originals]
            total = sample_pair(originals, sample["pairs"], catalog["seed"], outpaths)
        entry = dict(sample, source_pairs=total, source_bytes=source_bytes, source_sha256=source_hashes,
                     output_sha256=[digest(path) for path in outpaths], output_bytes=[path.stat().st_size for path in outpaths])
        report["samples"].append(entry)
        (ROOT / "provenance.json").write_text(json.dumps(report, indent=2) + "\n")
        print(f"  retained {sample['pairs']:,}/{total:,} pairs; temporary originals deleted", flush=True)
    for reference in catalog["references"]:
        path = ROOT / "references" / reference["path"]
        path.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix="hydra-reference-") as temporary:
            original = Path(temporary) / "reference"
            fetch(reference["url"], original)
            source_hash = digest(original)
            if reference["url"].endswith(".gz"):
                import shutil
                with gzip.open(original, "rb") as src, path.open("wb") as dst:
                    shutil.copyfileobj(src, dst)
            else:
                path.write_bytes(original.read_bytes())
        report["references"].append(dict(reference, source_sha256=source_hash, output_sha256=digest(path)))
    bed12_from_gtf()
    (ROOT / "provenance.json").write_text(json.dumps(report, indent=2) + "\n")
    make_configs(catalog)
    print("Ready: test_data/fastqs, references, configs and provenance.json", flush=True)


if __name__ == "__main__":
    main()
