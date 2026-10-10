#!/usr/bin/env python3
"""Prepare a bounded, paired reservoir sample of ENA ERR777695 for Hydra."""
import argparse
import copy
import csv
import json
import shutil
import tempfile
from pathlib import Path

import yaml

from prepare import digest, fetch, records, sample_pair

REPO = Path(__file__).resolve().parents[1]
RUN = "ERR777695"
READS = [
    ("f985d801c8e8adfd12cad5d3b7f31ffc", 61664840),
    ("b25b32d8f3e1c387c1f80bbc708ca747", 73798224),
]
REFERENCE = "rdp_19_toGenus_trainset.fa.gz"
REFERENCE_URL = "https://zenodo.org/records/14168771/files/" + REFERENCE
REFERENCE_MD5 = "390b8a359c45648adf538e72a1ee7e28"


def quality_profile(path):
    totals, q20, number = [], [], 0
    for record in records(path):
        quality = record[3].rstrip(b"\r\n")
        while len(totals) < len(quality):
            totals.append(0)
            q20.append(0)
        for index, value in enumerate(quality):
            score = value - 33
            totals[index] += score
            q20[index] += score >= 20
        number += 1
    return {"reads": number, "mean_q": [v / number for v in totals],
            "fraction_q20": [v / number for v in q20]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pairs", type=int, default=20000)
    parser.add_argument("--seed", type=int, default=20261011)
    parser.add_argument("--source-dir", type=Path,
                        help="Reuse downloaded R1.fastq.gz/R2.fastq.gz after MD5 verification")
    parser.add_argument("--verify", action="store_true", help="Verify completed public-fixture results without downloads")
    args = parser.parse_args()
    if not 1000 <= args.pairs <= 100000:
        parser.error("--pairs must be between 1000 and 100000")
    work = REPO / "test_data/work/amplicon_public"
    work.mkdir(parents=True, exist_ok=True)
    if args.verify:
        root = work / "results"
        with (root / "asv_counts.tsv").open() as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            ids = reader.fieldnames[1:]
            rows = list(reader)
        assert ids and len(rows) == 1 and rows[0]["sample"] == RUN
        retained = sum(int(rows[0][key]) for key in ids)
        tracking = list(csv.DictReader((root / "read_tracking.tsv").open(), delimiter="\t"))
        assert len(tracking) == 1 and int(tracking[0]["nonchimeric"]) == retained > 0
        source = json.loads((work / "sampling.json").read_text())
        assert int(tracking[0]["input"]) == source["pairs"]
        assert source["pairs"] >= int(tracking[0]["primer_retained"]) >= int(tracking[0]["filtered"]) >= retained
        assert [line[1:] for line in (root / "asv_sequences.fasta").read_text().splitlines() if line.startswith(">")] == ids
        taxonomy = list(csv.DictReader((root / "taxonomy.tsv").open(), delimiter="\t"))
        assert [row["asv_id"] for row in taxonomy] == ids
        assert any(row.get("Kingdom") == "Bacteria" for row in taxonomy)
        assert json.loads((root / "taxonomy_metadata.json").read_text())["status"] == "completed"
        metadata = json.loads((root / "run_metadata.json").read_text())
        assert len(metadata["error_learning_groups"]) == 2
        for group in metadata["error_learning_groups"]:
            assert (root / group["file"]).stat().st_size > 0
        assert len(list((root / "filtered").glob("*.fastq.gz"))) == 2
        assert (root / "dada2.rds").stat().st_size > 0
        qc = list(csv.DictReader((root / "multiqc/multiqc_report_data/multiqc_fastqc.txt").open(), delimiter="\t"))
        assert len(qc) == 4
        report = root / "multiqc/multiqc_report.html"
        assert report.stat().st_size > 1000 and "amplicon_read_retention" in report.read_text()
        for name, checksum in source["output_sha256"].items():
            assert digest(work / name) == checksum
        result = {"status": "passed", "asvs": len(ids), "tracking": tracking[0],
                  "tool_metadata": metadata, "taxonomy": "RDP 19 completed",
                  "report_sha256": digest(report)}
        (work / "validation.json").write_text(json.dumps(result, indent=2) + "\n")
        print(json.dumps(result, indent=2))
        return
    reference = REPO / "test_data/references/amplicon" / REFERENCE
    reference.parent.mkdir(parents=True, exist_ok=True)
    if not reference.exists() or digest(reference, "md5") != REFERENCE_MD5:
        fetch(REFERENCE_URL, reference)
    if digest(reference, "md5") != REFERENCE_MD5:
        raise ValueError("RDP reference MD5 mismatch")
    if shutil.disk_usage(work).free < 1024 ** 3:
        raise OSError("Need at least 1 GiB free for preparation and results")
    provenance = {"accession": RUN, "study": "PRJEB6244", "sample": "SAMEA3298272",
                  "paper_doi": "10.1093/nar/gku1341", "sample_alias": "M35",
                  "pairs": args.pairs, "seed": args.seed, "sources": []}
    outputs = [work / f"{RUN}_R{mate}.fastq.gz" for mate in (1, 2)]
    with tempfile.TemporaryDirectory(prefix="ena-", dir=work) as temporary:
        paths = []
        for mate, (md5, size) in enumerate(READS, 1):
            url = f"https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR777/{RUN}/{RUN}_{mate}.fastq.gz"
            source = (args.source_dir or Path(temporary)) / f"R{mate}.fastq.gz"
            if not source.exists():
                if args.source_dir:
                    raise FileNotFoundError(source)
                fetch(url, source)
            if source.stat().st_size != size or digest(source, "md5") != md5:
                raise ValueError(f"ENA size/MD5 mismatch: {source}")
            provenance["sources"].append({"url": url, "bytes": size, "md5": md5,
                                          "sha256": digest(source)})
            paths.append(source)
        provenance["source_pairs"] = sample_pair(paths, args.pairs, args.seed, outputs)
    provenance["output_sha256"] = {p.name: digest(p) for p in outputs}
    provenance["reference"] = {"url": REFERENCE_URL, "md5": REFERENCE_MD5,
                               "sha256": digest(reference), "release": "RDP 19 (2023-08-23)"}
    (work / "sampling.json").write_text(json.dumps(provenance, indent=2) + "\n")
    (work / "quality_profiles.json").write_text(json.dumps(
        {p.name: quality_profile(p) for p in outputs}, indent=2) + "\n")
    config = copy.deepcopy(yaml.safe_load((REPO / "config/config.yaml").read_text()))
    config["pipeline"] = "amplicon"
    a = config["amplicon"]
    a.update(output_dir=str(work / "results"), log_dir=str(work / "logs"), threads=2,
             mode="merge", samples_info={RUN: {"R1": str(outputs[0]), "R2": str(outputs[1]),
                                               "batch": RUN}})
    # Supplement Table S1: DI V4 F515A/805R target sequences (adapter tails excluded).
    # Table S3 labels M35 "515 & 805RA"; the deposited R1 starts match the DI
    # degenerate F515A sequence, including the B position, rather than fixed-C 515.
    a["primers"].update(mode="remove", forward="GTGBCAGCMGCCGCGGTAA",
                         reverse="GACTACHVGGGTATCTAATCC")
    # After removing 19/20 nt primers, retain raw cycles 20-219 and 21-200.
    # The measured R2 quality falls sharply beyond raw cycle 200. V4 retains
    # ample overlap with these lengths; merging must succeed without fallback.
    a["dada2"].update(trunc_len=[200, 180], seed=args.seed)
    a["taxonomy"].update(enabled=True, training_set=str(reference),
                          reference_version="RDP 19 (2023-08-23), DADA2 format DOI 10.5281/zenodo.14168771")
    configs = REPO / "test_data/configs"
    configs.mkdir(exist_ok=True)
    destination = configs / "amplicon_public.yaml"
    destination.write_text(yaml.safe_dump(config, sort_keys=False))
    print(f"Prepared {args.pairs} of {provenance['source_pairs']} paired reads; config: {destination}")


if __name__ == "__main__":
    main()
