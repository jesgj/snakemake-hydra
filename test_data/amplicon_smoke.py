#!/usr/bin/env python3
"""Exercise amplicon modes using the FASTQs bundled with pinned DADA2 (no downloads)."""

import argparse
import copy
import csv
import gzip
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from amplicon_config import normalize
sys.path.insert(0, str(REPO / "scripts"))
from amplicon_reads import audit


def run(config, cores, log, target=None):
    env = {**os.environ, "HYDRA_CONFIG": str(config)}
    command = ["pixi", "run", "--locked", "snakemake", "--cores", str(cores), "--jobs", "1"]
    if target:
        command.append(str(target))
    with log.open("w") as handle:
        return subprocess.run(command, cwd=REPO, env=env, stdout=handle, stderr=subprocess.STDOUT).returncode


def verify(root, taxonomy_enabled):
    with (root / "asv_counts.tsv").open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        ids = reader.fieldnames[1:]
        counts = {r["sample"]: sum(int(v) for k, v in r.items() if k != "sample") for r in reader}
    assert ids and sum(counts.values()) > 0
    tracking = list(csv.DictReader((root / "read_tracking.tsv").open(), delimiter="\t"))
    assert len(tracking) == 2
    for row in tracking:
        assert counts[row["sample"]] == int(row["nonchimeric"])
        assert int(row["input"]) >= int(row["primer_retained"]) >= int(row["filtered"]) >= int(row["nonchimeric"])
    fasta_ids = [line[1:] for line in (root / "asv_sequences.fasta").read_text().splitlines() if line.startswith(">")]
    assert fasta_ids == ids
    taxonomy = list(csv.DictReader((root / "taxonomy.tsv").open(), delimiter="\t"))
    assert [r["asv_id"] for r in taxonomy] == ids
    if taxonomy_enabled:
        assert any(r.get("Kingdom") == "Bacteria" for r in taxonomy)
    assert json.loads((root / "taxonomy_metadata.json").read_text())["status"] == ("completed" if taxonomy_enabled else "not_requested")
    report = root / "multiqc/multiqc_report.html"
    assert report.stat().st_size > 1000 and "amplicon_read_retention" in report.read_text()
    samples = json.loads((root / "resolved_config.json").read_text())["samples_info"]
    expected_qc = 2 * sum(2 if s.get("R2") else 1 for s in samples.values())
    fastqc = list(csv.DictReader((root / "multiqc/multiqc_report_data/multiqc_fastqc.txt").open(), delimiter="\t"))
    assert len(fastqc) == expected_qc, "Raw and primer-removed QC reports must remain distinct in MultiQC"
    return {"asvs": len(ids), "retained": counts, "report_sha256": hashlib.sha256(report.read_bytes()).hexdigest()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workdir", type=Path, default=REPO / "test_data/work/amplicon_smoke")
    args = parser.parse_args()
    work = args.workdir.resolve()
    work.mkdir(parents=True, exist_ok=True)
    ext = Path(subprocess.check_output(["Rscript", "-e", "cat(system.file('extdata', package='dada2'))"], text=True))
    base = yaml.safe_load((REPO / "config/config.yaml").read_text())
    base["pipeline"] = "amplicon"
    a = base["amplicon"]
    a["primers"]["mode"] = "already_removed"
    a["dada2"]["trunc_len"] = [240, 200]
    a["threads"] = 2
    a["samples_info"] = {f"sam{i}": {"R1": str(ext / f"sam{i}F.fastq.gz"), "R2": str(ext / f"sam{i}R.fastq.gz"), "batch": f"run{i}"} for i in (1, 2)}
    # Require the 5' primer on both mates, clip 3' read-through, reject orphan pairs.
    small = work / "primer_contract"
    small.mkdir(exist_ok=True)
    primers = ("ACGTTGCACAGTCTAGACGA", "GTCACGTAGATCGTACAGTC")
    rc = lambda p: p.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    expected_sequences = []
    for mate, suffix in ((1, "F"), (2, "R")):
        with gzip.open(ext / f"sam1{suffix}.fastq.gz", "rt") as src:
            header, sequence, plus, quality = [src.readline().strip() for _ in range(4)]
        expected_sequences.append(sequence)
        text = ""
        for index in (1, 2):
            prefix = primers[mate - 1] if index == 1 or mate == 1 else "TTTTTTTTTTTTTTTTTTTT"
            seq = prefix + sequence + rc(primers[2 - mate])
            text += f"@test{index}/{mate}\n{seq}\n+\n{'I' * len(prefix)}{quality}{'I' * len(primers[2 - mate])}\n"
        (small / f"R{mate}.fastq").write_text(text)
    cfg = copy.deepcopy(a)
    cfg["output_dir"] = str(small / "output")
    cfg["samples_info"] = {"contract": {"R1": str(small / "R1.fastq"), "R2": str(small / "R2.fastq")}}
    cfg["primers"].update(mode="remove", forward=primers[0], reverse=primers[1], min_overlap=10)
    resolved = small / "resolved.json"
    resolved.write_text(json.dumps(normalize(cfg)))
    with (small / "cutadapt.log").open("w") as log:
        subprocess.run([sys.executable, str(REPO / "scripts/amplicon_reads.py"), "prepare", "--config", str(resolved), "--sample", "contract"], check=True, stdout=log, stderr=subprocess.STDOUT)
    for mate in (1, 2):
        with gzip.open(small / "output/primers" / f"contract_R{mate}.fastq.gz", "rt") as handle:
            records = handle.read().splitlines()
            assert len(records) == 4 and records[1] == expected_sequences[mate - 1]
    broken = small / "broken_R2.fastq"
    broken.write_text((small / "R2.fastq").read_text().replace("@test1/2", "@different/2"))
    try:
        audit([small / "R1.fastq", broken])
    except ValueError as exc:
        assert "Unsynchronized" in str(exc)
    else:
        raise AssertionError("Mismatched mate IDs must fail")
    invalid = copy.deepcopy(a)
    invalid["primers"]["mode"] = "remove"
    try:
        normalize(invalid)
    except ValueError as exc:
        assert "primer" in str(exc)
    else:
        raise AssertionError("Missing primers must be rejected")
    invalid = copy.deepcopy(a)
    invalid["samples_info"]["sam2"].pop("R2")
    try:
        normalize(invalid)
    except ValueError as exc:
        assert "comparable region" in str(exc)
    else:
        raise AssertionError("Mixed merged/partial amplicons must be rejected")
    report = {"dada2_fixture_version": "1.38.0", "checks": {"invalid_configs": "passed", "primer_contract_and_mate_validation": "passed"}}
    (work / "validation.json").write_text(json.dumps(report, indent=2) + "\n")
    for variant in ("paired", "single", "reverse", "primers", "no_overlap"):
        cfg = copy.deepcopy(base)
        a = cfg["amplicon"]
        root = work / variant
        a["output_dir"] = str(root)
        a["log_dir"] = str(work / f"{variant}-logs")
        if variant == "paired":
            a["taxonomy"].update(enabled=True, training_set=str(ext / "example_train_set.fa.gz"), reference_version="DADA2 1.38.0 bundled test example")
            for i in (1, 2):
                a["samples_info"][f"sam{i}"]["batch"] = f"run{i}"
        elif variant == "single":
            # Test strict automatic discovery, including uncompressed single-end inputs.
            raw = work / "single_inputs"
            raw.mkdir(exist_ok=True)
            for i in (1, 2):
                with gzip.open(ext / f"sam{i}F.fastq.gz", "rb") as src:
                    (raw / f"sam{i}.fastq").write_bytes(src.read())
            a["samples_info"] = {}
            a["raw_fastqs_dir"] = str(raw)
        elif variant == "reverse":
            # Test TSV manifests and the explicit R2 fallback with reverse-complemented ASVs.
            sheet = work / "reverse_samples.tsv"
            sheet.write_text("sample\tR1\tR2\tmode\tbatch\n" + "".join(f"sam{i}\t{ext / f'sam{i}F.fastq.gz'}\t{ext / f'sam{i}R.fastq.gz'}\treverse\tdefault\n" for i in (1, 2)))
            a["samples_info"] = {}
            a["sample_sheet"] = str(sheet)
        elif variant == "primers":
            # Engineered 5' primer prefixes on real fixture reads test trimming;
            # these are software fixtures, not biological evidence of primer identity.
            raw = work / "primer_inputs"
            raw.mkdir(exist_ok=True)
            a["primers"].update(mode="remove", forward="ACGTTGCA", reverse="GTCACGTA", min_overlap=8)
            for i in (1, 2):
                for mate, suffix, primer in ((1, "F", "ACGTTGCA"), (2, "R", "GTCACGTA")):
                    target = raw / f"sam{i}_R{mate}.fastq.gz"
                    with gzip.open(ext / f"sam{i}{suffix}.fastq.gz", "rt") as src, gzip.open(target, "wt") as dst:
                        while header := src.readline():
                            seq, plus, qual = [src.readline() for _ in range(3)]
                            dst.writelines([header, primer + seq, plus, "I" * len(primer) + qual])
                    a["samples_info"][f"sam{i}"][f"R{mate}"] = str(target)
        else:
            a["dada2"]["trunc_len"] = [80, 80]
        config = work / f"{variant}.yaml"
        config.write_text(yaml.safe_dump(cfg, sort_keys=False))
        log = work / f"{variant}.log"
        print(f"Checking {variant}...", flush=True)
        code = run(config, 1 if variant == "single" else 2, log, root / "asv_counts.tsv" if variant == "no_overlap" else None)
        if variant == "no_overlap":
            assert code != 0 and "No merged reads" in (work / "no_overlap-logs/dada2.log.err").read_text()
            report["checks"][variant] = "expected failure (no silent fallback)"
        else:
            if code:
                raise RuntimeError(f"{variant} failed; inspect {log}")
            report["checks"][variant] = verify(root, variant == "paired")
            if variant == "primers":
                assert report["checks"][variant]["retained"] == report["checks"]["paired"]["retained"], "Primer removal altered retained reads"
                for i in (1, 2):
                    for mate, suffix in ((1, "F"), (2, "R")):
                        with gzip.open(root / "primers" / f"sam{i}_R{mate}.fastq.gz", "rb") as got, gzip.open(ext / f"sam{i}{suffix}.fastq.gz", "rb") as expected:
                            assert got.read() == expected.read()
        (work / "validation.json").write_text(json.dumps(report, indent=2) + "\n")
    # Saved DAG, rerun, and dependency tracking checks on the paired cohort.
    repeat = work / "repeat.log"
    assert run(work / "paired.yaml", 2, repeat) == 0
    assert "Nothing to be done" in repeat.read_text()
    print(f"Checks passed; current-run evidence: {work / 'validation.json'}", flush=True)


if __name__ == "__main__":
    main()
