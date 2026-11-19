#!/usr/bin/env python3
"""Focused shared FASTQ/read-group validator regressions; no Pixi or data needed."""

import gzip
import json
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))
from utils import discover_samples
from amplicon_reads import audit

PE_VALIDATOR = REPO / "scripts/validate_fastq_pairs.py"
BAM_VALIDATOR = REPO / "src/validate_bam_read_group.py"
MERGED_VALIDATOR = REPO / "src/validate_merged_bam_read_groups.py"


def write_fastq(path, records):
    with path.open("w") as handle:
        for header, sequence, quality in records:
            handle.write(f"{header}\n{sequence}\n+\n{quality}\n")


def run_validator(script, *args, stdin=None):
    return subprocess.run(
        [sys.executable, str(script), *map(str, args)],
        input=stdin,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )


def check_discovery():
    with tempfile.TemporaryDirectory(prefix="hydra-bare-mates-") as directory:
        root = Path(directory)
        for name in ("sample1.fastq.gz", "sample2.fastq.gz"):
            (root / name).touch()
        found = discover_samples(root)
        assert list(found) == ["sample1", "sample2"], found
        assert {sample: info["type"] for sample, info in found.items()} == {
            "sample1": "SE",
            "sample2": "SE",
        }
        (root / "pair_R1.fastq.gz").touch()
        (root / "pair_R2.fastq.gz").touch()
        found = discover_samples(root)
        assert found["pair"]["type"] == "PE"


def check_fastq_validators():
    with tempfile.TemporaryDirectory(prefix="hydra-fastq-") as directory:
        root = Path(directory)
        r1 = root / "sample_R1.fastq"
        r2 = root / "sample_R2.fastq"
        write_fastq(r1, [("@read1/1 A:B:C:D", "ACGT", "IIII")])
        write_fastq(r2, [("@read1/2 A:B:C:D", "TGCA", "IIII")])
        result = run_validator(PE_VALIDATOR, "--r1", r1, "--r2", r2, "--all")
        assert result.returncode == 0, result.stdout + result.stderr
        assert "status\tPASS" in result.stdout

        write_fastq(r2, [("@read2/2 A:B:C:D", "TGCA", "IIII")])
        result = run_validator(PE_VALIDATOR, "--r1", r1, "--r2", r2, "--all")
        assert result.returncode == 1, result.stdout + result.stderr
        assert "status\tFAIL" in result.stdout

        single = root / "single.fastq"
        write_fastq(single, [("@read1/1", "ACGT", "IIII")])
        result = run_validator(PE_VALIDATOR, "--r1", single, "--all")
        assert result.returncode == 0, result.stdout + result.stderr
        assert "mode\tsingle_end" in result.stdout

        compressed_single = root / "single.fq.gz"
        with gzip.open(compressed_single, "wt") as handle:
            handle.write("@read1/1\nACGT\n+\nIIII\n")
        result = run_validator(PE_VALIDATOR, "--r1", compressed_single, "--all")
        assert result.returncode == 0, result.stdout + result.stderr

        uncompressed_fq = root / "single.fq"
        write_fastq(uncompressed_fq, [("@read1/1", "ACGT", "IIII")])
        result = run_validator(PE_VALIDATOR, "--r1", uncompressed_fq, "--all")
        assert result.returncode == 0, result.stdout + result.stderr

        write_fastq(single, [("@read1/2", "ACGT", "IIII")])
        result = run_validator(PE_VALIDATOR, "--r1", single, "--all")
        assert result.returncode == 0, result.stdout + result.stderr
        assert "status\tPASS" in result.stdout
        assert audit([str(single)])[0] == 1

        swapped = run_validator(PE_VALIDATOR, "--r1", single, "--r2", r2, "--all")
        assert swapped.returncode != 0
        assert "expected read 1" in swapped.stdout

        unlabeled = root / "unlabeled.fastq"
        write_fastq(unlabeled, [("@read1", "ACGT", "IIII")])
        result = run_validator(PE_VALIDATOR, "--r1", unlabeled, "--r2", unlabeled, "--all")
        assert result.returncode == 2, result.stdout + result.stderr
        assert "same file" in result.stdout


def check_amplicon_same_file():
    with tempfile.TemporaryDirectory(prefix="hydra-amplicon-") as directory:
        single = Path(directory) / "sample.fastq"
        write_fastq(single, [("@read1", "ACGT", "IIII")])
        try:
            audit([str(single), str(single)])
        except ValueError as exc:
            assert "same file" in str(exc)
        else:
            raise AssertionError("Amplicon audit must reject identical R1/R2 paths")


def sam(*records, read_groups=()):
    header = ["@HD\tVN:1.6\tSO:unsorted", "@SQ\tSN:chr1\tLN:1000"]
    header.extend(f"@RG\t{fields}" for fields in read_groups)
    alignments = [
        f"{name}\t0\tchr1\t{position}\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:{group}"
        for name, position, group in records
    ]
    return "\n".join([*header, *alignments]) + "\n"


def check_bam_validators():
    expected = ["--expected-id", "rg1", "--expected-sample", "sample1",
                "--expected-library", "lib1", "--expected-platform", "ILLUMINA",
                "--expected-platform-unit", "unit1"]
    valid = sam(("read1", 1, "rg1"), read_groups=("ID:rg1\tSM:sample1\tLB:lib1\tPL:ILLUMINA\tPU:unit1",))
    result = run_validator(BAM_VALIDATOR, *expected, stdin=valid)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "status\tPASS_COMPLETE" in result.stdout

    missing = sam(("read1", 1, "other"), read_groups=("ID:rg1\tSM:sample1\tLB:lib1\tPL:ILLUMINA\tPU:unit1",))
    result = run_validator(BAM_VALIDATOR, *expected, stdin=missing)
    assert result.returncode == 1, result.stdout + result.stderr

    merged_args = [
        "--expected", "rg1", "sample1", "lib1", "ILLUMINA", "unit1",
        "--expected", "rg2", "sample1", "lib1", "ILLUMINA", "unit2",
    ]
    merged = sam(
        ("read1", 1, "rg1"), ("read2", 5, "rg2"),
        read_groups=(
            "ID:rg1\tSM:sample1\tLB:lib1\tPL:ILLUMINA\tPU:unit1",
            "ID:rg2\tSM:sample1\tLB:lib1\tPL:ILLUMINA\tPU:unit2",
        ),
    )
    result = run_validator(MERGED_VALIDATOR, *merged_args, stdin=merged)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "read_groups_checked\t2" in result.stdout

    missing_group = sam(("read1", 1, "rg1"), read_groups=("ID:rg1\tSM:sample1\tLB:lib1\tPL:ILLUMINA\tPU:unit1",))
    result = run_validator(MERGED_VALIDATOR, *merged_args, stdin=missing_group)
    assert result.returncode == 1, result.stdout + result.stderr


def main():
    check_discovery()
    check_fastq_validators()
    check_amplicon_same_file()
    check_bam_validators()
    print(json.dumps({"discovery": "passed", "fastq_validators": "passed",
                      "bam_read_group_validators": "passed"}, sort_keys=True))


if __name__ == "__main__":
    main()
