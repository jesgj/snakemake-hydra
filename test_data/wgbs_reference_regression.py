#!/usr/bin/env python3
"""Exercise WGBS reference staging, complete index outputs, and deletion repair."""
import argparse
import copy
import gzip
import json
import tempfile
import subprocess
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parents[1]
INDEX_PARTS = ("1", "2", "3", "4", "rev.1", "rev.2")
CONVERSIONS = (("CT_conversion", "BS_CT"), ("GA_conversion", "BS_GA"))


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def write_fasta(path, sequences):
    with path.open("w") as handle:
        for name, sequence in sequences.items():
            handle.write(f">{name}\n")
            for start in range(0, len(sequence), 60):
                handle.write(sequence[start:start + 60] + "\n")


def write_fastq(path, header, sequence):
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as handle:
            handle.write(f"{header}\n{sequence}\n+\n{'I' * len(sequence)}\n".encode())


def expected_index_files(staging, large_index):
    extension = "bt2l" if large_index else "bt2"
    index = staging / "Bisulfite_Genome"
    files = [
        index / conversion / f"{prefix}.{part}.{extension}"
        for conversion, prefix in CONVERSIONS
        for part in INDEX_PARTS
    ]
    files.extend(
        index / conversion / f"genome_mfa.{conversion}.fa"
        for conversion, _ in CONVERSIONS
    )
    return files


def make_fixture(workdir):
    source = workdir / "source_reference"
    staging = workdir / "reference"
    trimmed = workdir / "trimmed_fastqs"
    raw = workdir / "raw_fastqs"
    for directory in (source, trimmed, raw):
        directory.mkdir(parents=True)

    sequence = ("ACGTTGCAAGTCCTAGGATCCGTAACCGGTTAGCATGCA" * 100)[:4000]
    write_fasta(source / "configured.fa", {"chr1": sequence})
    write_fasta(source / "extra.fa", {"extra": "TGCA" * 500})
    write_fastq(trimmed / "sample_R1.trimmed.fq.gz", "@read1/1", sequence[100:180])
    write_fastq(trimmed / "sample_R2.trimmed.fq.gz", "@read1/2", reverse_complement(sequence[300:380]))
    return source, staging, trimmed, raw


def make_workflow(workdir, source, staging, trimmed, raw, large_index):
    base = yaml.safe_load((REPO / "config/config.yaml").read_text())["wgbs"]
    section = copy.deepcopy(base)
    for key in section:
        if key.endswith("_dir"):
            section[key] = str(workdir / "results" / key)
    section.update(
        pipeline="wgbs",
        ref_genome=str(source / "configured.fa"),
        reference_dir=str(staging),
        raw_fastqs_dir=str(raw),
        trimmed_dir=str(trimmed),
        samples_info={
            "sample": {
                "R1": str(raw / "sample_R1.fastq.gz"),
                "R2": str(raw / "sample_R2.fastq.gz"),
                "type": "PE",
            }
        },
    )
    section["bismark"] = {"threads": 1, "parallel": 1, "large_index": large_index, "extra_args": ""}
    config = workdir / "config.yaml"
    config.write_text(yaml.safe_dump(section, sort_keys=False))
    snakefile = workdir / "Snakefile"
    snakefile.write_text(
        f"configfile: {str(config)!r}\n"
        f"include: {str(REPO / 'workflows/wgbs.smk')!r}\n"
    )
    return config, snakefile, section


def run_snakemake(workdir, snakefile, targets, label):
    log = workdir / "logs" / f"{label}.log"
    log.parent.mkdir(parents=True, exist_ok=True)
    command = [
        "pixi", "run", "--locked", "snakemake",
        "--snakefile", str(snakefile),
        "--cores", "1", "--jobs", "1", "--printshellcmds",
        "--latency-wait", "5", *map(str, targets),
    ]
    with log.open("w") as handle:
        handle.write(" ".join(command) + "\n\n")
        handle.flush()
        result = subprocess.run(command, cwd=REPO, stdout=handle, stderr=subprocess.STDOUT, text=True)
    output = log.read_text()
    if result.returncode:
        raise RuntimeError(f"Snakemake failed ({label}); see {log}")
    return output


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workdir", type=Path,
                        default=REPO / "test_data/work/wgbs_reference_regression")
    parser.add_argument("--large-index", action="store_true", help="Exercise the .bt2l index family")
    args = parser.parse_args()
    workdir = args.workdir.resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    workdir = Path(tempfile.mkdtemp(prefix="run-", dir=workdir))

    source, staging, trimmed, raw = make_fixture(workdir)
    source_before = {path.name: path.read_bytes() for path in source.iterdir()}
    _, snakefile, section = make_workflow(workdir, source, staging, trimmed, raw, args.large_index)

    report = Path(section["alignment_dir"]) / "sample_pe_report.txt"
    staged_fasta = staging / "genome.fa"
    staged_fai = staging / "genome.fa.fai"
    index_files = expected_index_files(staging, args.large_index)
    # Mirror the default workflow target, which retains each index component.
    targets = [report, staged_fai, *index_files]

    run_snakemake(workdir, snakefile, targets, "initial")
    require(report.is_file() and report.stat().st_size, "Bismark alignment report was not produced")
    require(staged_fasta.is_file() and staged_fasta.stat().st_size, "Staged reference FASTA is missing")
    require(staged_fasta.read_bytes() == source_before["configured.fa"], "Staging did not copy the configured FASTA")
    require(staged_fai.is_file() and staged_fai.stat().st_size, "Staged reference .fai is missing")
    require({path.name for path in source.iterdir()} == set(source_before),
            "Source reference directory changed during indexing or alignment")
    require({path.name: path.read_bytes() for path in source.iterdir()} == source_before,
            "Source reference contents changed during indexing or alignment")
    require(not (source / "Bisulfite_Genome").exists(), "Bismark wrote its index into the source reference directory")
    require(not (source / "configured.fa.fai").exists(), "samtools wrote .fai into the source reference directory")

    staging_entries = sorted(path.name for path in staging.iterdir())
    require(staging_entries == ["Bisulfite_Genome", "genome.fa", "genome.fa.fai"],
            f"Unexpected entries in workflow-owned staging: {staging_entries}")
    require(all(path.is_file() and path.stat().st_size for path in index_files),
            "One or more declared Bismark index outputs is missing or empty")
    for conversion, _ in CONVERSIONS:
        converted = staging / "Bisulfite_Genome" / conversion / f"genome_mfa.{conversion}.fa"
        headers = [line.strip() for line in converted.read_text().splitlines() if line.startswith(">")]
        require(len(headers) == 1 and headers[0].startswith(">chr1"),
                f"Bismark indexed an extra FASTA in {converted}: {headers}")

    log = workdir / "logs" / "up-to-date-dry-run.log"
    command = [
        "pixi", "run", "--locked", "snakemake", "--snakefile", str(snakefile),
        "--cores", "1", "--jobs", "1", "--dry-run", *map(str, targets),
    ]
    with log.open("w") as handle:
        handle.write(" ".join(command) + "\n\n")
        handle.flush()
        subprocess.run(command, cwd=REPO, stdout=handle, stderr=subprocess.STDOUT, text=True, check=True)
    require("nothing to be done" in log.read_text().lower(),
            "A completed synthetic reference DAG was not reproducible without reruns")

    deleted_index = index_files[3]
    deleted_index.unlink()
    deletion_output = run_snakemake(workdir, snakefile, targets, "index-deletion-repair")
    require(deleted_index.is_file() and deleted_index.stat().st_size, "Deleted Bismark index output was not rebuilt")
    require("rule bismark_genome_preparation:" in deletion_output and "rule bismark_alignment:" in deletion_output,
            "Deleting an index output did not schedule index and alignment rebuild")

    staged_fai.unlink()
    fai_output = run_snakemake(workdir, snakefile, targets, "fai-deletion-repair")
    require(staged_fai.is_file() and staged_fai.stat().st_size, "Deleted staged .fai was not rebuilt")
    require("rule index_bismark_reference_fasta:" in fai_output,
            "Deleting the staged .fai did not schedule its producer")
    require("rule bismark_genome_preparation:" not in fai_output,
            "Deleting .fai unnecessarily rebuilt the complete Bismark index")

    validation = {
        "source_reference_unchanged": True,
        "staging_entries": staging_entries,
        "configured_fasta_staged": True,
        "extra_fasta_indexed": False,
        "declared_index_files": [str(path.relative_to(staging)) for path in index_files],
        "index_deletion_rebuilt": True,
        "fai_deletion_rebuilt_without_index_rebuild": True,
        "up_to_date_rerun": True,
    }
    (workdir / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    print(f"Passed WGBS reference staging/index regression; validation: {workdir / 'validation.json'}")


if __name__ == "__main__":
    main()
