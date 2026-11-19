#!/usr/bin/env python3
"""Regression checks for discovery and empty-peak quantification using real fixture BAMs."""
import argparse
import copy
import csv
import json
import os
import shutil
import subprocess
import sys
import tempfile
import warnings
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from utils import discover_samples


def check_discovery():
    with tempfile.TemporaryDirectory(prefix="hydra-discovery-") as directory:
        root = Path(directory)
        for name in ("d_R1.fastq.gz", "d_R1_001.fastq.gz", "d_R2.fastq.gz", "d_R2_001.fastq.gz"):
            (root / name).touch()
        try:
            discover_samples(root)
        except ValueError as exc:
            assert "Ambiguous" in str(exc)
        else:
            raise AssertionError("Duplicate mates must fail rather than create phantom samples")
        for path in root.iterdir():
            path.unlink()
        for name in ("a_R1_001.fq.gz", "a_R2_001.fq.gz", "b.fastq.gz", "c.1.fastq.gz", "c.2.fastq.gz"):
            (root / name).touch()
        (root / "ignored_R1.fastq.gz").mkdir()
        (root / "ignored_R1.fastq.gz.backup").touch()
        found = discover_samples(root)
        assert list(found) == ["a", "c", "b"]
        assert {s: i["type"] for s, i in found.items()} == {"a": "PE", "c": "PE", "b": "SE"}
        (root / "orphan_R2.fastq.gz").touch()
        try:
            discover_samples(root)
        except ValueError as exc:
            assert "Orphan R2" in str(exc)
        else:
            raise AssertionError("Orphan R2 must fail")
        (root / "orphan_R2.fastq.gz").unlink()
        for lane in ("L001", "L002"):
            for mate in (1, 2):
                (root / f"lane_{lane}_R{mate}.fastq.gz").touch()
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            found = discover_samples(root)
        assert "lane_L001" in found and "lane_L002" in found
        assert any("separate samples" in str(w.message) for w in caught)
    with tempfile.TemporaryDirectory(prefix="hydra-bare-mates-") as directory:
        root = Path(directory)
        for name in ("sample1.fastq.gz", "sample2.fastq.gz"):
            (root / name).touch()
        found = discover_samples(root)
        assert list(found) == ["sample1", "sample2"]
        assert {sample: info["type"] for sample, info in found.items()} == {
            "sample1": "SE", "sample2": "SE"
        }


def check_atac(config_path, work):
    original = yaml.safe_load(config_path.read_text())["atacseq"]
    samples = list(original["samples_info"])
    assert len(samples) == 2, "Use the two-sample ATAC fixture"
    # Synthetic peak intervals exercise edge cases on the real, already filtered BAMs.
    cases = {"all_empty": ("", "", 1, False),
             "one_empty": ("I\t100\t500\n", "", 1, True),
             "insufficient_support": ("I\t100\t500\n", "II\t100\t500\n", 2, False)}
    for name, (first, second, support, nonempty) in cases.items():
        root = work / name
        root.mkdir(parents=True, exist_ok=True)
        section = copy.deepcopy(original)
        # Spaces exercise the shared MultiQC path quoting as well.
        for key in ("peaks_dir", "consensus_peaks_dir", "featurecounts_dir", "qc_summary_dir", "multiqc_results_dir"):
            section[key] = str(root / key.replace("_", " "))
        section.update(pipeline="atacseq", pipeline_name="atacseq", consensus_peaks={"minimum_support": support})
        peakdir = Path(section["peaks_dir"])
        peakdir.mkdir(exist_ok=True)
        for sample, content in zip(samples, (first, second)):
            (peakdir / f"{sample}_peaks.narrowPeak").write_text(content)
        summary = Path(section["qc_summary_dir"]) / "sample_qc_summary.tsv"
        section["multiqc_input_files"] = [
            str(Path(section["qc_summary_dir"]) / "atac_qc_summary_mqc.tsv"),
            str(Path(section["featurecounts_dir"]) / "consensus_peak_counts.txt.summary"),
        ]
        cfg = root / "config.yaml"
        cfg.write_text(yaml.safe_dump(section))
        snakefile = root / "Snakefile"
        snakefile.write_text(
            f"configfile: {str(cfg)!r}\n"
            "SAMPLES = list(config['samples_info'])\n"
            f"include: {str(REPO / 'workflows/rules/atacseq/quantification.smk')!r}\n"
            f"include: {str(REPO / 'workflows/rules/multiqc.smk')!r}\n"
            "rule all:\n    default_target: True\n    input: rules.multiqc.output\n"
        )
        with (root / "run.log").open("w") as log:
            subprocess.run(["pixi", "run", "--locked", "snakemake", "--snakefile", str(snakefile),
                            "--cores", "2", "--jobs", "1", "--printshellcmds"], cwd=REPO,
                           env={**os.environ, "HYDRA_CONFIG": str(cfg)}, stdout=log, stderr=subprocess.STDOUT, check=True)
        bed = Path(section["consensus_peaks_dir"]) / "consensus_peaks.bed"
        assert bool(bed.read_text()) == nonempty
        matrix = Path(section["featurecounts_dir"]) / "consensus_peak_counts.txt"
        rows = [r for r in matrix.read_text().splitlines() if not r.startswith("#")]
        assert (len(rows) > 1) == nonempty
        assert len(rows[0].split("\t")) == 6 + len(samples)
        qc = list(csv.DictReader(summary.open(), delimiter="\t"))
        assert len(qc) == len(samples)
        if not nonempty:
            assert all(r["fragments_in_peaks"] == "0" and r["frip_pct"] == "NA" for r in qc)
        assert (Path(section["multiqc_results_dir"]) / "multiqc_report.html").stat().st_size > 1000
        print(f"Passed ATAC {name}", flush=True)


def _write_fastp_json(path, sample):
    payload = {
        "command": f"fastp --in1 {sample}.fastq.gz --out1 {sample}.trimmed.fastq.gz",
        "summary": {
            "before_filtering": {
                "total_reads": 1000,
                "total_bases": 100000,
                "q20_bases": 90000,
                "q30_bases": 80000,
                "q20_rate": 0.9,
                "q30_rate": 0.8,
                "read1_mean_length": 100,
                "gc_content": 0.5,
            },
            "after_filtering": {
                "total_reads": 900,
                "total_bases": 90000,
                "q20_bases": 85000,
                "q30_bases": 80000,
                "q20_rate": 0.94,
                "q30_rate": 0.88,
                "read1_mean_length": 100,
                "gc_content": 0.5,
            },
            "filtering_result": {
                "passed_filter_reads": 900,
                "low_quality_reads": 50,
                "too_many_N_reads": 0,
                "too_short_reads": 0,
                "too_long_reads": 0,
            },
            "duplication": {"rate": 0.1},
            "insert_size": {"peak": 200},
        },
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload))


def _run_multiqc_rule(root, sources, result_dir):
    log_dir = root / "logs" / "multiqc_regression"
    (log_dir / "multiqc").mkdir(parents=True, exist_ok=True)
    config_path = root / "config.yaml"
    config_path.write_text(yaml.safe_dump({
        "pipeline": "multiqc_regression",
        "pipeline_name": "multiqc_regression",
        "log_dir": str(log_dir),
        "multiqc_results_dir": str(result_dir),
        "multiqc_input_files": [str(source) for source in sources],
        # Deliberately broad: exact-source selection must ignore stale files here.
        "multiqc_analysis_dirs": [str(root / "inputs")],
    }))
    snakefile = root / "Snakefile"
    snakefile.write_text(
        f"configfile: {str(config_path)!r}\n"
        f"include: {str(REPO / 'workflows/rules/multiqc.smk')!r}\n"
        "rule all:\n    default_target: True\n    input: rules.multiqc.output\n"
    )
    with (root / f"{sources[0].name}.log").open("w") as log:
        subprocess.run(
            ["pixi", "run", "--locked", "snakemake", "--snakefile", str(snakefile),
             "--cores", "2", "--jobs", "1", "--printshellcmds"],
            cwd=REPO,
            env={**os.environ, "HYDRA_CONFIG": str(config_path)},
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )
    data_dir = result_dir / "multiqc_report_data"
    with (data_dir / "multiqc_fastp.txt").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    sources_text = (data_dir / "multiqc_sources.txt").read_text()
    return {row["Sample"] for row in rows}, sources_text


def check_multiqc_exact_sources(work):
    root = work / "multiqc_exact_sources"
    if root.exists():
        shutil.rmtree(root)
    inputs = root / "inputs"
    current = inputs / "current.fastp.json"
    stale = inputs / "stale.fastp.json"
    _write_fastp_json(current, "current")
    _write_fastp_json(stale, "stale")

    report_a = root / "report_current"
    samples, sources_text = _run_multiqc_rule(root, [current], report_a)
    assert samples == {"current"}, samples
    assert str(current) in sources_text and str(stale) not in sources_text

    # Repeating the same cohort must not pick up the sibling stale report.
    samples_again, sources_again = _run_multiqc_rule(root, [current], report_a)
    assert samples_again == {"current"}, samples_again
    assert sources_again == sources_text

    # Replacing the cohort in the same results directory must rebuild the
    # report and drop the previous sample even though its input still exists.
    samples, sources_text = _run_multiqc_rule(root, [stale], report_a)
    assert samples == {"stale"}, samples
    assert str(stale) in sources_text and str(current) not in sources_text

    shutil.rmtree(report_a / "multiqc_report_data")
    rebuilt_samples, rebuilt_sources = _run_multiqc_rule(root, [stale], report_a)
    assert rebuilt_samples == samples and rebuilt_sources == sources_text
    print("Passed MultiQC exact-source, stale-cohort and report-data rebuild regressions", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atac-config", type=Path)
    parser.add_argument("--workdir", type=Path, default=REPO / "test_data/work/review_regressions")
    args = parser.parse_args()
    check_discovery()
    print("Passed discovery regressions", flush=True)
    args.workdir.mkdir(parents=True, exist_ok=True)
    check_multiqc_exact_sources(args.workdir.resolve())
    if args.atac_config:
        check_atac(args.atac_config, args.workdir.resolve())
    (args.workdir / "validation.json").write_text(json.dumps(
        {
            "discovery": "passed",
            "multiqc_exact_sources": "passed",
            "atac_empty_peak_cases": "passed" if args.atac_config else "not_run",
        },
        indent=2,
    ) + "\n")


if __name__ == "__main__":
    main()
