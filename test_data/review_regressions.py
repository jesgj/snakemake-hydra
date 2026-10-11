#!/usr/bin/env python3
"""Regression checks for discovery and empty-peak quantification using real fixture BAMs."""
import argparse
import copy
import csv
import json
import os
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
        section["multiqc_input_files"] = [str(summary)]
        section["multiqc_analysis_dirs"] = [section["featurecounts_dir"], section["qc_summary_dir"]]
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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atac-config", type=Path)
    parser.add_argument("--workdir", type=Path, default=REPO / "test_data/work/review_regressions")
    args = parser.parse_args()
    check_discovery()
    print("Passed discovery regressions", flush=True)
    if args.atac_config:
        check_atac(args.atac_config, args.workdir.resolve())
    args.workdir.mkdir(parents=True, exist_ok=True)
    (args.workdir / "validation.json").write_text(json.dumps(
        {"discovery": "passed", "atac_empty_peak_cases": "passed" if args.atac_config else "not_run"}, indent=2) + "\n")


if __name__ == "__main__":
    main()
