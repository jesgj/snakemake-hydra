#!/usr/bin/env python3
"""Check successful smoke-run artifacts, save evidence, and remove their outputs."""
import argparse
import hashlib
import json
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent


def verify_final_outputs(config):
    """Check selected assay outputs in addition to the shared BAM/report checks."""
    mode = config["pipeline"]
    section = config[mode]
    required = []
    for sample, info in section["samples_info"].items():
        if mode == "rnaseq":
            required.append(Path(section["kallisto_output_dir"]) / sample / "abundance.tsv")
        elif mode == "wgbs":
            required.append(Path(section["methyldackel_dir"]) / f"{sample}_CpG.methylKit")
        else:
            read_type = info["type"].lower()
            required.append(Path(section["bigwig_dir"]) / f"{sample}_{read_type}.bw")
        if mode == "atacseq":
            # An empty peak set is legitimate, but the declared file must exist.
            peak = Path(section["peaks_dir"]) / f"{sample}_peaks.narrowPeak"
            if not peak.is_file():
                raise ValueError(f"Missing peak output: {peak}")
    if mode == "atacseq":
        required.extend([
            Path(section["featurecounts_dir"]) / "consensus_peak_counts.txt",
            Path(section["qc_summary_dir"]) / "sample_qc_summary.tsv",
        ])
    if mode == "rnaseq" and section["gene_body_coverage"]["enabled"]:
        prefix = Path(section["gene_body_coverage_dir"]) / "all_samples"
        fmt = section["gene_body_coverage"].get("format", "pdf")
        required.extend([Path(f"{prefix}.geneBodyCoverage.txt"), Path(f"{prefix}.geneBodyCoverage.curves.{fmt}")])
        if len(section["samples_info"]) >= 3:
            required.append(Path(f"{prefix}.geneBodyCoverage.heatMap.{fmt}"))
    for path in required:
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"Missing or empty assay output: {path}")
    return [str(path.relative_to(ROOT)) for path in required]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("variant", choices=[p.stem for p in (ROOT / "configs").glob("*.yaml")])
    parser.add_argument("--cores", type=int, default=2)
    parser.add_argument("--keep-outputs", action="store_true", help="Retain verified results for CI artifacts")
    parser.add_argument("--evidence", type=Path, default=ROOT / "execution.json",
                        help="Destination for execution checks and artifact hashes")
    args = parser.parse_args()
    config = yaml.safe_load((ROOT / "configs" / (args.variant + ".yaml")).read_text())
    section = config[config["pipeline"]]
    work = ROOT / "work" / args.variant
    report = Path(section["multiqc_results_dir"]) / "multiqc_report.html"
    if not report.is_file() or report.stat().st_size == 0:
        raise ValueError(f"Missing or empty final report: {report}")
    assay_outputs = verify_final_outputs(config)
    files = sorted(p for p in work.rglob("*") if p.is_file())
    bams = [p for p in files if p.suffix == ".bam"]
    if not bams:
        raise ValueError("No BAMs found for execution verification")
    subprocess.run(["samtools", "quickcheck", "-v", *map(str, bams)], check=True)
    mapped = {}
    for bam in bams:
        count = int(subprocess.check_output(["samtools", "view", "-c", "-F", "4", str(bam)], text=True))
        if not count:
            raise ValueError(f"No mapped reads in {bam}; retain outputs for diagnosis")
        mapped[str(bam.relative_to(ROOT))] = count
    artifacts = []
    peak_counts = {}
    for path in files:
        if path.suffix == ".narrowPeak":
            count = 0
            for line in path.read_text().splitlines():
                if not line or line.startswith(("#", "track")):
                    continue
                fields = line.split("\t")
                if len(fields) != 10 or not 0 <= int(fields[1]) < int(fields[2]):
                    raise ValueError(f"Invalid narrowPeak record in {path}")
                count += 1
            peak_counts[str(path.relative_to(ROOT))] = count
        h = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                h.update(chunk)
        artifacts.append({"path": str(path.relative_to(ROOT)), "bytes": path.stat().st_size, "sha256": h.hexdigest()})
    entry = {"variant": args.variant, "pipeline": config["pipeline"], "checked_at_utc": datetime.now(timezone.utc).isoformat(),
             "scheduler_cores": args.cores, "scheduler_jobs": 1,
             "multiqc_report_verified": True, "bam_quickcheck": True, "mapped_alignment_records": mapped,
             "peak_counts": peak_counts,
             "assay_outputs_verified": assay_outputs,
             "artifacts": artifacts, "output_bytes": sum(p["bytes"] for p in artifacts), "outputs_deleted": False}
    evidence = args.evidence
    evidence.parent.mkdir(parents=True, exist_ok=True)
    results = json.loads(evidence.read_text()) if evidence.exists() else {"runs": []}
    results["runs"].append(entry)
    evidence.write_text(json.dumps(results, indent=2) + "\n")
    if args.keep_outputs:
        print(f"{args.variant}: verified {len(bams)} BAMs, assay outputs and MultiQC; retained results")
        return
    # Work is selected from known configs and is always inside this fixture tree.
    shutil.rmtree(work)
    if config["pipeline"] == "wgbs":
        ref = Path(section["ref_genome"])
        ref.relative_to(ROOT / "references")
        index = ref.parent / "Bisulfite_Genome"
        if index.is_dir():
            entry["reference_index_bytes_deleted"] = sum(p.stat().st_size for p in index.rglob("*") if p.is_file())
            shutil.rmtree(index)
        Path(str(ref) + ".fai").unlink(missing_ok=True)
    entry["outputs_deleted"] = True
    evidence.write_text(json.dumps(results, indent=2) + "\n")
    print(f"{args.variant}: verified {len(bams)} BAMs and MultiQC; deleted {len(files)} output files ({entry['output_bytes']:,} bytes)")


if __name__ == "__main__":
    main()
