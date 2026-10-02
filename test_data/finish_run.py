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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("variant", choices=[p.stem for p in (ROOT / "configs").glob("*.yaml")])
    parser.add_argument("--cores", type=int, default=2)
    args = parser.parse_args()
    config = yaml.safe_load((ROOT / "configs" / (args.variant + ".yaml")).read_text())
    section = config[config["pipeline"]]
    work = ROOT / "work" / args.variant
    report = Path(section["multiqc_results_dir"]) / "multiqc_report.html"
    if not report.is_file() or report.stat().st_size == 0:
        raise ValueError(f"Missing or empty final report: {report}")
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
             "artifacts": artifacts, "output_bytes": sum(p["bytes"] for p in artifacts), "outputs_deleted": False}
    evidence = ROOT / "execution.json"
    results = json.loads(evidence.read_text()) if evidence.exists() else {"runs": []}
    results["runs"].append(entry)
    evidence.write_text(json.dumps(results, indent=2) + "\n")
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
