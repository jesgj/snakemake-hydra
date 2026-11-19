#!/usr/bin/env python3
"""Independently rerun the small RNA-seq alignment and quantification targets."""
import argparse
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parents[1]
BASE = REPO / "test_data/work/rnaseq"
sys.path.insert(0, str(REPO / "scripts"))
from replay_evidence import build_evidence, compare


def roles(section):
    return [(f"{sample}/{kind}", str(path))
            for sample in section["samples_info"]
            for kind, path in (
                ("sorted_alignment", Path(section["alignment_dir"]) / f"{sample}_pe.sorted.bam"),
                ("abundance", Path(section["kallisto_output_dir"]) / sample / "abundance.tsv"),
            )]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workdir", type=Path, required=True, help="Where to retain replay evidence and logs")
    parser.add_argument("--cores", type=int, default=2)
    args = parser.parse_args()
    evidence_dir = args.workdir.resolve()
    evidence_dir.mkdir(parents=True, exist_ok=True)
    source = REPO / "test_data/configs/rnaseq.yaml"
    config = yaml.safe_load(source.read_text())
    section = config["rnaseq"]
    baseline = build_evidence(roles(section), invocation="verified RNA-seq fixture")
    with tempfile.TemporaryDirectory(prefix="hydra-rnaseq-replay-") as tmp:
        root = Path(tmp)
        for key, value in section.items():
            if isinstance(value, str) and Path(value).is_relative_to(BASE):
                section[key] = str(root / Path(value).relative_to(BASE))
        candidate_roles = roles(section)
        if any(not Path(path).is_relative_to(root) for _, path in candidate_roles):
            raise ValueError("Replay requires the RNA-seq fixture config produced by test_data/prepare.py")
        replay_config = root / "config.yaml"
        replay_config.write_text(yaml.safe_dump(config, sort_keys=False))
        command = ["pixi", "run", "--locked", "snakemake", "--cores", str(args.cores), "--jobs", "1",
                   "--replace-workflow-config", "--configfile", str(replay_config), "--",
                   *(path for _, path in candidate_roles)]
        with (evidence_dir / "replay.log").open("w") as log:
            subprocess.run(command, cwd=REPO, env={**os.environ, "HYDRA_CONFIG": str(replay_config)},
                           stdout=log, stderr=subprocess.STDOUT, check=True)
        replay = build_evidence(candidate_roles, invocation=command)
    problems = compare(baseline, replay, identity=False)
    (evidence_dir / "replay.json").write_text(json.dumps({
        "baseline": baseline, "replay": replay, "differences": problems,
    }, indent=2) + "\n")
    if problems:
        raise SystemExit("RNA-seq replay differences:\n" + "\n".join(problems))
    print("RNA-seq clean replay: decoded alignment records and abundance tables match")


if __name__ == "__main__":
    main()
