#!/usr/bin/env python3
"""Validate amplicon FASTQs, remove configured primers, and run named FastQC outputs."""

import argparse
import gzip
import json
import shutil
import subprocess
import tempfile
from collections import Counter
from pathlib import Path

from validate_fastq_pairs import open_fastq, read_record, normalize_header


def audit(paths):
    lengths = [Counter() for _ in paths]
    count = 0
    from contextlib import ExitStack
    with ExitStack() as stack:
        handles = [stack.enter_context(open_fastq(Path(path))) for path in paths]
        while True:
            records = [read_record(h, Path(p), count + 1) for h, p in zip(handles, paths)]
            if all(record is None for record in records):
                break
            if any(record is None for record in records):
                raise ValueError(f"Unequal mate counts at record {count + 1}.")
            if len(records) == 2:
                names = [normalize_header(r[0], i + 1, count + 1, Path(p))[0] for i, (r, p) in enumerate(zip(records, paths))]
                if names[0] != names[1]:
                    raise ValueError(f"Unsynchronized mates at record {count + 1}.")
            for record, histogram in zip(records, lengths):
                if not record[1] or any(c not in "ACGTRYSWKMBDHVNacgtryswkmbdhvn" for c in record[1]):
                    raise ValueError(f"Empty or invalid nucleotide sequence at record {count + 1}.")
                if any(not 33 <= ord(c) <= 126 for c in record[3]):
                    raise ValueError(f"Invalid FASTQ quality character at record {count + 1}.")
                histogram[len(record[1])] += 1
            count += 1
    return count, [dict(sorted(h.items())) for h in lengths]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare", "qc"))
    parser.add_argument("--config", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()
    config = json.loads(Path(args.config).read_text())
    sample = config["samples_info"][args.sample]
    root = Path(config["output_dir"])
    raw = [sample[k] for k in ("R1", "R2") if sample.get(k)]
    prepared = [root / "primers" / f"{args.sample}_R{i}.fastq.gz" for i in range(1, len(raw) + 1)]
    if args.action == "qc":
        for stage, paths in (("raw_qc", raw), ("primer_qc", prepared)):
            outdir = root / stage / args.sample
            outdir.mkdir(parents=True, exist_ok=True)
            with tempfile.TemporaryDirectory() as tmp:
                links = []
                for i, path in enumerate(paths, 1):
                    link = Path(tmp) / f"{args.sample}_{stage}_R{i}.fastq{'.gz' if str(path).endswith('.gz') else ''}"
                    link.symlink_to(Path(path).resolve())
                    links.append(str(link))
                subprocess.run(["fastqc", "--outdir", str(outdir), "--threads", str(min(args.threads, len(links))), *links], check=True)
        return
    count, lengths = audit(raw)
    if not count:
        raise ValueError(f"Empty raw FASTQs for {args.sample}.")
    prepared[0].parent.mkdir(parents=True, exist_ok=True)
    primers = sample["primers"]
    if primers["mode"] == "remove":
        # Linked primers require the anchored 5' primer; optional 3' read-through
        # is clipped without requiring short reads to reach the opposite primer.
        def reverse_complement(primer):
            return primer.upper().translate(str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN"))[::-1]

        command = ["cutadapt", "--cores", str(args.threads), "--error-rate", str(primers["error_rate"]),
                   "--overlap", str(primers["min_overlap"]),
                   "--discard-untrimmed", "--json", str(root / "primers" / f"{args.sample}.cutadapt.json")]
        if primers.get("reverse"):
            command += ["-a", "^" + primers["forward"] + "..." + reverse_complement(primers["reverse"])]
        else:
            command += ["-g", "^" + primers["forward"]]
        command += ["-o", str(prepared[0])]
        if len(raw) == 2:
            command += ["-A", "^" + primers["reverse"] + "..." + reverse_complement(primers["forward"]), "--pair-filter=any", "-p", str(prepared[1])]
        with (root / "primers" / f"{args.sample}.cutadapt.txt").open("w") as report:
            subprocess.run(command + raw, check=True, stdout=report)
    else:
        (root / "primers" / f"{args.sample}.cutadapt.txt").write_text("Primer trimming not requested: inputs declared already_removed.\n")
        (root / "primers" / f"{args.sample}.cutadapt.json").write_text(
            json.dumps({"status": "not_requested", "primer_mode": "already_removed"}) + "\n")
        for source, destination in zip(raw, prepared):
            if str(source).endswith(".gz"):
                shutil.copyfile(source, destination)
            else:
                with open(source, "rb") as src, gzip.open(destination, "wb") as dst:
                    shutil.copyfileobj(src, dst)
    retained, retained_lengths = audit(prepared)
    report = {"sample": args.sample, "input": count, "primer_retained": retained,
              "raw_read_lengths": lengths, "primer_read_lengths": retained_lengths,
              "primer_mode": primers["mode"], "analysis_mode": sample["mode"], "batch": sample["batch"]}
    (root / "primers" / f"{args.sample}.reads.json").write_text(json.dumps(report, indent=2) + "\n")
    if not retained:
        raise ValueError(f"No primer-retained reads for {args.sample}; check primer sequences and orientation.")


if __name__ == "__main__":
    main()
