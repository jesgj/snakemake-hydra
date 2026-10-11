#!/usr/bin/env python3
"""Write a zero-row count matrix when no consensus peaks can be quantified."""
import argparse
import csv
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("bams", nargs="+")
    args = parser.parse_args()
    with open(args.matrix, "w", newline="") as handle:
        handle.write("# No consensus peaks; featureCounts skipped.\n")
        csv.writer(handle, delimiter="\t", lineterminator="\n").writerow(
            ["Geneid", "Chr", "Start", "End", "Strand", "Length", *args.bams]
        )
    # No fragments were processed by featureCounts, so FRiP is unavailable.
    with open(args.summary, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["Status", *args.bams])
        writer.writerow(["Assigned", *([0] * len(args.bams))])
    print("No consensus peaks: counting skipped; zero-row matrix written; FRiP is NA.", file=sys.stderr)


if __name__ == "__main__":
    main()
