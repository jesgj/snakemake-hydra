#!/usr/bin/env python3
"""Validate synchronization and provenance of one paired-end FASTQ pair."""

from __future__ import annotations

import argparse
import gzip
import re
import sys
from pathlib import Path
from typing import TextIO


class FastqValidationError(Exception):
    """Raised for malformed FASTQ input or invalid mate metadata."""


def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be at least 1")
    return parsed


def open_fastq(path: Path) -> TextIO:
    if not path.is_file():
        raise FastqValidationError(f"input does not exist or is not a file: {path}")
    try:
        if path.name.endswith(".gz"):
            return gzip.open(path, "rt", encoding="utf-8", errors="strict", newline="")
        return path.open("rt", encoding="utf-8", errors="strict", newline="")
    except OSError as exc:
        raise FastqValidationError(f"cannot open {path}: {exc}") from exc


def read_record(handle: TextIO, path: Path, record_number: int):
    lines: list[str] = []
    try:
        for line_index in range(4):
            line = handle.readline()
            if line == "":
                if line_index == 0:
                    return None
                raise FastqValidationError(
                    f"truncated FASTQ record {record_number} in {path}"
                )
            lines.append(line.rstrip("\r\n"))
    except (OSError, EOFError, UnicodeError) as exc:
        raise FastqValidationError(
            f"failed while reading record {record_number} from {path}: {exc}"
        ) from exc

    header, sequence, separator, quality = lines
    if not header.startswith("@"):
        raise FastqValidationError(
            f"record {record_number} in {path} has no '@' header"
        )
    if not separator.startswith("+"):
        raise FastqValidationError(
            f"record {record_number} in {path} has no '+' separator"
        )
    if len(sequence) != len(quality):
        raise FastqValidationError(
            f"record {record_number} in {path} has sequence/quality length mismatch"
        )
    return header, sequence, separator, quality


def normalize_header(header: str, expected_mate: int, record_number: int, path: Path):
    fields = header[1:].split()
    if not fields or not fields[0]:
        raise FastqValidationError(f"empty header at record {record_number} in {path}")

    cluster_id = fields[0]
    slash_match = re.search(r"/([12])$", cluster_id)
    if slash_match:
        observed = int(slash_match.group(1))
        if observed != expected_mate:
            raise FastqValidationError(
                f"record {record_number} in {path} is labelled read {observed}, "
                f"expected read {expected_mate}"
            )
        cluster_id = cluster_id[:-2]

    if len(fields) > 1:
        casava_match = re.match(r"^([12]):", fields[1])
        if casava_match:
            observed = int(casava_match.group(1))
            if observed != expected_mate:
                raise FastqValidationError(
                    f"record {record_number} in {path} is labelled read {observed}, "
                    f"expected read {expected_mate}"
                )

    parts = cluster_id.split(":")
    provenance = None
    if len(parts) >= 4:
        provenance = ":".join(parts[:4])
    return cluster_id, provenance


def format_values(values: set[str]) -> str:
    if not values:
        return "unparsed"
    return ",".join(sorted(values))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare normalized cluster identifiers in paired FASTQs and validate "
            "FASTQ structure, mate labels, provenance, and record counts."
        )
    )
    parser.add_argument("--r1", required=True, type=Path, help="R1 FASTQ[.gz]")
    parser.add_argument("--r2", required=True, type=Path, help="R2 FASTQ[.gz]")
    scope = parser.add_mutually_exclusive_group()
    scope.add_argument(
        "--records",
        type=positive_int,
        default=100,
        help="number of records to sample (default: 100)",
    )
    scope.add_argument(
        "--all",
        action="store_true",
        help="stream and validate every record in both files",
    )
    parser.add_argument(
        "--max-examples",
        type=positive_int,
        default=5,
        help="maximum mismatch examples to report (default: 5)",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    limit = None if args.all else args.records
    records_checked = 0
    mismatches = 0
    examples: list[str] = []
    r1_provenance: set[str] = set()
    r2_provenance: set[str] = set()
    complete = False

    try:
        with open_fastq(args.r1) as r1_handle, open_fastq(args.r2) as r2_handle:
            while limit is None or records_checked < limit:
                record_number = records_checked + 1
                r1_record = read_record(r1_handle, args.r1, record_number)
                r2_record = read_record(r2_handle, args.r2, record_number)

                if r1_record is None and r2_record is None:
                    complete = True
                    break
                if r1_record is None or r2_record is None:
                    mismatches += 1
                    missing = "R1" if r1_record is None else "R2"
                    examples.append(
                        f"record {record_number}: {missing} ended before its mate"
                    )
                    complete = True
                    break

                records_checked += 1
                r1_id, r1_run = normalize_header(
                    r1_record[0], 1, record_number, args.r1
                )
                r2_id, r2_run = normalize_header(
                    r2_record[0], 2, record_number, args.r2
                )
                if r1_run:
                    r1_provenance.add(r1_run)
                if r2_run:
                    r2_provenance.add(r2_run)

                if r1_id != r2_id:
                    mismatches += 1
                    if len(examples) < args.max_examples:
                        examples.append(
                            f"record {record_number}: R1={r1_id} R2={r2_id}"
                        )

            if limit is None and not complete:
                complete = True

    except FastqValidationError as exc:
        print("status\tERROR")
        print(f"error\t{exc}")
        return 2
    except (OSError, EOFError, UnicodeError) as exc:
        print("status\tERROR")
        print(f"error\t{exc}")
        return 2

    provenance_matches = r1_provenance == r2_provenance
    failed = mismatches > 0 or not provenance_matches
    if failed:
        status = "FAIL"
    elif complete:
        status = "PASS"
    else:
        status = "PASS_SAMPLED"

    print(f"status\t{status}")
    print(f"scope\t{'all_records' if args.all else f'first_{args.records}_records'}")
    print(f"records_checked\t{records_checked}")
    print(f"mismatches\t{mismatches}")
    print(f"r1_provenance\t{format_values(r1_provenance)}")
    print(f"r2_provenance\t{format_values(r2_provenance)}")
    print(f"provenance_match\t{'yes' if provenance_matches else 'no'}")
    for index, example in enumerate(examples, start=1):
        print(f"mismatch_example_{index}\t{example}")
    if status == "PASS_SAMPLED":
        print("note\tSampled validation does not prove equal full-file record counts")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
