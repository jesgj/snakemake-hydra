#!/usr/bin/env python3

import argparse
import sys


READ_GROUP_FIELDS = {
    "ID": "expected_id",
    "SM": "expected_sample",
    "LB": "expected_library",
    "PL": "expected_platform",
    "PU": "expected_platform_unit",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Validate one expected @RG header and every alignment RG tag in SAM "
            "data read from standard input."
        )
    )
    parser.add_argument("--expected-id", required=True)
    parser.add_argument("--expected-sample", required=True)
    parser.add_argument("--expected-library", required=True)
    parser.add_argument("--expected-platform", required=True)
    parser.add_argument("--expected-platform-unit", required=True)
    return parser.parse_args()


def fail(message):
    print(f"ERROR: {message}", file=sys.stderr)
    return 1


def main():
    args = parse_args()
    read_group_headers = []
    records_checked = 0
    missing_tags = 0
    mismatched_tags = 0

    for line_number, line in enumerate(sys.stdin, start=1):
        if line.startswith("@"):
            if line.startswith("@RG\t"):
                fields = {}
                for field in line.rstrip("\n").split("\t")[1:]:
                    if ":" in field:
                        key, value = field.split(":", 1)
                        fields[key] = value
                read_group_headers.append(fields)
            continue

        columns = line.rstrip("\n").split("\t")
        if len(columns) < 11:
            return fail(f"malformed SAM record at input line {line_number}")

        records_checked += 1
        record_read_groups = [
            field[5:] for field in columns[11:] if field.startswith("RG:Z:")
        ]
        if not record_read_groups:
            missing_tags += 1
        elif len(record_read_groups) != 1 or record_read_groups[0] != args.expected_id:
            mismatched_tags += 1

    if len(read_group_headers) != 1:
        return fail(
            f"expected exactly one @RG header, found {len(read_group_headers)}"
        )

    read_group_header = read_group_headers[0]
    mismatched_fields = []
    for sam_field, argument_name in READ_GROUP_FIELDS.items():
        expected_value = getattr(args, argument_name)
        actual_value = read_group_header.get(sam_field)
        if actual_value != expected_value:
            mismatched_fields.append(
                f"{sam_field} expected {expected_value!r}, found {actual_value!r}"
            )
    if mismatched_fields:
        return fail("read-group header mismatch: " + "; ".join(mismatched_fields))

    if missing_tags or mismatched_tags:
        return fail(
            "alignment read-group tags failed validation: "
            f"records={records_checked}, missing={missing_tags}, "
            f"mismatched={mismatched_tags}"
        )

    print("status\tPASS_COMPLETE")
    print(f"records_checked\t{records_checked}")
    print("missing_rg_tags\t0")
    print("mismatched_rg_tags\t0")
    for sam_field, argument_name in READ_GROUP_FIELDS.items():
        print(f"{sam_field}\t{getattr(args, argument_name)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
