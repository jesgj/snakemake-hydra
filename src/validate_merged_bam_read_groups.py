#!/usr/bin/env python3

import argparse
import sys


READ_GROUP_FIELDS = ["ID", "SM", "LB", "PL", "PU"]


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Validate the expected @RG headers and every alignment RG tag in "
            "merged SAM data read from standard input."
        )
    )
    parser.add_argument(
        "--expected",
        action="append",
        nargs=5,
        metavar=("ID", "SM", "LB", "PL", "PU"),
        required=True,
    )
    return parser.parse_args()


def fail(message):
    print(f"ERROR: {message}", file=sys.stderr)
    return 1


def main():
    args = parse_args()
    expected = {
        values[0]: dict(zip(READ_GROUP_FIELDS, values)) for values in args.expected
    }
    if len(expected) != len(args.expected):
        return fail("expected read-group IDs must be unique")

    observed_headers = {}
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
                read_group_id = fields.get("ID")
                if not read_group_id:
                    return fail("encountered an @RG header without an ID")
                if read_group_id in observed_headers:
                    return fail(f"duplicate @RG header ID {read_group_id!r}")
                observed_headers[read_group_id] = fields
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
        elif len(record_read_groups) != 1 or record_read_groups[0] not in expected:
            mismatched_tags += 1

    if set(observed_headers) != set(expected):
        missing = sorted(set(expected) - set(observed_headers))
        unexpected = sorted(set(observed_headers) - set(expected))
        return fail(
            f"merged @RG IDs differ; missing={missing}, unexpected={unexpected}"
        )

    mismatched_headers = []
    for read_group_id, expected_fields in expected.items():
        observed_fields = observed_headers[read_group_id]
        for field, expected_value in expected_fields.items():
            observed_value = observed_fields.get(field)
            if observed_value != expected_value:
                mismatched_headers.append(
                    f"{read_group_id}:{field} expected {expected_value!r}, "
                    f"found {observed_value!r}"
                )
    if mismatched_headers:
        return fail("read-group header mismatch: " + "; ".join(mismatched_headers))

    if missing_tags or mismatched_tags:
        return fail(
            "alignment read-group tags failed validation: "
            f"records={records_checked}, missing={missing_tags}, "
            f"mismatched={mismatched_tags}"
        )

    print("status\tPASS_COMPLETE")
    print(f"records_checked\t{records_checked}")
    print(f"read_groups_checked\t{len(expected)}")
    print(f"read_group_ids\t{','.join(sorted(expected))}")
    print("missing_rg_tags\t0")
    print("mismatched_rg_tags\t0")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
