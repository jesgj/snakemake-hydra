#!/usr/bin/env python3
"""Create or validate provenance for a read-only ATAC-seq Bowtie2 index.

A prebuilt index is only reused when a manifest next to its prefix binds the six
component SHA-256 digests to the configured source FASTA digest and the Bowtie2
builder version. Create the manifest once from known index provenance:

    pixi run python3 src/atac_bowtie2_manifest.py create \
        --prefix /data/indexes/yeast/genome \
        --reference ref/genome.fa \
        --builder-version 2.5.4 \
        --format standard

The ATAC workflow expects ``<prefix>.manifest.json`` and validates it before
alignment. Validate manually with:

    pixi run python3 src/atac_bowtie2_manifest.py validate \
        --prefix /data/indexes/yeast/genome \
        --reference ref/genome.fa
"""

import argparse
import gzip
import hashlib
import json
import os
import re
import subprocess
import sys
from pathlib import Path

SCHEMA = "hydra-atac-bowtie2-index-manifest"
SCHEMA_VERSION = 1
COMPONENT_SUFFIXES = {
    "standard": [
        ".1.bt2",
        ".2.bt2",
        ".3.bt2",
        ".4.bt2",
        ".rev.1.bt2",
        ".rev.2.bt2",
    ],
    "large": [
        ".1.bt2l",
        ".2.bt2l",
        ".3.bt2l",
        ".4.bt2l",
        ".rev.1.bt2l",
        ".rev.2.bt2l",
    ],
}
SHA256_PATTERN = re.compile(r"[0-9a-f]{64}")


class ManifestError(Exception):
    """A manifest or index provenance check failed."""


def default_manifest_path(prefix):
    return f"{os.fspath(prefix)}.manifest.json"


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def component_paths(prefix, index_format):
    return {
        suffix: f"{os.fspath(prefix)}{suffix}"
        for suffix in COMPONENT_SUFFIXES[index_format]
    }


def resolve_index_format(prefix, expected=None):
    complete = [
        index_format
        for index_format, suffixes in COMPONENT_SUFFIXES.items()
        if all(
            os.path.isfile(f"{os.fspath(prefix)}{suffix}") for suffix in suffixes
        )
    ]
    if not complete:
        details = "; ".join(
            f"{index_format} missing "
            + ", ".join(
                f"{os.fspath(prefix)}{suffix}"
                for suffix in suffixes
                if not os.path.isfile(f"{os.fspath(prefix)}{suffix}")
            )
            for index_format, suffixes in COMPONENT_SUFFIXES.items()
        )
        raise ManifestError(
            f"no complete Bowtie2 index for prefix {os.fspath(prefix)!r} ({details})"
        )
    if len(complete) > 1:
        raise ManifestError(
            "ambiguous Bowtie2 index: both .bt2 and .bt2l component families are "
            f"complete under prefix {os.fspath(prefix)!r}"
        )
    detected = complete[0]
    if expected is not None and expected != detected:
        raise ManifestError(
            f"index format {expected!r} was requested, but prefix "
            f"{os.fspath(prefix)!r} contains a complete {detected!r} index"
        )
    return detected


def _require_digest(manifest, label, value):
    if not isinstance(value, str) or not SHA256_PATTERN.fullmatch(value):
        raise ManifestError(f"{label} must be a lowercase SHA-256 hex digest")
    return value


def _builder_fields(manifest):
    builder = manifest.get("builder")
    if not isinstance(builder, dict):
        raise ManifestError("manifest builder must be an object with name and version")
    name = builder.get("name")
    version = builder.get("version")
    if not isinstance(name, str) or not name.strip():
        raise ManifestError("manifest builder.name must be a non-empty string")
    if not isinstance(version, str) or not version.strip():
        raise ManifestError("manifest builder.version must be a non-empty string")
    return name.strip(), version.strip()


def build_manifest(prefix, reference, builder, builder_version, index_format=None):
    if not isinstance(builder, str) or not builder.strip():
        raise ManifestError("builder must be a non-empty string")
    if not isinstance(builder_version, str) or not builder_version.strip():
        raise ManifestError("builder version must be a non-empty string")
    reference = os.fspath(reference)
    if not os.path.isfile(reference):
        raise ManifestError(f"source FASTA does not exist: {reference}")

    detected_format = resolve_index_format(prefix, index_format)
    paths = component_paths(prefix, detected_format)
    return {
        "schema": SCHEMA,
        "schema_version": SCHEMA_VERSION,
        "index_format": detected_format,
        "reference_sha256": sha256_file(reference),
        "components": {
            suffix: sha256_file(path) for suffix, path in paths.items()
        },
        "builder": {"name": builder.strip(), "version": builder_version.strip()},
    }


def write_manifest(manifest, output):
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + ".tmp")
    try:
        temporary.write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        os.replace(temporary, output)
    finally:
        if temporary.exists():
            temporary.unlink()


def read_manifest(manifest_path):
    try:
        with open(manifest_path, encoding="utf-8") as handle:
            manifest = json.load(handle)
    except FileNotFoundError as exc:
        raise ManifestError(f"manifest does not exist: {manifest_path}") from exc
    except (OSError, json.JSONDecodeError) as exc:
        raise ManifestError(f"cannot read manifest {manifest_path}: {exc}") from exc
    if not isinstance(manifest, dict):
        raise ManifestError("manifest must contain a JSON object")
    return manifest


def validate_manifest(manifest_path, prefix, reference, expected_format=None):
    manifest = read_manifest(manifest_path)
    if manifest.get("schema") != SCHEMA:
        raise ManifestError(
            f"unsupported manifest schema {manifest.get('schema')!r}; expected {SCHEMA!r}"
        )
    if manifest.get("schema_version") != SCHEMA_VERSION:
        raise ManifestError(
            "unsupported manifest schema_version "
            f"{manifest.get('schema_version')!r}; expected {SCHEMA_VERSION}"
        )

    index_format = manifest.get("index_format")
    if index_format not in COMPONENT_SUFFIXES:
        raise ManifestError(
            "manifest index_format must be one of "
            f"{sorted(COMPONENT_SUFFIXES)}, found {index_format!r}"
        )
    resolve_index_format(prefix, index_format)
    if expected_format is not None and expected_format != index_format:
        raise ManifestError(
            f"manifest index_format is {index_format!r}, but the workflow expects "
            f"{expected_format!r}"
        )

    components = manifest.get("components")
    if not isinstance(components, dict):
        raise ManifestError("manifest components must be an object")
    expected_suffixes = set(COMPONENT_SUFFIXES[index_format])
    missing = sorted(expected_suffixes - set(components))
    unexpected = sorted(set(components) - expected_suffixes)
    if missing or unexpected:
        raise ManifestError(
            f"manifest component suffixes do not match {index_format!r}; "
            f"missing={missing}, unexpected={unexpected}"
        )
    for suffix in sorted(expected_suffixes):
        _require_digest(
            manifest, f"components[{suffix!r}]", components.get(suffix)
        )

    reference = os.fspath(reference)
    if not os.path.isfile(reference):
        raise ManifestError(f"source FASTA does not exist: {reference}")
    reference_sha256 = _require_digest(
        manifest, "reference_sha256", manifest.get("reference_sha256")
    )
    actual_reference_sha256 = sha256_file(reference)
    if actual_reference_sha256 != reference_sha256:
        raise ManifestError(
            "configured source FASTA SHA-256 mismatch for "
            f"{reference}: manifest={reference_sha256}, actual={actual_reference_sha256}"
        )

    builder_name, builder_version = _builder_fields(manifest)

    missing_files = []
    mismatched = []
    for suffix, path in component_paths(prefix, index_format).items():
        if not os.path.isfile(path):
            missing_files.append(path)
            continue
        actual = sha256_file(path)
        if actual != components[suffix]:
            mismatched.append(
                f"{path} (manifest={components[suffix]}, actual={actual})"
            )
    if missing_files:
        raise ManifestError(
            "manifest lists components that are missing: " + ", ".join(missing_files)
        )
    if mismatched:
        raise ManifestError(
            "Bowtie2 index component SHA-256 mismatch: " + "; ".join(mismatched)
        )

    return {
        "index_format": index_format,
        "reference_sha256": reference_sha256,
        "component_count": len(expected_suffixes),
        "builder": f"{builder_name} {builder_version}",
    }


def validate_dictionary(prefix, reference):
    """Check names, order, and lengths independently of the supplied manifest."""
    observed = []
    output = subprocess.check_output(["bowtie2-inspect", "-s", prefix], text=True)
    for line in output.splitlines():
        if line.startswith("Sequence-"):
            fields = line.split("\t", 1)[1]
            name, length = fields.rsplit("\t", 1)
            observed.append((name.split()[0], int(length)))
    expected = []
    name, length = None, 0
    opener = gzip.open if str(reference).lower().endswith(".gz") else open
    with opener(reference, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    expected.append((name, length))
                name, length = line[1:].split()[0], 0
            else:
                length += len("".join(line.split()))
    if name is not None:
        expected.append((name, length))
    if not expected or expected != observed:
        raise ManifestError("Bowtie2 index contig names/order/lengths differ from the configured FASTA")


def print_summary(summary):
    print("status\tPASS_COMPLETE")
    print(f"index_format\t{summary['index_format']}")
    print(f"reference_sha256\t{summary['reference_sha256']}")
    print(f"components_checked\t{summary['component_count']}")
    print(f"builder\t{summary['builder']}")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    create = subparsers.add_parser(
        "create", help="write a manifest from a known prebuilt index and FASTA"
    )
    create.add_argument("--prefix", required=True)
    create.add_argument("--reference", "--ref-genome", dest="reference", required=True)
    create.add_argument("--format", dest="index_format", choices=COMPONENT_SUFFIXES)
    create.add_argument("--builder", default="bowtie2")
    create.add_argument("--builder-version", required=True)
    create.add_argument("--output")

    validate = subparsers.add_parser(
        "validate", help="verify a manifest against the index and configured FASTA"
    )
    validate.add_argument("--manifest")
    validate.add_argument("--prefix", required=True)
    validate.add_argument("--reference", "--ref-genome", dest="reference", required=True)
    validate.add_argument("--format", dest="index_format", choices=COMPONENT_SUFFIXES)
    validate.add_argument("--check-dictionary", action="store_true", help="also inspect contigs with Bowtie2")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    try:
        if args.command == "create":
            output = args.output or default_manifest_path(args.prefix)
            manifest = build_manifest(
                args.prefix,
                args.reference,
                args.builder,
                args.builder_version,
                args.index_format,
            )
            write_manifest(manifest, output)
            print(f"wrote_manifest\t{output}")
            print_summary(
                validate_manifest(
                    output, args.prefix, args.reference, manifest["index_format"]
                )
            )
        else:
            manifest_path = args.manifest or default_manifest_path(args.prefix)
            summary = validate_manifest(manifest_path, args.prefix, args.reference, args.index_format)
            if args.check_dictionary:
                validate_dictionary(args.prefix, args.reference)
            print_summary(summary)
    except (ManifestError, OSError, subprocess.CalledProcessError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
