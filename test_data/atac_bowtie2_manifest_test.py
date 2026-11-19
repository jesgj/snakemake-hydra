#!/usr/bin/env python3
"""No-download checks for the ATAC prebuilt Bowtie2 index manifest helper."""
import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
HELPER = REPO / "src" / "atac_bowtie2_manifest.py"
SUFFIXES = {
    "standard": [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"],
    "large": [".1.bt2l", ".2.bt2l", ".3.bt2l", ".4.bt2l", ".rev.1.bt2l", ".rev.2.bt2l"],
}


def run(*args, check=True):
    result = subprocess.run(
        [sys.executable, "-B", str(HELPER), *map(str, args)],
        capture_output=True,
        text=True,
    )
    if check and result.returncode:
        raise AssertionError(f"{args} failed: {result.stderr}")
    return result


def make_index(root, index_format):
    prefix = root / "index" / "genome"
    prefix.parent.mkdir(parents=True, exist_ok=True)
    contents = {}
    for suffix in SUFFIXES[index_format]:
        path = Path(f"{prefix}{suffix}")
        path.write_bytes(f"synthetic {index_format} index{suffix}".encode())
        contents[suffix] = path.read_bytes()
    reference = root / "genome.fa"
    reference.write_text(">chr1\nACGTACGT\n")
    return prefix, reference, contents


def test_format(index_format):
    with tempfile.TemporaryDirectory(prefix=f"hydra-atac-{index_format}-") as directory:
        root = Path(directory)
        prefix, reference, contents = make_index(root, index_format)
        manifest_path = Path(f"{prefix}.manifest.json")
        run(
            "create", "--prefix", prefix, "--reference", reference,
            "--format", index_format, "--builder-version", "2.5.4",
        )
        manifest = json.loads(manifest_path.read_text())
        assert manifest["index_format"] == index_format
        assert manifest["builder"] == {"name": "bowtie2", "version": "2.5.4"}
        assert manifest["reference_sha256"] == hashlib.sha256(reference.read_bytes()).hexdigest()
        assert manifest["components"] == {
            suffix: hashlib.sha256(data).hexdigest() for suffix, data in contents.items()
        }
        run("validate", "--prefix", prefix, "--reference", reference, "--format", index_format)
        assert run("validate", "--prefix", prefix, "--reference", reference,
                   "--format", "large" if index_format == "standard" else "standard",
                   check=False).returncode != 0

        component = Path(f"{prefix}{SUFFIXES[index_format][2]}")
        component.write_bytes(b"tampered")
        result = run("validate", "--prefix", prefix, "--reference", reference, check=False)
        assert result.returncode != 0 and "component SHA-256 mismatch" in result.stderr
        component.write_bytes(contents[SUFFIXES[index_format][2]])

        reference.write_text(">chr1\nchanged\n")
        result = run("validate", "--prefix", prefix, "--reference", reference, check=False)
        assert result.returncode != 0 and "source FASTA SHA-256 mismatch" in result.stderr
        reference.write_text(">chr1\nACGTACGT\n")

        manifest["builder"]["version"] = " "
        manifest_path.write_text(json.dumps(manifest))
        result = run("validate", "--prefix", prefix, "--reference", reference, check=False)
        assert result.returncode != 0 and "builder.version" in result.stderr


def test_missing_and_ambiguous():
    with tempfile.TemporaryDirectory(prefix="hydra-atac-invalid-") as directory:
        root = Path(directory)
        prefix, reference, _ = make_index(root, "standard")
        assert run("create", "--prefix", prefix, "--reference", reference,
                   "--builder-version", "2.5.4", check=False).returncode == 0
        Path(f"{prefix}.1.bt2").unlink()
        assert run("validate", "--prefix", prefix, "--reference", reference,
                   check=False).returncode != 0
        assert run("validate", "--prefix", prefix, "--reference", reference,
                   "--manifest", root / "absent.json", check=False).returncode != 0

        Path(f"{prefix}.1.bt2").write_bytes(b"restored")
        for suffix in SUFFIXES["large"]:
            Path(f"{prefix}{suffix}").write_bytes(b"mixed")
        assert run("create", "--prefix", prefix, "--reference", reference,
                   "--builder-version", "2.5.4", check=False).returncode != 0


def main():
    test_format("standard")
    test_format("large")
    test_missing_and_ambiguous()
    print("Passed ATAC Bowtie2 manifest checks")


if __name__ == "__main__":
    main()
