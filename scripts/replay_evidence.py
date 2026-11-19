#!/usr/bin/env python3
"""Exact, schema-aware evidence for fixture replay runs.

Identity: git HEAD/dirty, pixi.lock digest, supplied invocation. ``record
--stage pre`` keeps raw byte digests of declared inputs; ``--stage post`` keeps
canonical results for the explicit ordered ``replay_evidence.roles`` list
(YAML/JSON; a mapping keeps order). Alignments are decoded samtools output
without @PG, text drops trailing whitespace, JSON sorts keys, roots inventory
sorted files; HTML and logs are excluded as volatile.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SECTION, IDENTITY = "replay_evidence", ("git", "pixi_lock_sha256", "invocation")
EXCLUDED = {".html", ".htm", ".log", ".err", ".out"}
ALIGNMENTS = {".bam", ".sam", ".cram"}
TEXTS = {".bed", ".broadpeak", ".csv", ".methylkit", ".narrowpeak", ".saf", ".tab", ".tsv", ".txt"}

def digest(data) -> str:
    return hashlib.sha256(data.encode("utf-8", "surrogateescape") if isinstance(data, str) else data).hexdigest()

def _alignment(path: Path) -> dict:
    exe = shutil.which("samtools")
    if not exe:
        raise ValueError("samtools is required for decoded alignment comparisons; use locked Pixi")
    sha = hashlib.sha256()
    records = 0
    # Coordinate ties can arrive in a different order even with the same inputs
    # and effective thread budget. Compare the full multiset of decoded records.
    with tempfile.TemporaryDirectory(prefix="hydra-replay-") as scratch:
        sam_records = Path(scratch) / "alignments.sam"
        with sam_records.open("wb") as output, subprocess.Popen(
            [exe, "view", "--no-PG", "-h", str(path)], stdout=subprocess.PIPE
        ) as proc:
            for line in proc.stdout:
                if line.startswith(b"@PG\t"):
                    continue
                if line.startswith(b"@"):
                    sha.update(line)
                else:
                    output.write(line)
                    records += 1
            if proc.wait():
                raise ValueError(f"samtools could not decode {path}")
        subprocess.run(["sort", "-T", scratch, "-o", str(sam_records), str(sam_records)],
                       env={**os.environ, "LC_ALL": "C"}, check=True)
        with sam_records.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                sha.update(chunk)
    return {"kind": "alignment-decoded", "sha256": sha.hexdigest(), "records": records}

def _root(path: Path, leaf) -> dict:
    files = [{"path": child.relative_to(path).as_posix(), **leaf(child)}
             for child in sorted(path.rglob("*"))
             if child.is_file() and child.suffix.lower() not in EXCLUDED]
    return {"kind": "root", "files": len(files),
            "sha256": digest(json.dumps(files, sort_keys=True)), "entries": files}

def fingerprint(path, *, decode: bool = True) -> dict:
    """Canonical, schema-specific result; ``decode=False`` keeps raw input bytes."""
    entry = Path(path)
    if not entry.exists():
        return {"kind": "missing"}
    if entry.is_dir():
        return _root(entry, lambda child: fingerprint(child, decode=decode))
    suffix = entry.suffix.lower()
    if suffix in EXCLUDED:
        return {"kind": "excluded"}
    if decode and suffix in ALIGNMENTS:
        return _alignment(entry)
    if decode and suffix in TEXTS:
        text = "\n".join(line.rstrip() for line in
                         entry.read_bytes().decode("utf-8", "surrogateescape").splitlines())
        return {"kind": "text", "sha256": digest(text), "lines": len(text.splitlines())}
    if decode and suffix == ".json":
        canonical = json.dumps(json.loads(entry.read_text(encoding="utf-8")), sort_keys=True, separators=(",", ":"))
        return {"kind": "json", "sha256": digest(canonical)}
    with entry.open("rb") as handle:
        sha = hashlib.file_digest(handle, "sha256").hexdigest()
    return {"kind": "bytes", "sha256": sha, "bytes": entry.stat().st_size}

def load_roles(config, section: str = SECTION) -> list[tuple[str, str]]:
    """Read the explicit ordered role/path list from a fixture YAML (or JSON) file."""
    path = Path(config)
    try:
        import yaml
    except ModuleNotFoundError:
        data = json.loads(path.read_text(encoding="utf-8"))
    else:
        data = yaml.safe_load(path.read_text(encoding="utf-8"))
    block = data.get(section) if isinstance(data, dict) else None
    if block is None:
        raise ValueError(f"{path}: no explicit {section!r} section")
    raw = block.get("roles") if isinstance(block, dict) else block
    if isinstance(raw, dict):
        raw = [{"role": role, "path": value} for role, value in raw.items()]
    if not isinstance(raw, list) or not raw:
        raise ValueError(f"{path}: {section}.roles must be a nonempty list or mapping")
    bad = [item for item in raw if not isinstance(item, dict) or set(item) != {"role", "path"}
           or not all(isinstance(item.get(key), str) and item[key] for key in ("role", "path"))]
    if bad:
        raise ValueError(f"{path}: each role needs exactly nonempty 'role' and 'path': {bad[0]!r}")
    roles = [(item["role"], item["path"]) for item in raw]
    if len({role for role, _ in roles}) != len(roles):
        raise ValueError(f"{path}: duplicate role names")
    return roles

def _git(repo: Path) -> dict:
    def git(*args):
        done = subprocess.run(["git", "-C", str(repo), *args], capture_output=True, text=True)
        return done.stdout.strip() if done.returncode == 0 else None
    status = git("status", "--porcelain", "--untracked-files=no")
    return {"head": git("rev-parse", "HEAD"), "dirty": None if status is None else bool(status)}

def build_evidence(roles, *, repo_root=ROOT, invocation=None, stage="post") -> dict:
    """Identity plus declared roles; ``pre`` records byte digests, ``post`` results."""
    repo_root = Path(repo_root)
    lock = repo_root / "pixi.lock"
    evidence = {"version": 1, "stage": stage, "git": _git(repo_root), "invocation": invocation, "roles": [],
                "created_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                "pixi_lock_sha256": digest(lock.read_bytes()) if lock.is_file() else None}
    for role, declared in roles:
        resolved = repo_root / declared  # pathlib keeps an absolute operand
        evidence["roles"].append({"role": role, "path": declared,
                                  "result" if stage == "post" else "input":
                                      fingerprint(resolved, decode=stage == "post")})
    return evidence

def _diff(role: str, base: dict, current: dict) -> list[str]:
    if base.get("kind") != current.get("kind"):
        return [f"{role}: canonical kind {base.get('kind')!r} != {current.get('kind')!r}"]
    if base.get("kind") != "root":
        return [f"{role}: {base.get('kind')} sha256 {base.get('sha256')} != {current.get('sha256')}"]
    b = {entry["path"]: entry for entry in base.get("entries", [])}
    c = {entry["path"]: entry for entry in current.get("entries", [])}
    problems = []
    for rel in sorted(set(b) | set(c)):
        if rel not in c or rel not in b:
            problems.append(f"{role}/{rel}: only in {'baseline' if rel not in c else 'current'}")
        elif b[rel] != c[rel]:
            problems.extend(_diff(f"{role}/{rel}", b[rel], c[rel]))
    return problems

def compare(baseline: dict, current: dict, *, identity: bool = True) -> list[str]:
    """Exact comparison; returns human-readable mismatches (empty means equal)."""
    problems = []
    for key in IDENTITY if identity else ():
        if baseline.get(key) != current.get(key):
            problems.append(f"identity.{key}: {baseline.get(key)!r} != {current.get(key)!r}")
    base, cur = ({entry["role"]: entry for entry in doc.get("roles", [])} for doc in (baseline, current))
    if list(base) != list(cur):
        problems.append(f"role order: {list(base)} != {list(cur)}")
    for role in [name for name in base if name in cur]:
        pb, pc = base[role].get("result", base[role].get("input")), cur[role].get("result", cur[role].get("input"))
        missing = [side for side, payload in (("baseline", pb), ("current", pc))
                   if payload is None or payload.get("kind") == "missing"]
        if missing:
            problems.append(f"{role}: no digest in {' and '.join(missing)}")
        elif pb != pc:
            problems.extend(_diff(role, pb, pc))
    return problems

def load_evidence(path, role: str = "result") -> dict:
    """Evidence JSON, or a plain file/output root fingerprinted as one role."""
    entry = Path(path)
    try:
        data = (json.loads(entry.read_text(encoding="utf-8"))
                if entry.is_file() and entry.suffix.lower() == ".json" else None)
    except ValueError:
        data = None
    if not (isinstance(data, dict) and isinstance(data.get("roles"), list)):
        data = {"stage": "path", "roles": [{"role": role, "path": str(entry), "result": fingerprint(entry)}]}
    return data

def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)
    record = commands.add_parser("record", help="write pre- or post-run evidence JSON")
    for flag, kwargs in (("--config", {"type": Path, "required": True}),
                         ("--output", {"type": Path, "required": True}),
                         ("--stage", {"choices": ("pre", "post"), "default": "post"}),
                         ("--invocation", {"default": None}),
                         ("--repo-root", {"type": Path, "default": ROOT})):
        record.add_argument(flag, **kwargs)
    match = commands.add_parser("compare", help="compare evidence JSON or two files/roots")
    for flag, kwargs in (("--baseline", {"type": Path, "required": True}),
                         ("--current", {"type": Path, "required": True}),
                         ("--skip-identity", {"action": "store_true"})):
        match.add_argument(flag, **kwargs)
    args = parser.parse_args(argv)
    try:
        if args.command == "record":
            roles = load_roles(args.config)
            evidence = build_evidence(roles, repo_root=args.repo_root, invocation=args.invocation, stage=args.stage)
            key = "input" if args.stage == "pre" else "result"
            for entry in evidence["roles"]:
                if entry[key]["kind"] in {"missing", "excluded"}:
                    raise ValueError(f"Cannot record {entry['role']}: {entry['path']} is missing or excluded")
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(json.dumps(evidence, indent=2) + "\n")
            print(f"Wrote {args.stage} evidence for {len(roles)} roles: {args.output}")
            return 0
        baseline, current = load_evidence(args.baseline), load_evidence(args.current)
        identity = not args.skip_identity and "git" in baseline and "git" in current
        problems = compare(baseline, current, identity=identity)
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2
    print("\n".join(problems) or "equal")
    return 0 if not problems else 1

if __name__ == "__main__":
    sys.exit(main())
