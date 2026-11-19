#!/usr/bin/env python3
"""Check the change range and select the workflow fixtures CI needs to run."""

import argparse
import ast
import json
import os
import re
import runpy
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
REGISTRY = Path(__file__).with_name("pipelines.json")
DOC_FILES = {"README.md", "AGENTS.md", "backlog.md", "LICENSE", "test_data/README.md",
             "docs/pipeline_explainer.html", "logo-snake.svg", "logo.png", "pipeline_v2.svg"}


def git(*args, root=ROOT, **kwargs):
    return subprocess.check_output(["git", *args], cwd=root, **kwargs)


def owners(registry, path):
    explicit = registry["dependencies"].get(path, [])
    return set(explicit) | {
        pipeline for pipeline in registry["pipelines"]
        if path == f"workflows/{pipeline}.smk" or path.startswith(f"workflows/rules/{pipeline}/")
    }


def validate_registry(registry, root=ROOT):
    pipelines = registry["pipelines"]
    modules = set(re.findall(r"^module (\w+):", (root / "Snakefile").read_text(), re.MULTILINE))
    if not modules or modules != set(pipelines):
        raise ValueError(f"Registry and Snakefile modules differ: {sorted(modules ^ set(pipelines))}")
    for path, consumers in registry["dependencies"].items():
        if not consumers or set(consumers) - set(pipelines):
            raise ValueError(f"Invalid pipeline dependencies for {path}: {consumers}")

    fixtures = {}
    amplicon_jobs = []
    for pipeline, spec in pipelines.items():
        if spec.get("bundled_16s") is True and not spec.get("fixtures"):
            amplicon_jobs.append(pipeline)
        elif spec.get("fixtures") and "bundled_16s" not in spec:
            for fixture in spec["fixtures"]:
                variant, cores = fixture["variant"], fixture["cores"]
                if variant in fixtures or type(cores) is not int or cores < 1:
                    raise ValueError(f"Duplicate fixture or invalid core count: {fixture}")
                fixtures[variant] = pipeline
        else:
            raise ValueError(f"Pipeline {pipeline} must have at least one CI job")
        # Check literal includes recursively; nonliteral includes and script dependencies
        # still need explicit registry entries if they cross pipeline boundaries.
        pending = [root / "workflows" / f"{pipeline}.smk"]
        seen = set()
        while pending:
            workflow = pending.pop()
            if workflow in seen:
                continue
            seen.add(workflow)
            for include in re.findall(r'''^\s*include:\s*["']([^"']+)["']''', workflow.read_text(), re.MULTILINE):
                rule = (workflow.parent / include).resolve()
                path = rule.relative_to(root).as_posix()
                if pipeline not in owners(registry, path):
                    raise ValueError(f"Register {pipeline} as a consumer of {path}")
                pending.append(rule)

    if amplicon_jobs != ["amplicon"]:
        raise ValueError("The bundled 16S job must cover the amplicon module")
    if fixtures != runpy.run_path(str(root / "test_data/prepare.py"))["VARIANTS"]:
        raise ValueError("Registry fixture variants and pipeline owners must match test_data/prepare.py")


def select(registry, paths, full_reason=None):
    pipelines = registry["pipelines"]
    affected = set(pipelines) if full_reason else set()
    reasons = []
    for path in sorted(set(paths)):
        matches = owners(registry, path)
        if matches:
            reasons.append((path, ", ".join(sorted(matches))))
            affected.update(matches)
        elif path in DOC_FILES or (path.startswith("docs/") and path.endswith(".md")):
            reasons.append((path, "lightweight checks only"))
        else:
            reasons.append((path, "full suite (global or unmapped)"))
            affected.update(pipelines)

    matrix = {"include": [fixture for name, spec in pipelines.items() if name in affected
                          for fixture in spec.get("fixtures", [])]}
    outputs = {"matrix": matrix, "has_fixtures": bool(matrix["include"]),
               "run_amplicon": any(pipelines[name].get("bundled_16s") for name in affected)}
    jobs = [f"{fixture['variant']} ({fixture['cores']} cores)" for fixture in matrix["include"]]
    if outputs["run_amplicon"]:
        jobs.append("bundled 16S smoke")
    summary = ["## CI selection", "", f"Pipelines: {', '.join(sorted(affected)) or 'none'}",
               f"Fixture jobs: {', '.join(jobs) or 'none (lightweight checks only)'}"]
    if full_reason:
        summary.extend(("", full_reason))
    if reasons:
        # Indented code block keeps arbitrary Git filenames inert in the summary.
        summary.extend(("", *(f"    {json.dumps(path)} -> {reason}" for path, reason in reasons)))
    return outputs, "\n".join(summary) + "\n"


def diff_range(event, base, head, root=ROOT):
    if event == "pull_request":
        if not base:
            raise ValueError("A pull request needs a base commit")
        base = git("merge-base", base, head, root=root).decode().strip()
    elif event == "push" and not base:
        raise ValueError("A push needs its before SHA")
    elif event == "workflow_dispatch":
        try:
            base = git("rev-parse", "--verify", "--quiet", f"{head}^", root=root).decode().strip()
        except subprocess.CalledProcessError as exc:
            if exc.returncode != 1:
                raise
            base = "0" * 40
    if base == "0" * 40:
        base = git("hash-object", "-w", "-t", "tree", "--stdin", root=root, input=b"").decode().strip()
    return base, head


def changed_paths(base, head, root=ROOT):
    # Disable rename detection so moving a file counts as both a deletion and an addition.
    raw = git("diff", "--no-renames", "--name-only", "-z", base, head, "--", root=root)
    return [os.fsdecode(path) for path in raw.split(b"\0") if path]


def check_python(root=ROOT):
    for path in git("ls-files", "-z", "--", "*.py", root=root).split(b"\0"):
        if path:
            name = os.fsdecode(path)
            ast.parse((root / name).read_bytes(), filename=name)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--event", choices=("pull_request", "push", "workflow_dispatch"), required=True)
    parser.add_argument("--base", default="")
    parser.add_argument("--head", default="HEAD")
    args = parser.parse_args()
    registry = json.loads(REGISTRY.read_text())
    validate_registry(registry)
    check_python()
    base, head = diff_range(args.event, args.base, args.head)
    subprocess.run(["git", "diff", "--check", base, head, "--"], cwd=ROOT, check=True)
    full_reason = None if args.event == "pull_request" else f"Full suite for {args.event}."
    outputs, summary = select(registry, changed_paths(base, head), full_reason)
    print(summary, end="")
    if path := os.environ.get("GITHUB_STEP_SUMMARY"):
        with open(path, "a") as handle:
            handle.write(summary)
    if path := os.environ.get("GITHUB_OUTPUT"):
        with open(path, "a") as handle:
            for key, value in outputs.items():
                handle.write(f"{key}={json.dumps(value, separators=(',', ':'))}\n")


if __name__ == "__main__":
    main()
