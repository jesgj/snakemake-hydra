#!/usr/bin/env python3
"""Focused tests for the dependency map and CI planner; no Pixi or data needed."""

import copy
import json
import subprocess
import tempfile
import unittest
from pathlib import Path

import plan


class SelectionTests(unittest.TestCase):
    def setUp(self):
        self.registry = json.loads(plan.REGISTRY.read_text())

    def check_selection(self, paths, expected, *, amplicon=False, full_reason=None):
        outputs, summary = plan.select(self.registry, paths, full_reason)
        self.assertEqual(expected, [job["variant"] for job in outputs["matrix"]["include"]])
        self.assertEqual(bool(expected), outputs["has_fixtures"])
        self.assertEqual(amplicon, outputs["run_amplicon"])
        return summary

    def test_registry_covers_workflows_and_fixtures(self):
        plan.validate_registry(self.registry)
        stale = copy.deepcopy(self.registry)
        del stale["pipelines"]["amplicon"]
        with self.assertRaisesRegex(ValueError, "Snakefile modules differ"):
            plan.validate_registry(stale)
        stale = copy.deepcopy(self.registry)
        stale["dependencies"]["workflows/rules/bam_qc.smk"].remove("wgbs")
        with self.assertRaisesRegex(ValueError, "Register wgbs as a consumer"):
            plan.validate_registry(stale)
        self.registry["pipelines"]["wgbs"]["fixtures"].pop()
        with self.assertRaisesRegex(ValueError, "must match test_data/prepare.py"):
            plan.validate_registry(self.registry)

    def test_each_pipeline_selects_its_variants(self):
        for pipeline, spec in self.registry["pipelines"].items():
            with self.subTest(pipeline=pipeline):
                self.check_selection([f"workflows/{pipeline}.smk"],
                                     [entry["variant"] for entry in spec.get("fixtures", [])],
                                     amplicon=spec.get("bundled_16s", False))

    def test_multi_pipeline_union_no_duplicates(self):
        summary = self.check_selection(
            ["README.md", "workflows/rnaseq.smk", "workflows/rules/rnaseq/alignment.smk", "workflows/wgbs.smk"],
            ["rnaseq", "wgbs", "wgbs_full"],
        )
        self.assertIn('"workflows/rules/rnaseq/alignment.smk" -> rnaseq', summary)

    def test_shared_rules_and_helpers(self):
        self.check_selection(["workflows/rules/bam_qc.smk"],
                             ["wgbs", "wgbs_full", "chipseq", "chipseq_se", "cutrun", "atacseq_macs3", "atacseq_genrich"])
        self.check_selection(["scripts/validate_fastq_pairs.py"], ["atacseq_macs3", "atacseq_genrich"], amplicon=True)
        self.check_selection(["test_data/prepare.py"], self.all_real_fixtures())

    def all_real_fixtures(self):
        return [job["variant"] for spec in self.registry["pipelines"].values() for job in spec.get("fixtures", [])]

    def test_global_unknown_and_full_runs(self):
        for path in ("Snakefile", "pixi.lock", "config/config.yaml", ".github/workflows/ci.yml",
                     "scripts/future.py", "workflows/unregistered.smk", "docs/future.py"):
            with self.subTest(path=path):
                self.check_selection([path], self.all_real_fixtures(), amplicon=True)
        self.check_selection(["README.md"], self.all_real_fixtures(), amplicon=True, full_reason="Full run")

    def test_docs_or_empty_diff(self):
        self.check_selection(["README.md", "docs/amplicon.md", "docs/pipeline_explainer.html"], [])
        self.check_selection([], [])


class DiffTests(unittest.TestCase):
    def test_ranges_include_earlier_commits_and_both_sides_of_renames(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)

            def git(*args):
                return plan.git(*args, root=root, stderr=subprocess.PIPE).decode().strip()

            def write(name):
                path = root / name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(f"fixture for {name}\n")

            def commit():
                git("add", "-A")
                git("commit", "-qm", "fixture change")
                return git("rev-parse", "HEAD")

            git("init", "-q", "-b", "main")
            git("config", "user.email", "planner@example.invalid")
            git("config", "user.name", "Planner")
            old = "workflows/rules/rnaseq/old rule.smk"
            new = "workflows/rules/chip_cr/new rule.smk"
            removed = "workflows/rules/atacseq/obsolete.smk"
            for path in (old, removed):
                write(path)
            start = commit()
            git("switch", "-qc", "topic")
            write("workflows/wgbs.smk")
            intermediate = commit()
            (root / new).parent.mkdir(parents=True)
            git("mv", old, new)
            (root / removed).unlink()
            head = commit()
            git("switch", "-q", "main")
            write("workflows/amplicon.smk")  # Base-only changes should not influence the PR plan.
            base = commit()

            diff_base, diff_head = plan.diff_range("pull_request", base, head, root)
            self.assertEqual((start, head), (diff_base, diff_head))
            self.assertCountEqual(plan.changed_paths(diff_base, diff_head, root),
                                  [old, new, removed, "workflows/wgbs.smk"])
            self.assertEqual((start, head), plan.diff_range("push", start, head, root))
            self.assertEqual((intermediate, head), plan.diff_range("workflow_dispatch", "", head, root))
            empty, _ = plan.diff_range("push", "0" * 40, start, root)
            self.assertCountEqual(plan.changed_paths(empty, start, root), [old, removed])


if __name__ == "__main__":
    unittest.main()
