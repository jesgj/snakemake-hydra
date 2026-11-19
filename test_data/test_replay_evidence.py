#!/usr/bin/env python3
"""Focused stdlib tests for scripts/replay_evidence.py; no Pixi or data needed.

Run from the repository root:

    python3 -B test_data/test_replay_evidence.py
"""
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
import replay_evidence as evidence


class FingerprintTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def write(self, name, content, mode="w"):
        path = self.root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        if mode == "w":
            path.write_text(content)
        else:
            path.write_bytes(content)
        return path

    def test_text_drops_trailing_whitespace_only(self):
        a = self.write("a.tsv", "gene\tcount  \nA\t3\n")
        b = self.write("b.tsv", "gene\tcount\nA\t3\r\n")
        self.assertEqual(evidence.fingerprint(a), evidence.fingerprint(b))
        c = self.write("c.tsv", "gene\tcount\nA\t4\n")
        self.assertNotEqual(evidence.fingerprint(a), evidence.fingerprint(c))

    def test_html_and_logs_are_excluded(self):
        html = self.write("report.html", "<html>generated 12:00</html>")
        log = self.write("run.log", "2026-01-01 12:00 start\n")
        for path in (html, log):
            self.assertEqual({"kind": "excluded"}, evidence.fingerprint(path))

    def test_json_ignores_key_order_but_not_values(self):
        a = self.write("a.json", '{"x": 1, "y": [2]}')
        b = self.write("b.json", '{"y": [2], "x": 1}')
        self.assertEqual(evidence.fingerprint(a), evidence.fingerprint(b))
        c = self.write("c.json", '{"x": 2, "y": [2]}')
        self.assertNotEqual(evidence.fingerprint(a), evidence.fingerprint(c))

    def test_missing_and_raw_input_digest(self):
        self.assertEqual({"kind": "missing"}, evidence.fingerprint(self.root / "absent.bam"))
        raw = self.write("in.fastq", b"@r\nACGT\n+\nIIII\n", mode="b")
        decoded = evidence.fingerprint(raw, decode=False)
        self.assertEqual("bytes", decoded["kind"])
        self.assertEqual(evidence.digest(raw.read_bytes()), decoded["sha256"])

    def test_root_is_sorted_inventory_and_skips_html(self):
        self.write("out/b.tsv", "b\n1\n")
        self.write("out/a.tsv", "a\n1\n")
        self.write("out/multiqc_report.html", "<html>12:00</html>")
        self.write("out/run.log", "12:00 start\n")
        first = evidence.fingerprint(self.root / "out")
        self.assertEqual(2, first["files"])
        self.assertEqual(["a.tsv", "b.tsv"], [entry["path"] for entry in first["entries"]])
        self.write("out/run.log", "13:00 start\n")
        self.assertEqual(first, evidence.fingerprint(self.root / "out"))
        self.write("out/a.tsv", "a\n2\n")
        self.assertNotEqual(first, evidence.fingerprint(self.root / "out"))

    @unittest.skipUnless(evidence.shutil.which("samtools"), "samtools not on PATH")
    def test_alignment_ignores_pg_lines_and_decodes(self):
        header = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100\n"
        record = "r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:x\n"
        sam = self.write("a.sam", header + "@PG\tID:bwa\tCL:run --out /abs/a.sam\n" + record)
        pg = self.write("b.sam", header + "@PG\tID:bwa\tCL:run --out /abs/b.sam\n" + record)
        decoded = evidence.fingerprint(sam)
        self.assertEqual("alignment-decoded", decoded["kind"])
        self.assertEqual(decoded, evidence.fingerprint(pg))
        changed = self.write("c.sam", header + record.replace("ACGT", "ACGA"))
        self.assertNotEqual(decoded, evidence.fingerprint(changed))


class EvidenceTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        (self.root / "pixi.lock").write_text("lock-v1\n")
        (self.root / "out").mkdir()
        (self.root / "out" / "counts.tsv").write_text("gene\tcount\nA\t3\n")

    def tearDown(self):
        self.tmp.cleanup()

    def evidence(self, *, stage="post", invocation="run --cores 2", roles=None):
        roles = roles if roles is not None else [("counts", "out/counts.tsv")]
        return evidence.build_evidence(roles, repo_root=self.root, invocation=invocation, stage=stage)

    def test_pre_records_byte_digests_and_identity(self):
        doc = self.evidence(stage="pre")
        entry = doc["roles"][0]
        self.assertEqual("counts", entry["role"])
        self.assertEqual(evidence.digest((self.root / "out/counts.tsv").read_bytes()), entry["input"]["sha256"])
        self.assertEqual(evidence.digest(b"lock-v1\n"), doc["pixi_lock_sha256"])
        self.assertEqual("run --cores 2", doc["invocation"])
        self.assertNotIn("result", entry)

    def test_compare_identity_is_exact_and_ignores_created_at(self):
        self.assertEqual([], evidence.compare(self.evidence(), self.evidence()))
        self.assertEqual(["identity.invocation: 'run --cores 2' != 'run --cores 4'"],
                         evidence.compare(self.evidence(), self.evidence(invocation="run --cores 4")))
        self.assertEqual([], evidence.compare(self.evidence(), self.evidence(), identity=False))

    def test_compare_results_is_exact_and_ignores_declared_paths(self):
        before = self.evidence()
        self.assertEqual([], evidence.compare(before, self.evidence(), identity=False))
        (self.root / "out" / "counts.tsv").write_text("gene\tcount\nA\t4\n")
        problems = evidence.compare(before, self.evidence(), identity=False)
        self.assertEqual(1, len(problems))
        self.assertIn("counts: text sha256", problems[0])
        moved = [("counts", "elsewhere/counts.tsv")]
        (self.root / "elsewhere").mkdir()
        (self.root / "elsewhere" / "counts.tsv").write_text("gene\tcount\nA\t3\n")
        self.assertEqual([], evidence.compare(before, self.evidence(roles=moved), identity=False))

    def test_compare_rejects_missing_on_either_side(self):
        missing = [("counts", "out/absent.tsv")]
        problems = evidence.compare(self.evidence(roles=missing), self.evidence(roles=missing), identity=False)
        self.assertIn("counts: no digest in baseline and current", problems)
        problems = evidence.compare(self.evidence(roles=missing), self.evidence(), identity=False)
        self.assertIn("counts: no digest in baseline", problems)

    def test_compare_flags_changed_root_entries(self):
        before = evidence.build_evidence([("root", "out")], repo_root=self.root, stage="post")
        (self.root / "out" / "extra.tsv").write_text("y\n")
        after = evidence.build_evidence([("root", "out")], repo_root=self.root, stage="post")
        problems = evidence.compare(before, after, identity=False)
        self.assertEqual(1, len(problems))
        self.assertIn("root/extra.tsv: only in current", problems[0])


class ConfigAndCliTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        (self.root / "counts.tsv").write_text("gene\tcount\nA\t3\n")

    def tearDown(self):
        self.tmp.cleanup()

    def config(self, body):
        path = self.root / "fixture.json"  # JSON is valid YAML; avoids a hard PyYAML need
        path.write_text(json.dumps(body))
        return path

    def test_roles_from_list_and_mapping_keep_explicit_order(self):
        listed = self.config({"replay_evidence": {"roles": [
            {"role": "second", "path": "b.tsv"}, {"role": "first", "path": "a.tsv"}]}})
        self.assertEqual([("second", "b.tsv"), ("first", "a.tsv")], evidence.load_roles(listed))
        mapped = self.config({"replay_evidence": {"roles": {"z": "z.tsv", "a": "a.tsv"}}})
        self.assertEqual([("z", "z.tsv"), ("a", "a.tsv")], evidence.load_roles(mapped))

    def test_roles_reject_implicit_or_duplicate_specs(self):
        cases = [{"pipeline": "rnaseq"}, {"replay_evidence": {"roles": []}},
                 {"replay_evidence": {"roles": [{"role": "x"}]}},
                 {"replay_evidence": {"roles": [{"role": "x", "path": "a"}, {"role": "x", "path": "b"}]}}]
        for body in cases:
            with self.subTest(body=body), self.assertRaises(ValueError):
                evidence.load_roles(self.config(body))

    def cli(self, *args):
        return subprocess.run([sys.executable, str(ROOT / "scripts" / "replay_evidence.py"), *map(str, args)],
                              capture_output=True, text=True)

    def test_cli_record_and_compare(self):
        config = self.config({"replay_evidence": {"roles": [{"role": "counts", "path": "counts.tsv"}]}})
        pre, post = self.root / "pre.json", self.root / "post.json"
        for stage, output in (("pre", pre), ("post", post)):
            done = self.cli("record", "--config", config, "--output", output, "--stage", stage,
                            "--invocation", "run --cores 2", "--repo-root", self.root)
            self.assertEqual(0, done.returncode, done.stderr)
        self.assertEqual("bytes", json.loads(pre.read_text())["roles"][0]["input"]["kind"])
        self.assertEqual("text", json.loads(post.read_text())["roles"][0]["result"]["kind"])
        done = self.cli("compare", "--baseline", post, "--current", post)
        self.assertEqual(0, done.returncode, done.stdout)
        self.assertIn("equal", done.stdout)
        (self.root / "counts.tsv").write_text("gene\tcount\nA\t9\n")
        changed = self.root / "changed.json"
        self.cli("record", "--config", config, "--output", changed, "--repo-root", self.root)
        done = self.cli("compare", "--baseline", post, "--current", changed, "--skip-identity")
        self.assertEqual(1, done.returncode)
        self.assertIn("counts: text sha256", done.stdout)


if __name__ == "__main__":
    unittest.main()
