#!/usr/bin/env python3
"""Exercise prepared fixture manifests without downloading or running a workflow."""
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import prepare

CATALOG = json.loads((prepare.ROOT / "datasets.json").read_text())


class PreparedInputsTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        root_patch = patch.object(prepare, "ROOT", self.root)
        root_patch.start()
        self.addCleanup(root_patch.stop)
        self.provenance = self.root / "provenance.json"
        self.variants = ["rnaseq", "atacseq_macs3", "wgbs_full"]
        self.paths = prepare.prepared_paths(CATALOG, self.variants)
        for path in self.paths:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(path.name)
        self.report = {"seed": CATALOG["seed"], "prepared_sha256": prepare.prepared_checksums(CATALOG, self.variants)}
        self.write_report()

    def write_report(self):
        self.provenance.write_text(json.dumps(self.report))

    def verify(self, variants=None):
        prepare.verify_prepared(CATALOG, self.variants if variants is None else variants, self.provenance)

    def test_round_trip_covers_derived_inputs(self):
        paths = {path.relative_to(self.root).as_posix() for path in self.paths}
        self.assertIn("references/rnaseq/genes.bed12", paths)
        self.assertIn("references/tair10/genome.fa", paths)
        self.assertNotIn("references/wgbs/genome.fa", paths)
        splits = [path for path in paths if path.startswith("fastqs/atacseq/split_runs/")]
        self.assertEqual(len(splits), 8)
        self.verify()

    def test_corruption_is_rejected(self):
        self.paths[0].write_text("unexpected content")
        with self.assertRaisesRegex(ValueError, "SHA-256 mismatch"):
            self.verify()

    def test_missing_input_is_rejected(self):
        self.paths[-1].unlink()
        with self.assertRaisesRegex(ValueError, "Missing or invalid prepared fixture"):
            self.verify()

    def test_wrong_fixture_or_incomplete_manifest_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "does not match selected inputs"):
            self.verify(["wgbs"])
        self.report["prepared_sha256"].pop(next(iter(self.report["prepared_sha256"])))
        self.write_report()
        with self.assertRaisesRegex(ValueError, "does not match selected inputs"):
            self.verify()


if __name__ == "__main__":
    unittest.main()
