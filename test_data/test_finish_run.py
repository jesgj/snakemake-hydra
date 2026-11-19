#!/usr/bin/env python3
"""Dependency-light negative tests for the real-data fixture verification contracts.

Run from the repository root:

    python3 -B test_data/test_finish_run.py

These tests exercise the pure inventory and table checks in finish_run.py. They
intentionally do not need Pixi, samtools, or any generated pipeline outputs.
"""
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import finish_run


def write(path, content):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content)
    return path


def touch_nonempty(path, content="fixture"):
    return write(path, content)


def make_section(root, pipeline):
    return {
        "pipeline": pipeline,
        "samples_info": {"sample_a": {"type": "PE"}},
        "alignment_dir": str(root / "alignment_dir"),
        "marked_bam_dir": str(root / "marked_bam_dir"),
        "dedup_dir": str(root / "dedup_dir"),
        "sorted_filtered_bam_dir": str(root / "sorted_filtered_bam_dir"),
        "filtered_bam_dir": str(root / "filtered_bam_dir"),
        "kallisto_output_dir": str(root / "kallisto_output_dir"),
        "methyldackel_dir": str(root / "methyldackel_dir"),
        "mbias_dir": str(root / "mbias_dir"),
        "consensus_peaks_dir": str(root / "consensus_peaks_dir"),
        "featurecounts_dir": str(root / "featurecounts_dir"),
        "qc_summary_dir": str(root / "qc_summary_dir"),
    }


class BamInventoryTests(unittest.TestCase):
    def test_missing_expected_bam_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            section = make_section(root, "rnaseq")
            root.mkdir(exist_ok=True)
            with self.assertRaisesRegex(ValueError, "Missing or empty expected BAM output"):
                finish_run.verify_bam_inventory(section, root)

    def test_missing_expected_bai_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            section = make_section(root, "chip_cr")
            section["samples_info"] = {"treat_rep1": {"type": "SE"}}
            for bam in finish_run.expected_bam_outputs(section):
                # Keep only the BAM so the missing-index branch is exercised.
                touch_nonempty(bam)
            with self.assertRaisesRegex(ValueError, "Missing or empty expected BAM index"):
                finish_run.verify_bam_inventory(section, root)

    def test_stale_sample_bam_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            section = make_section(root, "chip_cr")
            section["samples_info"] = {"treat_rep1": {"type": "SE"}}
            for bam, index in finish_run.expected_bam_outputs(section).items():
                touch_nonempty(bam)
                if index is not None:
                    touch_nonempty(index)
            stale = Path(section["alignment_dir"]) / "removed_sample_se.sorted.bam"
            touch_nonempty(stale)
            with self.assertRaisesRegex(ValueError, "Unexpected BAM output.*stale sample artifact"):
                finish_run.verify_bam_inventory(section, root)

    def test_expected_bams_and_indexes_are_accepted(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            section = make_section(root, "wgbs")
            for bam, index in finish_run.expected_bam_outputs(section).items():
                touch_nonempty(bam)
                if index is not None:
                    touch_nonempty(index)
            expected = finish_run.verify_bam_inventory(section, root)
            self.assertEqual(len(expected), 3)


class ScientificTableTests(unittest.TestCase):
    def test_missing_abundance_table_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "sample_a" / "abundance.tsv"
            with self.assertRaisesRegex(ValueError, "Missing or empty RNA-seq abundance table"):
                finish_run.check_abundance_table(path, "sample_a")

    def test_malformed_rnaseq_abundance_table_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "abundance.tsv"
            write(path, "target_id\tlength\teff_length\test_counts\ttpm\n"
                        "tx1\t100\t90\t12\tnot_a_number\n")
            with self.assertRaisesRegex(ValueError, "Malformed RNA-seq abundance"):
                finish_run.check_abundance_table(path, "sample_a")

    def test_malformed_wgbs_methylkit_table_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "sample_a_CpG.methylKit"
            write(path, "\t".join(finish_run.METHYLKIT_HEADER) + "\n"
                        "I.10\tI\t10\tF\t3\t50.0\n")
            with self.assertRaisesRegex(ValueError, "Malformed WGBS methylKit record"):
                finish_run.check_methylkit_table(path, "sample_a")

    def test_wgbs_mbias_table_rejects_stale_sample_identity(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "all_samples_mbias_options.tsv"
            write(path, "sample\toptions\nsample_a\t--OT 0\nremoved_sample\t--OT 1\n")
            with self.assertRaisesRegex(ValueError, "WGBS mbias option sample mismatch"):
                finish_run.check_mbias_options_table(path, ["sample_a"])

    def test_atac_featurecounts_rejects_stale_sample_column(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "consensus_peak_counts.txt"
            write(path, "Geneid\tChr\tStart\tEnd\tStrand\tLength\t"
                        "/scratch/sample_a_pe.filtered.sorted.bam\t"
                        "/scratch/removed_sample_pe.filtered.sorted.bam\n"
                        "peak_1\tI\t1\t10\t.\t9\t4\t2\n")
            with self.assertRaisesRegex(ValueError, "Unexpected ATAC featureCounts matrix sample columns"):
                finish_run.check_featurecounts_matrix(path, ["sample_a", "sample_b"], expected_rows=1)

    def test_malformed_atac_featurecounts_count_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "consensus_peak_counts.txt"
            write(path, "Geneid\tChr\tStart\tEnd\tStrand\tLength\t"
                        "/scratch/sample_a_pe.filtered.sorted.bam\n"
                        "peak_1\tI\t1\t10\t.\t9\tnot_a_count\n")
            with self.assertRaisesRegex(ValueError, "Malformed ATAC featureCounts count"):
                finish_run.check_featurecounts_matrix(path, ["sample_a"], expected_rows=1)

    def test_atac_qc_summary_rejects_missing_or_stale_samples(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "sample_qc_summary.tsv"
            rows = ["\t".join(finish_run.ATAC_QC_COLUMNS)]
            rows.append("sample_a\t10\t9\t1\t10.0\t0\t0.0\t2\t9\t22.22")
            rows.append("removed_sample\t8\t7\t1\t12.5\t0\t0.0\t1\t7\t14.28")
            write(path, "\n".join(rows) + "\n")
            with self.assertRaisesRegex(ValueError, "ATAC QC summary sample mismatch"):
                finish_run.check_atac_qc_summary(path, ["sample_a"])

    def test_atac_cross_file_row_count_is_enforced(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            section = make_section(root, "atacseq")
            write(Path(section["consensus_peaks_dir"]) / "consensus_peaks.bed",
                  "I\t10\t20\n")
            write(Path(section["consensus_peaks_dir"]) / "consensus_peaks.saf",
                  "GeneID\tChr\tStart\tEnd\tStrand\npeak_1\tI\t11\t20\t.\n")
            write(Path(section["featurecounts_dir"]) / "consensus_peak_counts.txt",
                  "# fixture\nGeneid\tChr\tStart\tEnd\tStrand\tLength\t"
                  "/scratch/sample_a_pe.filtered.sorted.bam\n")
            write(Path(section["featurecounts_dir"]) / "consensus_peak_counts.txt.summary",
                  "Status\t/scratch/sample_a_pe.filtered.sorted.bam\n"
                  "Assigned\t0\n")
            with self.assertRaisesRegex(ValueError, "row count differs from consensus peaks"):
                finish_run.check_consensus_cross_files(section, ["sample_a"])

    def test_empty_atac_consensus_files_are_valid(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            section = make_section(root, "atacseq")
            write(Path(section["consensus_peaks_dir"]) / "consensus_peaks.bed", "")
            write(Path(section["consensus_peaks_dir"]) / "consensus_peaks.saf",
                  "GeneID\tChr\tStart\tEnd\tStrand\n")
            write(Path(section["featurecounts_dir"]) / "consensus_peak_counts.txt",
                  "# No consensus peaks; featureCounts skipped.\n"
                  "Geneid\tChr\tStart\tEnd\tStrand\tLength\t"
                  "/scratch/sample_a_pe.filtered.sorted.bam\n")
            write(Path(section["featurecounts_dir"]) / "consensus_peak_counts.txt.summary",
                  "Status\t/scratch/sample_a_pe.filtered.sorted.bam\n"
                  "Assigned\t0\n")
            observed = finish_run.check_consensus_cross_files(section, ["sample_a"])
            self.assertEqual(observed, {"consensus_intervals": 0, "count_matrix_rows": 0})


if __name__ == "__main__":
    unittest.main()
