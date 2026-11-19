#!/usr/bin/env python3
"""Check successful smoke-run artifacts, save evidence, and remove their outputs."""
import argparse
import hashlib
import json
import math
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent

ABUNDANCE_HEADER = ["target_id", "length", "eff_length", "est_counts", "tpm"]
METHYLKIT_HEADER = ["chrBase", "chr", "base", "strand", "coverage", "freqC", "freqT"]
MBIAS_OPTIONS_HEADER = ["sample", "options"]
FEATURECOUNTS_COLUMNS = ["Geneid", "Chr", "Start", "End", "Strand", "Length"]
ATAC_QC_COLUMNS = [
    "sample",
    "aligned_reads",
    "final_filtered_reads",
    "duplicate_reads",
    "duplication_rate_pct",
    "mitochondrial_reads",
    "mitochondrial_fraction_pct",
    "fragments_in_peaks",
    "featurecounts_total_fragments",
    "frip_pct",
]
ATAC_QC_SAMPLE_COLUMN = 0


def _split_line(line):
    return line.rstrip("\r\n").split("\t")


def _data_lines(path, label, skip_comments=False):
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"Missing or empty {label}: {path}")
    lines = path.read_text().splitlines()
    if skip_comments:
        lines = [line for line in lines if line and not line.startswith("#")]
    else:
        lines = [line for line in lines if line]
    if not lines:
        raise ValueError(f"Missing table header in {label}: {path}")
    return lines


def _number(value, label, path, line_number, minimum=None, integral=False):
    try:
        observed = int(value) if integral else float(value)
    except (TypeError, ValueError):
        raise ValueError(f"Malformed {label} value {value!r} at {path}:{line_number}") from None
    if not integral and not math.isfinite(observed):
        raise ValueError(f"Malformed {label} value {value!r} at {path}:{line_number}")
    if minimum is not None and observed < minimum:
        raise ValueError(f"Malformed {label} value {value!r} at {path}:{line_number}")
    return observed


def _validate_header(path, label, expected, observed):
    if observed != expected:
        raise ValueError(
            f"Unexpected {label} header in {path}: expected {expected}, found {observed}"
        )


def expected_bam_outputs(section, mode=None):
    """Return the explicit BAM -> required BAI contract for the selected assay."""
    mode = mode or section.get("pipeline")
    if not mode:
        raise ValueError("A pipeline mode is required for BAM contracts")
    samples = list(section["samples_info"])
    if not samples:
        raise ValueError(f"No configured {mode} samples to verify")
    expected = {}

    def add(bam, index=None):
        bam_path = Path(bam)
        index_path = None if index is None else Path(index)
        if bam_path in expected:
            raise ValueError(f"Duplicate expected BAM contract: {bam_path}")
        expected[bam_path] = index_path

    if mode == "rnaseq":
        for sample in samples:
            aligned = Path(section["alignment_dir"]) / f"{sample}_pe.sorted.bam"
            marked = Path(section["marked_bam_dir"]) / f"{sample}_pe.markdup.bam"
            add(aligned, f"{aligned}.bai")
            # Picard is instructed to create an index with its own suffix.
            add(marked, Path(section["marked_bam_dir"]) / f"{sample}_pe.markdup.bai")
    elif mode == "wgbs":
        for sample in samples:
            add(Path(section["alignment_dir"]) / f"{sample}_pe.bam")
            add(Path(section["dedup_dir"]) / f"{sample}_pe.deduplicated.bam")
            filtered = Path(section["sorted_filtered_bam_dir"]) / f"{sample}_pe.filtered.sorted.bam"
            add(filtered, f"{filtered}.bai")
    elif mode == "chip_cr":
        for sample, info in section["samples_info"].items():
            read_type = str(info.get("type", "")).lower()
            if read_type not in {"pe", "se"}:
                raise ValueError(f"ChIP/CUT&RUN sample {sample!r} has invalid type {info.get('type')!r}")
            add(Path(section["alignment_dir"]) / f"{sample}_{read_type}.sorted.bam")
            filtered = Path(section["filtered_bam_dir"]) / f"{sample}_{read_type}.filtered.sorted.bam"
            add(filtered, f"{filtered}.bai")
    elif mode == "atacseq":
        for sample, info in section["samples_info"].items():
            runs = info.get("runs") or {"run1": info}
            for run in runs:
                add(Path(section["alignment_dir"]) / "runs" / f"{sample}__{run}.sorted.bam")
            merged = Path(section["alignment_dir"]) / f"{sample}_pe.sorted.bam"
            add(merged, f"{merged}.bai")
            filtered = Path(section["filtered_bam_dir"]) / f"{sample}_pe.filtered.sorted.bam"
            add(filtered, f"{filtered}.bai")
    else:
        raise ValueError(f"Unsupported pipeline for BAM contracts: {mode}")
    return expected


def verify_bam_inventory(section, work, mode=None):
    """Check explicit BAM/index presence and reject stale or unexpected BAM artifacts."""
    expected = expected_bam_outputs(section, mode)
    missing = [path for path in expected if not path.is_file() or path.stat().st_size == 0]
    if missing:
        raise ValueError("Missing or empty expected BAM output: " + ", ".join(map(str, missing)))
    missing_indexes = [
        index for index in expected.values()
        if index is not None and (not index.is_file() or index.stat().st_size == 0)
    ]
    if missing_indexes:
        raise ValueError("Missing or empty expected BAM index: " + ", ".join(map(str, missing_indexes)))
    discovered = {path.resolve() for path in Path(work).rglob("*.bam") if path.is_file()}
    declared = {path.resolve() for path in expected}
    unexpected = sorted(discovered - declared)
    if unexpected:
        raise ValueError(
            "Unexpected BAM output (possible stale sample artifact): "
            + ", ".join(map(str, unexpected))
        )
    return expected


def check_abundance_table(path, sample):
    """Validate kallisto abundance schema and numeric columns for one sample."""
    lines = _data_lines(path, f"RNA-seq abundance table for {sample}")
    _validate_header(path, "RNA-seq abundance", ABUNDANCE_HEADER, _split_line(lines[0]))
    for line_number, line in enumerate(lines[1:], 2):
        fields = _split_line(line)
        if len(fields) != len(ABUNDANCE_HEADER):
            raise ValueError(f"Malformed RNA-seq abundance record in {path}:{line_number}")
        if not fields[0]:
            raise ValueError(f"Missing transcript ID in {path}:{line_number}")
        for column, value in zip(ABUNDANCE_HEADER[1:], fields[1:]):
            _number(value, f"RNA-seq abundance {column}", path, line_number, minimum=0)
    return True


def check_methylkit_table(path, sample):
    """Validate MethylDackel methylKit schema and value ranges for one sample."""
    lines = _data_lines(path, f"WGBS methylKit table for {sample}")
    _validate_header(path, "WGBS methylKit", METHYLKIT_HEADER, _split_line(lines[0]))
    for line_number, line in enumerate(lines[1:], 2):
        fields = _split_line(line)
        if len(fields) != len(METHYLKIT_HEADER):
            raise ValueError(f"Malformed WGBS methylKit record in {path}:{line_number}")
        chr_base, chrom, base, strand, coverage, freq_c, freq_t = fields
        if not chr_base or not chrom:
            raise ValueError(f"Missing WGBS methylKit locus in {path}:{line_number}")
        if strand not in {"F", "R"}:
            raise ValueError(f"Malformed WGBS methylKit strand in {path}:{line_number}: {strand!r}")
        _number(base, "WGBS methylKit base", path, line_number, minimum=1, integral=True)
        _number(coverage, "WGBS methylKit coverage", path, line_number, minimum=0, integral=True)
        # MethylDackel emits fixed-point percentages, sometimes padded with spaces.
        for ratio in (freq_c.strip(), freq_t.strip()):
            _number(ratio, "WGBS methylKit frequency", path, line_number, minimum=0)
    return True


def check_mbias_options_table(path, samples):
    """Validate the aggregate MethylDackel options table and its sample identity."""
    lines = _data_lines(path, "WGBS mbias options table")
    _validate_header(path, "WGBS mbias options", MBIAS_OPTIONS_HEADER, _split_line(lines[0]))
    observed = []
    for line_number, line in enumerate(lines[1:], 2):
        fields = _split_line(line)
        if len(fields) != 2 or not fields[0]:
            raise ValueError(f"Malformed WGBS mbias options record in {path}:{line_number}")
        observed.append(fields[0])
    if len(observed) != len(set(observed)):
        raise ValueError(f"Duplicate WGBS mbias option samples in {path}")
    missing = sorted(set(samples) - set(observed))
    unexpected = sorted(set(observed) - set(samples))
    if missing or unexpected:
        raise ValueError(
            f"WGBS mbias option sample mismatch in {path}: missing {missing}, unexpected {unexpected}"
        )
    return True


def _featurecounts_bam_names(samples):
    return [f"{sample}_pe.filtered.sorted.bam" for sample in samples]


def _check_bam_columns(path, label, header, samples, prefix_columns):
    expected_names = _featurecounts_bam_names(samples)
    if len(header) != prefix_columns + len(samples):
        raise ValueError(
            f"Unexpected {label} column count in {path}: "
            f"expected {prefix_columns + len(samples)}, found {len(header)}"
        )
    observed_names = [Path(column).name for column in header[prefix_columns:]]
    if observed_names != expected_names:
        raise ValueError(
            f"Unexpected {label} sample columns in {path}: "
            f"expected {expected_names}, found {observed_names}"
        )


def check_featurecounts_matrix(path, samples, expected_rows=None):
    """Validate the ATAC consensus count matrix schema, sample identity, and row count."""
    lines = _data_lines(path, "ATAC featureCounts matrix", skip_comments=True)
    header = _split_line(lines[0])
    if header[:len(FEATURECOUNTS_COLUMNS)] != FEATURECOUNTS_COLUMNS:
        raise ValueError(f"Unexpected ATAC featureCounts header in {path}: {header[:6]}")
    _check_bam_columns(path, "ATAC featureCounts matrix", header, samples, len(FEATURECOUNTS_COLUMNS))
    rows = 0
    for line_number, line in enumerate(lines[1:], 2):
        fields = _split_line(line)
        if len(fields) != len(header):
            raise ValueError(f"Malformed ATAC featureCounts record in {path}:{line_number}")
        if not fields[0]:
            raise ValueError(f"Missing ATAC consensus peak ID in {path}:{line_number}")
        _number(fields[2], "ATAC featureCounts start", path, line_number, minimum=1, integral=True)
        _number(fields[3], "ATAC featureCounts end", path, line_number, minimum=1, integral=True)
        _number(fields[5], "ATAC featureCounts length", path, line_number, minimum=1, integral=True)
        for value in fields[6:]:
            _number(value, "ATAC featureCounts count", path, line_number, minimum=0, integral=True)
        rows += 1
    if expected_rows is not None and rows != expected_rows:
        raise ValueError(
            f"ATAC featureCounts row count differs from consensus peaks in {path}: "
            f"{rows} != {expected_rows}"
        )
    return rows


def check_featurecounts_summary(path, samples):
    """Validate the featureCounts summary columns and nonnegative counters."""
    lines = _data_lines(path, "ATAC featureCounts summary")
    header = _split_line(lines[0])
    if not header or header[0] != "Status":
        raise ValueError(f"Unexpected ATAC featureCounts summary header in {path}: {header}")
    _check_bam_columns(path, "ATAC featureCounts summary", header, samples, 1)
    statuses = []
    for line_number, line in enumerate(lines[1:], 2):
        fields = _split_line(line)
        if len(fields) != len(header) or not fields[0]:
            raise ValueError(f"Malformed ATAC featureCounts summary record in {path}:{line_number}")
        for value in fields[1:]:
            _number(value, "ATAC featureCounts summary count", path, line_number, minimum=0, integral=True)
        statuses.append(fields[0])
    if "Assigned" not in statuses:
        raise ValueError(f"ATAC featureCounts summary in {path} has no Assigned row")
    return True


def check_atac_qc_summary(path, samples):
    """Validate the ATAC QC table schema and one row per configured sample."""
    lines = _data_lines(path, "ATAC QC summary")
    _validate_header(path, "ATAC QC summary", ATAC_QC_COLUMNS, _split_line(lines[0]))
    observed = []
    for line_number, line in enumerate(lines[1:], 2):
        fields = _split_line(line)
        if len(fields) != len(ATAC_QC_COLUMNS) or not fields[ATAC_QC_SAMPLE_COLUMN]:
            raise ValueError(f"Malformed ATAC QC summary record in {path}:{line_number}")
        observed.append(fields[ATAC_QC_SAMPLE_COLUMN])
        for column, value in zip(ATAC_QC_COLUMNS[1:], fields[1:]):
            if value != "NA":
                _number(value, f"ATAC QC {column}", path, line_number, minimum=0)
    if len(observed) != len(set(observed)):
        raise ValueError(f"Duplicate samples in ATAC QC summary {path}")
    missing = sorted(set(samples) - set(observed))
    unexpected = sorted(set(observed) - set(samples))
    if missing or unexpected:
        raise ValueError(
            f"ATAC QC summary sample mismatch in {path}: missing {missing}, unexpected {unexpected}"
        )
    return True


def check_consensus_cross_files(section, samples):
    """Validate BED/SAF/count-table correspondence, including empty peak sets."""
    bed_path = Path(section["consensus_peaks_dir"]) / "consensus_peaks.bed"
    saf_path = Path(section["consensus_peaks_dir"]) / "consensus_peaks.saf"
    matrix_path = Path(section["featurecounts_dir"]) / "consensus_peak_counts.txt"
    summary_path = Path(section["featurecounts_dir"]) / "consensus_peak_counts.txt.summary"
    if not bed_path.is_file():
        raise ValueError(f"Missing ATAC consensus BED: {bed_path}")
    bed_intervals = []
    # An empty consensus peak set is a valid terminal result.
    for line_number, line in enumerate(bed_path.read_text().splitlines(), 1):
        if not line or line.startswith("#"):
            continue
        fields = _split_line(line)
        if len(fields) < 3:
            raise ValueError(f"Malformed ATAC consensus BED record in {bed_path}:{line_number}")
        start = _number(fields[1], "ATAC consensus start", bed_path, line_number, minimum=0, integral=True)
        end = _number(fields[2], "ATAC consensus end", bed_path, line_number, minimum=0, integral=True)
        if start >= end:
            raise ValueError(f"Invalid ATAC consensus interval in {bed_path}:{line_number}")
        bed_intervals.append(tuple(fields[:3]))
    saf_lines = _data_lines(saf_path, "ATAC consensus SAF")
    _validate_header(
        saf_path, "ATAC consensus SAF", ["GeneID", "Chr", "Start", "End", "Strand"], _split_line(saf_lines[0])
    )
    for line_number, line in enumerate(saf_lines[1:], 2):
        fields = _split_line(line)
        if len(fields) != 5 or fields[4] != ".":
            raise ValueError(f"Malformed ATAC consensus SAF record in {saf_path}:{line_number}")
        _number(fields[2], "ATAC consensus SAF start", saf_path, line_number, minimum=1, integral=True)
        _number(fields[3], "ATAC consensus SAF end", saf_path, line_number, minimum=1, integral=True)
    if len(saf_lines) - 1 != len(bed_intervals):
        raise ValueError(
            f"ATAC consensus SAF/BED row mismatch: {saf_path} has {len(saf_lines) - 1}, "
            f"{bed_path} has {len(bed_intervals)}"
        )
    matrix_rows = check_featurecounts_matrix(matrix_path, samples, expected_rows=len(bed_intervals))
    check_featurecounts_summary(summary_path, samples)
    return {"consensus_intervals": len(bed_intervals), "count_matrix_rows": matrix_rows}


def verify_scientific_tables(config):
    """Validate selected-assay tables and report the checked paths."""
    section = config[config["pipeline"]]
    mode = config["pipeline"]
    samples = list(section["samples_info"])
    checked = {}
    if mode == "rnaseq":
        for sample in samples:
            path = Path(section["kallisto_output_dir"]) / sample / "abundance.tsv"
            check_abundance_table(path, sample)
            checked[str(path.relative_to(ROOT))] = "kallisto abundance"
            # The committed fixture deliberately exercises seeded bootstrapping.
            if "--bootstrap-samples=2" in section.get("kallisto", {}).get("extra_args", "").split():
                info_path = path.with_name("run_info.json")
                if json.loads(info_path.read_text()).get("n_bootstraps") != 2:
                    raise ValueError(f"Expected two kallisto bootstraps in {info_path}")
                hdf5 = path.with_name("abundance.h5")
                with hdf5.open("rb") as handle:
                    if handle.read(8) != b"\x89HDF\r\n\x1a\n":
                        raise ValueError(f"Expected kallisto HDF5 output at {hdf5}")
                checked[str(info_path.relative_to(ROOT))] = "two seeded kallisto bootstraps"
    elif mode == "wgbs":
        for sample in samples:
            path = Path(section["methyldackel_dir"]) / f"{sample}_CpG.methylKit"
            check_methylkit_table(path, sample)
            checked[str(path.relative_to(ROOT))] = "methylKit CpG"
        options_path = Path(section["mbias_dir"]) / "all_samples_mbias_options.tsv"
        check_mbias_options_table(options_path, samples)
        checked[str(options_path.relative_to(ROOT))] = "mbias options"
    elif mode == "atacseq":
        path = Path(section["qc_summary_dir"]) / "sample_qc_summary.tsv"
        check_atac_qc_summary(path, samples)
        checked[str(path.relative_to(ROOT))] = "ATAC QC summary"
        cross_files = check_consensus_cross_files(section, samples)
        checked[str((Path(section["featurecounts_dir"]) / "consensus_peak_counts.txt").relative_to(ROOT))] = (
            f"featureCounts matrix ({cross_files['count_matrix_rows']} rows)"
        )
    return checked


def verify_bam_files(expected):
    """Check BAM integrity, mapped records, and that every required index is usable."""
    bams = list(expected)
    subprocess.run(["samtools", "quickcheck", "-v", *map(str, bams)], check=True)
    mapped = {}
    for bam in bams:
        count = int(subprocess.check_output(["samtools", "view", "-c", "-F", "4", str(bam)], text=True))
        if not count:
            raise ValueError(f"No mapped reads in {bam}; retain outputs for diagnosis")
        mapped[str(bam.relative_to(ROOT))] = count
    indexed_bams = [(bam, index) for bam, index in expected.items() if index is not None]
    for bam, index in indexed_bams:
        # idxstats may fall back to a full scan when the index is corrupt.
        # An indexed region query requires the exact index declared above.
        header = subprocess.check_output(["samtools", "view", "-H", str(bam)], text=True)
        contig = next((field[3:] for line in header.splitlines() if line.startswith("@SQ\t")
                       for field in line.split("\t")[1:] if field.startswith("SN:")), None)
        if not contig:
            raise ValueError(f"No reference contig in indexed BAM: {bam}")
        subprocess.run(
            ["samtools", "view", "-c", "-X", str(bam), str(index), contig],
            check=True,
            stdout=subprocess.DEVNULL,
        )
    return mapped, [str(index.relative_to(ROOT)) for _, index in indexed_bams]


def verify_atac_run_merges(section):
    """Check that every run contributes all records with its original RG metadata."""
    evidence = {}
    for sample, info in section["samples_info"].items():
        if "runs" not in info:
            continue
        merged = Path(section["alignment_dir"]) / f"{sample}_pe.sorted.bam"
        header = subprocess.check_output(["samtools", "view", "-H", str(merged)], text=True)
        groups = [dict(field.split(":", 1) for field in line.split("\t")[1:])
                  for line in header.splitlines() if line.startswith("@RG\t")]
        expected = [run["read_group"] for run in info["runs"].values()]
        if len(groups) != len(expected) or {g["ID"] for g in groups} != {g["id"] for g in expected}:
            raise ValueError(f"Merged read-group IDs differ for {sample}")
        counts = {}
        for run, run_info in info["runs"].items():
            rg = run_info["read_group"]
            observed = next(g for g in groups if g["ID"] == rg["id"])
            for field, key in (("SM", "sample"), ("LB", "library"), ("PL", "platform"), ("PU", "platform_unit")):
                if observed.get(field) != rg[key]:
                    raise ValueError(f"Merged {sample}/{run} changed RG {field}")
            bam = Path(section["alignment_dir"]) / "runs" / f"{sample}__{run}.sorted.bam"
            original = int(subprocess.check_output(["samtools", "view", "-c", str(bam)], text=True))
            retained = int(subprocess.check_output(["samtools", "view", "-c", "-r", rg["id"], str(merged)], text=True))
            if original == 0 or retained != original:
                raise ValueError(f"Merged {sample}/{run} lost or duplicated records: {original} -> {retained}")
            counts[rg["id"]] = retained
        total = int(subprocess.check_output(["samtools", "view", "-c", str(merged)], text=True))
        if total != sum(counts.values()):
            raise ValueError(f"Merged {sample} has missing or unexpected RG tags")
        evidence[sample] = {"records_by_read_group": counts, "merged_records": total}
    return evidence


def verify_final_outputs(config):
    """Check selected assay outputs in addition to the shared BAM/report checks."""
    mode = config["pipeline"]
    section = config[mode]
    required = []
    for sample, info in section["samples_info"].items():
        if mode == "rnaseq":
            required.append(Path(section["kallisto_output_dir"]) / sample / "abundance.tsv")
        elif mode == "wgbs":
            required.append(Path(section["methyldackel_dir"]) / f"{sample}_CpG.methylKit")
        else:
            read_type = info["type"].lower()
            required.append(Path(section["bigwig_dir"]) / f"{sample}_{read_type}.bw")
        if mode == "atacseq":
            # An empty peak set is legitimate, but the declared file must exist.
            peak = Path(section["peaks_dir"]) / f"{sample}_peaks.narrowPeak"
            if not peak.is_file():
                raise ValueError(f"Missing peak output: {peak}")
    if mode == "atacseq":
        required.extend([
            Path(section["consensus_peaks_dir"]) / "consensus_peaks.bed",
            Path(section["consensus_peaks_dir"]) / "consensus_peaks.saf",
            Path(section["featurecounts_dir"]) / "consensus_peak_counts.txt",
            Path(section["featurecounts_dir"]) / "consensus_peak_counts.txt.summary",
            Path(section["qc_summary_dir"]) / "sample_qc_summary.tsv",
        ])
    if mode == "rnaseq" and section["gene_body_coverage"]["enabled"]:
        prefix = Path(section["gene_body_coverage_dir"]) / "all_samples"
        fmt = section["gene_body_coverage"].get("format", "pdf")
        required.extend([Path(f"{prefix}.geneBodyCoverage.txt"), Path(f"{prefix}.geneBodyCoverage.curves.{fmt}")])
        if len(section["samples_info"]) >= 3:
            required.append(Path(f"{prefix}.geneBodyCoverage.heatMap.{fmt}"))
    for path in required:
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"Missing or empty assay output: {path}")
    return [str(path.relative_to(ROOT)) for path in required]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("variant", choices=[p.stem for p in (ROOT / "configs").glob("*.yaml")])
    parser.add_argument("--cores", type=int, default=2)
    parser.add_argument("--keep-outputs", action="store_true", help="Retain verified results for CI artifacts")
    parser.add_argument("--evidence", type=Path,
                        help="Destination for execution checks and artifact hashes")
    args = parser.parse_args()
    import yaml  # Keep the validation helpers importable without third-party modules.

    config = yaml.safe_load((ROOT / "configs" / (args.variant + ".yaml")).read_text())
    section = config[config["pipeline"]]
    work = ROOT / "work" / args.variant
    report = Path(section["multiqc_results_dir"]) / "multiqc_report.html"
    sources = report.parent / "multiqc_report_data" / "multiqc_sources.txt"
    for path in (report, sources):
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"Missing or empty MultiQC artifact: {path}")
    assay_outputs = verify_final_outputs(config)
    table_contracts = verify_scientific_tables(config)
    expected_bams = verify_bam_inventory(section, work, config["pipeline"])
    run_merges = verify_atac_run_merges(section) if config["pipeline"] == "atacseq" else {}
    files = sorted(p for p in work.rglob("*") if p.is_file())
    mapped, bam_indexes = verify_bam_files(expected_bams)
    artifacts = []
    peak_counts = {}
    for path in files:
        if path.suffix == ".narrowPeak":
            count = 0
            for line in path.read_text().splitlines():
                if not line or line.startswith(("#", "track")):
                    continue
                fields = line.split("\t")
                if len(fields) != 10 or not 0 <= int(fields[1]) < int(fields[2]):
                    raise ValueError(f"Invalid narrowPeak record in {path}")
                count += 1
            peak_counts[str(path.relative_to(ROOT))] = count
        h = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                h.update(chunk)
        artifacts.append({"path": str(path.relative_to(ROOT)), "bytes": path.stat().st_size, "sha256": h.hexdigest()})
    entry = {"variant": args.variant, "pipeline": config["pipeline"], "checked_at_utc": datetime.now(timezone.utc).isoformat(),
             "scheduler_cores": args.cores, "scheduler_jobs": 1,
             "multiqc_report_verified": True, "bam_quickcheck": True, "mapped_alignment_records": mapped,
             "expected_bams_verified": [str(path.relative_to(ROOT)) for path in expected_bams],
             "bam_indexes_verified": bam_indexes,
             "table_contracts_verified": table_contracts,
             "peak_counts": peak_counts,
             "assay_outputs_verified": assay_outputs,
             "atac_run_merges_verified": run_merges,
             "artifacts": artifacts, "output_bytes": sum(p["bytes"] for p in artifacts), "outputs_deleted": False}
    evidence = args.evidence or ROOT / "work" / "evidence" / f"{args.variant}.json"
    evidence.parent.mkdir(parents=True, exist_ok=True)
    results = json.loads(evidence.read_text()) if evidence.exists() else {"runs": []}
    results["runs"].append(entry)
    evidence.write_text(json.dumps(results, indent=2) + "\n")
    if args.keep_outputs:
        print(f"{args.variant}: verified {len(expected_bams)} BAMs, assay outputs and MultiQC; retained results")
        return
    # Work is selected from known configs and is always inside this fixture tree.
    # Staged WGBS references and indexes are inside the variant work directory.
    # Never delete files next to the source FASTA.
    shutil.rmtree(work)
    entry["outputs_deleted"] = True
    evidence.write_text(json.dumps(results, indent=2) + "\n")
    print(f"{args.variant}: verified {len(expected_bams)} BAMs and MultiQC; deleted {len(files)} output files ({entry['output_bytes']:,} bytes)")


if __name__ == "__main__":
    main()
