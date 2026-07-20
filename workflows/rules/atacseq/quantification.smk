import csv
import os
import re


FILTERED_BAM_DIR = config["filtered_bam_dir"]
PEAKS_DIR = config["peaks_dir"]
CONSENSUS_PEAKS_DIR = config["consensus_peaks_dir"]
FEATURECOUNTS_DIR = config["featurecounts_dir"]
QC_SUMMARY_DIR = config["qc_summary_dir"]
BAM_QC_DIR = config["bam_qc_dir"]
FILTERED_BAM_QC_DIR = config["filtered_bam_qc_dir"]
DUPLICATION_QC_DIR = config["duplication_qc_dir"]

CONSENSUS_CONFIG = config.get("consensus_peaks", {})
CONSENSUS_MINIMUM_SUPPORT = int(CONSENSUS_CONFIG.get("minimum_support", 1))
if not 1 <= CONSENSUS_MINIMUM_SUPPORT <= len(SAMPLES):
    raise ValueError(
        "ATAC-seq consensus_peaks.minimum_support must be between 1 and "
        f"the number of samples ({len(SAMPLES)})."
    )

FEATURECOUNTS_CONFIG = config.get("featurecounts", {})
FEATURECOUNTS_THREADS = int(FEATURECOUNTS_CONFIG.get("threads", 4))
FEATURECOUNTS_MIN_MAPQ = int(FEATURECOUNTS_CONFIG.get("minimum_mapping_quality", 0))
FEATURECOUNTS_EXTRA_ARGS = FEATURECOUNTS_CONFIG.get("extra_args", "")
if FEATURECOUNTS_THREADS < 1:
    raise ValueError("ATAC-seq featurecounts.threads must be at least 1.")
if FEATURECOUNTS_MIN_MAPQ < 0:
    raise ValueError("ATAC-seq featurecounts.minimum_mapping_quality cannot be negative.")

MITOCHONDRIAL_CONTIG_CONFIG = config.get("mitochondrial_contigs", ["MT", "chrM"])
if (
    not isinstance(MITOCHONDRIAL_CONTIG_CONFIG, (list, tuple))
    or not MITOCHONDRIAL_CONTIG_CONFIG
    or not all(isinstance(contig, str) and contig for contig in MITOCHONDRIAL_CONTIG_CONFIG)
):
    raise ValueError(
        "ATAC-seq mitochondrial_contigs must be a non-empty list of contig names."
    )
MITOCHONDRIAL_CONTIGS = set(MITOCHONDRIAL_CONTIG_CONFIG)

CONSENSUS_BED = os.path.join(CONSENSUS_PEAKS_DIR, "consensus_peaks.bed")
CONSENSUS_SAF = os.path.join(CONSENSUS_PEAKS_DIR, "consensus_peaks.saf")
FEATURECOUNTS_MATRIX = os.path.join(FEATURECOUNTS_DIR, "consensus_peak_counts.txt")
FEATURECOUNTS_SUMMARY = f"{FEATURECOUNTS_MATRIX}.summary"
ATAC_QC_SUMMARY = os.path.join(QC_SUMMARY_DIR, "sample_qc_summary.tsv")
ATAC_QC_MULTIQC = os.path.join(QC_SUMMARY_DIR, "atac_qc_summary_mqc.tsv")
SORTED_PEAK_PATTERN = os.path.join(CONSENSUS_PEAKS_DIR, "{sample}.peaks.sorted.bed")


rule sort_peak_bed:
    """
    Normalizes caller-specific narrowPeak output to sorted BED3 intervals.
    """
    input:
        peak = os.path.join(PEAKS_DIR, "{sample}_peaks.narrowPeak")
    output:
        bed = temp(SORTED_PEAK_PATTERN)
    threads: 1
    log:
        os.path.join("logs", config["pipeline"], "consensus_peaks", "{sample}.sort.log")
    shell:
        r"""
        set -euo pipefail
        (cut -f 1-3 {input.peak:q} | \
        pixi run bedtools sort -i - > {output.bed}) > {log}.out 2> {log}.err
        test -s {output.bed}
        """


rule consensus_peaks:
    """
    Builds a configurable sample-supported union of per-sample peak calls.
    """
    input:
        peaks = expand(SORTED_PEAK_PATTERN, sample=SAMPLES)
    output:
        bed = CONSENSUS_BED
    params:
        minimum_support = CONSENSUS_MINIMUM_SUPPORT,
        sample_count = len(SAMPLES)
    threads: 1
    log:
        os.path.join("logs", config["pipeline"], "consensus_peaks", "consensus_peaks.log")
    shell:
        r"""
        set -euo pipefail
        if [ {params.sample_count} -eq 1 ]; then
            (pixi run bedtools merge -i {input.peaks:q} > {output.bed}) > {log}.out 2> {log}.err
        else
            (pixi run bedtools multiinter -i {input.peaks:q} | \
            awk -v minimum_support={params.minimum_support} 'BEGIN {{OFS="\t"}} $4 >= minimum_support {{print $1, $2, $3}}' | \
            pixi run bedtools merge -i - > {output.bed}) > {log}.out 2> {log}.err
        fi
        test -s {output.bed}
        """


rule consensus_peaks_to_saf:
    """
    Converts zero-based BED consensus peaks to one-based inclusive SAF.
    """
    input:
        bed = CONSENSUS_BED
    output:
        saf = CONSENSUS_SAF
    run:
        with open(input.bed) as source, open(output.saf, "w") as destination:
            destination.write("GeneID\tChr\tStart\tEnd\tStrand\n")
            peak_number = 0
            for line in source:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 3:
                    continue
                peak_number += 1
                destination.write(
                    f"peak_{peak_number}\t{fields[0]}\t{int(fields[1]) + 1}\t{fields[2]}\t.\n"
                )
        if peak_number == 0:
            raise ValueError("The ATAC-seq consensus peak BED contains no intervals.")


rule featurecounts:
    """
    Counts paired-end fragments from all filtered BAMs over consensus peaks.
    """
    input:
        saf = CONSENSUS_SAF,
        bams = expand(
            os.path.join(FILTERED_BAM_DIR, "{sample}_pe.filtered.sorted.bam"),
            sample=SAMPLES,
        ),
        bais = expand(
            os.path.join(FILTERED_BAM_DIR, "{sample}_pe.filtered.sorted.bam.bai"),
            sample=SAMPLES,
        )
    output:
        matrix = FEATURECOUNTS_MATRIX,
        summary = FEATURECOUNTS_SUMMARY
    params:
        minimum_mapping_quality = FEATURECOUNTS_MIN_MAPQ,
        extra = FEATURECOUNTS_EXTRA_ARGS
    threads: FEATURECOUNTS_THREADS
    log:
        os.path.join("logs", config["pipeline"], "featurecounts", "featurecounts.log")
    shell:
        r"""
        pixi run featureCounts \
            -p --countReadPairs -B -C \
            -Q {params.minimum_mapping_quality} \
            -F SAF \
            -T {threads} \
            -a {input.saf} \
            -o {output.matrix} \
            {params.extra} \
            {input.bams:q} > {log}.out 2> {log}.err
        """


def _atac_flagstat_total(path):
    with open(path) as handle:
        for line in handle:
            if " in total " in line:
                match = re.match(r"\s*([0-9,]+)\s+\+", line)
                if match:
                    return int(match.group(1).replace(",", ""))
    return None


def _atac_sambamba_duplication_metrics(path):
    with open(path) as handle:
        contents = handle.read()
    patterns = {
        "sorted_end_pairs": r"sorted (\d+) end pairs",
        "single_ends": r"and (\d+) single ends",
        "unmatched_pairs": r"among them (\d+) unmatched",
        "duplicate_reads": r"found (\d+) duplicates",
    }
    metrics = {}
    for name, pattern in patterns.items():
        match = re.search(pattern, contents)
        if not match:
            return None, None
        metrics[name] = int(match.group(1))

    examined_reads = (
        2 * metrics["sorted_end_pairs"]
        + metrics["single_ends"]
        - metrics["unmatched_pairs"]
    )
    duplication_rate = (
        100 * metrics["duplicate_reads"] / examined_reads if examined_reads else None
    )
    return metrics["duplicate_reads"], duplication_rate


def _atac_idxstats_metrics(path):
    total_mapped = 0
    mitochondrial_mapped = 0
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4 or fields[0] == "*":
                continue
            mapped = int(fields[2])
            total_mapped += mapped
            if fields[0] in MITOCHONDRIAL_CONTIGS:
                mitochondrial_mapped += mapped
    mitochondrial_fraction = (
        100 * mitochondrial_mapped / total_mapped if total_mapped else None
    )
    return mitochondrial_mapped, mitochondrial_fraction


def _atac_featurecounts_metrics(path):
    with open(path, newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        rows = list(reader)
    if not rows or len(rows[0]) - 1 != len(SAMPLES):
        raise ValueError("featureCounts summary columns do not match the configured ATAC-seq samples.")
    per_sample = {sample: {"assigned": 0, "total": 0} for sample in SAMPLES}
    for row in rows[1:]:
        if not row:
            continue
        for index, sample in enumerate(SAMPLES, start=1):
            count = int(row[index])
            per_sample[sample]["total"] += count
            if row[0] == "Assigned":
                per_sample[sample]["assigned"] = count
    return per_sample


def _atac_format_metric(value, digits=4):
    if value is None:
        return "NA"
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    return str(value)


rule atac_qc_summary:
    """
    Consolidates core ATAC-seq QC metrics and emits MultiQC custom content.
    """
    input:
        aligned_flagstats = expand(
            os.path.join(BAM_QC_DIR, "{sample}_pe.flagstat.txt"), sample=SAMPLES
        ),
        filtered_flagstats = expand(
            os.path.join(FILTERED_BAM_QC_DIR, "{sample}_pe.filtered.flagstat.txt"),
            sample=SAMPLES,
        ),
        duplication_metrics = expand(
            os.path.join(DUPLICATION_QC_DIR, "{sample}_pe.markdup.metrics.txt"),
            sample=SAMPLES,
        ),
        idxstats = expand(
            os.path.join(BAM_QC_DIR, "{sample}_pe.idxstats.txt"), sample=SAMPLES
        ),
        featurecounts_summary = FEATURECOUNTS_SUMMARY
    output:
        table = ATAC_QC_SUMMARY,
        multiqc = ATAC_QC_MULTIQC
    run:
        featurecounts_metrics = _atac_featurecounts_metrics(input.featurecounts_summary)
        fieldnames = [
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
        rows = []
        for index, sample in enumerate(SAMPLES):
            duplicate_reads, duplication_rate = _atac_sambamba_duplication_metrics(
                input.duplication_metrics[index]
            )
            mitochondrial_reads, mitochondrial_fraction = _atac_idxstats_metrics(
                input.idxstats[index]
            )
            assigned = featurecounts_metrics[sample]["assigned"]
            counted_total = featurecounts_metrics[sample]["total"]
            frip = 100 * assigned / counted_total if counted_total else None
            rows.append(
                {
                    "sample": sample,
                    "aligned_reads": _atac_format_metric(
                        _atac_flagstat_total(input.aligned_flagstats[index])
                    ),
                    "final_filtered_reads": _atac_format_metric(
                        _atac_flagstat_total(input.filtered_flagstats[index])
                    ),
                    "duplicate_reads": _atac_format_metric(duplicate_reads),
                    "duplication_rate_pct": _atac_format_metric(duplication_rate),
                    "mitochondrial_reads": _atac_format_metric(mitochondrial_reads),
                    "mitochondrial_fraction_pct": _atac_format_metric(
                        mitochondrial_fraction
                    ),
                    "fragments_in_peaks": _atac_format_metric(assigned),
                    "featurecounts_total_fragments": _atac_format_metric(counted_total),
                    "frip_pct": _atac_format_metric(frip),
                }
            )

        with open(output.table, "w", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(rows)

        with open(output.multiqc, "w", newline="") as handle:
            handle.write("# id: 'atac_qc_summary'\n")
            handle.write("# section_name: 'ATAC-seq QC Summary'\n")
            handle.write(
                "# description: 'Core alignment, filtering, duplication, mitochondrial, and FRiP metrics.'\n"
            )
            handle.write("# format: 'tsv'\n")
            handle.write("# plot_type: 'table'\n")
            handle.write("# pconfig:\n")
            handle.write("#   id: 'atac_qc_summary_table'\n")
            handle.write("#   title: 'ATAC-seq QC Summary'\n")
            writer = csv.DictWriter(
                handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(rows)
