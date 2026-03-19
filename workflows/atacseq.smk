import os
import sys

sys.path.insert(0, os.path.abspath("src"))
from utils import prepare_sample_data, validate_existing_file, validate_required_path

# --- CONFIGURATION ---
RAW_DIR = config["raw_fastqs_dir"]
QC_DIR = config.get("qc_dir", os.path.join("results", "atacseq", "qc_raw"))
TRIMMED_DIR = config.get("trimmed_dir", os.path.join("results", "atacseq", "trimmed_fastqs"))
QC_TRIMMED_DIR = config.get("qc_trimmed_dir", os.path.join("results", "atacseq", "qc_trimmed"))
REF_GENOME = config.get("ref_genome")
BOWTIE2_INDEX_DIR = config.get("bowtie2_index_dir")
ALIGNMENT_DIR = config.get("alignment_dir", os.path.join("results", "atacseq", "aligned_bams"))
MARKED_BAM_DIR = config.get("marked_bam_dir", os.path.join("results", "atacseq", "marked_bams"))
DUPLICATION_QC_DIR = config.get("duplication_qc_dir", os.path.join("results", "atacseq", "duplication_qc"))
BAM_QC_DIR = config.get("bam_qc_dir", os.path.join("results", "atacseq", "bam_qc"))
FILTERED_BAM_DIR = config.get("filtered_bam_dir", os.path.join("results", "atacseq", "filtered_bams"))
FILTERED_BAM_QC_DIR = config.get("filtered_bam_qc_dir", os.path.join("results", "atacseq", "filtered_bam_qc"))
DEEPTOOLS_DIR = config.get("deeptools_dir", os.path.join("results", "atacseq", "deeptools"))
FRAGMENT_QC_DIR = config.get("fragment_qc_dir", os.path.join("results", "atacseq", "fragment_qc"))
BIGWIG_DIR = config.get("bigwig_dir", os.path.join("results", "atacseq", "bigwigs"))
PEAKS_DIR = config.get("peaks_dir", os.path.join("results", "atacseq", "peaks"))
PEAK_CALLER = config.get("peak_caller", "macs3")

# Reference inputs are user-managed. Only the configured Bowtie2 index directory is created here.
validate_existing_file(REF_GENOME, "ref_genome")
validate_required_path(BOWTIE2_INDEX_DIR, "bowtie2_index_dir")

if PEAK_CALLER not in ["macs3", "genrich"]:
    raise ValueError("ATAC-seq 'peak_caller' must be either 'macs3' or 'genrich'.")

# Make config available to included rules
config["raw_fastqs_dir"] = RAW_DIR
config["qc_dir"] = QC_DIR
config["trimmed_dir"] = TRIMMED_DIR
config["qc_trimmed_dir"] = QC_TRIMMED_DIR
config["ref_genome"] = REF_GENOME
config["bowtie2_index_dir"] = BOWTIE2_INDEX_DIR
config["alignment_dir"] = ALIGNMENT_DIR
config["marked_bam_dir"] = MARKED_BAM_DIR
config["duplication_qc_dir"] = DUPLICATION_QC_DIR
config["bam_qc_dir"] = BAM_QC_DIR
config["filtered_bam_dir"] = FILTERED_BAM_DIR
config["filtered_bam_qc_dir"] = FILTERED_BAM_QC_DIR
config["deeptools_dir"] = DEEPTOOLS_DIR
config["fragment_qc_dir"] = FRAGMENT_QC_DIR
config["bigwig_dir"] = BIGWIG_DIR
config["peaks_dir"] = PEAKS_DIR
config["peak_caller"] = PEAK_CALLER

# Ensure output directories exist
os.makedirs(QC_DIR, exist_ok=True)
os.makedirs(TRIMMED_DIR, exist_ok=True)
os.makedirs(QC_TRIMMED_DIR, exist_ok=True)
os.makedirs(BOWTIE2_INDEX_DIR, exist_ok=True)
os.makedirs(ALIGNMENT_DIR, exist_ok=True)
os.makedirs(MARKED_BAM_DIR, exist_ok=True)
os.makedirs(DUPLICATION_QC_DIR, exist_ok=True)
os.makedirs(BAM_QC_DIR, exist_ok=True)
os.makedirs(FILTERED_BAM_DIR, exist_ok=True)
os.makedirs(FILTERED_BAM_QC_DIR, exist_ok=True)
os.makedirs(DEEPTOOLS_DIR, exist_ok=True)
os.makedirs(FRAGMENT_QC_DIR, exist_ok=True)
os.makedirs(BIGWIG_DIR, exist_ok=True)
os.makedirs(PEAKS_DIR, exist_ok=True)

# Ensure log directories exist
os.makedirs(os.path.join("logs", config["pipeline"], "fastqc_raw"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "fastp"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "fastqc_trimmed"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "bowtie2_build"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "bowtie2_align"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "markduplicates"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "bam_qc"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "samtools_index_filtered"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "sambamba_filter"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "deeptools"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "fragment_qc"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "bamCoverage"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], PEAK_CALLER), exist_ok=True)
if PEAK_CALLER == "genrich":
    os.makedirs(os.path.join("logs", config["pipeline"], "genrich_qname_sort"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "multiqc"), exist_ok=True)

# --- SAMPLE DISCOVERY ---
SAMPLES_INFO, SAMPLES = prepare_sample_data(config)

if not SAMPLES:
    raise ValueError("No ATAC-seq samples were found. Provide 'samples_info' or valid FASTQs in 'raw_fastqs_dir'.")

invalid_samples = []
for sample, info in SAMPLES_INFO.items():
    has_pair = "R1" in info and "R2" in info
    declared_type = info.get("type")
    if has_pair and declared_type in [None, "PE"]:
        info["type"] = "PE"
        continue
    invalid_samples.append(sample)

config["samples_info"] = SAMPLES_INFO
if invalid_samples:
    raise ValueError(
        "ATAC-seq currently supports paired-end samples only (type=PE with R1 and R2). "
        f"Invalid samples: {', '.join(invalid_samples)}"
    )


def get_all_filtered_bams(samples):
    return [os.path.join(FILTERED_BAM_DIR, f"{sample}_pe.filtered.sorted.bam") for sample in samples]


def get_all_filtered_bais(samples):
    return [f"{bam}.bai" for bam in get_all_filtered_bams(samples)]


ALL_FILTERED_BAMS = get_all_filtered_bams(SAMPLES)
ALL_FILTERED_BAIS = get_all_filtered_bais(SAMPLES)
HAS_MULTI_BAM_DEEPTOOLS = len(ALL_FILTERED_BAMS) >= 2


def get_atacseq_outputs(samples):
    outputs = []
    outputs.extend(expand(os.path.join(QC_DIR, "{sample}_R1_raw_fastqc.html"), sample=samples))
    outputs.extend(expand(os.path.join(QC_DIR, "{sample}_R2_raw_fastqc.html"), sample=samples))
    outputs.extend(expand(os.path.join(QC_TRIMMED_DIR, "{sample}_R1_trimmed_fastqc.html"), sample=samples))
    outputs.extend(expand(os.path.join(QC_TRIMMED_DIR, "{sample}_R2_trimmed_fastqc.html"), sample=samples))
    outputs.extend(expand(os.path.join(ALIGNMENT_DIR, "{sample}_pe.sorted.bam"), sample=samples))
    outputs.extend(expand(os.path.join(DUPLICATION_QC_DIR, "{sample}_pe.markdup.metrics.txt"), sample=samples))
    outputs.extend(expand(os.path.join(BAM_QC_DIR, "{sample}_pe.stats.txt"), sample=samples))
    outputs.extend(expand(os.path.join(BAM_QC_DIR, "{sample}_pe.flagstat.txt"), sample=samples))
    outputs.extend(expand(os.path.join(FILTERED_BAM_DIR, "{sample}_pe.filtered.sorted.bam"), sample=samples))
    outputs.extend(expand(os.path.join(FILTERED_BAM_QC_DIR, "{sample}_pe.filtered.stats.txt"), sample=samples))
    outputs.extend(expand(os.path.join(FILTERED_BAM_QC_DIR, "{sample}_pe.filtered.flagstat.txt"), sample=samples))
    outputs.append(os.path.join(DEEPTOOLS_DIR, "fingerprints.png"))
    outputs.append(os.path.join(DEEPTOOLS_DIR, "fingerprints.metrics.tab"))
    if len(samples) >= 2:
        outputs.append(os.path.join(DEEPTOOLS_DIR, "bam_correlation_heatmap.png"))
        outputs.append(os.path.join(DEEPTOOLS_DIR, "bam_correlation_matrix.tab"))
    outputs.extend(expand(os.path.join(FRAGMENT_QC_DIR, "{sample}_pe.insert_size_metrics.txt"), sample=samples))
    outputs.extend(expand(os.path.join(FRAGMENT_QC_DIR, "{sample}_pe.insert_size_histogram.pdf"), sample=samples))
    outputs.extend(expand(os.path.join(BIGWIG_DIR, "{sample}_pe.bw"), sample=samples))
    outputs.extend(expand(os.path.join(PEAKS_DIR, "{sample}_peaks.narrowPeak"), sample=samples))
    if PEAK_CALLER == "macs3":
        outputs.extend(expand(os.path.join(PEAKS_DIR, "{sample}_summits.bed"), sample=samples))
        outputs.extend(expand(os.path.join(PEAKS_DIR, "{sample}_treat_pileup.bdg"), sample=samples))
        outputs.extend(expand(os.path.join(PEAKS_DIR, "{sample}_peaks.xls"), sample=samples))
    return outputs


def get_atacseq_multiqc_analysis_dirs():
    return list(
        dict.fromkeys(
            [
                QC_DIR,
                TRIMMED_DIR,
                QC_TRIMMED_DIR,
                ALIGNMENT_DIR,
                DUPLICATION_QC_DIR,
                BAM_QC_DIR,
                FILTERED_BAM_DIR,
                FILTERED_BAM_QC_DIR,
                DEEPTOOLS_DIR,
                FRAGMENT_QC_DIR,
                BIGWIG_DIR,
                PEAKS_DIR,
                os.path.join("logs", config["pipeline"]),
            ]
        )
    )


# --- MultiQC Configuration ---
config["pipeline_name"] = "atacseq"
config["multiqc_results_dir"] = config.get("multiqc_results_dir")
if not config["multiqc_results_dir"]:
    raise ValueError("ATAC-seq requires 'multiqc_results_dir' in the config.")
config["multiqc_input_files"] = get_atacseq_outputs(SAMPLES)
config["multiqc_analysis_dirs"] = get_atacseq_multiqc_analysis_dirs()


# --- MODULE INCLUSION ---
include: "rules/qc.smk"
include: "rules/trimming_and_qc.smk"
include: "rules/atacseq/alignment.smk"
include: "rules/bam_qc.smk"
include: "rules/deeptools_qc.smk"
include: "rules/atacseq/markduplicates.smk"
include: "rules/atacseq/filter_bam.smk"
include: "rules/atacseq/fragment_qc.smk"
include: "rules/atacseq/bigwig.smk"
include: "rules/atacseq/peak_calling.smk"
include: "rules/multiqc.smk"


# --- BAM QC INSTANTIATION ---
use rule samtools_stats_generic as samtools_stats_aligned with:
    input:
        bam = os.path.join(ALIGNMENT_DIR, "{sample}_pe.sorted.bam")
    output:
        stats = os.path.join(BAM_QC_DIR, "{sample}_pe.stats.txt")
    log:
        os.path.join("logs", config["pipeline"], "bam_qc", "{sample}_pe_stats.log")

use rule samtools_flagstat_generic as samtools_flagstat_aligned with:
    input:
        bam = os.path.join(ALIGNMENT_DIR, "{sample}_pe.sorted.bam")
    output:
        flagstat = os.path.join(BAM_QC_DIR, "{sample}_pe.flagstat.txt")
    log:
        os.path.join("logs", config["pipeline"], "bam_qc", "{sample}_pe_flagstat.log")

use rule samtools_stats_generic as samtools_stats_filtered with:
    input:
        bam = os.path.join(FILTERED_BAM_DIR, "{sample}_pe.filtered.sorted.bam")
    output:
        stats = os.path.join(FILTERED_BAM_QC_DIR, "{sample}_pe.filtered.stats.txt")
    log:
        os.path.join("logs", config["pipeline"], "bam_qc", "{sample}_pe_filtered_stats.log")

use rule samtools_flagstat_generic as samtools_flagstat_filtered with:
    input:
        bam = os.path.join(FILTERED_BAM_DIR, "{sample}_pe.filtered.sorted.bam")
    output:
        flagstat = os.path.join(FILTERED_BAM_QC_DIR, "{sample}_pe.filtered.flagstat.txt")
    log:
        os.path.join("logs", config["pipeline"], "bam_qc", "{sample}_pe_filtered_flagstat.log")

use rule samtools_index_bam_generic as samtools_index_filtered_bam with:
    input:
        bam = os.path.join(FILTERED_BAM_DIR, "{sample}_pe.filtered.sorted.bam")
    output:
        bai = os.path.join(FILTERED_BAM_DIR, "{sample}_pe.filtered.sorted.bam.bai")
    log:
        os.path.join("logs", config["pipeline"], "samtools_index_filtered", "{sample}_pe.log")

MBS_ARGS = config.get("deeptools", {}).get("multiBamSummary", {}).get("extra_args", "--binSize 10000")
PC_ARGS = config.get("deeptools", {}).get("plotCorrelation", {}).get("extra_args", "-p heatmap --corMethod spearman --skipZeros")
PF_ARGS = config.get("deeptools", {}).get("plotFingerprint", {}).get("extra_args", "")

use rule plotFingerprint_generic as plotFingerprint with:
    input:
        bams = ALL_FILTERED_BAMS,
        bais = ALL_FILTERED_BAIS
    output:
        plot = os.path.join(DEEPTOOLS_DIR, "fingerprints.png"),
        metrics = os.path.join(DEEPTOOLS_DIR, "fingerprints.metrics.tab")
    params:
        extra = PF_ARGS
    log:
        os.path.join("logs", config["pipeline"], "deeptools", "plotFingerprint.log")

if HAS_MULTI_BAM_DEEPTOOLS:
    use rule multiBamSummary_generic as multiBamSummary with:
        input:
            bams = ALL_FILTERED_BAMS,
            bais = ALL_FILTERED_BAIS
        output:
            npz = os.path.join(DEEPTOOLS_DIR, "read_coverage.npz")
        params:
            extra = MBS_ARGS
        log:
            os.path.join("logs", config["pipeline"], "deeptools", "multiBamSummary.log")

    use rule plotCorrelation_generic as plotCorrelation with:
        input:
            npz = os.path.join(DEEPTOOLS_DIR, "read_coverage.npz")
        output:
            heatmap = os.path.join(DEEPTOOLS_DIR, "bam_correlation_heatmap.png"),
            matrix = os.path.join(DEEPTOOLS_DIR, "bam_correlation_matrix.tab")
        params:
            extra = PC_ARGS
        log:
            os.path.join("logs", config["pipeline"], "deeptools", "plotCorrelation.log")


# --- FINAL TARGETS ---
rule all:
    input:
        get_atacseq_outputs(SAMPLES),
        os.path.join(config["multiqc_results_dir"], "multiqc_report.html")
    default_target: True
