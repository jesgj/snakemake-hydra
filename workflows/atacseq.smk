import os
import sys

sys.path.insert(0, os.path.abspath("src"))
from utils import prepare_sample_data

# --- CONFIGURATION ---
RAW_DIR = config["raw_fastqs_dir"]
QC_DIR = config["qc_dir"]
TRIMMED_DIR = config["trimmed_dir"]
QC_TRIMMED_DIR = config["qc_trimmed_dir"]
REF_GENOME = config["ref_genome"]
BOWTIE2_INDEX_DIR = config["bowtie2_index_dir"]
ALIGNMENT_DIR = config["alignment_dir"]
BAM_QC_DIR = config["bam_qc_dir"]
FILTERED_BAM_DIR = config["filtered_bam_dir"]
FILTERED_BAM_QC_DIR = config["filtered_bam_qc_dir"]
FRAGMENT_QC_DIR = config["fragment_qc_dir"]
BIGWIG_DIR = config["bigwig_dir"]
PEAKS_DIR = config["peaks_dir"]

# Make config available to included rules
config["raw_fastqs_dir"] = RAW_DIR
config["qc_dir"] = QC_DIR
config["trimmed_dir"] = TRIMMED_DIR
config["qc_trimmed_dir"] = QC_TRIMMED_DIR
config["ref_genome"] = REF_GENOME
config["bowtie2_index_dir"] = BOWTIE2_INDEX_DIR
config["alignment_dir"] = ALIGNMENT_DIR
config["bam_qc_dir"] = BAM_QC_DIR
config["filtered_bam_dir"] = FILTERED_BAM_DIR
config["filtered_bam_qc_dir"] = FILTERED_BAM_QC_DIR
config["fragment_qc_dir"] = FRAGMENT_QC_DIR
config["bigwig_dir"] = BIGWIG_DIR
config["peaks_dir"] = PEAKS_DIR

# Ensure output directories exist
os.makedirs(QC_DIR, exist_ok=True)
os.makedirs(TRIMMED_DIR, exist_ok=True)
os.makedirs(QC_TRIMMED_DIR, exist_ok=True)
os.makedirs(BOWTIE2_INDEX_DIR, exist_ok=True)
os.makedirs(ALIGNMENT_DIR, exist_ok=True)
os.makedirs(BAM_QC_DIR, exist_ok=True)
os.makedirs(FILTERED_BAM_DIR, exist_ok=True)
os.makedirs(FILTERED_BAM_QC_DIR, exist_ok=True)
os.makedirs(FRAGMENT_QC_DIR, exist_ok=True)
os.makedirs(BIGWIG_DIR, exist_ok=True)
os.makedirs(PEAKS_DIR, exist_ok=True)

# Ensure log directories exist
os.makedirs(os.path.join("logs", config["pipeline"], "fastqc_raw"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "fastp"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "fastqc_trimmed"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "bowtie2_build"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "bowtie2_align"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "bam_qc"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "sambamba_filter"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "fragment_qc"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "bamCoverage"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "macs3"), exist_ok=True)
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


def get_atacseq_outputs(samples):
    outputs = []
    outputs.extend(expand(os.path.join(QC_DIR, "{sample}_R1_raw_fastqc.html"), sample=samples))
    outputs.extend(expand(os.path.join(QC_DIR, "{sample}_R2_raw_fastqc.html"), sample=samples))
    outputs.extend(expand(os.path.join(QC_TRIMMED_DIR, "{sample}_R1_trimmed_fastqc.html"), sample=samples))
    outputs.extend(expand(os.path.join(QC_TRIMMED_DIR, "{sample}_R2_trimmed_fastqc.html"), sample=samples))
    outputs.extend(expand(os.path.join(ALIGNMENT_DIR, "{sample}_pe.sorted.bam"), sample=samples))
    outputs.extend(expand(os.path.join(BAM_QC_DIR, "{sample}_pe.stats.txt"), sample=samples))
    outputs.extend(expand(os.path.join(BAM_QC_DIR, "{sample}_pe.flagstat.txt"), sample=samples))
    outputs.extend(expand(os.path.join(FILTERED_BAM_DIR, "{sample}_pe.filtered.sorted.bam"), sample=samples))
    outputs.extend(expand(os.path.join(FILTERED_BAM_QC_DIR, "{sample}_pe.filtered.stats.txt"), sample=samples))
    outputs.extend(expand(os.path.join(FILTERED_BAM_QC_DIR, "{sample}_pe.filtered.flagstat.txt"), sample=samples))
    outputs.extend(expand(os.path.join(FRAGMENT_QC_DIR, "{sample}_pe.insert_size_metrics.txt"), sample=samples))
    outputs.extend(expand(os.path.join(FRAGMENT_QC_DIR, "{sample}_pe.insert_size_histogram.pdf"), sample=samples))
    outputs.extend(expand(os.path.join(BIGWIG_DIR, "{sample}_pe.bw"), sample=samples))
    outputs.extend(expand(os.path.join(PEAKS_DIR, "{sample}_peaks.narrowPeak"), sample=samples))
    outputs.extend(expand(os.path.join(PEAKS_DIR, "{sample}_summits.bed"), sample=samples))
    outputs.extend(expand(os.path.join(PEAKS_DIR, "{sample}_treat_pileup.bdg"), sample=samples))
    outputs.extend(expand(os.path.join(PEAKS_DIR, "{sample}_peaks.xls"), sample=samples))
    return outputs


# --- MultiQC Configuration ---
config["pipeline_name"] = "atacseq"
config["multiqc_results_dir"] = "results/atacseq"
config["multiqc_input_files"] = get_atacseq_outputs(SAMPLES)


# --- MODULE INCLUSION ---
include: "rules/qc.smk"
include: "rules/trimming_and_qc.smk"
include: "rules/atacseq/alignment.smk"
include: "rules/bam_qc.smk"
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


# --- FINAL TARGETS ---
rule all:
    input:
        get_atacseq_outputs(SAMPLES),
        os.path.join(config["multiqc_results_dir"], "multiqc_report.html")
    default_target: True
