# WGBS module workflow
import os
import sys

sys.path.insert(0, os.path.abspath("src"))
from utils import (
    fastq_validation_scope,
    normalize_read_group_metadata,
    prepare_sample_data,
    validate_required_path,
)

# --- CONFIGURATION ---
# Define directories from config
RAW_DIR = config["raw_fastqs_dir"]
QC_DIR = config["qc_dir"]
TRIMMED_DIR = config["trimmed_dir"]
QC_TRIMMED_DIR = config["qc_trimmed_dir"]
FASTQ_PAIR_VALIDATION_DIR = config.get(
    "fastq_pair_validation_dir",
    os.path.join("results", "wgbs", "fastq_pair_validation"),
)
READ_GROUP_QC_DIR = config.get(
    "read_group_qc_dir", os.path.join("results", "wgbs", "read_group_qc")
)
FASTQ_VALIDATION_SCOPE = fastq_validation_scope(config, "wgbs.fastq_validation")
ALIGN_DIR = config["alignment_dir"]
DEDUP_DIR = config["dedup_dir"]
DEDUP_BAM_QC_DIR = config["dedup_bam_qc_dir"]
FILTERED_BAM_DIR = config["filtered_bam_dir"]
SORTED_FILTERED_BAM_DIR = config["sorted_filtered_bam_dir"]
FILTERED_BAM_QC_DIR = config["filtered_bam_qc_dir"]
MBIAS_DIR = config["mbias_dir"]
METHYLDACKEL_DIR = config["methyldackel_dir"]
METHYLDACKEL_MERGECONTEXT_DIR = config["methyldackel_mergecontext_dir"]
REF_GENOME = config.get("ref_genome")
REFERENCE_DIR = config.get("reference_dir", os.path.join("results", "wgbs", "reference"))
LOG_DIR = config.get("log_dir", os.path.join("logs", config["pipeline"]))

# The configured reference FASTA is user-managed and read-only. Bismark and the
# methylation tools use a workflow-owned staged copy so they cannot index extra
# FASTAs or write .fai/index files into the source reference directory.
validate_required_path(REF_GENOME, "ref_genome")
STAGED_REF_GENOME = os.path.join(REFERENCE_DIR, "genome.fa")
SOURCE_REF_PATH = os.path.realpath(REF_GENOME)
INDEX_BUILD_PATH = os.path.realpath(os.path.join(REFERENCE_DIR, "Bisulfite_Genome"))
if (
    os.path.realpath(REFERENCE_DIR) == os.path.dirname(SOURCE_REF_PATH)
    or os.path.realpath(STAGED_REF_GENOME) == SOURCE_REF_PATH
    or os.path.commonpath([SOURCE_REF_PATH, INDEX_BUILD_PATH]) == INDEX_BUILD_PATH
):
    raise ValueError(
        "WGBS reference_dir must be a workflow-owned directory separate from the "
        f"source FASTA directory: {REF_GENOME}"
    )

# Make config available to included rules
config["raw_fastqs_dir"] = RAW_DIR
config["qc_dir"] = QC_DIR
config["trimmed_dir"] = TRIMMED_DIR
config["qc_trimmed_dir"] = QC_TRIMMED_DIR
config["fastq_pair_validation_dir"] = FASTQ_PAIR_VALIDATION_DIR
config["read_group_qc_dir"] = READ_GROUP_QC_DIR
config["fastq_validation_scope"] = FASTQ_VALIDATION_SCOPE
config["alignment_dir"] = ALIGN_DIR
config["dedup_dir"] = DEDUP_DIR
config["dedup_bam_qc_dir"] = DEDUP_BAM_QC_DIR
config["filtered_bam_dir"] = FILTERED_BAM_DIR
config["sorted_filtered_bam_dir"] = SORTED_FILTERED_BAM_DIR
config["filtered_bam_qc_dir"] = FILTERED_BAM_QC_DIR
config["mbias_dir"] = MBIAS_DIR
config["methyldackel_dir"] = METHYLDACKEL_DIR
config["methyldackel_mergecontext_dir"] = METHYLDACKEL_MERGECONTEXT_DIR
config["ref_genome"] = REF_GENOME
config["reference_dir"] = REFERENCE_DIR
config["reference_fasta"] = STAGED_REF_GENOME
config["reference_fasta_fai"] = STAGED_REF_GENOME + ".fai"
config["log_dir"] = LOG_DIR


# Ensure output directories exist
os.makedirs(config["qc_dir"], exist_ok=True)
os.makedirs(config["trimmed_dir"], exist_ok=True)
os.makedirs(config["qc_trimmed_dir"], exist_ok=True)
os.makedirs(FASTQ_PAIR_VALIDATION_DIR, exist_ok=True)
os.makedirs(READ_GROUP_QC_DIR, exist_ok=True)
os.makedirs(config["alignment_dir"], exist_ok=True)
os.makedirs(config["dedup_dir"], exist_ok=True)
os.makedirs(config["dedup_bam_qc_dir"], exist_ok=True)
os.makedirs(config["filtered_bam_dir"], exist_ok=True)
os.makedirs(config["sorted_filtered_bam_dir"], exist_ok=True)
os.makedirs(config["filtered_bam_qc_dir"], exist_ok=True)
os.makedirs(config["mbias_dir"], exist_ok=True)
os.makedirs(config["methyldackel_dir"], exist_ok=True)
os.makedirs(config["methyldackel_mergecontext_dir"], exist_ok=True)

# Ensure log directories exist
os.makedirs(os.path.join(LOG_DIR, "fastqc_raw"), exist_ok=True)
os.makedirs(os.path.join(LOG_DIR, "fastp"), exist_ok=True)
os.makedirs(os.path.join(LOG_DIR, "fastqc_trimmed"), exist_ok=True)
os.makedirs(os.path.join(LOG_DIR, "fastq_pair_validation"), exist_ok=True)
os.makedirs(os.path.join(LOG_DIR, "read_group_validation"), exist_ok=True)
os.makedirs(os.path.join(LOG_DIR, "multiqc"), exist_ok=True)
os.makedirs(os.path.join(LOG_DIR, "bismark_reference"), exist_ok=True)
os.makedirs(os.path.join(LOG_DIR, "bismark_genome_preparation"), exist_ok=True)
os.makedirs(os.path.join(LOG_DIR, "bismark"), exist_ok=True)
os.makedirs(os.path.join(LOG_DIR, "deduplicate_bismark"), exist_ok=True)
os.makedirs(os.path.join(LOG_DIR, "sambamba_filter"), exist_ok=True)
os.makedirs(os.path.join(LOG_DIR, "sambamba_sort"), exist_ok=True)
os.makedirs(os.path.join(LOG_DIR, "methyldackel_mbias"), exist_ok=True)
os.makedirs(os.path.join(LOG_DIR, "methyldackel_extract_methylkit"), exist_ok=True)
os.makedirs(os.path.join(LOG_DIR, "methyldackel_extract_mergecontext"), exist_ok=True)
os.makedirs(os.path.join(LOG_DIR, "samtools_stats_dedup"), exist_ok=True)
os.makedirs(os.path.join(LOG_DIR, "samtools_flagstat_dedup"), exist_ok=True)
os.makedirs(os.path.join(LOG_DIR, "picard_collect_alignment_metrics_dedup"), exist_ok=True)
os.makedirs(os.path.join(LOG_DIR, "samtools_stats_filtered"), exist_ok=True)
os.makedirs(os.path.join(LOG_DIR, "samtools_flagstat_filtered"), exist_ok=True)
os.makedirs(os.path.join(LOG_DIR, "picard_collect_alignment_metrics_filtered"), exist_ok=True)
os.makedirs(os.path.join(LOG_DIR, "samtools_index_filtered"), exist_ok=True)


# --- SAMPLE DISCOVERY ---
SAMPLES_INFO, SAMPLES = prepare_sample_data(config)

if not SAMPLES:
    raise ValueError("No WGBS samples were found. Provide 'samples_info' or valid FASTQs in 'raw_fastqs_dir'.")

invalid_samples = []
for sample, info in SAMPLES_INFO.items():
    has_pair = "R1" in info and "R2" in info
    declared_type = info.get("type")
    if has_pair and declared_type in [None, "PE"]:
        info["type"] = "PE"
        continue
    invalid_samples.append(sample)

SAMPLES_INFO = normalize_read_group_metadata(SAMPLES_INFO, "wgbs.samples_info")
config["samples_info"] = SAMPLES_INFO
if invalid_samples:
    raise ValueError(
        "WGBS currently supports paired-end samples only (type=PE with R1 and R2). "
        f"Invalid samples: {', '.join(invalid_samples)}"
    )


# --- HELPER FUNCTION FOR OUTPUTS ---
def get_wgbs_outputs(samples):
    outputs = []
    outputs.extend([os.path.join(config["qc_dir"], f"{sample}_{read}_raw_fastqc.html") for sample in samples for read in ["R1", "R2"]])
    outputs.extend([os.path.join(config["qc_trimmed_dir"], f"{sample}_{read}_trimmed_fastqc.html") for sample in samples for read in ["R1", "R2"]])
    outputs.extend(expand(os.path.join(config["alignment_dir"], "{sample}_pe_report.txt"), sample=samples))
    outputs.extend(expand(os.path.join(config["dedup_dir"], "{sample}_pe.deduplicated.bam"), sample=samples))
    outputs.extend(expand(os.path.join(config["dedup_dir"], "{sample}_pe.deduplication_report.txt"), sample=samples))
    outputs.extend(expand(os.path.join(config["sorted_filtered_bam_dir"], "{sample}_pe.filtered.sorted.bam"), sample=samples))
    outputs.extend(expand(os.path.join(config["sorted_filtered_bam_dir"], "{sample}_pe.filtered.sorted.bam.bai"), sample=samples))
    outputs.extend(expand(os.path.join(config["dedup_bam_qc_dir"], "{sample}.dedup.stats.txt"), sample=samples))
    outputs.extend(expand(os.path.join(config["dedup_bam_qc_dir"], "{sample}.dedup.flagstat.txt"), sample=samples))
    #outputs.extend(expand(os.path.join(config["dedup_bam_qc_dir"], "{sample}.dedup.alignment_summary_metrics.txt"), sample=samples))
    outputs.extend(expand(os.path.join(config["filtered_bam_qc_dir"], "{sample}.filtered.stats.txt"), sample=samples))
    outputs.extend(expand(os.path.join(config["filtered_bam_qc_dir"], "{sample}.filtered.flagstat.txt"), sample=samples))
    outputs.extend(expand(os.path.join(config["filtered_bam_qc_dir"], "{sample}.filtered.alignment_summary_metrics.txt"), sample=samples))
    outputs.extend(expand(os.path.join(config["mbias_dir"], "{sample}.mbias.txt"), sample=samples))
    outputs.extend(expand(os.path.join(config["mbias_dir"], "{sample}.options.txt"), sample=samples))
    outputs.append(os.path.join(config["mbias_dir"], "all_samples_mbias_options.tsv"))
    outputs.extend(expand(os.path.join(config["methyldackel_dir"], "{sample}_CpG.methylKit"), sample=samples))
    outputs.extend(expand(os.path.join(config["methyldackel_mergecontext_dir"], "{sample}_CpG.bedGraph"), sample=samples))
    return outputs


def get_wgbs_validation_outputs(samples):
    outputs = []
    outputs.extend(expand(os.path.join(FASTQ_PAIR_VALIDATION_DIR, "{sample}.fastq_pair_validation.tsv"), sample=samples))
    outputs.extend(expand(os.path.join(READ_GROUP_QC_DIR, "{sample}_pe.read_group_validation.tsv"), sample=samples))
    return outputs


def get_wgbs_multiqc_inputs(samples):
    inputs = []
    inputs.extend([os.path.join(config["qc_dir"], f"{sample}_{read}_raw_fastqc.zip") for sample in samples for read in ["R1", "R2"]])
    inputs.extend(expand(os.path.join(config["trimmed_dir"], "{sample}.fastp.json"), sample=samples))
    inputs.extend([os.path.join(config["qc_trimmed_dir"], f"{sample}_{read}_trimmed_fastqc.zip") for sample in samples for read in ["R1", "R2"]])
    inputs.extend(expand(os.path.join(config["alignment_dir"], "{sample}_pe_report.txt"), sample=samples))
    inputs.extend(expand(os.path.join(config["dedup_dir"], "{sample}_pe.deduplication_report.txt"), sample=samples))
    inputs.extend(expand(os.path.join(config["dedup_bam_qc_dir"], "{sample}.dedup.stats.txt"), sample=samples))
    inputs.extend(expand(os.path.join(config["dedup_bam_qc_dir"], "{sample}.dedup.flagstat.txt"), sample=samples))
    inputs.extend(expand(os.path.join(config["filtered_bam_qc_dir"], "{sample}.filtered.stats.txt"), sample=samples))
    inputs.extend(expand(os.path.join(config["filtered_bam_qc_dir"], "{sample}.filtered.flagstat.txt"), sample=samples))
    inputs.extend(expand(os.path.join(config["filtered_bam_qc_dir"], "{sample}.filtered.alignment_summary_metrics.txt"), sample=samples))
    return inputs

# --- MultiQC Configuration ---
config["pipeline_name"] = "wgbs"
config["multiqc_results_dir"] = config.get("multiqc_results_dir", os.path.join("results", "wgbs"))
config["multiqc_input_files"] = get_wgbs_multiqc_inputs(SAMPLES)


# --- MODULE INCLUSION ---

# 1. Shared FASTQ and read-group validation
include: "rules/fastq_validation.smk"
include: "rules/read_group_validation.smk"

# 2. QC on raw files
include: "rules/qc.smk"

# 2. Trimming and QC on trimmed files
include: "rules/trimming_and_qc.smk"

# 3. Alignment
include: "rules/wgbs/alignment.smk"

# 4. Generic BAM QC
include: "rules/bam_qc.smk"

# 5. Filtering and Sorting of BAMs
include: "rules/wgbs/filter_bam.smk"

# 6. Methylation Bias
include: "rules/wgbs/mbias.smk"

# 7. Methylation Extraction
include: "rules/wgbs/methyldackel.smk"

# 8. MultiQC report
include: "rules/multiqc.smk"


use rule validate_fastq_pair_generic as validate_fastq_pair with:
    input:
        r1=lambda wildcards: SAMPLES_INFO[wildcards.sample]["R1"],
        r2=lambda wildcards: SAMPLES_INFO[wildcards.sample]["R2"],
        validator="scripts/validate_fastq_pairs.py"
    output:
        report=os.path.join(
            FASTQ_PAIR_VALIDATION_DIR, "{sample}.fastq_pair_validation.tsv"
        )
    params:
        scope=FASTQ_VALIDATION_SCOPE
    log:
        os.path.join(LOG_DIR, "fastq_pair_validation", "{sample}.log")

use rule validate_bam_read_group_generic as validate_bam_read_group with:
    input:
        bam=os.path.join(ALIGN_DIR, "{sample}_pe.bam"),
        validator="src/validate_bam_read_group.py"
    output:
        report=os.path.join(
            READ_GROUP_QC_DIR, "{sample}_pe.read_group_validation.tsv"
        )
    params:
        rg_id=lambda wildcards: SAMPLES_INFO[wildcards.sample]["read_group"]["id"],
        rg_sm=lambda wildcards: SAMPLES_INFO[wildcards.sample]["read_group"]["sample"],
        rg_lb=lambda wildcards: SAMPLES_INFO[wildcards.sample]["read_group"]["library"],
        rg_pl=lambda wildcards: SAMPLES_INFO[wildcards.sample]["read_group"]["platform"],
        rg_pu=lambda wildcards: SAMPLES_INFO[wildcards.sample]["read_group"]["platform_unit"]
    log:
        os.path.join(LOG_DIR, "read_group_validation", "{sample}_pe.log")


# --- BAM QC INSTANTIATION ---

# Deduplicated BAM QC stays limited to samtools because deduplicate_bismark does
# not guarantee coordinate-sorted output for Picard alignment metrics.
use rule samtools_stats_generic as samtools_stats_dedup with:
    input:
        bam = os.path.join(DEDUP_DIR, "{sample}_pe.deduplicated.bam")
    output:
        stats = os.path.join(DEDUP_BAM_QC_DIR, "{sample}.dedup.stats.txt")
    log:
        os.path.join(LOG_DIR, "samtools_stats_dedup", "{sample}.log")

use rule samtools_flagstat_generic as samtools_flagstat_dedup with:
    input:
        bam = os.path.join(DEDUP_DIR, "{sample}_pe.deduplicated.bam")
    output:
        flagstat = os.path.join(DEDUP_BAM_QC_DIR, "{sample}.dedup.flagstat.txt")
    log:
        os.path.join(LOG_DIR, "samtools_flagstat_dedup", "{sample}.log")

#use rule picard_collect_alignment_metrics_generic as picard_collect_alignment_metrics_dedup with:
#    input:
#        bam = os.path.join(DEDUP_DIR, "{sample}_pe.deduplicated.bam"),
#        ref = config["ref_genome"]
#    output:
#        metrics = os.path.join(DEDUP_BAM_QC_DIR, "{sample}.dedup.alignment_summary_metrics.txt")
#    log:
#        os.path.join(LOG_DIR, "picard_collect_alignment_metrics_dedup", "{sample}.log")

# Filtered BAMs
use rule samtools_stats_generic as samtools_stats_filtered with:
    input:
        bam = os.path.join(SORTED_FILTERED_BAM_DIR, "{sample}_pe.filtered.sorted.bam")
    output:
        stats = os.path.join(FILTERED_BAM_QC_DIR, "{sample}.filtered.stats.txt")
    log:
        os.path.join(LOG_DIR, "samtools_stats_filtered", "{sample}.log")

use rule samtools_flagstat_generic as samtools_flagstat_filtered with:
    input:
        bam = os.path.join(SORTED_FILTERED_BAM_DIR, "{sample}_pe.filtered.sorted.bam")
    output:
        flagstat = os.path.join(FILTERED_BAM_QC_DIR, "{sample}.filtered.flagstat.txt")
    log:
        os.path.join(LOG_DIR, "samtools_flagstat_filtered", "{sample}.log")

use rule picard_collect_alignment_metrics_generic as picard_collect_alignment_metrics_filtered with:
    input:
        bam = os.path.join(SORTED_FILTERED_BAM_DIR, "{sample}_pe.filtered.sorted.bam"),
        ref = config["reference_fasta"]
    output:
        metrics = os.path.join(FILTERED_BAM_QC_DIR, "{sample}.filtered.alignment_summary_metrics.txt")
    log:
        os.path.join(LOG_DIR, "picard_collect_alignment_metrics_filtered", "{sample}.log")


# --- FINAL TARGETS ---

rule all:
    input:
        get_wgbs_outputs(SAMPLES),
        get_wgbs_validation_outputs(SAMPLES),
        BISMARK_INDEX_FILES,
        config["reference_fasta_fai"],
        rules.multiqc.output
    default_target: True
