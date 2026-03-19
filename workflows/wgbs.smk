# WGBS module workflow
import os
import sys

sys.path.insert(0, os.path.abspath("src"))
from utils import prepare_sample_data, validate_existing_file

# --- CONFIGURATION ---
# Define directories from config
RAW_DIR = config["raw_fastqs_dir"]
QC_DIR = config["qc_dir"]
TRIMMED_DIR = config["trimmed_dir"]
QC_TRIMMED_DIR = config["qc_trimmed_dir"]
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

# Reference inputs are user-managed. Bismark will build its index next to ref_genome.
validate_existing_file(REF_GENOME, "ref_genome")

# Make config available to included rules
config["raw_fastqs_dir"] = RAW_DIR
config["qc_dir"] = QC_DIR
config["trimmed_dir"] = TRIMMED_DIR
config["qc_trimmed_dir"] = QC_TRIMMED_DIR
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


# Ensure output directories exist
os.makedirs(config["qc_dir"], exist_ok=True)
os.makedirs(config["trimmed_dir"], exist_ok=True)
os.makedirs(config["qc_trimmed_dir"], exist_ok=True)
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
os.makedirs(os.path.join("logs", config["pipeline"], "fastqc_raw"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "fastp"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "fastqc_trimmed"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "multiqc"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "bismark_genome_preparation"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "bismark"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "deduplicate_bismark"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "sambamba_filter"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "sambamba_sort"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "methyldackel_mbias"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "methyldackel_extract_methylkit"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "methyldackel_extract_mergecontext"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "samtools_stats_dedup"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "samtools_flagstat_dedup"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "picard_collect_alignment_metrics_dedup"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "samtools_stats_filtered"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "samtools_flagstat_filtered"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "picard_collect_alignment_metrics_filtered"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "samtools_index_filtered"), exist_ok=True)


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

# --- MultiQC Configuration ---
config["pipeline_name"] = "wgbs"
config["multiqc_results_dir"] = config.get("multiqc_results_dir")
if not config["multiqc_results_dir"]:
    raise ValueError("WGBS requires 'multiqc_results_dir' in the config.")
config["multiqc_analysis_dirs"] = list(
    dict.fromkeys(
        [
            config["multiqc_results_dir"],
            os.path.join("logs", config["pipeline"]),
            ALIGN_DIR,
            DEDUP_DIR,
        ]
    )
)
config["multiqc_input_files"] = get_wgbs_outputs(SAMPLES)


# --- MODULE INCLUSION ---

# 1. QC on raw files
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


# --- BAM QC INSTANTIATION ---

# Deduplicated BAM QC stays limited to samtools because deduplicate_bismark does
# not guarantee coordinate-sorted output for Picard alignment metrics.
use rule samtools_stats_generic as samtools_stats_dedup with:
    input:
        bam = os.path.join(DEDUP_DIR, "{sample}_pe.deduplicated.bam")
    output:
        stats = os.path.join(DEDUP_BAM_QC_DIR, "{sample}.dedup.stats.txt")
    log:
        os.path.join("logs", config["pipeline"], "samtools_stats_dedup", "{sample}.log")

use rule samtools_flagstat_generic as samtools_flagstat_dedup with:
    input:
        bam = os.path.join(DEDUP_DIR, "{sample}_pe.deduplicated.bam")
    output:
        flagstat = os.path.join(DEDUP_BAM_QC_DIR, "{sample}.dedup.flagstat.txt")
    log:
        os.path.join("logs", config["pipeline"], "samtools_flagstat_dedup", "{sample}.log")

#use rule picard_collect_alignment_metrics_generic as picard_collect_alignment_metrics_dedup with:
#    input:
#        bam = os.path.join(DEDUP_DIR, "{sample}_pe.deduplicated.bam"),
#        ref = config["ref_genome"]
#    output:
#        metrics = os.path.join(DEDUP_BAM_QC_DIR, "{sample}.dedup.alignment_summary_metrics.txt")
#    log:
#        os.path.join("logs", config["pipeline"], "picard_collect_alignment_metrics_dedup", "{sample}.log")

# Filtered BAMs
use rule samtools_stats_generic as samtools_stats_filtered with:
    input:
        bam = os.path.join(SORTED_FILTERED_BAM_DIR, "{sample}_pe.filtered.sorted.bam")
    output:
        stats = os.path.join(FILTERED_BAM_QC_DIR, "{sample}.filtered.stats.txt")
    log:
        os.path.join("logs", config["pipeline"], "samtools_stats_filtered", "{sample}.log")

use rule samtools_flagstat_generic as samtools_flagstat_filtered with:
    input:
        bam = os.path.join(SORTED_FILTERED_BAM_DIR, "{sample}_pe.filtered.sorted.bam")
    output:
        flagstat = os.path.join(FILTERED_BAM_QC_DIR, "{sample}.filtered.flagstat.txt")
    log:
        os.path.join("logs", config["pipeline"], "samtools_flagstat_filtered", "{sample}.log")

use rule picard_collect_alignment_metrics_generic as picard_collect_alignment_metrics_filtered with:
    input:
        bam = os.path.join(SORTED_FILTERED_BAM_DIR, "{sample}_pe.filtered.sorted.bam"),
        ref = config["ref_genome"]
    output:
        metrics = os.path.join(FILTERED_BAM_QC_DIR, "{sample}.filtered.alignment_summary_metrics.txt")
    log:
        os.path.join("logs", config["pipeline"], "picard_collect_alignment_metrics_filtered", "{sample}.log")


# --- FINAL TARGETS ---

rule all:
    input:
        get_wgbs_outputs(SAMPLES),
        os.path.join(config["multiqc_results_dir"], "multiqc_report.html")
    default_target: True
