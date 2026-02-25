import os
import re
from collections import defaultdict
import sys

sys.path.insert(0, os.path.abspath("src"))
from utils import prepare_sample_data

# --- CONFIGURATION ---
# Define directories from config
RAW_DIR = config["raw_fastqs_dir"]
QC_DIR = config["qc_dir"]
TRIMMED_DIR = config["trimmed_dir"]
QC_TRIMMED_DIR = config["qc_trimmed_dir"]
ALIGNMENT_DIR = config["alignment_dir"]
BOWTIE2_INDEX_DIR = config["bowtie2_index_dir"]
BAM_QC_DIR = config["bam_qc_dir"]
FILTERED_BAM_DIR = config["filtered_bam_dir"]
FILTERED_BAM_QC_DIR = config["filtered_bam_qc_dir"]
DEEPTOOLS_DIR = config["deeptools_dir"]
BIGWIG_DIR = config["bigwig_dir"]
SUBTRACTED_BIGWIG_DIR = config["subtracted_bigwig_dir"]

# Make config available to included rules
config["raw_fastqs_dir"] = RAW_DIR
config["qc_dir"] = QC_DIR
config["trimmed_dir"] = TRIMMED_DIR
config["qc_trimmed_dir"] = QC_TRIMMED_DIR
config["alignment_dir"] = ALIGNMENT_DIR
config["bowtie2_index_dir"] = BOWTIE2_INDEX_DIR
config["bam_qc_dir"] = BAM_QC_DIR
config["filtered_bam_dir"] = FILTERED_BAM_DIR
config["filtered_bam_qc_dir"] = FILTERED_BAM_QC_DIR
config["deeptools_dir"] = DEEPTOOLS_DIR
config["bigwig_dir"] = BIGWIG_DIR
config["subtracted_bigwig_dir"] = SUBTRACTED_BIGWIG_DIR


# Ensure output directories exist
os.makedirs(QC_DIR, exist_ok=True)
os.makedirs(TRIMMED_DIR, exist_ok=True)
os.makedirs(QC_TRIMMED_DIR, exist_ok=True)
os.makedirs(ALIGNMENT_DIR, exist_ok=True)
os.makedirs(BOWTIE2_INDEX_DIR, exist_ok=True)
os.makedirs(BAM_QC_DIR, exist_ok=True)
os.makedirs(FILTERED_BAM_DIR, exist_ok=True)
os.makedirs(FILTERED_BAM_QC_DIR, exist_ok=True)
os.makedirs(DEEPTOOLS_DIR, exist_ok=True)
os.makedirs(BIGWIG_DIR, exist_ok=True)
os.makedirs(SUBTRACTED_BIGWIG_DIR, exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "fastqc_raw"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "fastp"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "fastqc_trimmed"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "bowtie2_align"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "bam_qc"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "sambamba_filter"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "deeptools"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "bamCoverage"), exist_ok=True)


# --- SAMPLE DISCOVERY ---
SAMPLES_INFO, SAMPLES = prepare_sample_data(config)

def get_subtraction_pairs(samples_info, bigwig_dir):
    pairs = []
    
    chip_samples = [s for s in samples_info if '_input' not in s]
    input_samples = [s for s in samples_info if '_input' in s]

    for chip_sample in chip_samples:
        # Construct the expected input name, e.g., "liver_rep1" -> "liver_input_rep1"
        parts = chip_sample.split('_')
        if len(parts) == 2: # Expects "tissue_repX" format
            base, rep = parts
            expected_input = f"{base}_input_{rep}"
            
            if expected_input in input_samples:
                chip_info = samples_info[chip_sample]
                input_info = samples_info[expected_input]
                
                chip_read_type = 'pe' if chip_info['type'] == 'PE' else 'se'
                # Input and chip might have different read types, but for pairing we assume they match
                # The bw file will have the correct read_type suffix from its own rule
                input_read_type = 'pe' if input_info['type'] == 'PE' else 'se'

                chip_bw = os.path.join(bigwig_dir, f"{chip_sample}_{chip_read_type}.bw")
                input_bw = os.path.join(bigwig_dir, f"{expected_input}_{input_read_type}.bw")
                
                pairs.append({
                    'chip_bw': chip_bw, 
                    'input_bw': input_bw, 
                    'sample': chip_sample, 
                    'read_type': chip_read_type
                })
                
    return pairs

SUBTRACTION_PAIRS = get_subtraction_pairs(SAMPLES_INFO, BIGWIG_DIR)

# Pass SUBTRACTION_PAIRS to config so it's accessible in included rules
config['subtraction_pairs'] = SUBTRACTION_PAIRS


# --- MultiQC Configuration ---
config["pipeline_name"] = "chipseq_cutrun"
config["multiqc_results_dir"] = "results/chipseq_cutrun"

pe_samples = [s for s, i in SAMPLES_INFO.items() if i['type'] == 'PE']
se_samples = [s for s, i in SAMPLES_INFO.items() if i['type'] == 'SE']

final_outputs = []
# Raw QC
final_outputs.extend(expand(os.path.join(QC_DIR, "{sample}_R1_raw_fastqc.html"), sample=pe_samples))
final_outputs.extend(expand(os.path.join(QC_DIR, "{sample}_R2_raw_fastqc.html"), sample=pe_samples))
final_outputs.extend(expand(os.path.join(QC_DIR, "{sample}_raw_fastqc.html"), sample=se_samples))
# Trimmed QC
final_outputs.extend(expand(os.path.join(QC_TRIMMED_DIR, "{sample}_R1_trimmed_fastqc.html"), sample=pe_samples))
final_outputs.extend(expand(os.path.join(QC_TRIMMED_DIR, "{sample}_R2_trimmed_fastqc.html"), sample=pe_samples))
final_outputs.extend(expand(os.path.join(QC_TRIMMED_DIR, "{sample}_SE_trimmed_fastqc.html"), sample=se_samples))
# Aligned BAM QC
final_outputs.extend(expand(os.path.join(BAM_QC_DIR, "{sample}_pe.stats.txt"), sample=pe_samples))
final_outputs.extend(expand(os.path.join(BAM_QC_DIR, "{sample}_pe.flagstat.txt"), sample=pe_samples))
final_outputs.extend(expand(os.path.join(BAM_QC_DIR, "{sample}_pe.alignment_summary_metrics.txt"), sample=pe_samples))
final_outputs.extend(expand(os.path.join(BAM_QC_DIR, "{sample}_se.stats.txt"), sample=se_samples))
final_outputs.extend(expand(os.path.join(BAM_QC_DIR, "{sample}_se.flagstat.txt"), sample=se_samples))
final_outputs.extend(expand(os.path.join(BAM_QC_DIR, "{sample}_se.alignment_summary_metrics.txt"), sample=se_samples))
# Filtered BAMs
final_outputs.extend(expand(os.path.join(FILTERED_BAM_DIR, "{sample}_pe.filtered.sorted.bam"), sample=pe_samples))
final_outputs.extend(expand(os.path.join(FILTERED_BAM_DIR, "{sample}_se.filtered.sorted.bam"), sample=se_samples))
# Filtered BAM QC
final_outputs.extend(expand(os.path.join(FILTERED_BAM_QC_DIR, "{sample}_pe.stats.txt"), sample=pe_samples))
final_outputs.extend(expand(os.path.join(FILTERED_BAM_QC_DIR, "{sample}_pe.flagstat.txt"), sample=pe_samples))
final_outputs.extend(expand(os.path.join(FILTERED_BAM_QC_DIR, "{sample}_pe.alignment_summary_metrics.txt"), sample=pe_samples))
final_outputs.extend(expand(os.path.join(FILTERED_BAM_QC_DIR, "{sample}_se.stats.txt"), sample=se_samples))
final_outputs.extend(expand(os.path.join(FILTERED_BAM_QC_DIR, "{sample}_se.flagstat.txt"), sample=se_samples))
final_outputs.extend(expand(os.path.join(FILTERED_BAM_QC_DIR, "{sample}_se.alignment_summary_metrics.txt"), sample=se_samples))
# Deeptools
final_outputs.extend([
    os.path.join(DEEPTOOLS_DIR, "bam_correlation_heatmap.png"),
    os.path.join(DEEPTOOLS_DIR, "fingerprints.png"),
    os.path.join(DEEPTOOLS_DIR, "heatmap.png"),
])
# Bigwigs
final_outputs.extend(expand(os.path.join(BIGWIG_DIR, "{sample}_pe.bw"), sample=pe_samples))
final_outputs.extend(expand(os.path.join(BIGWIG_DIR, "{sample}_se.bw"), sample=se_samples))
# Subtracted Bigwigs
final_outputs.extend([os.path.join(SUBTRACTED_BIGWIG_DIR, f"{pair['sample']}_{pair['read_type']}.subtracted.bw") for pair in config.get('subtraction_pairs', [])])

config["multiqc_input_files"] = final_outputs


# --- MODULE INCLUSION ---

# 1. QC on raw files
include: "rules/qc.smk"

# 2. Trimming and QC on trimmed files
include: "rules/trimming_and_qc.smk"

# 3. Alignment
include: "rules/chip_cr/alignment.smk"

# 4. Generic BAM QC
include: "rules/bam_qc.smk"

# 5. Filtering and Deduplication
include: "rules/chip_cr/filter_bam.smk"

# 6. BigWig generation and subtraction
include: "rules/chip_cr/bigwig.smk"

# 7. Heatmap computeMatrix + plotHeatmap
include: "rules/chip_cr/heatmap.smk"

# 8. MultiQC report
include: "rules/multiqc.smk"


# --- BAM QC INSTANTIATION ---

# Aligned BAMs
use rule samtools_stats_generic as samtools_stats_aligned with:
    input:
        bam = os.path.join(ALIGNMENT_DIR, "{sample}_{read_type}.sorted.bam")
    output:
        stats = os.path.join(BAM_QC_DIR, "{sample}_{read_type}.stats.txt")
    log:
        os.path.join("logs", config["pipeline"], "bam_qc", "{sample}_{read_type}_stats.log")

use rule samtools_flagstat_generic as samtools_flagstat_aligned with:
    input:
        bam = os.path.join(ALIGNMENT_DIR, "{sample}_{read_type}.sorted.bam")
    output:
        flagstat = os.path.join(BAM_QC_DIR, "{sample}_{read_type}.flagstat.txt")
    log:
        os.path.join("logs", config["pipeline"], "bam_qc", "{sample}_{read_type}_flagstat.log")

use rule picard_collect_alignment_metrics_generic as picard_collect_alignment_metrics_aligned with:
    input:
        bam = os.path.join(ALIGNMENT_DIR, "{sample}_{read_type}.sorted.bam"),
        ref = config["ref_genome"]
    output:
        metrics = os.path.join(BAM_QC_DIR, "{sample}_{read_type}.alignment_summary_metrics.txt")
    log:
        os.path.join("logs", config["pipeline"], "bam_qc", "{sample}_{read_type}_picard_metrics.log")

# Filtered BAMs
use rule samtools_stats_generic as samtools_stats_filtered with:
    input:
        bam = os.path.join(FILTERED_BAM_DIR, "{sample}_{read_type}.filtered.sorted.bam")
    output:
        stats = os.path.join(FILTERED_BAM_QC_DIR, "{sample}_{read_type}.stats.txt")
    log:
        os.path.join("logs", config["pipeline"], "bam_qc", "{sample}_{read_type}_filtered_stats.log")

use rule samtools_flagstat_generic as samtools_flagstat_filtered with:
    input:
        bam = os.path.join(FILTERED_BAM_DIR, "{sample}_{read_type}.filtered.sorted.bam")
    output:
        flagstat = os.path.join(FILTERED_BAM_QC_DIR, "{sample}_{read_type}.flagstat.txt")
    log:
        os.path.join("logs", config["pipeline"], "bam_qc", "{sample}_{read_type}_filtered_flagstat.log")

use rule picard_collect_alignment_metrics_generic as picard_collect_alignment_metrics_filtered with:
    input:
        bam = os.path.join(FILTERED_BAM_DIR, "{sample}_{read_type}.filtered.sorted.bam"),
        ref = config["ref_genome"]
    output:
        metrics = os.path.join(FILTERED_BAM_QC_DIR, "{sample}_{read_type}.alignment_summary_metrics.txt")
    log:
        os.path.join("logs", config["pipeline"], "bam_qc", "{sample}_{read_type}_filtered_picard_metrics.log")


# --- DEEPTOOLS RULES (from former filtered_bam_qc.smk) ---

def get_all_filtered_bams(samples_info):
    bams = []
    for sample, info in samples_info.items():
        read_type = 'pe' if info['type'] == 'PE' else 'se'
        bams.append(os.path.join(FILTERED_BAM_DIR, f"{sample}_{read_type}.filtered.sorted.bam"))
    return bams

ALL_FILTERED_BAMS = get_all_filtered_bams(SAMPLES_INFO)

MBS_ARGS = config.get("deeptools", {}).get("multiBamSummary", {}).get("extra_args", "")
PC_ARGS = config.get("deeptools", {}).get("plotCorrelation", {}).get("extra_args", "")
PF_ARGS = config.get("deeptools", {}).get("plotFingerprint", {}).get("extra_args", "")

rule multiBamSummary:
    """
    Computes read coverages for multiple BAM files.
    """
    input:
        bams = ALL_FILTERED_BAMS
    output:
        npz = os.path.join(DEEPTOOLS_DIR, "read_coverage.npz")
    params:
        extra = MBS_ARGS
    threads: 8
    log:
        os.path.join("logs", config["pipeline"], "deeptools", "multiBamSummary.log")
    shell:
        "pixi run multiBamSummary bins -b {input.bams} -o {output.npz} -p {threads} {params.extra} > {log}.out 2> {log}.err"

rule plotCorrelation:
    """
    Creates a heatmap of correlations based on multiBamSummary output.
    """
    input:
        npz = os.path.join(DEEPTOOLS_DIR, "read_coverage.npz")
    output:
        heatmap = os.path.join(DEEPTOOLS_DIR, "bam_correlation_heatmap.png"),
        matrix = os.path.join(DEEPTOOLS_DIR, "bam_correlation_matrix.tab")
    params:
        extra = PC_ARGS
    threads: 1
    log:
        os.path.join("logs", config["pipeline"], "deeptools", "plotCorrelation.log")
    shell:
        "pixi run plotCorrelation -in {input.npz} -o {output.heatmap} --outFileCorMatrix {output.matrix} {params.extra} > {log}.out 2> {log}.err"

rule plotFingerprint:
    """
    Generates fingerprints for each BAM file to assess ChIP-seq quality.
    """
    input:
        bams = ALL_FILTERED_BAMS
    output:
        plot = os.path.join(DEEPTOOLS_DIR, "fingerprints.png"),
        metrics = os.path.join(DEEPTOOLS_DIR, "fingerprints.metrics.tab")
    params:
        extra = PF_ARGS
    threads: 8
    log:
        os.path.join("logs", config["pipeline"], "deeptools", "plotFingerprint.log")
    shell:
        "pixi run plotFingerprint -b {input.bams} -o {output.plot} --outRawCounts {output.metrics} -p {threads} {params.extra} > {log}.out 2> {log}.err"


# --- FINAL TARGETS ---
#print("Final outputs:", final_outputs)
rule all:
    input:
        final_outputs,
        os.path.join(config["multiqc_results_dir"], "multiqc_report.html")
    default_target: True
