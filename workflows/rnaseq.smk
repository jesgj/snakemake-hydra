# workflows/rnaseq.smk
import os
import re
from collections import defaultdict
import sys

sys.path.insert(0, os.path.abspath("src"))
from utils import prepare_sample_data

# --- CONFIGURATION ---
RAW_DIR = config["raw_fastqs_dir"]
QC_DIR = config["qc_dir"]
TRIMMED_DIR = config["trimmed_dir"]
QC_TRIMMED_DIR = config["qc_trimmed_dir"]
# Kallisto
KALLISTO_OUTPUT_DIR = config["kallisto_output_dir"]
TRANSCRIPTOME_FASTA = config["transcriptome_fasta"]
KALLISTO_INDEX = config["kallisto_index"]
# HISAT2
REF_GENOME = config["ref_genome"]
HISAT2_INDEX_DIR = config["hisat2_index_dir"]
ALIGNMENT_DIR = config["alignment_dir"]
MARKED_BAM_DIR = config["marked_bam_dir"]
DUPLICATION_QC_DIR = config["duplication_qc_dir"]
BIGWIG_DIR = config["bigwig_dir"]
DEEPTOOLS_DIR = config.get("deeptools_dir", os.path.join("results", "rnaseq", "deeptools"))
GENE_BODY_COVERAGE_CONFIG = config.get("gene_body_coverage", {})
GENE_BODY_COVERAGE_ENABLED = GENE_BODY_COVERAGE_CONFIG.get("enabled", True)
GENE_BODY_COVERAGE_DIR = config.get(
    "gene_body_coverage_dir", os.path.join("results", "rnaseq", "gene_body_coverage")
)
GENE_BODY_COVERAGE_FORMAT = GENE_BODY_COVERAGE_CONFIG.get("format", "pdf")


# Make config available to included rules
config["raw_fastqs_dir"] = RAW_DIR
config["qc_dir"] = QC_DIR
config["trimmed_dir"] = TRIMMED_DIR
config["qc_trimmed_dir"] = QC_TRIMMED_DIR
config["kallisto_output_dir"] = KALLISTO_OUTPUT_DIR
config["transcriptome_fasta"] = TRANSCRIPTOME_FASTA
config["kallisto_index"] = KALLISTO_INDEX
config["ref_genome"] = REF_GENOME
config["hisat2_index_dir"] = HISAT2_INDEX_DIR
config["alignment_dir"] = ALIGNMENT_DIR
config["marked_bam_dir"] = MARKED_BAM_DIR
config["duplication_qc_dir"] = DUPLICATION_QC_DIR
config["bigwig_dir"] = BIGWIG_DIR
config["deeptools_dir"] = DEEPTOOLS_DIR
config["gene_body_coverage_dir"] = GENE_BODY_COVERAGE_DIR
config["gene_body_coverage"] = GENE_BODY_COVERAGE_CONFIG


# Ensure output directories exist
os.makedirs(QC_DIR, exist_ok=True)
os.makedirs(TRIMMED_DIR, exist_ok=True)
os.makedirs(QC_TRIMMED_DIR, exist_ok=True)
os.makedirs(KALLISTO_OUTPUT_DIR, exist_ok=True)
os.makedirs(HISAT2_INDEX_DIR, exist_ok=True)
os.makedirs(ALIGNMENT_DIR, exist_ok=True)
os.makedirs(MARKED_BAM_DIR, exist_ok=True)
os.makedirs(DUPLICATION_QC_DIR, exist_ok=True)
os.makedirs(BIGWIG_DIR, exist_ok=True)
os.makedirs(DEEPTOOLS_DIR, exist_ok=True)
if GENE_BODY_COVERAGE_ENABLED:
    os.makedirs(GENE_BODY_COVERAGE_DIR, exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "fastqc_raw"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "fastp"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "fastqc_trimmed"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "kallisto_quant"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "hisat2_align"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "hisat2_build"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "kallisto_index"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "samtools_index"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "markduplicates"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "deeptools"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "bamCoverage"), exist_ok=True)
if GENE_BODY_COVERAGE_ENABLED:
    os.makedirs(os.path.join("logs", config["pipeline"], "gene_body_coverage"), exist_ok=True)


# --- SAMPLE DISCOVERY ---
SAMPLES_INFO, SAMPLES = prepare_sample_data(config)

if not SAMPLES:
    raise ValueError("No RNA-seq samples were found. Provide 'samples_info' or valid FASTQs in 'raw_fastqs_dir'.")

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
        "RNA-seq currently supports paired-end samples only (type=PE with R1 and R2). "
        f"Invalid samples: {', '.join(invalid_samples)}"
    )


def get_all_aligned_bams(samples):
    return [os.path.join(ALIGNMENT_DIR, f"{sample}_pe.sorted.bam") for sample in samples]


def get_all_aligned_bais(samples):
    return [f"{bam}.bai" for bam in get_all_aligned_bams(samples)]


ALL_ALIGNED_BAMS = get_all_aligned_bams(SAMPLES)
ALL_ALIGNED_BAIS = get_all_aligned_bais(SAMPLES)
HAS_MULTI_BAM_DEEPTOOLS = len(ALL_ALIGNED_BAMS) >= 2


# --- HELPER FUNCTION FOR OUTPUTS ---
def get_markduplicates_outputs(samples):
    outputs = []
    outputs.extend(expand(os.path.join(MARKED_BAM_DIR, "{sample}_pe.markdup.bam"), sample=samples))
    outputs.extend(expand(os.path.join(MARKED_BAM_DIR, "{sample}_pe.markdup.bam.bai"), sample=samples))
    outputs.extend(expand(os.path.join(DUPLICATION_QC_DIR, "{sample}_pe.markdup.metrics.txt"), sample=samples))
    return outputs


def get_gene_body_coverage_outputs(samples):
    if not GENE_BODY_COVERAGE_ENABLED:
        return []

    output_prefix = os.path.join(GENE_BODY_COVERAGE_DIR, "all_samples")
    outputs = [
        f"{output_prefix}.geneBodyCoverage.txt",
        f"{output_prefix}.geneBodyCoverage.curves.{GENE_BODY_COVERAGE_FORMAT}",
    ]
    if len(samples) >= 3:
        outputs.append(f"{output_prefix}.geneBodyCoverage.heatMap.{GENE_BODY_COVERAGE_FORMAT}")
    return outputs


def get_rnaseq_outputs(samples):
    outputs = []
    outputs.extend(expand(os.path.join(QC_DIR, "{sample}_R1_raw_fastqc.html"), sample=samples))
    outputs.extend(expand(os.path.join(QC_DIR, "{sample}_R2_raw_fastqc.html"), sample=samples))
    outputs.extend(expand(os.path.join(QC_TRIMMED_DIR, "{sample}_R1_trimmed_fastqc.html"), sample=samples))
    outputs.extend(expand(os.path.join(QC_TRIMMED_DIR, "{sample}_R2_trimmed_fastqc.html"), sample=samples))
    outputs.extend(expand(os.path.join(KALLISTO_OUTPUT_DIR, "{sample}", "abundance.tsv"), sample=samples))
    outputs.extend(expand(os.path.join(ALIGNMENT_DIR, "{sample}_pe.sorted.bam"), sample=samples))
    outputs.extend(expand(os.path.join(ALIGNMENT_DIR, "{sample}_pe.hisat2.summary.txt"), sample=samples))
    outputs.extend(get_markduplicates_outputs(samples))
    outputs.extend(get_gene_body_coverage_outputs(samples))
    if len(samples) >= 2:
        outputs.append(os.path.join(DEEPTOOLS_DIR, "bam_correlation_heatmap.png"))
        outputs.append(os.path.join(DEEPTOOLS_DIR, "bam_correlation_matrix.tab"))
    outputs.extend(expand(os.path.join(BIGWIG_DIR, "{sample}_pe.bw"), sample=samples))
    return outputs


def get_rnaseq_multiqc_inputs(samples):
    inputs = []
    inputs.extend(expand(os.path.join(QC_DIR, "{sample}_R1_raw_fastqc.zip"), sample=samples))
    inputs.extend(expand(os.path.join(QC_DIR, "{sample}_R2_raw_fastqc.zip"), sample=samples))
    inputs.extend(expand(os.path.join(TRIMMED_DIR, "{sample}.fastp.html"), sample=samples))
    inputs.extend(expand(os.path.join(TRIMMED_DIR, "{sample}.fastp.json"), sample=samples))
    inputs.extend(expand(os.path.join(QC_TRIMMED_DIR, "{sample}_R1_trimmed_fastqc.zip"), sample=samples))
    inputs.extend(expand(os.path.join(QC_TRIMMED_DIR, "{sample}_R2_trimmed_fastqc.zip"), sample=samples))
    inputs.extend(expand(os.path.join(KALLISTO_OUTPUT_DIR, "{sample}", "abundance.tsv"), sample=samples))
    inputs.extend(expand(os.path.join(ALIGNMENT_DIR, "{sample}_pe.hisat2.summary.txt"), sample=samples))
    inputs.extend(expand(os.path.join(DUPLICATION_QC_DIR, "{sample}_pe.markdup.metrics.txt"), sample=samples))
    inputs.extend(get_gene_body_coverage_outputs(samples))
    return inputs


def get_rnaseq_multiqc_analysis_dirs():
    analysis_dirs = [
        QC_DIR,
        TRIMMED_DIR,
        QC_TRIMMED_DIR,
        KALLISTO_OUTPUT_DIR,
        ALIGNMENT_DIR,
        MARKED_BAM_DIR,
        DUPLICATION_QC_DIR,
        BIGWIG_DIR,
        DEEPTOOLS_DIR,
        os.path.join("logs", config["pipeline"]),
    ]
    if GENE_BODY_COVERAGE_ENABLED:
        analysis_dirs.append(GENE_BODY_COVERAGE_DIR)
    return list(dict.fromkeys(analysis_dirs))

# --- MultiQC Configuration ---
config["pipeline_name"] = "rnaseq"
config["multiqc_results_dir"] = "results/rnaseq"
config["multiqc_input_files"] = get_rnaseq_multiqc_inputs(SAMPLES)
config["multiqc_analysis_dirs"] = get_rnaseq_multiqc_analysis_dirs()


# --- MODULE INCLUSION ---

# 1. QC on raw files (using the generic rule)
include: "rules/qc.smk"

# 2. Trimming and QC on trimmed files (using the generic rule)
include: "rules/trimming_and_qc.smk"

# 3. Pseudo-alignment
include: "rules/rnaseq/pseudoalignment.smk"

# 4. Alignment
include: "rules/rnaseq/alignment.smk"

# 5. Duplicate marking QC
include: "rules/rnaseq/markduplicates.smk"

# 6. Shared deepTools QC
include: "rules/deeptools_qc.smk"

# 7. Optional gene body coverage QC
if GENE_BODY_COVERAGE_ENABLED:
    include: "rules/rnaseq/gene_body_coverage.smk"

# 8. BigWig generation
include: "rules/rnaseq/bigwig.smk"

# 9. MultiQC report
include: "rules/multiqc.smk"


MBS_ARGS = config.get("deeptools", {}).get("multiBamSummary", {}).get("extra_args", "--binSize 10000")
PC_ARGS = config.get("deeptools", {}).get("plotCorrelation", {}).get("extra_args", "-p heatmap --corMethod spearman --skipZeros")

if HAS_MULTI_BAM_DEEPTOOLS:
    use rule multiBamSummary_generic as multiBamSummary with:
        input:
            bams = ALL_ALIGNED_BAMS,
            bais = ALL_ALIGNED_BAIS
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
        get_rnaseq_outputs(SAMPLES),
        # MultiQC report
        os.path.join(config["multiqc_results_dir"], "multiqc_report.html")
