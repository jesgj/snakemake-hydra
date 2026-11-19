import os
import re
from collections import defaultdict
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
    os.path.join("results", "chipseq_cutrun", "fastq_pair_validation"),
)
READ_GROUP_QC_DIR = config.get(
    "read_group_qc_dir", os.path.join("results", "chipseq_cutrun", "read_group_qc")
)
FASTQ_VALIDATION_SCOPE = fastq_validation_scope(
    config, "chip_cr.fastq_validation"
)
REF_GENOME = config.get("ref_genome")
GENE_BED = config.get("gene_bed")
ALIGNMENT_DIR = config["alignment_dir"]
BOWTIE2_INDEX_DIR = config.get("bowtie2_index_dir")
BAM_QC_DIR = config["bam_qc_dir"]
FILTERED_BAM_DIR = config["filtered_bam_dir"]
FILTERED_BAM_QC_DIR = config["filtered_bam_qc_dir"]
DEEPTOOLS_DIR = config["deeptools_dir"]
BIGWIG_DIR = config["bigwig_dir"]
SUBTRACTED_BIGWIG_DIR = config["subtracted_bigwig_dir"]

# Reference inputs are user-managed. Only the configured Bowtie2 index directory is created here.
validate_required_path(REF_GENOME, "ref_genome")
validate_required_path(GENE_BED, "gene_bed")
validate_required_path(BOWTIE2_INDEX_DIR, "bowtie2_index_dir")

# Make config available to included rules
config["raw_fastqs_dir"] = RAW_DIR
config["qc_dir"] = QC_DIR
config["trimmed_dir"] = TRIMMED_DIR
config["qc_trimmed_dir"] = QC_TRIMMED_DIR
config["fastq_pair_validation_dir"] = FASTQ_PAIR_VALIDATION_DIR
config["read_group_qc_dir"] = READ_GROUP_QC_DIR
config["fastq_validation_scope"] = FASTQ_VALIDATION_SCOPE
config["ref_genome"] = REF_GENOME
config["gene_bed"] = GENE_BED
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
os.makedirs(FASTQ_PAIR_VALIDATION_DIR, exist_ok=True)
os.makedirs(READ_GROUP_QC_DIR, exist_ok=True)
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
os.makedirs(os.path.join("logs", config["pipeline"], "fastq_pair_validation"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "read_group_validation"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "bowtie2_align"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "bam_qc"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "samtools_index_filtered"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "sambamba_filter"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "deeptools"), exist_ok=True)
os.makedirs(os.path.join("logs", config["pipeline"], "bamCoverage"), exist_ok=True)


# --- SAMPLE DISCOVERY ---
SAMPLES_INFO, SAMPLES = prepare_sample_data(config)

if not SAMPLES:
    raise ValueError("No ChIP-seq/CUT&RUN samples were found. Provide 'samples_info' or valid FASTQs in 'raw_fastqs_dir'.")


def normalize_chip_cr_samples(samples_info):
    normalized_samples = {}
    invalid_samples = []

    for sample, info in samples_info.items():
        normalized_info = dict(info)
        has_r1 = "R1" in normalized_info
        has_r2 = "R2" in normalized_info
        declared_type = normalized_info.get("type")

        if has_r1 and has_r2:
            inferred_type = "PE"
        elif has_r1:
            inferred_type = "SE"
        else:
            invalid_samples.append(f"{sample} (missing R1)")
            continue

        if declared_type is None:
            normalized_info["type"] = inferred_type
        elif declared_type not in ["PE", "SE"]:
            invalid_samples.append(f"{sample} (type must be PE or SE)")
            continue
        elif declared_type != inferred_type:
            invalid_samples.append(
                f"{sample} (declared type {declared_type} does not match provided reads)"
            )
            continue

        normalized_samples[sample] = normalized_info

    if invalid_samples:
        raise ValueError(
            "Invalid ChIP-seq/CUT&RUN samples. Define R1 for all samples, add R2 only for paired-end data, "
            "and use type 'PE' or 'SE' when declared. Invalid samples: " + ", ".join(invalid_samples)
        )

    return normalized_samples


SAMPLES_INFO = normalize_chip_cr_samples(SAMPLES_INFO)
SAMPLES_INFO = normalize_read_group_metadata(SAMPLES_INFO, "chip_cr.samples_info")
SAMPLES = list(SAMPLES_INFO.keys())
config["samples_info"] = SAMPLES_INFO


INPUT_SAMPLE_PATTERN = re.compile(r"^(?P<base>.+)_input_rep(?P<rep>\d+)$")
CHIP_SAMPLE_PATTERN = re.compile(r"^(?P<base>.+)_rep(?P<rep>\d+)$")

def get_subtraction_pairs(samples_info, bigwig_dir):
    pairs = []

    input_samples = {
        sample: info for sample, info in samples_info.items() if INPUT_SAMPLE_PATTERN.match(sample)
    }

    for chip_sample, chip_info in samples_info.items():
        chip_match = CHIP_SAMPLE_PATTERN.match(chip_sample)
        if not chip_match or INPUT_SAMPLE_PATTERN.match(chip_sample):
            continue

        expected_input = f"{chip_match.group('base')}_input_rep{chip_match.group('rep')}"
        if expected_input in input_samples:
            input_info = input_samples[expected_input]

            chip_read_type = 'pe' if chip_info['type'] == 'PE' else 'se'
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
config["multiqc_results_dir"] = config.get("multiqc_results_dir", os.path.join("results", "chipseq_cutrun"))

pe_samples = [s for s, i in SAMPLES_INFO.items() if i['type'] == 'PE']
se_samples = [s for s, i in SAMPLES_INFO.items() if i['type'] == 'SE']
PE_SAMPLE_PATTERN = "|".join(re.escape(sample) for sample in pe_samples) or r"(?!)"


def get_all_filtered_bams(samples_info):
    bams = []
    for sample, info in samples_info.items():
        read_type = 'pe' if info['type'] == 'PE' else 'se'
        bams.append(os.path.join(FILTERED_BAM_DIR, f"{sample}_{read_type}.filtered.sorted.bam"))
    return bams


def get_all_filtered_bais(samples_info):
    return [f"{bam}.bai" for bam in get_all_filtered_bams(samples_info)]


ALL_FILTERED_BAMS = get_all_filtered_bams(SAMPLES_INFO)
ALL_FILTERED_BAIS = get_all_filtered_bais(SAMPLES_INFO)
HAS_MULTI_BAM_DEEPTOOLS = len(ALL_FILTERED_BAMS) >= 2


def get_chip_cr_multiqc_inputs(samples_info):
    inputs = []
    for sample, info in samples_info.items():
        read_type = "pe" if info["type"] == "PE" else "se"
        if read_type == "pe":
            inputs.extend([
                os.path.join(QC_DIR, f"{sample}_R1_raw_fastqc.zip"),
                os.path.join(QC_DIR, f"{sample}_R2_raw_fastqc.zip"),
                os.path.join(TRIMMED_DIR, f"{sample}.fastp.json"),
                os.path.join(QC_TRIMMED_DIR, f"{sample}_R1_trimmed_fastqc.zip"),
                os.path.join(QC_TRIMMED_DIR, f"{sample}_R2_trimmed_fastqc.zip"),
            ])
        else:
            inputs.extend([
                os.path.join(QC_DIR, f"{sample}_raw_fastqc.zip"),
                os.path.join(TRIMMED_DIR, f"{sample}_se.fastp.json"),
                os.path.join(QC_TRIMMED_DIR, f"{sample}_SE_trimmed_fastqc.zip"),
            ])
        inputs.extend([
            os.path.join("logs", config["pipeline"], "bowtie2_align", f"{sample}_{read_type}.log.err"),
            os.path.join("logs", config["pipeline"], "sambamba_markdup", f"{sample}_{read_type}.log.err"),
            os.path.join(BAM_QC_DIR, f"{sample}_{read_type}.stats.txt"),
            os.path.join(BAM_QC_DIR, f"{sample}_{read_type}.flagstat.txt"),
            os.path.join(BAM_QC_DIR, f"{sample}_{read_type}.alignment_summary_metrics.txt"),
            os.path.join(FILTERED_BAM_QC_DIR, f"{sample}_{read_type}.stats.txt"),
            os.path.join(FILTERED_BAM_QC_DIR, f"{sample}_{read_type}.flagstat.txt"),
            os.path.join(FILTERED_BAM_QC_DIR, f"{sample}_{read_type}.alignment_summary_metrics.txt"),
        ])

    inputs.append(os.path.join(DEEPTOOLS_DIR, "fingerprints.metrics.tab"))
    if HAS_MULTI_BAM_DEEPTOOLS:
        inputs.append(os.path.join(DEEPTOOLS_DIR, "bam_correlation_matrix.tab"))
    return inputs

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
    os.path.join(DEEPTOOLS_DIR, "fingerprints.png"),
    os.path.join(DEEPTOOLS_DIR, "fingerprints.metrics.tab"),
    os.path.join(DEEPTOOLS_DIR, "heatmap.png"),
])
if HAS_MULTI_BAM_DEEPTOOLS:
    final_outputs.extend([
        os.path.join(DEEPTOOLS_DIR, "bam_correlation_heatmap.png"),
        os.path.join(DEEPTOOLS_DIR, "bam_correlation_matrix.tab"),
    ])
# Bigwigs
final_outputs.extend(expand(os.path.join(BIGWIG_DIR, "{sample}_pe.bw"), sample=pe_samples))
final_outputs.extend(expand(os.path.join(BIGWIG_DIR, "{sample}_se.bw"), sample=se_samples))
# Subtracted Bigwigs
final_outputs.extend([os.path.join(SUBTRACTED_BIGWIG_DIR, f"{pair['sample']}_{pair['read_type']}.subtracted.bw") for pair in config.get('subtraction_pairs', [])])

validation_outputs = []
validation_outputs.extend(expand(os.path.join(FASTQ_PAIR_VALIDATION_DIR, "{sample}.fastq_pair_validation.tsv"), sample=pe_samples))
validation_outputs.extend(expand(os.path.join(FASTQ_PAIR_VALIDATION_DIR, "{sample}.fastq_single_validation.tsv"), sample=se_samples))
validation_outputs.extend(expand(os.path.join(READ_GROUP_QC_DIR, "{sample}_pe.read_group_validation.tsv"), sample=pe_samples))
validation_outputs.extend(expand(os.path.join(READ_GROUP_QC_DIR, "{sample}_se.read_group_validation.tsv"), sample=se_samples))

config["multiqc_input_files"] = get_chip_cr_multiqc_inputs(SAMPLES_INFO)


# --- MODULE INCLUSION ---

# 1. Shared FASTQ and read-group validation
include: "rules/fastq_validation.smk"
include: "rules/read_group_validation.smk"

# 2. QC on raw files
include: "rules/qc.smk"

# 2. Trimming and QC on trimmed files
include: "rules/trimming_and_qc.smk"

# 3. Alignment
include: "rules/chip_cr/alignment.smk"

# 4. Generic BAM QC
include: "rules/bam_qc.smk"

# 5. Shared deepTools QC
include: "rules/deeptools_qc.smk"

# 6. Filtering and Deduplication
include: "rules/chip_cr/filter_bam.smk"

# 7. BigWig generation and subtraction
include: "rules/chip_cr/bigwig.smk"

# 8. Heatmap computeMatrix + plotHeatmap
include: "rules/chip_cr/heatmap.smk"

# 9. MultiQC report
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
    wildcard_constraints:
        sample=PE_SAMPLE_PATTERN
    params:
        scope=FASTQ_VALIDATION_SCOPE
    log:
        os.path.join(
            "logs", config["pipeline"], "fastq_pair_validation", "{sample}.log"
        )

use rule validate_fastq_single_generic as validate_fastq_single with:
    input:
        r1=lambda wildcards: SAMPLES_INFO[wildcards.sample]["R1"],
        validator="scripts/validate_fastq_pairs.py"
    output:
        report=os.path.join(
            FASTQ_PAIR_VALIDATION_DIR, "{sample}.fastq_single_validation.tsv"
        )
    wildcard_constraints:
        sample=SE_SAMPLE_PATTERN
    params:
        scope=FASTQ_VALIDATION_SCOPE
    log:
        os.path.join(
            "logs", config["pipeline"], "fastq_pair_validation", "{sample}.log"
        )

use rule validate_bam_read_group_generic as validate_bam_read_group with:
    input:
        bam=os.path.join(ALIGNMENT_DIR, "{sample}_{read_type}.sorted.bam"),
        validator="src/validate_bam_read_group.py"
    output:
        report=os.path.join(
            READ_GROUP_QC_DIR, "{sample}_{read_type}.read_group_validation.tsv"
        )
    wildcard_constraints:
        read_type="pe|se"
    params:
        rg_id=lambda wildcards: SAMPLES_INFO[wildcards.sample]["read_group"]["id"],
        rg_sm=lambda wildcards: SAMPLES_INFO[wildcards.sample]["read_group"]["sample"],
        rg_lb=lambda wildcards: SAMPLES_INFO[wildcards.sample]["read_group"]["library"],
        rg_pl=lambda wildcards: SAMPLES_INFO[wildcards.sample]["read_group"]["platform"],
        rg_pu=lambda wildcards: SAMPLES_INFO[wildcards.sample]["read_group"]["platform_unit"]
    log:
        os.path.join(
            "logs", config["pipeline"], "read_group_validation", "{sample}_{read_type}.log"
        )


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

use rule samtools_index_bam_generic as samtools_index_filtered_bam with:
    input:
        bam = os.path.join(FILTERED_BAM_DIR, "{sample}_{read_type}.filtered.sorted.bam")
    output:
        bai = os.path.join(FILTERED_BAM_DIR, "{sample}_{read_type}.filtered.sorted.bam.bai")
    log:
        os.path.join("logs", config["pipeline"], "samtools_index_filtered", "{sample}_{read_type}.log")

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
#print("Final outputs:", final_outputs)
rule all:
    input:
        final_outputs,
        validation_outputs,
        rules.multiqc.output
    default_target: True
