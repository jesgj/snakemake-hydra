# workflows/rules/alignment.smk
import os

# --- CONFIGURATION ---
SOURCE_REF_GENOME = config["ref_genome"]
REF_DIR = config["reference_dir"]
STAGED_REF_GENOME = config["reference_fasta"]
STAGED_REF_FAI = config["reference_fasta_fai"]
TRIMMED_DIR = config["trimmed_dir"]
ALIGN_DIR = config["alignment_dir"]
DEDUP_DIR = config["dedup_dir"]
READ_GROUP_QC_DIR = config.get(
    "read_group_qc_dir", os.path.join("results", "wgbs", "read_group_qc")
)
SAMPLES_INFO = config["samples_info"]
BISMARK_EXTRA_ARGS = config.get("bismark", {}).get("extra_args", "")
BISMARK_THREADS = int(config.get("bismark", {}).get("threads", 66))
BISMARK_PARALLEL = int(config.get("bismark", {}).get("parallel", 8))
BISMARK_LARGE_INDEX = config.get("bismark", {}).get("large_index", False)
if not isinstance(BISMARK_LARGE_INDEX, bool):
    raise ValueError("WGBS bismark.large_index must be true or false.")
LOG_DIR = config.get("log_dir", os.path.join("logs", config["pipeline"]))

if BISMARK_THREADS < 1 or BISMARK_PARALLEL < 1:
    raise ValueError("WGBS bismark.threads and bismark.parallel must be positive integers.")
if BISMARK_THREADS < BISMARK_PARALLEL:
    raise ValueError(
        "WGBS bismark.threads must be greater than or equal to bismark.parallel so each "
        "parallel Bismark worker receives at least one Bowtie2 thread."
    )


def _read_group_line(wildcards):
    read_group = SAMPLES_INFO[wildcards.sample]["read_group"]
    fields = [
        "@RG",
        f"ID:{read_group['id']}",
        f"SM:{read_group['sample']}",
        f"LB:{read_group['library']}",
        f"PL:{read_group['platform']}",
        f"PU:{read_group['platform_unit']}",
    ]
    return "\t".join(fields)


# Pinned Bismark 3.1.0 emits both converted genomes and six Bowtie2 files per
# conversion. These are the complete, real index outputs; the staging directory
# contains only the staged reference FASTA besides Bisulfite_Genome.
BISMARK_INDEX_DIR = os.path.join(REF_DIR, "Bisulfite_Genome")
BISMARK_INDEX_EXTENSION = "bt2l" if BISMARK_LARGE_INDEX else "bt2"
BISMARK_INDEX_SUFFIXES = [f".{part}.{BISMARK_INDEX_EXTENSION}" for part in ("1", "2", "3", "4", "rev.1", "rev.2")]
BISMARK_CONVERTED_STEM = f"{os.path.splitext(os.path.basename(STAGED_REF_GENOME))[0]}_mfa"
BISMARK_CONVERSIONS = (("CT_conversion", "BS_CT"), ("GA_conversion", "BS_GA"))
BISMARK_INDEX_FILES = [
    os.path.join(BISMARK_INDEX_DIR, conversion, f"{prefix}{suffix}")
    for conversion, prefix in BISMARK_CONVERSIONS
    for suffix in BISMARK_INDEX_SUFFIXES
] + [
    os.path.join(BISMARK_INDEX_DIR, conversion, f"{BISMARK_CONVERTED_STEM}.{conversion}.fa")
    for conversion, _ in BISMARK_CONVERSIONS
]

# --- RULES ---

rule stage_bismark_reference:
    """
    Copies the configured reference into workflow-owned staging.
    """
    input:
        ref = SOURCE_REF_GENOME
    output:
        fasta = STAGED_REF_GENOME
    params:
        staging_dir = REF_DIR
    threads: 1
    log:
        os.path.join(LOG_DIR, "bismark_reference", "stage_reference.log")
    shell:
        """
        set -euo pipefail
        mkdir -p -- {params.staging_dir:q}
        rm -f -- {output.fasta:q}
        cp -- {input.ref:q} {output.fasta:q} > {log:q}.out 2> {log:q}.err
        """

rule index_bismark_reference_fasta:
    """
    Creates the staged reference .fai used by downstream methylation tools.
    """
    input:
        fasta = STAGED_REF_GENOME
    output:
        fai = STAGED_REF_FAI
    threads: 1
    log:
        os.path.join(LOG_DIR, "bismark_reference", "samtools_faidx.log")
    shell:
        """
        set -euo pipefail
        rm -f {output.fai:q}
        pixi run samtools faidx {input.fasta:q} > {log:q}.out 2> {log:q}.err
        """

rule bismark_genome_preparation:
    """
    Builds and tracks every Bismark index component for the staged reference.
    """
    input:
        ref = STAGED_REF_GENOME
    output:
        index = BISMARK_INDEX_FILES
    params:
        ref_dir = REF_DIR,
        index_dir = BISMARK_INDEX_DIR,
        large_index = "--large-index" if BISMARK_LARGE_INDEX else ""
    threads: 1
    log:
        os.path.join(LOG_DIR, "bismark_genome_preparation", "bismark_genome_preparation.log")
    shell:
        """
        set -euo pipefail
        # Bismark scans all FASTAs directly under the genome directory.
        unexpected=$(find -L {params.ref_dir:q} -maxdepth 1 -type f \
            \\( -iname '*.fa' -o -iname '*.fasta' -o -iname '*.fna' \\) \
            ! -name genome.fa -print -quit)
        if [[ -n "$unexpected" ]]; then
            printf 'Unexpected FASTA in WGBS staging directory: %s\\n' "$unexpected" >&2
            exit 1
        fi
        # Bismark cannot safely resume a partial Bisulfite_Genome directory.
        rm -rf -- {params.index_dir:q}
        pixi run bismark_genome_preparation --bowtie2 --verbose {params.large_index} {params.ref_dir:q} > {log:q}.out 2> {log:q}.err
        """

rule bismark_alignment:
    """
    Aligns trimmed paired-end reads to the bisulfite genome using Bismark.
    """
    input:
        r1 = os.path.join(TRIMMED_DIR, "{sample}_R1.trimmed.fq.gz"),
        r2 = os.path.join(TRIMMED_DIR, "{sample}_R2.trimmed.fq.gz"),
        genome = STAGED_REF_GENOME,
        index = BISMARK_INDEX_FILES
    output:
        bam = os.path.join(ALIGN_DIR, "{sample}_pe.bam"),
        report = os.path.join(ALIGN_DIR, "{sample}_pe_report.txt")
    params:
        ref_dir = REF_DIR,
        parallel = lambda wildcards, threads: min(BISMARK_PARALLEL, threads),
        # Pinned Bismark rejects explicit -p 1; omission uses one thread.
        bowtie2_args = lambda wildcards, threads: (
            f"-p {threads // min(BISMARK_PARALLEL, threads)}"
            if threads // min(BISMARK_PARALLEL, threads) > 1 else ""
        ),
        align_dir = ALIGN_DIR,
        extra = BISMARK_EXTRA_ARGS,
        rg_line = _read_group_line
    threads: BISMARK_THREADS
    log:
        os.path.join(LOG_DIR, "bismark", "{sample}.log")
    shell:
        """
        set -euo pipefail
        R1_BASE=$(basename {input.r1:q})
        R1_STEM=${{R1_BASE%%.fastq.gz}}
        R1_STEM=${{R1_STEM%%.fq.gz}}
        pixi run bismark --bowtie2 --parallel {params.parallel} {params.bowtie2_args} {params.extra} \
        --genome {params.ref_dir:q} \
        -1 {input.r1:q} -2 {input.r2:q} \
        -o {params.align_dir:q} > {log:q}.out 2> {log:q}.err
        # Multicore Bismark does not support --basename, so normalize its default
        # paired-end filenames back to the workflow's stable sample-based outputs.
        mv {params.align_dir:q}/"${{R1_STEM}}_bismark_bt2_pe.bam" {output.bam:q}
        mv {params.align_dir:q}/"${{R1_STEM}}_bismark_bt2_PE_report.txt" {output.report:q}
        # Perl Bismark 3.1 can set only ID/SM; add the complete deterministic
        # read-group header and tag for downstream Picard-compatible BAMs.
        TMP_RG={output.bam:q}.rg.tmp.bam
        trap 'rm -f "$TMP_RG"' EXIT
        pixi run samtools addreplacerg -m overwrite_all \
            -r {params.rg_line:q} -o "$TMP_RG" {output.bam:q} >> {log:q}.out 2>> {log:q}.err
        mv "$TMP_RG" {output.bam:q}
        """

rule deduplicate_bismark:
    """
    Deduplicates a Bismark BAM file and moves it to the deduplication directory.
    """
    input:
        bam = os.path.join(ALIGN_DIR, "{sample}_pe.bam"),
        read_group_validation = os.path.join(
            READ_GROUP_QC_DIR, "{sample}_pe.read_group_validation.tsv"
        )
    output:
        dedup_bam = os.path.join(DEDUP_DIR, "{sample}_pe.deduplicated.bam"),
        report = os.path.join(DEDUP_DIR, "{sample}_pe.deduplication_report.txt")
    params:
        outdir = DEDUP_DIR
    log:
        os.path.join(LOG_DIR, "deduplicate_bismark", "{sample}.log")
    shell:
        """
        pixi run deduplicate_bismark -p --bam {input.bam} --output_dir {params.outdir} > {log}.out 2> {log}.err
        """
