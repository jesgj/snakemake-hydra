# workflows/rules/rnaseq/alignment.smk
import os

# --- CONFIGURATION ---
REF_GENOME = config["ref_genome"]
HISAT2_INDEX_DIR = config["hisat2_index_dir"]
TRIMMED_DIR = config["trimmed_dir"]
ALIGNMENT_DIR = config["alignment_dir"]
SAMPLES_INFO = config["samples_info"]
HISAT2_EXTRA_ARGS = config.get("hisat2", {}).get("extra_args", "")

# The basename for the index, derived from the genome file name
REF_BASENAME = os.path.splitext(os.path.basename(REF_GENOME))[0]
HISAT2_INDEX_PREFIX = os.path.join(HISAT2_INDEX_DIR, REF_BASENAME)
HISAT2_INDEX_FILES = expand(f"{HISAT2_INDEX_PREFIX}.{{idx}}.ht2", idx=range(1, 9))

# --- RULES ---

rule hisat2_build:
    """
    Builds a HISAT2 index for the reference genome.
    """
    input:
        ref=REF_GENOME
    output:
        index=HISAT2_INDEX_FILES
    params:
        prefix=HISAT2_INDEX_PREFIX
    threads: 1
    log:
        os.path.join("logs", config["pipeline"], "hisat2_build", "hisat2_build.log")
    shell:
        "pixi run hisat2-build {input.ref} {params.prefix} > {log}.out 2> {log}.err"

rule hisat2_align_pe:
    """
    Aligns trimmed paired-end reads to the reference genome using HISAT2,
    and sorts the output into a BAM file.
    """
    input:
        r1=os.path.join(TRIMMED_DIR, "{sample}_R1.trimmed.fq.gz"),
        r2=os.path.join(TRIMMED_DIR, "{sample}_R2.trimmed.fq.gz"),
        index=HISAT2_INDEX_FILES
    output:
        bam=os.path.join(ALIGNMENT_DIR, "{sample}_pe.sorted.bam"),
        summary=os.path.join(ALIGNMENT_DIR, "{sample}_pe.hisat2.summary.txt")
    params:
        extra=HISAT2_EXTRA_ARGS,
        index_prefix=HISAT2_INDEX_PREFIX,
        rg_id=lambda wildcards: SAMPLES_INFO[wildcards.sample]["read_group"]["id"],
        rg_sm=lambda wildcards: f"SM:{SAMPLES_INFO[wildcards.sample]['read_group']['sample']}",
        rg_lb=lambda wildcards: f"LB:{SAMPLES_INFO[wildcards.sample]['read_group']['library']}",
        rg_pl=lambda wildcards: f"PL:{SAMPLES_INFO[wildcards.sample]['read_group']['platform']}",
        rg_pu=lambda wildcards: f"PU:{SAMPLES_INFO[wildcards.sample]['read_group']['platform_unit']}"
    threads: 8
    log:
        os.path.join("logs", config["pipeline"], "hisat2_align", "{sample}_pe.log")
    shell:
        """
        set -euo pipefail
        (pixi run hisat2 -p {threads} {params.extra} \
        --rg-id {params.rg_id:q} \
        --rg {params.rg_sm:q} \
        --rg {params.rg_lb:q} \
        --rg {params.rg_pl:q} \
        --rg {params.rg_pu:q} \
        --summary-file {output.summary:q} -x {params.index_prefix:q} \
        -1 {input.r1:q} -2 {input.r2:q} | \
        pixi run samtools view -bS - | \
        pixi run samtools sort -@ {threads} - -o {output.bam:q}) > {log:q}.out 2> {log:q}.err
        """

rule samtools_index_aligned_bam:
    """
    Builds an index for a sorted RNA-seq BAM file.
    """
    input:
        bam = os.path.join(ALIGNMENT_DIR, "{sample}_pe.sorted.bam")
    output:
        bai = os.path.join(ALIGNMENT_DIR, "{sample}_pe.sorted.bam.bai")
    threads: 4
    log:
        os.path.join("logs", config["pipeline"], "samtools_index", "{sample}_pe.log")
    shell:
        "pixi run samtools index -@ {threads} {input.bam} {output.bai} > {log}.out 2> {log}.err"
