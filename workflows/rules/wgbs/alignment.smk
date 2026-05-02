# workflows/rules/alignment.smk
import os
from snakemake.io import touch

# --- CONFIGURATION ---
REF_GENOME = config["ref_genome"]
REF_DIR = os.path.dirname(REF_GENOME)
TRIMMED_DIR = config["trimmed_dir"]
ALIGN_DIR = config["alignment_dir"]
DEDUP_DIR = config["dedup_dir"]
BISMARK_EXTRA_ARGS = config.get("bismark", {}).get("extra_args", "")
BISMARK_THREADS = int(config.get("bismark", {}).get("threads", 66))
BISMARK_PARALLEL = int(config.get("bismark", {}).get("parallel", 8))
LOG_DIR = config.get("log_dir", os.path.join("logs", config["pipeline"]))

if BISMARK_THREADS < 1 or BISMARK_PARALLEL < 1:
    raise ValueError("WGBS bismark.threads and bismark.parallel must be positive integers.")
if BISMARK_THREADS < BISMARK_PARALLEL:
    raise ValueError(
        "WGBS bismark.threads must be greater than or equal to bismark.parallel so each "
        "parallel Bismark worker receives at least one Bowtie2 thread."
    )

BISMARK_BOWTIE2_THREADS = max(1, BISMARK_THREADS // BISMARK_PARALLEL)


# --- RULES ---

rule bismark_genome_preparation:
    """
    Creates a Bismark index for the reference genome if it doesn't exist.
    """
    input:
        REF_GENOME
    output:
        touch(os.path.join(REF_DIR, "Bisulfite_Genome", "bismark_index_created.OK"))
    params:
        ref_dir = REF_DIR
    threads: 1
    log:
        os.path.join(LOG_DIR, "bismark_genome_preparation", "bismark_genome_preparation.log")
    shell:
        """
        pixi run bismark_genome_preparation --bowtie2 --verbose {params.ref_dir} > {log}.out 2> {log}.err
        """

rule bismark_alignment:
    """
    Aligns trimmed paired-end reads to the bisulfite genome using Bismark.
    """
    input:
        r1 = os.path.join(TRIMMED_DIR, "{sample}_R1.trimmed.fq.gz"),
        r2 = os.path.join(TRIMMED_DIR, "{sample}_R2.trimmed.fq.gz"),
        index = os.path.join(REF_DIR, "Bisulfite_Genome", "bismark_index_created.OK")
    output:
        bam = os.path.join(ALIGN_DIR, "{sample}_pe.bam"),
        report = os.path.join(ALIGN_DIR, "{sample}_pe_report.txt")
    params:
        ref_dir = REF_DIR,
        parallel = BISMARK_PARALLEL,
        bowtie2_threads = BISMARK_BOWTIE2_THREADS,
        align_dir = ALIGN_DIR,
        extra = BISMARK_EXTRA_ARGS
    threads: BISMARK_THREADS
    log:
        os.path.join(LOG_DIR, "bismark", "{sample}.log")
    shell:
        """
        set -euo pipefail
        R1_BASE=$(basename {input.r1})
        R1_STEM=${{R1_BASE%%.fastq.gz}}
        R1_STEM=${{R1_STEM%%.fq.gz}}
        pixi run bismark --bowtie2 --parallel {params.parallel} -p {params.bowtie2_threads} {params.extra} \
        --genome {params.ref_dir} \
        -1 {input.r1} -2 {input.r2} \
        -o {params.align_dir} > {log}.out 2> {log}.err
        # Multicore Bismark does not support --basename, so normalize its default
        # paired-end filenames back to the workflow's stable sample-based outputs.
        mv {params.align_dir}/${{R1_STEM}}_bismark_bt2_pe.bam {output.bam}
        mv {params.align_dir}/${{R1_STEM}}_bismark_bt2_PE_report.txt {output.report}
        """

rule deduplicate_bismark:
    """
    Deduplicates a Bismark BAM file and moves it to the deduplication directory.
    """
    input:
        bam = os.path.join(ALIGN_DIR, "{sample}_pe.bam")
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
