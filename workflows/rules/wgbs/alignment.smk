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
        os.path.join("logs", config["pipeline"], "bismark_genome_preparation", "bismark_genome_preparation.log")
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
        parallel = 8,
        align_dir = ALIGN_DIR,
        extra = BISMARK_EXTRA_ARGS
    threads: 66
    log:
        os.path.join("logs", config["pipeline"], "bismark", "{sample}.log")
    shell:
        """
        pixi run bismark --bowtie2 -p {params.parallel} {params.extra} \
        --genome {params.ref_dir} \
        -1 {input.r1} -2 {input.r2} \
        -o {params.align_dir} --basename {wildcards.sample} > {log}.out 2> {log}.err
        # Bismark creates a report file with `_PE_` instead of `_pe_`.
        # We rename it to have a consistent naming convention.
        mv {params.align_dir}/{wildcards.sample}_PE_report.txt {output.report}
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
        os.path.join("logs", config["pipeline"], "deduplicate_bismark", "{sample}.log")
    shell:
        """
        pixi run deduplicate_bismark -p --bam {input.bam} --output_dir {params.outdir} > {log}.out 2> {log}.err
        """
