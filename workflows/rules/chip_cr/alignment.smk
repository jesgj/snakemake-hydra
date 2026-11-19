# workflows/rules/chip_cr/alignment.smk
import os

# --- CONFIGURATION ---
REF_GENOME = config["ref_genome"]
BOWTIE2_INDEX_DIR = config["bowtie2_index_dir"]
TRIMMED_DIR = config["trimmed_dir"]
ALIGNMENT_DIR = config["alignment_dir"]
SAMPLES_INFO = config["samples_info"]
BOWTIE2_EXTRA_ARGS = config.get("bowtie2", {}).get("extra_args", "")

# The basename for the index, derived from the genome file name
REF_BASENAME = os.path.splitext(os.path.basename(REF_GENOME))[0]
BOWTIE2_INDEX_PREFIX = os.path.join(BOWTIE2_INDEX_DIR, REF_BASENAME)
BOWTIE2_INDEX_FILES = [
    f"{BOWTIE2_INDEX_PREFIX}.1.bt2",
    f"{BOWTIE2_INDEX_PREFIX}.2.bt2",
    f"{BOWTIE2_INDEX_PREFIX}.3.bt2",
    f"{BOWTIE2_INDEX_PREFIX}.4.bt2",
    f"{BOWTIE2_INDEX_PREFIX}.rev.1.bt2",
    f"{BOWTIE2_INDEX_PREFIX}.rev.2.bt2",
]

# --- RULES ---

rule bowtie2_build:
    """
    Builds a Bowtie2 index for the reference genome if it doesn't exist.
    """
    input:
        ref = REF_GENOME
    output:
        index = BOWTIE2_INDEX_FILES
    params:
        prefix = BOWTIE2_INDEX_PREFIX
    threads: 1
    log:
        os.path.join("logs", config["pipeline"], "bowtie2_build.log")
    shell:
        "pixi run bowtie2-build \"{input.ref}\" \"{params.prefix}\" > {log}.out 2> {log}.err"

rule bowtie2_align_pe:
    """
    Aligns trimmed paired-end reads to the reference genome using Bowtie2,
    and sorts the output into a BAM file.
    """
    input:
        r1 = os.path.join(TRIMMED_DIR, "{sample}_R1.trimmed.fq.gz"),
        r2 = os.path.join(TRIMMED_DIR, "{sample}_R2.trimmed.fq.gz"),
        index = BOWTIE2_INDEX_FILES
    output:
        bam = os.path.join(ALIGNMENT_DIR, "{sample}_pe.sorted.bam")
    params:
        extra = BOWTIE2_EXTRA_ARGS,
        index_prefix = BOWTIE2_INDEX_PREFIX,
        rg_id = lambda wildcards: SAMPLES_INFO[wildcards.sample]["read_group"]["id"],
        rg_sm = lambda wildcards: f"SM:{SAMPLES_INFO[wildcards.sample]['read_group']['sample']}",
        rg_lb = lambda wildcards: f"LB:{SAMPLES_INFO[wildcards.sample]['read_group']['library']}",
        rg_pl = lambda wildcards: f"PL:{SAMPLES_INFO[wildcards.sample]['read_group']['platform']}",
        rg_pu = lambda wildcards: f"PU:{SAMPLES_INFO[wildcards.sample]['read_group']['platform_unit']}"
    threads: 8
    log:
        out=os.path.join("logs", config["pipeline"], "bowtie2_align", "{sample}_pe.log.out"),
        err=os.path.join("logs", config["pipeline"], "bowtie2_align", "{sample}_pe.log.err")
    shell:
        """
        set -euo pipefail
        (pixi run bowtie2 -p {threads} {params.extra} -x {params.index_prefix:q} \
        --rg-id {params.rg_id:q} \
        --rg {params.rg_sm:q} --rg {params.rg_lb:q} \
        --rg {params.rg_pl:q} --rg {params.rg_pu:q} \
        -1 {input.r1:q} -2 {input.r2:q} | \
        pixi run samtools view -bS - | \
        pixi run samtools sort -@ {threads} - -o {output.bam:q}) > {log.out:q} 2> {log.err:q}
        """

rule bowtie2_align_se:
    """
    Aligns trimmed single-end reads to the reference genome using Bowtie2,
    and sorts the output into a BAM file.
    """
    input:
        r1 = os.path.join(TRIMMED_DIR, "{sample}_SE.trimmed.fq.gz"),
        index = BOWTIE2_INDEX_FILES
    output:
        bam = os.path.join(ALIGNMENT_DIR, "{sample}_se.sorted.bam")
    params:
        extra = BOWTIE2_EXTRA_ARGS,
        index_prefix = BOWTIE2_INDEX_PREFIX,
        rg_id = lambda wildcards: SAMPLES_INFO[wildcards.sample]["read_group"]["id"],
        rg_sm = lambda wildcards: f"SM:{SAMPLES_INFO[wildcards.sample]['read_group']['sample']}",
        rg_lb = lambda wildcards: f"LB:{SAMPLES_INFO[wildcards.sample]['read_group']['library']}",
        rg_pl = lambda wildcards: f"PL:{SAMPLES_INFO[wildcards.sample]['read_group']['platform']}",
        rg_pu = lambda wildcards: f"PU:{SAMPLES_INFO[wildcards.sample]['read_group']['platform_unit']}"
    threads: 8
    log:
        out=os.path.join("logs", config["pipeline"], "bowtie2_align", "{sample}_se.log.out"),
        err=os.path.join("logs", config["pipeline"], "bowtie2_align", "{sample}_se.log.err")
    shell:
        """
        set -euo pipefail
        (pixi run bowtie2 -p {threads} {params.extra} -x {params.index_prefix:q} \
        --rg-id {params.rg_id:q} \
        --rg {params.rg_sm:q} --rg {params.rg_lb:q} \
        --rg {params.rg_pl:q} --rg {params.rg_pu:q} \
        -U {input.r1:q} | \
        pixi run samtools view -bS - | \
        pixi run samtools sort -@ {threads} - -o {output.bam:q}) > {log.out:q} 2> {log.err:q}
        """
