# workflows/rules/chip_cr/alignment.smk
import os

# --- CONFIGURATION ---
REF_GENOME = config["ref_genome"]
BOWTIE2_INDEX_DIR = config["bowtie2_index_dir"]
TRIMMED_DIR = config["trimmed_dir"]
ALIGNMENT_DIR = config["alignment_dir"]
BOWTIE2_EXTRA_ARGS = config.get("bowtie2", {}).get("extra_args", "")

# The basename for the index, derived from the genome file name
REF_BASENAME = os.path.splitext(os.path.basename(REF_GENOME))[0]
BOWTIE2_INDEX_PREFIX = os.path.join(BOWTIE2_INDEX_DIR, REF_BASENAME)

# --- RULES ---

rule bowtie2_build:
    """
    Builds a Bowtie2 index for the reference genome if it doesn't exist.
    """
    input:
        ref = REF_GENOME
    output:
        # Bowtie2 build creates multiple files with extensions like .1.bt2, .2.bt2, etc.
        # We use a sentinel file to mark completion.
        sentinel = os.path.join(BOWTIE2_INDEX_DIR, "index_built.OK")
    params:
        prefix = BOWTIE2_INDEX_PREFIX
    threads: 1
    log:
        os.path.join("logs", config["pipeline"], "bowtie2_build.log")
    shell:
        "pixi run bowtie2-build {input.ref} {params.prefix} > {log}.out 2> {log}.err && touch {output.sentinel}"

rule bowtie2_align_pe:
    """
    Aligns trimmed paired-end reads to the reference genome using Bowtie2,
    and sorts the output into a BAM file.
    """
    input:
        r1 = os.path.join(TRIMMED_DIR, "{sample}_R1.trimmed.fq.gz"),
        r2 = os.path.join(TRIMMED_DIR, "{sample}_R2.trimmed.fq.gz"),
        index_sentinel = os.path.join(BOWTIE2_INDEX_DIR, "index_built.OK")
    output:
        bam = os.path.join(ALIGNMENT_DIR, "{sample}_pe.sorted.bam")
    params:
        extra = BOWTIE2_EXTRA_ARGS,
        index_prefix = BOWTIE2_INDEX_PREFIX
    threads: 8
    log:
        os.path.join("logs", config["pipeline"], "bowtie2_align", "{sample}_pe.log")
    shell:
        """
        pixi run bowtie2 -p {threads} {params.extra} -x {params.index_prefix} -1 {input.r1} -2 {input.r2} | \
        pixi run samtools view -bS - | \
        pixi run samtools sort -@ {threads} - -o {output.bam}) 2> {log}
        """

rule bowtie2_align_se:
    """
    Aligns trimmed single-end reads to the reference genome using Bowtie2,
    and sorts the output into a BAM file.
    """
    input:
        r1 = os.path.join(TRIMMED_DIR, "{sample}_SE.trimmed.fq.gz"),
        index_sentinel = os.path.join(BOWTIE2_INDEX_DIR, "index_built.OK")
    output:
        bam = os.path.join(ALIGNMENT_DIR, "{sample}_se.sorted.bam")
    params:
        extra = BOWTIE2_EXTRA_ARGS,
        index_prefix = BOWTIE2_INDEX_PREFIX
    threads: 8
    log:
        os.path.join("logs", config["pipeline"], "bowtie2_align", "{sample}_se.log")
    shell:
        """
        pixi run bowtie2 -p {threads} {params.extra} -x {params.index_prefix} -U {input.r1} | \
        pixi run samtools view -bS - | \
        pixi run samtools sort -@ {threads} - -o {output.bam}) 2> {log}
        """
