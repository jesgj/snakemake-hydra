# workflows/rules/rnaseq/alignment.smk
import os

# --- CONFIGURATION ---
REF_GENOME = config["ref_genome"]
HISAT2_INDEX_DIR = config["hisat2_index_dir"]
TRIMMED_DIR = config["trimmed_dir"]
ALIGNMENT_DIR = config["alignment_dir"]
HISAT2_EXTRA_ARGS = config.get("hisat2", {}).get("extra_args", "")

# The basename for the index, derived from the genome file name
REF_BASENAME = os.path.splitext(os.path.basename(REF_GENOME))[0]
HISAT2_INDEX_PREFIX = os.path.join(HISAT2_INDEX_DIR, REF_BASENAME)

# --- RULES ---

rule hisat2_build:
    """
    Builds a HISAT2 index for the reference genome.
    """
    input:
        ref=REF_GENOME
    output:
        # hisat2-build creates multiple files, we use a sentinel file to track completion
        sentinel=os.path.join(HISAT2_INDEX_DIR, "index_built.OK")
    params:
        prefix=HISAT2_INDEX_PREFIX
    threads: 1
    log:
        os.path.join("logs", config["pipeline"], "hisat2_build", "hisat2_build.log")
    shell:
        "pixi run hisat2-build {input.ref} {params.prefix} > {log}.out 2> {log}.err && touch {output.sentinel}"

rule hisat2_align_pe:
    """
    Aligns trimmed paired-end reads to the reference genome using HISAT2,
    and sorts the output into a BAM file.
    """
    input:
        r1=os.path.join(TRIMMED_DIR, "{sample}_R1.trimmed.fq.gz"),
        r2=os.path.join(TRIMMED_DIR, "{sample}_R2.trimmed.fq.gz"),
        index_sentinel=os.path.join(HISAT2_INDEX_DIR, "index_built.OK")
    output:
        bam=os.path.join(ALIGNMENT_DIR, "{sample}_pe.sorted.bam")
    params:
        extra=HISAT2_EXTRA_ARGS,
        index_prefix=HISAT2_INDEX_PREFIX
    threads: 8
    log:
        os.path.join("logs", config["pipeline"], "hisat2_align", "{sample}_pe.log")
    shell:
        """
        (pixi run hisat2 -p {threads} {params.extra} -x {params.index_prefix} -1 {input.r1} -2 {input.r2} | \
        pixi run samtools view -bS - | \
        pixi run samtools sort -@ {threads} - -o {output.bam}) 2> {log}
        """
