import os

# --- CONFIGURATION ---
REF_GENOME = config["ref_genome"]
BOWTIE2_INDEX_DIR = config["bowtie2_index_dir"]
TRIMMED_DIR = config["trimmed_dir"]
ALIGNMENT_DIR = config["alignment_dir"]
BOWTIE2_EXTRA_ARGS = config.get("bowtie2", {}).get("extra_args", "")

REF_BASENAME = os.path.basename(REF_GENOME)
for suffix in [".fa.gz", ".fasta.gz", ".fna.gz", ".fa", ".fasta", ".fna"]:
    if REF_BASENAME.endswith(suffix):
        REF_BASENAME = REF_BASENAME[: -len(suffix)]
        break
if not REF_BASENAME:
    REF_BASENAME = "genome"
BOWTIE2_INDEX_PREFIX = os.path.join(BOWTIE2_INDEX_DIR, REF_BASENAME)


rule bowtie2_build:
    """
    Builds a Bowtie2 index for the reference genome if it does not exist.
    """
    input:
        ref = REF_GENOME
    output:
        sentinel = os.path.join(BOWTIE2_INDEX_DIR, "index_built.OK")
    params:
        prefix = BOWTIE2_INDEX_PREFIX
    threads: 1
    log:
        os.path.join("logs", config["pipeline"], "bowtie2_build", "bowtie2_build.log")
    shell:
        "pixi run bowtie2-build {input.ref} {params.prefix} > {log}.out 2> {log}.err && "
        "pixi run bowtie2-inspect -n {params.prefix} > /dev/null 2>> {log}.err && "
        "touch {output.sentinel}"


rule bowtie2_align_pe:
    """
    Aligns trimmed paired-end reads with Bowtie2 and sorts the resulting BAM.
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
    threads: 4
    log:
        os.path.join("logs", config["pipeline"], "bowtie2_align", "{sample}_pe.log")
    shell:
        """
        set -euo pipefail
        pixi run bowtie2-inspect -n {params.index_prefix} > /dev/null 2>> {log}.err
        (pixi run bowtie2 -p {threads} {params.extra} -x {params.index_prefix} -1 {input.r1} -2 {input.r2} | \
        pixi run samtools view -bS - | \
        pixi run samtools sort -@ {threads} -o {output.bam} -) > {log}.out 2> {log}.err
        """
