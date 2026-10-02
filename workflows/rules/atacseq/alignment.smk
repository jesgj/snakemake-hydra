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
        "pixi run bowtie2-build {input.ref} {params.prefix} > {log}.out 2> {log}.err && touch {output.sentinel}"


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
        index_prefix = BOWTIE2_INDEX_PREFIX,
        sample = lambda wildcards: wildcards.sample,
        bowtie_threads = lambda wildcards, threads: max(1, threads // 2),
        sort_extra_threads = lambda wildcards, threads: max(0, threads - max(1, threads // 2) - 1)
    threads: 4
    log:
        os.path.join("logs", config["pipeline"], "bowtie2_align", "{sample}_pe.log")
    shell:
        """
        (
          BOWTIE_CMD=(pixi run bowtie2 -p {params.bowtie_threads} {params.extra}
            -x {params.index_prefix:q} --rg-id {params.sample:q}
            --rg "SM:{params.sample}" --rg "LB:{params.sample}"
            --rg "PL:ILLUMINA" --rg "PU:{params.sample}"
            -1 {input.r1:q} -2 {input.r2:q})
          if [ {threads} -eq 1 ]; then
            TMP_SAM={output.bam:q}.alignment.tmp.sam
            trap 'rm -f "$TMP_SAM"' EXIT
            "${{BOWTIE_CMD[@]}}" -S "$TMP_SAM"
            pixi run samtools sort -@ 0 -o {output.bam:q} "$TMP_SAM"
          else
            "${{BOWTIE_CMD[@]}}" | \
            pixi run samtools sort -@ {params.sort_extra_threads} -o {output.bam:q} -
          fi
        ) > {log}.out 2> {log}.err
        """
