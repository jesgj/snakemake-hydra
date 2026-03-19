import os

ALIGNMENT_DIR = config["alignment_dir"]
MARKED_BAM_DIR = config["marked_bam_dir"]
DUPLICATION_QC_DIR = config["duplication_qc_dir"]
PICARD_CONFIG = config.get("picard", {})
PICARD_JAVA_OPTS = PICARD_CONFIG.get("java_opts", "-Xmx4g")
PICARD_MARKDUP_EXTRA_ARGS = PICARD_CONFIG.get("markduplicates", {}).get("extra_args", "")


rule picard_markduplicates:
    """
    Marks duplicate reads in RNA-seq BAM files and writes duplication metrics.
    """
    input:
        bam=os.path.join(ALIGNMENT_DIR, "{sample}_pe.sorted.bam")
    output:
        bam=os.path.join(MARKED_BAM_DIR, "{sample}_pe.markdup.bam"),
        bai=os.path.join(MARKED_BAM_DIR, "{sample}_pe.markdup.bai"),
        metrics=os.path.join(DUPLICATION_QC_DIR, "{sample}_pe.markdup.metrics.txt")
    params:
        java_opts=PICARD_JAVA_OPTS,
        extra=PICARD_MARKDUP_EXTRA_ARGS
    threads: 4
    log:
        os.path.join("logs", config["pipeline"], "markduplicates", "{sample}_pe.log")
    shell:
        """
        pixi run picard {params.java_opts} MarkDuplicates \
            -I "{input.bam}" \
            -O "{output.bam}" \
            -M "{output.metrics}" \
            --CREATE_INDEX true \
            --REMOVE_DUPLICATES false \
            {params.extra} > {log}.out 2> {log}.err
        """
