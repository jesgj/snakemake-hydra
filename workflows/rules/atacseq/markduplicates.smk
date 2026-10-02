import os

ALIGNMENT_DIR = config["alignment_dir"]
READ_GROUP_QC_DIR = config["read_group_qc_dir"]
MARKED_BAM_DIR = config["marked_bam_dir"]
DUPLICATION_QC_DIR = config["duplication_qc_dir"]
PICARD_CONFIG = config.get("picard", {})
PICARD_JAVA_OPTS = PICARD_CONFIG.get("java_opts", "-Xmx2g")
PICARD_MARKDUP_EXTRA_ARGS = PICARD_CONFIG.get("markduplicates", {}).get(
    "extra_args", ""
)


rule picard_markduplicates:
    """
    Marks duplicate reads in ATAC-seq BAM files and writes Picard metrics.
    """
    input:
        bam=os.path.join(ALIGNMENT_DIR, "{sample}_pe.sorted.bam"),
        merge_validation=os.path.join(
            READ_GROUP_QC_DIR, "{sample}_merged.read_group_validation.tsv"
        )
    output:
        bam=temp(os.path.join(MARKED_BAM_DIR, "{sample}_pe.markdup.bam")),
        metrics=os.path.join(DUPLICATION_QC_DIR, "{sample}_pe.markdup.metrics.txt")
    params:
        java_opts=PICARD_JAVA_OPTS,
        extra=PICARD_MARKDUP_EXTRA_ARGS
    threads: 1
    log:
        os.path.join("logs", config["pipeline"], "markduplicates", "{sample}_pe.log")
    shell:
        """
        pixi run picard {params.java_opts} MarkDuplicates \
            --INPUT {input.bam:q} \
            --OUTPUT {output.bam:q} \
            --METRICS_FILE {output.metrics:q} \
            --REMOVE_DUPLICATES false \
            --ASSUME_SORT_ORDER coordinate \
            {params.extra} \
            > {log:q}.out 2> {log:q}.err
        """
