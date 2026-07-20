import os

ALIGNMENT_DIR = config["alignment_dir"]
MARKED_BAM_DIR = config["marked_bam_dir"]
DUPLICATION_QC_DIR = config["duplication_qc_dir"]
SAMBAMBA_CONFIG = config.get("sambamba", {})
SAMBAMBA_MARKDUP_THREADS = max(1, int(SAMBAMBA_CONFIG.get("markdup_threads", 4)))
SAMBAMBA_MARKDUP_EXTRA_ARGS = SAMBAMBA_CONFIG.get("markdup_extra_args", "")


rule sambamba_markduplicates:
    """
    Marks duplicate reads in ATAC-seq BAM files and writes the Sambamba report.
    """
    input:
        bam=os.path.join(ALIGNMENT_DIR, "{sample}_pe.sorted.bam")
    output:
        bam=temp(os.path.join(MARKED_BAM_DIR, "{sample}_pe.markdup.bam")),
        metrics=os.path.join(DUPLICATION_QC_DIR, "{sample}_pe.markdup.metrics.txt")
    params:
        extra=SAMBAMBA_MARKDUP_EXTRA_ARGS
    threads: SAMBAMBA_MARKDUP_THREADS
    shell:
        """
        pixi run sambamba markdup \
            -t {threads} \
            {params.extra} \
            "{input.bam}" \
            "{output.bam}" > "{output.metrics}" 2>&1
        """
