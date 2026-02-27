import os

# --- CONFIGURATION ---
FILTERED_BAM_DIR = config["filtered_bam_dir"]
FRAGMENT_QC_DIR = config["fragment_qc_dir"]
PICARD_JAVA_OPTS = config.get("picard", {}).get("java_opts", "-Xmx2g")


rule picard_collect_insert_size_metrics:
    """
    Generates insert-size metrics and histogram for ATAC-seq fragment QC.
    """
    input:
        bam = os.path.join(FILTERED_BAM_DIR, "{sample}_pe.filtered.sorted.bam")
    output:
        metrics = os.path.join(FRAGMENT_QC_DIR, "{sample}_pe.insert_size_metrics.txt"),
        histogram = os.path.join(FRAGMENT_QC_DIR, "{sample}_pe.insert_size_histogram.pdf")
    params:
        java_opts = PICARD_JAVA_OPTS
    threads: 1
    log:
        os.path.join("logs", config["pipeline"], "fragment_qc", "{sample}.log")
    shell:
        """
        pixi run picard {params.java_opts} CollectInsertSizeMetrics \
            -I {input.bam} \
            -O {output.metrics} \
            -H {output.histogram} \
            -M 0.5 > {log}.out 2> {log}.err
        """
