# workflows/rules/bam_qc.smk
import os

# This file defines a generic set of BAM QC rules.
# They are intended to be instantiated using `use rule ... as ... with:`,
# which will override the input, output, log, and other parameters.

rule samtools_stats_generic:
    """
    Generic rule for samtools stats. Override I/O in parent workflow.
    """
    input:
        bam = "path/to/input.bam"
    output:
        stats = "path/to/output.stats.txt"
    threads: 4
    log:
        "logs/samtools_stats.log"
    shell:
        "pixi run samtools stats -@ {threads} {input.bam} > {output.stats} 2> {log}"

rule samtools_flagstat_generic:
    """
    Generic rule for samtools flagstat. Override I/O in parent workflow.
    """
    input:
        bam = "path/to/input.bam"
    output:
        flagstat = "path/to/output.flagstat.txt"
    threads: 1
    log:
        "logs/samtools_flagstat.log"
    shell:
        "pixi run samtools flagstat {input.bam} > {output.flagstat} 2> {log}"

rule picard_collect_alignment_metrics_generic:
    """
    Generic rule for Picard CollectAlignmentSummaryMetrics. Override I/O in parent workflow.
    """
    input:
        bam = "path/to/input.bam",
        ref = config["ref_genome"]
    output:
        metrics = "path/to/output.metrics.txt"
    threads: 4
    params:
        java_opts=config.get("picard", {}).get("java_opts", "-Xmx4g")
    log:
        "logs/picard_metrics.log"
    shell:
        """
        mkdir -p "$(dirname "{output.metrics}")"
        pixi run picard {params.java_opts} CollectAlignmentSummaryMetrics \
            -I {input.bam} \
            -O {output.metrics} --IS_BISULFITE_SEQUENCED true > {log}.out 2> {log}.err
        """

