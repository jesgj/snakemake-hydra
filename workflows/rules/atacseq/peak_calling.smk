import os

# --- CONFIGURATION ---
FILTERED_BAM_DIR = config["filtered_bam_dir"]
PEAKS_DIR = config["peaks_dir"]
MACS3_GENOME_SIZE = str(config.get("macs3", {}).get("genome_size", "1.4e9"))
MACS3_QVALUE = config.get("macs3", {}).get("qvalue", 0.05)
MACS3_EXTRA_ARGS = config.get("macs3", {}).get("extra_args", "")


rule macs3_callpeak:
    """
    Calls ATAC-seq peaks with MACS3 using BAMPE and ATAC-specific settings.
    """
    input:
        bam = os.path.join(FILTERED_BAM_DIR, "{sample}_pe.filtered.sorted.bam")
    output:
        narrowpeak = os.path.join(PEAKS_DIR, "{sample}_peaks.narrowPeak"),
        summits = os.path.join(PEAKS_DIR, "{sample}_summits.bed"),
        pileup_bdg = os.path.join(PEAKS_DIR, "{sample}_treat_pileup.bdg"),
        peaks_xls = os.path.join(PEAKS_DIR, "{sample}_peaks.xls")
    params:
        sample_name = lambda wildcards: wildcards.sample,
        outdir = PEAKS_DIR,
        genome_size = MACS3_GENOME_SIZE,
        qvalue = MACS3_QVALUE,
        extra = MACS3_EXTRA_ARGS
    threads: 1
    log:
        os.path.join("logs", config["pipeline"], "macs3", "{sample}.log")
    shell:
        """
        pixi run macs3 callpeak \
            -t {input.bam} \
            -f BAMPE \
            -g {params.genome_size} \
            -n {params.sample_name} \
            --outdir {params.outdir} \
            -q {params.qvalue} \
            --nomodel \
            --shift -75 \
            --extsize 150 \
            --keep-dup all \
            -B \
            --call-summits \
            {params.extra} > {log}.out 2> {log}.err
        """
