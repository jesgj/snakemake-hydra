import os
import math

# --- CONFIGURATION ---
FILTERED_BAM_DIR = config["filtered_bam_dir"]
PEAKS_DIR = config["peaks_dir"]
PEAK_CALLER = config.get("peak_caller", "macs3")
MACS3_GENOME_SIZE = str(config.get("macs3", {}).get("genome_size", "1.4e9"))
MACS3_QVALUE = float(config.get("macs3", {}).get("qvalue", 0.05))
MACS3_EXTRA_ARGS = config.get("macs3", {}).get("extra_args", "")
GENRICH_QVALUE = float(config.get("genrich", {}).get("qvalue", 0.05))
GENRICH_MIN_AUC = float(config.get("genrich", {}).get("min_auc", 20.0))
GENRICH_MIN_LENGTH = int(config.get("genrich", {}).get("min_length", 0))
GENRICH_MAX_GAP = int(config.get("genrich", {}).get("max_gap", 100))
GENRICH_EXTRA_ARGS = config.get("genrich", {}).get("extra_args", "")

if not 0 < MACS3_QVALUE <= 1:
    raise ValueError("ATAC-seq macs3.qvalue must be greater than 0 and at most 1.")
if MACS3_GENOME_SIZE not in {"hs", "mm", "ce", "dm"}:
    try:
        genome_size_value = float(MACS3_GENOME_SIZE)
        if not math.isfinite(genome_size_value) or genome_size_value <= 0:
            raise ValueError
    except ValueError:
        raise ValueError(
            "ATAC-seq macs3.genome_size must be positive or one of hs, mm, ce, or dm."
        )
if not 0 < GENRICH_QVALUE <= 1:
    raise ValueError("ATAC-seq genrich.qvalue must be greater than 0 and at most 1.")
if (
    not math.isfinite(GENRICH_MIN_AUC)
    or GENRICH_MIN_AUC < 0
    or GENRICH_MIN_LENGTH < 0
    or GENRICH_MAX_GAP < 0
):
    raise ValueError(
        "ATAC-seq Genrich min_auc, min_length, and max_gap values cannot be negative."
    )


if PEAK_CALLER == "macs3":
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
                --keep-dup all \
                -B \
                --call-summits \
                {params.extra} > {log}.out 2> {log}.err
            """
elif PEAK_CALLER == "genrich":
    rule samtools_sort_name_for_genrich:
        """
        Creates a queryname-sorted BAM intermediate for Genrich peak calling.
        """
        input:
            bam = os.path.join(FILTERED_BAM_DIR, "{sample}_pe.filtered.sorted.bam")
        output:
            bam = temp(os.path.join(FILTERED_BAM_DIR, "{sample}_pe.filtered.qname.bam"))
        threads: 4
        log:
            os.path.join("logs", config["pipeline"], "genrich_qname_sort", "{sample}.log")
        shell:
            """
            pixi run samtools sort -n -@ {threads} -o "{output.bam}" "{input.bam}" > {log}.out 2> {log}.err
            """

    rule genrich_callpeak:
        """
        Calls ATAC-seq peaks with Genrich in ATAC-seq mode.
        """
        input:
            bam = os.path.join(FILTERED_BAM_DIR, "{sample}_pe.filtered.qname.bam")
        output:
            narrowpeak = os.path.join(PEAKS_DIR, "{sample}_peaks.narrowPeak")
        params:
            qvalue = GENRICH_QVALUE,
            min_auc = GENRICH_MIN_AUC,
            min_length = GENRICH_MIN_LENGTH,
            max_gap = GENRICH_MAX_GAP,
            extra = GENRICH_EXTRA_ARGS
        threads: 1
        log:
            os.path.join("logs", config["pipeline"], "genrich", "{sample}.log")
        shell:
            """
            pixi run Genrich \
                -t "{input.bam}" \
                -o "{output.narrowpeak}" \
                -j \
                -q {params.qvalue} \
                -a {params.min_auc} \
                -l {params.min_length} \
                -g {params.max_gap} \
                -v \
                {params.extra} > {log}.out 2> {log}.err
            """
