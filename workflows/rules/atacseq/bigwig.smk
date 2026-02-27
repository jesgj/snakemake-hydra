import os

# --- CONFIGURATION ---
FILTERED_BAM_DIR = config["filtered_bam_dir"]
BIGWIG_DIR = config["bigwig_dir"]
BC_NORM = config.get("deeptools", {}).get("bamCoverage", {}).get("normalize_using", "CPM")
BC_ARGS = config.get("deeptools", {}).get("bamCoverage", {}).get("extra_args", "")


rule bamCoverage:
    """
    Generates normalized bigWig tracks from filtered BAMs.
    """
    input:
        bam = os.path.join(FILTERED_BAM_DIR, "{sample}_pe.filtered.sorted.bam")
    output:
        bigwig = os.path.join(BIGWIG_DIR, "{sample}_pe.bw")
    params:
        normalize = BC_NORM,
        extra = BC_ARGS
    threads: 4
    log:
        os.path.join("logs", config["pipeline"], "bamCoverage", "{sample}_pe.log")
    shell:
        """
        pixi run bamCoverage -b {input.bam} -o {output.bigwig} \
            --normalizeUsing {params.normalize} \
            {params.extra} \
            -p {threads} > {log}.out 2> {log}.err
        """
