# workflows/rules/chip_cr/bigwig.smk
import os

# --- CONFIGURATION ---
FILTERED_BAM_DIR = config["filtered_bam_dir"]
BIGWIG_DIR = config["bigwig_dir"]
SUBTRACTED_BIGWIG_DIR = config["subtracted_bigwig_dir"]

# Deeptools parameters
BC_NORM = config.get("deeptools", {}).get("bamCoverage", {}).get("normalize_using", "RPKM")
BC_ARGS = config.get("deeptools", {}).get("bamCoverage", {}).get("extra_args", "")
BWC_METHOD = config.get("deeptools", {}).get("bigwigCompare", {}).get("method", "subtract")
BWC_ARGS = config.get("deeptools", {}).get("bigwigCompare", {}).get("extra_args", "")

# --- BIGWIG GENERATION ---

rule bamCoverage:
    """
    Generates a normalized bigWig file from a filtered BAM file.
    """
    input:
        bam = os.path.join(FILTERED_BAM_DIR, "{sample}_{read_type}.filtered.sorted.bam")
    output:
        bigwig = os.path.join(BIGWIG_DIR, "{sample}_{read_type}.bw")
    params:
        normalize = BC_NORM,
        extra = BC_ARGS
    threads: 8
    log:
        os.path.join("logs", config["pipeline"], "bamCoverage", "{sample}_{read_type}.log")
    shell:
        """
pixi run bamCoverage -b {input.bam} -o {output.bigwig} \
            --normalizeUsing {params.normalize} \
            {params.extra} \
            -p {threads} > {log}.out 2> {log}.err
        """


# --- BIGWIG SUBTRACTION ---

def get_input_bigwig_for_chip(wildcards):
    """
    Finds the corresponding input bigWig for a given ChIP sample.
    Expects pairs in config['subtraction_pairs'] filled by the workflow.
    """
    for pair in config['subtraction_pairs']:
        if pair['sample'] == wildcards.sample and pair['read_type'] == wildcards.read_type:
            return pair['input_bw']

    from snakemake.exceptions import MissingInputException
    raise MissingInputException(f"No input bigwig found for {wildcards.sample}_{wildcards.read_type}")


rule bigwigCompare:
    """
    Subtracts the input signal from the ChIP signal using bigwigCompare.
    """
    input:
        chip_bw = os.path.join(BIGWIG_DIR, "{sample}_{read_type}.bw"),
        input_bw = get_input_bigwig_for_chip
    output:
        subtracted_bw = os.path.join(SUBTRACTED_BIGWIG_DIR, "{sample}_{read_type}.subtracted.bw")
    params:
        method = BWC_METHOD,
        extra = BWC_ARGS
    threads: 1
    log:
        os.path.join("logs", config["pipeline"], "deeptools", "{sample}_{read_type}_bigwigCompare.log")
    shell:
        """
pixi run bigwigCompare \
            -b1 {input.chip_bw} \
            -b2 {input.input_bw} \
            -o {output.subtracted_bw} \
            --operation {params.method} \
            {params.extra} > {log}.out 2> {log}.err
        """

