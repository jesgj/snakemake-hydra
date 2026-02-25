# workflows/rules/wgbs/methyldackel.smk
import os

# --- CONFIGURATION ---
SORTED_FILTERED_BAM_DIR = config["sorted_filtered_bam_dir"]
REF_GENOME = config["ref_genome"]
MBIAS_DIR = config["mbias_dir"]
METHYLDACKEL_DIR = config["methyldackel_dir"]
METHYLDACKEL_MERGECONTEXT_DIR = config["methyldackel_mergecontext_dir"]
METHYLDACKEL_EXTRACT_EXTRA_ARGS = config.get("methyldackel_extract", {}).get("extra_args", "")


# --- RULES ---

rule methyldackel_extract_methylkit:
    """
    Extracts methylation calls for methylKit.
    """
    input:
        bam=os.path.join(SORTED_FILTERED_BAM_DIR, "{sample}_pe.filtered.sorted.bam"),
        ref=REF_GENOME,
        options=os.path.join(MBIAS_DIR, "{sample}.options.txt")
    output:
        methylkit_file=os.path.join(METHYLDACKEL_DIR, "{sample}_CpG.methylKit")
    params:
        output_prefix=os.path.join(METHYLDACKEL_DIR, "{sample}"),
        extra_opts=config.get("methyldackel_extract", {}).get("extra_opts", "--minOppositeDepth 10 --maxVariantFrac 0.5"),
        extra=METHYLDACKEL_EXTRACT_EXTRA_ARGS
    threads: 1
    log:
        os.path.join("logs", config["pipeline"], "methyldackel_extract_methylkit", "{sample}.log")
    shell:
        """
        options=$(cat {input.options})
        echo ${{options}}
        pixi run MethylDackel extract --methylKit ${{options}} {params.extra_opts} {params.extra} -o {params.output_prefix} {input.ref} {input.bam} > {log}.out 2> {log}.err
        """

rule methyldackel_extract_mergecontext:
    """
    Extracts methylation calls with merged contexts for DSS.
    """
    input:
        bam=os.path.join(SORTED_FILTERED_BAM_DIR, "{sample}_pe.filtered.sorted.bam"),
        ref=REF_GENOME,
        options=os.path.join(MBIAS_DIR, "{sample}.options.txt")
    output:
        bedgraph=os.path.join(METHYLDACKEL_MERGECONTEXT_DIR, "{sample}_CpG.bedGraph")
    params:
        output_prefix=lambda wildcards: os.path.join(METHYLDACKEL_MERGECONTEXT_DIR, wildcards.sample),
        extra_opts=config.get("methyldackel_extract", {}).get("extra_opts", "--minOppositeDepth 10 --maxVariantFrac 0.5"),
        extra=METHYLDACKEL_EXTRACT_EXTRA_ARGS
    threads: 1
    log:
        os.path.join("logs", config["pipeline"], "methyldackel_extract_mergecontext", "{sample}.log")
    shell:
        """
        options=$(cat {input.options})
        pixi run MethylDackel extract --mergeContext ${{options}} {params.extra_opts} {params.extra} -o {params.output_prefix} {input.ref} {input.bam} > {log}.out 2> {log}.err
        """
