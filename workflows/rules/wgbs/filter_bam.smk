# workflows/rules/filter_bam.smk
import os

# --- CONFIGURATION ---
DEDUP_DIR = config["dedup_dir"]
FILTERED_BAM_DIR = config["filtered_bam_dir"]
SORTED_FILTERED_BAM_DIR = config["sorted_filtered_bam_dir"]
SAMBAMBA_EXTRA_ARGS = config.get("sambamba", {}).get("extra_args", "")

# --- RULES ---

rule sambamba_filter:
    """
    Filters a BAM file using sambamba view.
    """
    input:
        bam=os.path.join(DEDUP_DIR, "{sample}_pe.deduplicated.bam")
    output:
        filtered_bam=os.path.join(FILTERED_BAM_DIR, "{sample}_pe.filtered.bam")
    params:
        extra=SAMBAMBA_EXTRA_ARGS
    threads: 8
    log:
        os.path.join("logs", config["pipeline"], "sambamba_filter", "{sample}.log")
    shell:
        """
        pixi run sambamba view -t {threads} -f bam -h -F "{params.extra}" {input.bam} -o {output.filtered_bam} 2> {log}
        """

rule sambamba_sort:
    """
    Sorts a filtered BAM file using sambamba sort.
    """
    input:
        bam=os.path.join(FILTERED_BAM_DIR, "{sample}_pe.filtered.bam")
    output:
        sorted_bam=os.path.join(SORTED_FILTERED_BAM_DIR, "{sample}_pe.filtered.sorted.bam")
    threads: 8
    log:
        os.path.join("logs", config["pipeline"], "sambamba_sort", "{sample}.log")
    shell:
        """
        pixi run sambamba sort -t {threads} -o {output.sorted_bam} {input.bam} 2> {log}
        """

