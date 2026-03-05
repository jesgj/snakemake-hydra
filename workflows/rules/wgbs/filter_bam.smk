# workflows/rules/filter_bam.smk
import os

# --- CONFIGURATION ---
DEDUP_DIR = config["dedup_dir"]
FILTERED_BAM_DIR = config["filtered_bam_dir"]
SORTED_FILTERED_BAM_DIR = config["sorted_filtered_bam_dir"]
SAMBAMBA_EXTRA_ARGS = config.get("sambamba", {}).get("extra_args", "")
SAMBAMBA_FILTER_THREADS = int(config.get("sambamba", {}).get("filter_threads", 8))
SAMBAMBA_SORT_THREADS = int(config.get("sambamba", {}).get("sort_threads", 8))

if SAMBAMBA_FILTER_THREADS < 1 or SAMBAMBA_SORT_THREADS < 1:
    raise ValueError("WGBS sambamba.filter_threads and sambamba.sort_threads must be positive integers.")

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
    threads: SAMBAMBA_FILTER_THREADS
    log:
        os.path.join("logs", config["pipeline"], "sambamba_filter", "{sample}.log")
    shell:
        """
        set -euo pipefail
        pixi run sambamba view -t {threads} -f bam -h -F "{params.extra}" {input.bam} -o {output.filtered_bam} > {log}.out 2> {log}.err
        """

rule sambamba_sort:
    """
    Sorts a filtered BAM file using sambamba sort.
    """
    input:
        bam=os.path.join(FILTERED_BAM_DIR, "{sample}_pe.filtered.bam")
    output:
        sorted_bam=os.path.join(SORTED_FILTERED_BAM_DIR, "{sample}_pe.filtered.sorted.bam")
    threads: SAMBAMBA_SORT_THREADS
    log:
        os.path.join("logs", config["pipeline"], "sambamba_sort", "{sample}.log")
    shell:
        """
        set -euo pipefail
        pixi run sambamba sort -t {threads} -o {output.sorted_bam} {input.bam} > {log}.out 2> {log}.err
        """

rule samtools_index_filtered_bam:
    """
    Builds an index for a sorted, filtered WGBS BAM file.
    """
    input:
        bam=os.path.join(SORTED_FILTERED_BAM_DIR, "{sample}_pe.filtered.sorted.bam")
    output:
        bai=os.path.join(SORTED_FILTERED_BAM_DIR, "{sample}_pe.filtered.sorted.bam.bai")
    threads: 4
    log:
        os.path.join("logs", config["pipeline"], "samtools_index_filtered", "{sample}.log")
    shell:
        """
        pixi run samtools index -@ {threads} {input.bam} {output.bai} > {log}.out 2> {log}.err
        """
