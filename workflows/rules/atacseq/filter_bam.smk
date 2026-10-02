import os

# --- CONFIGURATION ---
MARKED_BAM_DIR = config["marked_bam_dir"]
FILTERED_BAM_DIR = config["filtered_bam_dir"]
SAMBAMBA_CONFIG = config.get("sambamba", {})
SAMBAMBA_FILTER_THREADS = max(1, int(SAMBAMBA_CONFIG.get("filter_threads", 2)))
SAMBAMBA_SORT_THREADS = max(1, int(SAMBAMBA_CONFIG.get("sort_threads", 2)))
SAMBAMBA_VIEW_EXTRA_ARGS = config.get(
    "sambamba",
    {},
).get(
    "view_extra_args",
    "not unmapped and proper_pair and mapping_quality >= 10 and not (ref_name =~ /^(MT|chrM)/)",
)

if SAMBAMBA_VIEW_EXTRA_ARGS:
    SAMBAMBA_FILTER_EXPRESSION = f"not duplicate and ({SAMBAMBA_VIEW_EXTRA_ARGS})"
else:
    SAMBAMBA_FILTER_EXPRESSION = "not duplicate"


rule sambamba_filter_dedup_sort:
    """
    Removes duplicate-marked reads, filters reads, and sorts BAM files for ATAC-seq.
    """
    input:
        bam = os.path.join(MARKED_BAM_DIR, "{sample}_pe.markdup.bam")
    output:
        filtered_sorted_bam = os.path.join(FILTERED_BAM_DIR, "{sample}_pe.filtered.sorted.bam")
    params:
        filter_expression = SAMBAMBA_FILTER_EXPRESSION,
        filter_threads = lambda wildcards, threads: max(1, min(SAMBAMBA_FILTER_THREADS, threads // 2)),
        sort_threads = lambda wildcards, threads: max(1, min(SAMBAMBA_SORT_THREADS, threads - max(1, min(SAMBAMBA_FILTER_THREADS, threads // 2))))
    threads: SAMBAMBA_FILTER_THREADS + SAMBAMBA_SORT_THREADS
    log:
        os.path.join("logs", config["pipeline"], "sambamba_filter", "{sample}_pe.log")
    shell:
        """
        (
          if [ {threads} -eq 1 ]; then
            TMP_BAM={output.filtered_sorted_bam:q}.filter.tmp.bam
            trap 'rm -f "$TMP_BAM"' EXIT
            pixi run sambamba view -t 1 -f bam -F {params.filter_expression:q} -o "$TMP_BAM" {input.bam:q}
            pixi run sambamba sort -t 1 -o {output.filtered_sorted_bam:q} "$TMP_BAM"
          else
            pixi run sambamba view -t {params.filter_threads} -f bam -F {params.filter_expression:q} {input.bam:q} | \
            pixi run sambamba sort -t {params.sort_threads} -o {output.filtered_sorted_bam:q} /dev/stdin
          fi
        ) > {log}.out 2> {log}.err
        """
