import os

from utils import validate_sambamba_filter

# --- CONFIGURATION ---
MARKED_BAM_DIR = config["marked_bam_dir"]
FILTERED_BAM_DIR = config["filtered_bam_dir"]
SAMBAMBA_CONFIG = config.get("sambamba", {})
SAMBAMBA_FILTER_THREADS = max(1, int(SAMBAMBA_CONFIG.get("filter_threads", 2)))
SAMBAMBA_SORT_THREADS = max(1, int(SAMBAMBA_CONFIG.get("sort_threads", 2)))
SAMBAMBA_VIEW_EXTRA_ARGS = validate_sambamba_filter(
    SAMBAMBA_CONFIG.get(
        "view_extra_args",
        "not unmapped and proper_pair and mapping_quality >= 10 and not (ref_name =~ /^(MT|chrM)/)",
    ),
    "atacseq.sambamba.view_extra_args",
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
        filter_threads = SAMBAMBA_FILTER_THREADS,
        sort_threads = SAMBAMBA_SORT_THREADS
    threads: SAMBAMBA_FILTER_THREADS + SAMBAMBA_SORT_THREADS
    log:
        os.path.join("logs", config["pipeline"], "sambamba_filter", "{sample}_pe.log")
    shell:
        """
        set -euo pipefail
        (pixi run sambamba view -t {params.filter_threads} -f bam -F {params.filter_expression:q} {input.bam:q} | \
        pixi run sambamba sort -t {params.sort_threads} -o {output.filtered_sorted_bam:q} /dev/stdin) > {log:q}.out 2> {log:q}.err
        """
