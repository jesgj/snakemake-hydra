import os

# --- CONFIGURATION ---
MARKED_BAM_DIR = config["marked_bam_dir"]
FILTERED_BAM_DIR = config["filtered_bam_dir"]
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
        filter_expression = SAMBAMBA_FILTER_EXPRESSION
    threads: 4
    log:
        os.path.join("logs", config["pipeline"], "sambamba_filter", "{sample}_pe.log")
    shell:
        """
        (pixi run sambamba view -t {threads} -f bam -F '{params.filter_expression}' "{input.bam}" | \
        pixi run sambamba sort -t {threads} -o "{output.filtered_sorted_bam}" /dev/stdin) > {log}.out 2> {log}.err
        """
