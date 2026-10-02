# workflows/rules/chip_cr/filter_bam.smk
import os
import shlex

from utils import validate_sambamba_filter

# --- CONFIGURATION ---
ALIGNMENT_DIR = config["alignment_dir"]
FILTERED_BAM_DIR = config["filtered_bam_dir"]
SAMBAMBA_MARKDUP_EXTRA_ARGS = config.get("sambamba", {}).get("markdup_extra_args", "")
SAMBAMBA_VIEW_EXTRA_ARGS = validate_sambamba_filter(
    config.get("sambamba", {}).get("view_extra_args", ""),
    "chip_cr.sambamba.view_extra_args",
)

# --- HELPER FUNCTION ---
def get_aligned_bam(wildcards):
    """
    Returns the path to the aligned BAM file for a sample,
    using the {read_type} wildcard ('pe' or 'se').
    """
    return os.path.join(ALIGNMENT_DIR, f"{wildcards.sample}_{wildcards.read_type}.sorted.bam")

# --- RULE ---

rule sambamba_markdup:
    """Mark or remove duplicates using a regular BAM output."""
    input:
        bam = get_aligned_bam
    output:
        bam = temp(os.path.join(FILTERED_BAM_DIR, "{sample}_{read_type}.dedup.bam"))
    params:
        extra = SAMBAMBA_MARKDUP_EXTRA_ARGS
    threads: 8
    log:
        os.path.join("logs", config["pipeline"], "sambamba_markdup", "{sample}_{read_type}.log")
    shell:
        "pixi run sambamba markdup -t {threads} {params.extra} {input.bam:q} {output.bam:q} > {log}.out 2> {log}.err"


rule sambamba_filter_dedup_sort:
    """
    Removes duplicates, applies filters, and sorts the BAM file.
    """
    input:
        bam = os.path.join(FILTERED_BAM_DIR, "{sample}_{read_type}.dedup.bam")
    output:
        filtered_sorted_bam = os.path.join(FILTERED_BAM_DIR, "{sample}_{read_type}.filtered.sorted.bam")
    params:
        filter_args = "-F " + shlex.quote(SAMBAMBA_VIEW_EXTRA_ARGS) if SAMBAMBA_VIEW_EXTRA_ARGS else "",
        view_threads = lambda wildcards, threads: max(1, threads // 2),
        sort_threads = lambda wildcards, threads: max(1, threads - threads // 2)
    threads: 8
    log:
        os.path.join("logs", config["pipeline"], "sambamba_filter", "{sample}_{read_type}.log")
    shell:
        """
        set -euo pipefail
        (
          if [ {threads} -eq 1 ]; then
            # A one-core job runs the tools sequentially to respect its CPU budget.
            TMP_BAM={output.filtered_sorted_bam:q}.filter.tmp.bam
            trap 'rm -f "$TMP_BAM"' EXIT
            pixi run sambamba view -t 1 -f bam {params.filter_args} -o "$TMP_BAM" {input.bam:q}
            pixi run sambamba sort -t 1 -o {output.filtered_sorted_bam:q} "$TMP_BAM"
          else
            pixi run sambamba view -t {params.view_threads} -f bam {params.filter_args} {input.bam:q} | \
            pixi run sambamba sort -t {params.sort_threads} -o {output.filtered_sorted_bam:q} /dev/stdin
          fi
        ) > {log:q}.out 2> {log:q}.err
        """
