# workflows/rules/chip_cr/filter_bam.smk
import os

# --- CONFIGURATION ---
ALIGNMENT_DIR = config["alignment_dir"]
FILTERED_BAM_DIR = config["filtered_bam_dir"]
SAMBAMBA_MARKDUP_EXTRA_ARGS = config.get("sambamba", {}).get("markdup_extra_args", "")
SAMBAMBA_VIEW_EXTRA_ARGS = config.get("sambamba", {}).get("view_extra_args", "")

# --- HELPER FUNCTION ---
def get_aligned_bam(wildcards):
    """
    Returns the path to the aligned BAM file for a sample,
    using the {read_type} wildcard ('pe' or 'se').
    """
    return os.path.join(ALIGNMENT_DIR, f"{wildcards.sample}_{wildcards.read_type}.sorted.bam")

# --- RULE ---

rule sambamba_filter_dedup_sort:
    """
    Removes duplicates, applies filters, and sorts the BAM file.
    """
    input:
        bam = get_aligned_bam
    output:
        filtered_sorted_bam = os.path.join(FILTERED_BAM_DIR, "{sample}_{read_type}.filtered.sorted.bam")
    params:
        markdup_extra = SAMBAMBA_MARKDUP_EXTRA_ARGS,
        view_extra = SAMBAMBA_VIEW_EXTRA_ARGS
    threads: 8
    log:
        os.path.join("logs", config["pipeline"], "sambamba_filter", "{sample}_{read_type}.log")
    shell:
        """
        pixi run bash -c "
          sambamba markdup -t {threads} {params.markdup_extra} {input.bam} /dev/stdout |
          sambamba view -t {threads} -f bam -F '{params.view_extra}' /dev/stdin |
          sambamba sort -t {threads} -o {output.filtered_sorted_bam} /dev/stdin
        " 2> {log}
        """
