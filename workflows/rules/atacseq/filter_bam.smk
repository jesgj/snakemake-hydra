import os

# --- CONFIGURATION ---
ALIGNMENT_DIR = config["alignment_dir"]
FILTERED_BAM_DIR = config["filtered_bam_dir"]
SAMBAMBA_MARKDUP_EXTRA_ARGS = config.get("sambamba", {}).get("markdup_extra_args", "-r")
SAMBAMBA_VIEW_EXTRA_ARGS = config.get(
    "sambamba",
    {},
).get(
    "view_extra_args",
    "not unmapped and proper_pair and mapping_quality >= 10 and not (ref_name =~ /^(MT|chrM)/)",
)


rule sambamba_filter_dedup_sort:
    """
    Removes duplicates, filters reads, and sorts BAM files for ATAC-seq.
    """
    input:
        bam = os.path.join(ALIGNMENT_DIR, "{sample}_pe.sorted.bam")
    output:
        filtered_sorted_bam = os.path.join(FILTERED_BAM_DIR, "{sample}_pe.filtered.sorted.bam")
    params:
        markdup_extra = SAMBAMBA_MARKDUP_EXTRA_ARGS,
        view_extra = SAMBAMBA_VIEW_EXTRA_ARGS
    threads: 4
    log:
        os.path.join("logs", config["pipeline"], "sambamba_filter", "{sample}_pe.log")
    shell:
        """
        pixi run bash -c "
          sambamba markdup -t {threads} {params.markdup_extra} {input.bam} /dev/stdout |
          sambamba view -t {threads} -f bam -F '{params.view_extra}' /dev/stdin |
          sambamba sort -t {threads} -o {output.filtered_sorted_bam} /dev/stdin
        " > {log}.out 2> {log}.err
        """
