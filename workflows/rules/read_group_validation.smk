# Shared per-BAM and merged-BAM read-group rules.
#
# Parent workflows override the placeholder input/output/log values with
# `use rule ... as ... with:`. Merged validation is only used by pipelines that
# actually merge multiple run or lane BAMs.

READ_GROUP_VALIDATOR = "src/validate_bam_read_group.py"
MERGED_READ_GROUP_VALIDATOR = "src/validate_merged_bam_read_groups.py"

rule validate_bam_read_group_generic:
    """Validate one BAM's complete @RG header and every alignment RG tag."""
    input:
        bam="path/to/input.bam",
        validator=READ_GROUP_VALIDATOR
    output:
        report="path/to/report.read_group_validation.tsv"
    params:
        rg_id="expected-id",
        rg_sm="expected-sample",
        rg_lb="expected-library",
        rg_pl="ILLUMINA",
        rg_pu="expected-platform-unit"
    threads: 1
    log:
        "logs/read_group_validation/bam.log"
    shell:
        """
        set -euo pipefail
        pixi run samtools view -h {input.bam:q} | \
        pixi run python3 {input.validator:q} \
            --expected-id {params.rg_id:q} \
            --expected-sample {params.rg_sm:q} \
            --expected-library {params.rg_lb:q} \
            --expected-platform {params.rg_pl:q} \
            --expected-platform-unit {params.rg_pu:q} \
            > {output.report:q} 2> {log:q}.err
        """

rule validate_merged_bam_read_groups_generic:
    """Validate expected @RG headers and tags in a merged BAM."""
    input:
        bam="path/to/merged.bam",
        validator=MERGED_READ_GROUP_VALIDATOR
    output:
        report="path/to/report.merged.read_group_validation.tsv"
    params:
        expected=""
    threads: 1
    log:
        "logs/read_group_validation/merged.log"
    shell:
        """
        set -euo pipefail
        pixi run samtools view -h {input.bam:q} | \
        pixi run python3 {input.validator:q} {params.expected} \
            > {output.report:q} 2> {log:q}.err
        """
