# Shared FASTQ structure and mate-synchronization rules.
#
# Parent workflows override the placeholder input/output/log values with
# `use rule ... as ... with:`. The validator itself is a declared input so an
# edited script invalidates existing validation reports.

FASTQ_VALIDATION_SCRIPT = "scripts/validate_fastq_pairs.py"

rule validate_fastq_pair_generic:
    """Validate one paired-end FASTQ pair before trimming."""
    input:
        r1="path/to/R1.fastq.gz",
        r2="path/to/R2.fastq.gz",
        validator=FASTQ_VALIDATION_SCRIPT
    output:
        report="path/to/report.fastq_pair_validation.tsv"
    params:
        scope="--all"
    threads: 1
    log:
        "logs/fastq_validation/pair.log"
    shell:
        """
        set -euo pipefail
        pixi run python3 {input.validator:q} \
            --r1 {input.r1:q} \
            --r2 {input.r2:q} \
            {params.scope} > {output.report:q} 2> {log:q}.err
        """

rule validate_fastq_single_generic:
    """Validate one single-end FASTQ with the shared reader/audit."""
    input:
        r1="path/to/R1.fastq.gz",
        validator=FASTQ_VALIDATION_SCRIPT
    output:
        report="path/to/report.fastq_single_validation.tsv"
    params:
        scope="--all"
    threads: 1
    log:
        "logs/fastq_validation/single.log"
    shell:
        """
        set -euo pipefail
        pixi run python3 {input.validator:q} \
            --r1 {input.r1:q} \
            {params.scope} > {output.report:q} 2> {log:q}.err
        """
