import os


SAMPLES_INFO = config["samples_info"]
FASTQ_PAIR_VALIDATION_DIR = config["fastq_pair_validation_dir"]
FASTQ_VALIDATION_CONFIG = config.get("fastq_validation", {})
FASTQ_VALIDATION_ALL_RECORDS = FASTQ_VALIDATION_CONFIG.get("all_records", True)
FASTQ_VALIDATION_RECORDS = int(FASTQ_VALIDATION_CONFIG.get("records", 100))

if not isinstance(FASTQ_VALIDATION_ALL_RECORDS, bool):
    raise ValueError("ATAC-seq fastq_validation.all_records must be true or false.")
if FASTQ_VALIDATION_RECORDS < 1:
    raise ValueError("ATAC-seq fastq_validation.records must be a positive integer.")


rule validate_fastq_pair:
    """
    Validates normalized cluster IDs, mate labels, FASTQ structure, and counts.
    """
    input:
        r1=lambda wildcards: SAMPLES_INFO[wildcards.sample]["runs"][wildcards.run]["R1"],
        r2=lambda wildcards: SAMPLES_INFO[wildcards.sample]["runs"][wildcards.run]["R2"]
    output:
        report=os.path.join(
            FASTQ_PAIR_VALIDATION_DIR,
            "{sample}__{run}.fastq_pair_validation.tsv",
        )
    params:
        scope="--all" if FASTQ_VALIDATION_ALL_RECORDS else f"--records {FASTQ_VALIDATION_RECORDS}"
    log:
        os.path.join(
            "logs", config["pipeline"], "fastq_pair_validation", "{sample}__{run}.log"
        )
    shell:
        """
        pixi run python3 scripts/validate_fastq_pairs.py \
            --r1 {input.r1:q} \
            --r2 {input.r2:q} \
            {params.scope} > {output.report:q} 2> {log:q}.err
        """
