import os

ALIGNMENT_DIR = config["alignment_dir"]
GENE_BODY_COVERAGE_DIR = config["gene_body_coverage_dir"]
SAMPLES_INFO = config["samples_info"]

GENE_BODY_COVERAGE_CONFIG = config.get("gene_body_coverage", {})
REFGENE_BED = GENE_BODY_COVERAGE_CONFIG.get("refgene_bed")
MIN_MRNA_LENGTH = GENE_BODY_COVERAGE_CONFIG.get("minimum_length", 100)
OUTPUT_FORMAT = GENE_BODY_COVERAGE_CONFIG.get("format", "pdf")

if not REFGENE_BED:
    raise ValueError(
        "RNA-seq gene body coverage requires 'gene_body_coverage.refgene_bed' in the config."
    )

if MIN_MRNA_LENGTH < 100:
    raise ValueError(
        "RNA-seq gene body coverage requires 'gene_body_coverage.minimum_length' to be at least 100."
    )

if OUTPUT_FORMAT not in ["pdf", "png", "jpeg"]:
    raise ValueError(
        "RNA-seq gene body coverage requires 'gene_body_coverage.format' to be one of: pdf, png, jpeg."
    )


def get_aligned_bams(samples_info):
    return [os.path.join(ALIGNMENT_DIR, f"{sample}_pe.sorted.bam") for sample in samples_info]


def get_aligned_bais(samples_info):
    return [os.path.join(ALIGNMENT_DIR, f"{sample}_pe.sorted.bam.bai") for sample in samples_info]


GENE_BODY_COVERAGE_PREFIX = os.path.join(GENE_BODY_COVERAGE_DIR, "all_samples")
GENE_BODY_COVERAGE_OUTPUTS = [
    f"{GENE_BODY_COVERAGE_PREFIX}.geneBodyCoverage.txt",
    f"{GENE_BODY_COVERAGE_PREFIX}.geneBodyCoverage.curves.{OUTPUT_FORMAT}",
]

if len(SAMPLES_INFO) >= 3:
    GENE_BODY_COVERAGE_OUTPUTS.append(
        f"{GENE_BODY_COVERAGE_PREFIX}.geneBodyCoverage.heatMap.{OUTPUT_FORMAT}"
    )


rule gene_body_coverage:
    """
    Calculates RNA-seq read coverage across gene bodies with RSeQC.
    """
    input:
        bams=get_aligned_bams(SAMPLES_INFO),
        bais=get_aligned_bais(SAMPLES_INFO),
        refgene=REFGENE_BED
    output:
        GENE_BODY_COVERAGE_OUTPUTS
    params:
        bam_input=",".join(get_aligned_bams(SAMPLES_INFO)),
        min_mrna_length=MIN_MRNA_LENGTH,
        output_format=OUTPUT_FORMAT,
        prefix=GENE_BODY_COVERAGE_PREFIX
    threads: 1
    log:
        os.path.join("logs", config["pipeline"], "gene_body_coverage", "all_samples.log")
    shell:
        """
        pixi run geneBody_coverage.py \
            -i "{params.bam_input}" \
            -r "{input.refgene}" \
            -l {params.min_mrna_length} \
            -f {params.output_format} \
            -o "{params.prefix}" > {log}.out 2> {log}.err
        """
