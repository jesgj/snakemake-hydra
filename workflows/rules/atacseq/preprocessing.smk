import os


SAMPLES_INFO = config["samples_info"]
QC_DIR = config["qc_dir"]
TRIMMED_DIR = config["trimmed_dir"]
QC_TRIMMED_DIR = config["qc_trimmed_dir"]
FASTQ_PAIR_VALIDATION_DIR = config["fastq_pair_validation_dir"]
FASTP_EXTRA_ARGS = config.get("fastp", {}).get("extra_args", "")


def _atac_run_info(wildcards):
    return SAMPLES_INFO[wildcards.sample]["runs"][wildcards.run]


rule fastqc_raw_run:
    """
    Runs FastQC on one original ATAC-seq run-level FASTQ pair.
    """
    input:
        r1=lambda wildcards: _atac_run_info(wildcards)["R1"],
        r2=lambda wildcards: _atac_run_info(wildcards)["R2"]
    output:
        html_r1=os.path.join(QC_DIR, "{sample}__{run}_R1_raw_fastqc.html"),
        html_r2=os.path.join(QC_DIR, "{sample}__{run}_R2_raw_fastqc.html"),
        zip_r1=os.path.join(QC_DIR, "{sample}__{run}_R1_raw_fastqc.zip"),
        zip_r2=os.path.join(QC_DIR, "{sample}__{run}_R2_raw_fastqc.zip")
    params:
        outdir=QC_DIR
    threads: 2
    log:
        os.path.join("logs", config["pipeline"], "fastqc_raw", "{sample}__{run}.log")
    shell:
        """
        set -euo pipefail
        pixi run fastqc -o {params.outdir:q} -t {threads} {input.r1:q} {input.r2:q} \
            > {log:q}.out 2> {log:q}.err
        R1_BASE=$(basename {input.r1:q})
        R1_STEM=${{R1_BASE%%.fastq.gz}}
        R1_STEM=${{R1_STEM%%.fq.gz}}
        R2_BASE=$(basename {input.r2:q})
        R2_STEM=${{R2_BASE%%.fastq.gz}}
        R2_STEM=${{R2_STEM%%.fq.gz}}
        mv {params.outdir:q}/${{R1_STEM}}_fastqc.html {output.html_r1:q}
        mv {params.outdir:q}/${{R1_STEM}}_fastqc.zip {output.zip_r1:q}
        mv {params.outdir:q}/${{R2_STEM}}_fastqc.html {output.html_r2:q}
        mv {params.outdir:q}/${{R2_STEM}}_fastqc.zip {output.zip_r2:q}
        """


rule fastp_trim_run:
    """
    Trims one validated ATAC-seq run-level FASTQ pair.
    """
    input:
        r1=lambda wildcards: _atac_run_info(wildcards)["R1"],
        r2=lambda wildcards: _atac_run_info(wildcards)["R2"],
        pair_validation=os.path.join(
            FASTQ_PAIR_VALIDATION_DIR, "{sample}__{run}.fastq_pair_validation.tsv"
        )
    output:
        trimmed_r1=os.path.join(TRIMMED_DIR, "{sample}__{run}_R1.trimmed.fq.gz"),
        trimmed_r2=os.path.join(TRIMMED_DIR, "{sample}__{run}_R2.trimmed.fq.gz"),
        html=os.path.join(TRIMMED_DIR, "{sample}__{run}.fastp.html"),
        json=os.path.join(TRIMMED_DIR, "{sample}__{run}.fastp.json")
    params:
        extra=FASTP_EXTRA_ARGS
    threads: 4
    log:
        os.path.join("logs", config["pipeline"], "fastp", "{sample}__{run}.log")
    shell:
        """
        pixi run fastp -i {input.r1:q} -I {input.r2:q} \
            -o {output.trimmed_r1:q} -O {output.trimmed_r2:q} \
            -h {output.html:q} -j {output.json:q} \
            {params.extra} -w {threads} > {log:q}.out 2> {log:q}.err
        """


rule fastqc_trimmed_run:
    """
    Runs FastQC on one trimmed ATAC-seq run-level FASTQ pair.
    """
    input:
        r1=os.path.join(TRIMMED_DIR, "{sample}__{run}_R1.trimmed.fq.gz"),
        r2=os.path.join(TRIMMED_DIR, "{sample}__{run}_R2.trimmed.fq.gz")
    output:
        html_r1=os.path.join(QC_TRIMMED_DIR, "{sample}__{run}_R1_trimmed_fastqc.html"),
        html_r2=os.path.join(QC_TRIMMED_DIR, "{sample}__{run}_R2_trimmed_fastqc.html"),
        zip_r1=os.path.join(QC_TRIMMED_DIR, "{sample}__{run}_R1_trimmed_fastqc.zip"),
        zip_r2=os.path.join(QC_TRIMMED_DIR, "{sample}__{run}_R2_trimmed_fastqc.zip")
    params:
        outdir=QC_TRIMMED_DIR
    threads: 2
    log:
        os.path.join("logs", config["pipeline"], "fastqc_trimmed", "{sample}__{run}.log")
    shell:
        """
        set -euo pipefail
        pixi run fastqc -o {params.outdir:q} -t {threads} {input.r1:q} {input.r2:q} \
            > {log:q}.out 2> {log:q}.err
        R1_BASE=$(basename {input.r1:q})
        R1_STEM=${{R1_BASE%%.fq.gz}}
        R2_BASE=$(basename {input.r2:q})
        R2_STEM=${{R2_BASE%%.fq.gz}}
        mv {params.outdir:q}/${{R1_STEM}}_fastqc.html {output.html_r1:q}
        mv {params.outdir:q}/${{R1_STEM}}_fastqc.zip {output.zip_r1:q}
        mv {params.outdir:q}/${{R2_STEM}}_fastqc.html {output.html_r2:q}
        mv {params.outdir:q}/${{R2_STEM}}_fastqc.zip {output.zip_r2:q}
        """
