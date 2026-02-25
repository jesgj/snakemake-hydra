# workflows/rules/qc.smk
import os

# This module now expects SAMPLES_INFO to be in the config
SAMPLES_INFO = config['samples_info']
QC_DIR = config["qc_dir"]

rule fastqc_raw_pe:
    """
    Runs FastQC on raw paired-end reads for a single sample.
    """
    input:
        r1 = lambda wildcards: SAMPLES_INFO[wildcards.sample]['R1'],
        r2 = lambda wildcards: SAMPLES_INFO[wildcards.sample]['R2']
    output:
        html_r1 = os.path.join(QC_DIR, "{sample}_R1_raw_fastqc.html"),
        html_r2 = os.path.join(QC_DIR, "{sample}_R2_raw_fastqc.html"),
        zip_r1 = os.path.join(QC_DIR, "{sample}_R1_raw_fastqc.zip"),
        zip_r2 = os.path.join(QC_DIR, "{sample}_R2_raw_fastqc.zip")
    params:
        outdir = QC_DIR
    threads: 2
    log:
        os.path.join("logs", config["pipeline"], "fastqc_raw", "{sample}_pe.log")
    shell:
        """
        # 1. Run FastQC
        pixi run fastqc -o {params.outdir} -t {threads} {input.r1} {input.r2} > {log}.out 2> {log}.err

        # --- Rename logic for R1 ---
        # Get base filename
        R1_BASE=$(basename {input.r1})
        # Remove .fastq.gz or .fq.gz to get the 'stem'
        R1_STEM=${{R1_BASE%%.fastq.gz}}
        R1_STEM=${{R1_STEM%%.fq.gz}}
        
        # Move using the real name that FastQC generated
        mv {params.outdir}/${{R1_STEM}}_fastqc.html {output.html_r1}
        mv {params.outdir}/${{R1_STEM}}_fastqc.zip {output.zip_r1}

        # --- Rename logic for R2 ---
        R2_BASE=$(basename {input.r2})
        R2_STEM=${{R2_BASE%%.fastq.gz}}
        R2_STEM=${{R2_STEM%%.fq.gz}}

        mv {params.outdir}/${{R2_STEM}}_fastqc.html {output.html_r2}
        mv {params.outdir}/${{R2_STEM}}_fastqc.zip {output.zip_r2}
        """

rule fastqc_raw_se:
    """
    Runs FastQC on raw single-end reads for a single sample.
    """
    input:
        r1 = lambda wildcards: SAMPLES_INFO[wildcards.sample]['R1']
    output:
        html = os.path.join(QC_DIR, "{sample}_raw_fastqc.html"),
        zip = os.path.join(QC_DIR, "{sample}_raw_fastqc.zip")
    params:
        outdir = QC_DIR
    threads: 1
    log:
        os.path.join("logs", config["pipeline"], "fastqc_raw", "{sample}_se.log")
    shell:
        """
        # 1. Run FastQC
        pixi run fastqc -o {params.outdir} -t {threads} {input.r1} > {log}.out 2> {log}.err

        # --- Rename logic ---
        # Get base filename
        R1_BASE=$(basename {input.r1})
        # Remove .fastq.gz or .fq.gz to get the 'stem'
        R1_STEM=${{R1_BASE%%.fastq.gz}}
        R1_STEM=${{R1_STEM%%.fq.gz}}
        
        # Move using the real name that FastQC generated
        mv {params.outdir}/${{R1_STEM}}_fastqc.html {output.html}
        mv {params.outdir}/${{R1_STEM}}_fastqc.zip {output.zip}
        """
