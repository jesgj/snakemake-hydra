# workflows/rules/rnaseq/pseudoalignment.smk
import os

# --- CONFIGURATION ---
TRIMMED_DIR = config["trimmed_dir"]
KALLISTO_INDEX = config["kallisto_index"]
KALLISTO_OUTPUT_DIR = config["kallisto_output_dir"]
TRANSCRIPTOME_FASTA = config["transcriptome_fasta"]
KALLISTO_EXTRA_ARGS = config.get("kallisto", {}).get("extra_args", "")

# --- RULES ---

rule kallisto_index:
    """
    Builds a kallisto index from a transcriptome FASTA file.
    """
    input:
        fasta=TRANSCRIPTOME_FASTA
    output:
        index=KALLISTO_INDEX
    threads: 8
    log:
        os.path.join("logs", config["pipeline"], "kallisto_index", "kallisto_index.log")
    shell:
        "pixi run kallisto index -i {output.index} {input.fasta} > {log}.out 2> {log}.err"

rule kallisto_quant:
    """
    Runs kallisto quant for pseudo-alignment and quantification.
    """
    input:
        index=KALLISTO_INDEX,
        r1=os.path.join(TRIMMED_DIR, "{sample}_R1.trimmed.fq.gz"),
        r2=os.path.join(TRIMMED_DIR, "{sample}_R2.trimmed.fq.gz")
    output:
        # Kallisto outputs a directory, so we target a specific file within it
        tsv=os.path.join(KALLISTO_OUTPUT_DIR, "{sample}", "abundance.tsv")
    params:
        extra=KALLISTO_EXTRA_ARGS,
        output_dir=os.path.join(KALLISTO_OUTPUT_DIR, "{sample}")
    threads: 8
    log:
        os.path.join("logs", config["pipeline"], "kallisto_quant", "{sample}.log")
    shell:
        """
        pixi run kallisto quant -i {input.index} -o {params.output_dir} \
        -t {threads} {params.extra} {input.r1} {input.r2} > {log}.out 2> {log}.err
        """
