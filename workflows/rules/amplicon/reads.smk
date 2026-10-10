import json
import re

PE_PATTERN = "|".join(re.escape(s) for s in SAMPLES if SAMPLES_INFO[s].get("R2")) or r"(?!)"
SE_PATTERN = "|".join(re.escape(s) for s in SAMPLES if not SAMPLES_INFO[s].get("R2")) or r"(?!)"

rule resolved_config:
    output:
        config=CONFIG_JSON
    params:
        content=json.dumps(config, sort_keys=True, indent=2)
    log:
        os.path.join(LOG_DIR, "resolved_config.log")
    run:
        with open(output.config, "w") as handle:
            handle.write(f"{params.content}\n")


rule remove_primers_se:
    input:
        config=CONFIG_JSON,
        reads=lambda w: [SAMPLES_INFO[w.sample][m] for m in ("R1", "R2") if SAMPLES_INFO[w.sample].get(m)],
        script="scripts/amplicon_reads.py",
        validator="scripts/validate_fastq_pairs.py"
    output:
        reads=os.path.join(OUT, "primers", "{sample}_R1.fastq.gz"),
        stats=os.path.join(OUT, "primers", "{sample}.reads.json"),
        cutadapt_report=os.path.join(OUT, "primers", "{sample}.cutadapt.txt")
    wildcard_constraints:
        sample=SE_PATTERN
    params:
        sample=lambda w: w.sample
    threads: THREADS
    log:
        os.path.join(LOG_DIR, "primers", "{sample}.log")
    shell:
        "pixi run python {input.script:q} prepare --config {input.config:q} --sample {params.sample:q} --threads {threads} > {log:q}.out 2> {log:q}.err"


use rule remove_primers_se as remove_primers_pe with:
    output:
        reads=os.path.join(OUT, "primers", "{sample}_R1.fastq.gz"),
        stats=os.path.join(OUT, "primers", "{sample}.reads.json"),
        cutadapt_report=os.path.join(OUT, "primers", "{sample}.cutadapt.txt"),
        r2=os.path.join(OUT, "primers", "{sample}_R2.fastq.gz")
    wildcard_constraints:
        sample=PE_PATTERN


rule fastqc:
    input:
        config=CONFIG_JSON,
        prepared=lambda w: [os.path.join(OUT, "primers", f"{w.sample}_{m}.fastq.gz") for m in ("R1", "R2") if SAMPLES_INFO[w.sample].get(m)],
        raw=lambda w: [SAMPLES_INFO[w.sample][m] for m in ("R1", "R2") if SAMPLES_INFO[w.sample].get(m)],
        script="scripts/amplicon_reads.py",
        validator="scripts/validate_fastq_pairs.py"
    output:
        raw=directory(os.path.join(OUT, "raw_qc", "{sample}")),
        prepared=directory(os.path.join(OUT, "primer_qc", "{sample}"))
    params:
        sample=lambda w: w.sample
    threads: min(2, THREADS)
    log:
        os.path.join(LOG_DIR, "fastqc", "{sample}.log")
    shell:
        "pixi run python {input.script:q} qc --config {input.config:q} --sample {params.sample:q} --threads {threads} > {log:q}.out 2> {log:q}.err"
