PREPARED_READS = [os.path.join(OUT, "primers", f"{s}_{m}.fastq.gz") for s in SAMPLES for m in ("R1", "R2") if SAMPLES_INFO[s].get(m)]
PREPARED_STATS = expand(os.path.join(OUT, "primers", "{sample}.reads.json"), sample=SAMPLES)

rule denoise:
    input:
        config=CONFIG_JSON,
        reads=PREPARED_READS,
        stats=PREPARED_STATS,
        script="scripts/amplicon_dada2.R"
    output:
        counts=os.path.join(OUT, "asv_counts.tsv"),
        sequences=os.path.join(OUT, "asv_sequences.fasta"),
        tracking=os.path.join(OUT, "read_tracking.tsv"),
        rds=os.path.join(OUT, "dada2.rds"),
        metadata=os.path.join(OUT, "run_metadata.json"),
        multiqc=os.path.join(OUT, "amplicon_read_retention_mqc.tsv"),
        filtered=directory(os.path.join(OUT, "filtered")),
        errors=directory(os.path.join(OUT, "errors"))
    threads: THREADS
    log:
        os.path.join(LOG_DIR, "dada2.log")
    shell:
        "pixi run Rscript {input.script:q} {input.config:q} {threads} > {log:q}.out 2> {log:q}.err"


rule taxonomy:
    input:
        config=CONFIG_JSON,
        rds=os.path.join(OUT, "dada2.rds"),
        reference=[config["taxonomy"]["training_set"]] if config["taxonomy"]["enabled"] else [],
        script="scripts/amplicon_taxonomy.R"
    output:
        table=os.path.join(OUT, "taxonomy.tsv"),
        metadata=os.path.join(OUT, "taxonomy_metadata.json")
    threads: THREADS
    log:
        os.path.join(LOG_DIR, "taxonomy.log")
    shell:
        "pixi run Rscript {input.script:q} {input.config:q} {threads} > {log:q}.out 2> {log:q}.err"


rule multiqc:
    input:
        raw=expand(os.path.join(OUT, "raw_qc", "{sample}"), sample=SAMPLES),
        primer=expand(os.path.join(OUT, "primer_qc", "{sample}"), sample=SAMPLES),
        retention=os.path.join(OUT, "amplicon_read_retention_mqc.tsv"),
        stats=PREPARED_STATS,
        cutadapt=expand(os.path.join(OUT, "primers", "{sample}.cutadapt.txt"), sample=SAMPLES),
        cutadapt_json=expand(os.path.join(OUT, "primers", "{sample}.cutadapt.json"), sample=SAMPLES),
        taxonomy=os.path.join(OUT, "taxonomy.tsv")
    output:
        report=os.path.join(OUT, "multiqc", "multiqc_report.html")
    params:
        scans=lambda w, input: [*input.raw, *input.primer, *input.cutadapt, input.retention],
        outdir=lambda w, output: os.path.dirname(output.report)
    log:
        os.path.join(LOG_DIR, "multiqc.log")
    shell:
        "pixi run multiqc {params.scans:q} --outdir {params.outdir:q} --filename multiqc_report.html --force > {log:q}.out 2> {log:q}.err"
