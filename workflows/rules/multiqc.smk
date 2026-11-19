# workflows/rules/multiqc.smk
import os

# This rule assumes the following config variables are set by the parent workflow:
# - pipeline_name: e.g., "chip_cr", "rnaseq", "wgbs", "atacseq"
# - multiqc_input_files: Exact files consumed by MultiQC. Each path must also be
#   a declared output or log of another rule so report rebuilds track this cohort.
# - multiqc_results_dir: The top-level results directory for the pipeline.

_MULTIQC_RESULTS_DIR = config["multiqc_results_dir"]
_PIPELINE_NAME = config["pipeline_name"]
_LOG_PIPELINE = config.get("pipeline", _PIPELINE_NAME)
_LOG_DIR = config.get("log_dir", os.path.join("logs", _LOG_PIPELINE))
_MULTIQC_SOURCES = list(dict.fromkeys(config["multiqc_input_files"]))

if not _MULTIQC_SOURCES:
    raise ValueError("MultiQC requires at least one exact source in multiqc_input_files.")


rule multiqc:
    input:
        _MULTIQC_SOURCES
    output:
        report = os.path.join(_MULTIQC_RESULTS_DIR, "multiqc_report.html"),
        data = directory(os.path.join(_MULTIQC_RESULTS_DIR, "multiqc_report_data"))
    params:
        output_dir = _MULTIQC_RESULTS_DIR,
        filename = "multiqc_report.html"
    log:
        run = os.path.join(_LOG_DIR, "multiqc", "multiqc.log"),
        sources = os.path.join(_LOG_DIR, "multiqc", "multiqc.sources.txt")
    shell:
        """
        printf '%s\\n' {input:q} > {log.sources:q}
        pixi run multiqc --file-list {log.sources:q} --data-dir \
            -o {params.output_dir:q} --filename {params.filename:q} --force \
            > {log.run:q}.out 2> {log.run:q}.err
        """
