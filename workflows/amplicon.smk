import os
import sys

sys.path.insert(0, os.path.abspath("src"))
from amplicon_config import normalize

config = normalize(config)
SAMPLES_INFO = config["samples_info"]
SAMPLES = list(SAMPLES_INFO)
OUT = config["output_dir"]
LOG_DIR = config["log_dir"]
THREADS = config["threads"]
CONFIG_JSON = os.path.join(OUT, "resolved_config.json")

include: "rules/amplicon/reads.smk"
include: "rules/amplicon/dada2.smk"

rule all:
    input:
        os.path.join(OUT, "asv_counts.tsv"),
        os.path.join(OUT, "asv_sequences.fasta"),
        os.path.join(OUT, "read_tracking.tsv"),
        os.path.join(OUT, "taxonomy.tsv"),
        os.path.join(OUT, "taxonomy_metadata.json"),
        os.path.join(OUT, "run_metadata.json"),
        os.path.join(OUT, "multiqc", "multiqc_report.html")
