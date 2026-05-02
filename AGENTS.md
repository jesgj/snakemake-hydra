# Agent Guidance for `snakemake-hydra`

## Source Of Truth
- This is a Snakemake + Pixi NGS pipeline, not a Python package; workflow behavior lives in `Snakefile`, `workflows/*.smk`, `workflows/rules/**/*.smk`, and `src/utils.py`.
- `Snakefile` reads `config/config.yaml`, builds config objects for all four modules, then imports only the selected pipeline with prefixed rule names like `rnaseq_*`, `wgbs_*`, `chip_cr_*`, or `atacseq_*`.
- Supported `config["pipeline"]` values are `rnaseq`, `wgbs`, `chip_cr`, and `atacseq`; focused rule/target runs must match the selected pipeline or the rule will not be imported.
- `pixi.toml` has dependencies but no tasks; use explicit `pixi run ...` commands.
- No CI, pre-commit, formatter, linter, typechecker, or unit-test config is present; do not invent repo-wide `pytest`, `ruff`, `black`, etc.
- `config/config.yaml` contains machine-specific absolute paths; preserve them unless the task is explicitly about configuration.

## Commands
- Install/update the managed environment: `pixi install`.
- Dry-run the selected pipeline: `pixi run snakemake -n --cores 1`.
- Run the selected pipeline: `pixi run snakemake --cores 8 --printshellcmds --latency-wait 60`.
- Re-run incomplete jobs: `pixi run snakemake --cores 8 --rerun-incomplete`.
- Dry-run one selected-pipeline rule path: `pixi run snakemake -n --cores 1 --until <prefixed_rule_name>`.
- Dry-run one concrete output: `pixi run snakemake -n --cores 1 <path/to/output>`.
- DAG only: `pixi run snakemake --dag --cores 1`; render with Graphviz using `pixi run snakemake --dag --cores 1 | dot -Tpng > dag.png`.

## Validation Strategy
- For config-only edits, a selected-pipeline dry-run is usually the best verification.
- For rule edits, prefer `--until <prefixed_rule_name>` first, then a concrete output if paths or filenames changed.
- Useful concrete outputs: `results/rnaseq/kallisto/sample/abundance.tsv`, `results/wgbs/methyldackel/sample_CpG.methylKit`, `results/atacseq/peaks/sample_peaks.narrowPeak`.
- When changing MultiQC dependencies, dry-run the final report at `<multiqc_results_dir>/multiqc_report.html`.
- Keep verification scoped to the active pipeline; switching `config["pipeline"]` is a config change, not a harmless test setup detail.

## Workflow Conventions
- Each top-level workflow reads config, validates required keys, normalizes values back into `config[...]`, builds output/MultiQC lists, then includes rule files.
- Included rules expect parent workflows to set keys such as `samples_info`, output directories, `multiqc_input_files`, and `multiqc_analysis_dirs`.
- Workflow modules import helpers via `sys.path.insert(0, os.path.abspath("src"))` followed by direct imports from `utils`.
- Add user-facing config keys to the owning workflow, `config/config.yaml`, and `README.md`; add notable fixes/features to `backlog.md`.
- Use `os.path.join(...)` for paths and preserve established suffixes such as `_pe`, `_se`, `_R1`, `_R2`, `.filtered.sorted.bam`, and `.markdup.metrics.txt`.
- Shell blocks should run tools through `pixi run` and redirect stdout/stderr in the existing `{log}.out` / `{log}.err` style.

## Samples And References
- `src/utils.py` falls back from `<pipeline>.samples_info` to FASTQ discovery in `raw_fastqs_dir`; discovery supports `_R1/_R2`, `_1/_2`, `.R1/.R2`, and single-end `.fastq.gz`/`.fq.gz` files.
- `rnaseq`, `wgbs`, and `atacseq` reject non-paired samples; `chip_cr` accepts PE and SE and infers `type` from presence of `R2` when omitted.
- ChIP/CUT&RUN input subtraction is name-driven: signal `<base>_repN` pairs with input `<base>_input_repN`; `<base>` may contain underscores.
- Reference values are required config keys but are only presence-checked at workflow parse time; actual files/permissions fail later in tool rules.
- The pipeline does not download references. RNA-seq can build `kallisto_index` and HISAT2 indexes, WGBS creates `Bisulfite_Genome/` next to `ref_genome`, and ChIP/ATAC create configured Bowtie2 index directories.

## Pipeline Gotchas
- RNA-seq gene body coverage is optional via `rnaseq.gene_body_coverage.enabled`; `refgene_bed` is required only when enabled, and its heatmap output exists only with at least 3 samples.
- RNA-seq, ChIP/CUT&RUN, and ATAC-seq correlation heatmaps are generated only when at least 2 BAMs are available; do not add unconditional MultiQC/report dependencies for those files.
- WGBS uses `wgbs.log_dir` when set; other pipelines default logs under `logs/<pipeline>`.
- WGBS Bismark thread allocation derives Bowtie2 threads from `wgbs.bismark.threads / wgbs.bismark.parallel`; keep `threads >= parallel`.
- ATAC-seq supports only `peak_caller: "macs3"` or `"genrich"`; Genrich uses a queryname-sorted intermediate but both paths expose the shared `.narrowPeak` output.
