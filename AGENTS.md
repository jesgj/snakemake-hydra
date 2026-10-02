# Agent Guidance for `snakemake-hydra`

## Source Of Truth

- This is a Snakemake + Pixi NGS pipeline, not a Python package; workflow behavior lives in `Snakefile`, `workflows/*.smk`, `workflows/rules/**/*.smk`, and `src/utils.py`.
- `Snakefile` reads `config/config.yaml`, builds config objects for all four modules, then imports only the selected pipeline with prefixed rule names like `rnaseq_*`, `wgbs_*`, `chip_cr_*`, or `atacseq_*`.
- Keep all four top-level pipeline config mappings present: `Snakefile` accesses each before selecting rules. Inside a module, `config` is that pipeline's subsection plus `pipeline`; use `config["samples_info"]`, not `config["rnaseq"]["samples_info"]`.
- Supported `config["pipeline"]` values are `rnaseq`, `wgbs`, `chip_cr`, and `atacseq`; focused rule/target runs must match the selected pipeline or the rule will not be imported.
- `pixi.toml` targets `linux-64`, has dependencies but no tasks, and is paired with committed `pixi.lock`; use explicit `pixi run ...` commands. Keep manifest and lockfile consistent when changing dependencies.
- No CI, pre-commit, formatter, linter, typechecker, or unit-test config is present; do not invent repo-wide `pytest`, `ruff`, `black`, etc. Snakemake's built-in `--lint` is available without a repository linter config.
- Inspect `git status --short` before editing and preserve existing user changes. Treat sample, reference, and output paths in `config/config.yaml` as user configuration; do not rewrite them to make validation pass.
- `config/config.local.yaml` is ignored by Git but is not automatically loaded by `Snakefile`.

## Where To Make Changes

- Pipeline setup, defaults, sample validation, output lists, and rule instantiation: `workflows/<pipeline>.smk`.
- Tool commands and their declared inputs/outputs: `workflows/rules/<pipeline>/*.smk`.
- Shared FASTQ QC/trimming: `workflows/rules/qc.smk` and `workflows/rules/trimming_and_qc.smk`.
- Shared BAM/deepTools QC and reporting: `workflows/rules/bam_qc.smk`, `workflows/rules/deeptools_qc.smk`, and `workflows/rules/multiqc.smk`. Check the including workflows when changing these shared interfaces.
- Sample discovery and path validation helpers: `src/utils.py`. The shell scripts under `src/` are Slurm launch examples with fixed resources/targets, not the workflow entrypoint.

## Commands

Run from the repository root; config paths and the `src` import depend on the working directory.

- Install the managed environment: `pixi install`.
- Dry-run the selected pipeline: `pixi run snakemake -n --cores 1`.
- Run the selected pipeline: `pixi run snakemake --cores 8 --printshellcmds --latency-wait 60`.
- Re-run incomplete jobs: `pixi run snakemake --cores 8 --rerun-incomplete`.
- Dry-run one selected-pipeline rule path: `pixi run snakemake -n --cores 1 --until <prefixed_rule_name>`.
- Dry-run one concrete output: `pixi run snakemake -n --cores 1 <path/to/output>`.
- DAG only: `pixi run snakemake --dag --cores 1`; render with Graphviz using `pixi run snakemake --dag --cores 1 | dot -Tpng > dag.png`. Graphviz `dot` is optional and is not declared in `pixi.toml`.

## Validation Strategy

- Documentation-only changes: verify claims against source and run `git diff --check`; a pipeline run is unnecessary.
- For config-only edits, a selected-pipeline dry-run is usually the best verification.
- For rule edits, prefer `--until <prefixed_rule_name>` first, then a concrete output if paths or filenames changed.
- For workflow/rule edits, also use `pixi run snakemake --lint` when the selected workflow can be parsed. Address findings introduced by the change; distinguish existing warnings and recommendations that do not fit the repository's Pixi environment strategy.
- Useful output templates: `<kallisto_output_dir>/<sample>/abundance.tsv`, `<methyldackel_dir>/<sample>_CpG.methylKit`, `<peaks_dir>/<sample>_peaks.narrowPeak`. Substitute actual selected-pipeline config values and sample IDs; do not copy placeholders literally.
- When changing MultiQC dependencies, dry-run the final report at `<multiqc_results_dir>/multiqc_report.html`.
- Keep verification scoped to the active pipeline; switching `config["pipeline"]` is a config change, not a harmless test setup detail.
- Dry-runs and DAG generation still execute module-level Python, including sample discovery and `os.makedirs` for output, log, and index directories. Inspect configured paths first.
- A successful dry-run verifies DAG construction, not tool execution or biological results. If missing inputs or permissions block validation, report the command and blocker; do not create fake production inputs. Run real analysis only when it is part of the task.

## Workflow Conventions

- Each top-level workflow reads config, validates required keys, normalizes values back into `config[...]`, builds output/MultiQC lists, then includes rule files.
- Included rules expect parent workflows to set keys such as `samples_info`, output directories, `multiqc_input_files`, and `multiqc_analysis_dirs`.
- `multiqc_input_files` establishes dependencies that must finish before reporting; `multiqc_analysis_dirs` specifies where MultiQC searches. Update both as needed when adding reports or moving outputs; the report output directory alone is not a substitute for the analysis directories.
- Preserve explicit BAM-index dependencies for tools that need indexed BAMs. Match each producer's filename: RNA-seq Picard emits `.markdup.bai`, while samtools-indexed BAMs here use `.bam.bai`.
- Workflow modules import helpers via `sys.path.insert(0, os.path.abspath("src"))` followed by direct imports from `utils`.
- Add user-facing config keys to the owning workflow, `config/config.yaml`, and `README.md`; add notable fixes/features to `backlog.md`.
- Use `os.path.join(...)` for paths and preserve established suffixes such as `_pe`, `_se`, `_R1`, `_R2`, `.filtered.sorted.bam`, and `.markdup.metrics.txt`.
- Shell blocks should run tools through `pixi run` and redirect stdout/stderr in the existing `{log}.out` / `{log}.err` style.

## Snakemake Best Practices

Apply these to new or changed rules within the existing module/Pixi architecture; they are not a request to migrate unrelated workflows. Check feature availability against the Snakemake version in `pixi.lock` before adopting newer syntax.

### Design The DAG

- Give each rule a clear transformation; connect stages through files, not assumed execution order. Declare every consumed data/reference/index file in `input`, including files otherwise hidden in command arguments.
- Prefer named `input`, `output`, and `params`. Keep non-file options in `params`; give every output one unambiguous producer.
- Use wildcards for individual jobs and `expand()` for aggregate targets. Constrain ambiguous wildcards; keep input functions deterministic and free of side effects.
- Keep sample selection explicit through `SAMPLES_INFO`. Do not glob generated output folders to decide which jobs should run: a clean checkout and a partial run must request the same intended results.
- Keep the default `all` target explicit and config-aware. Optional analyses must gate rules, final targets, and reporting dependencies consistently.
- Reuse shared rules through the existing includes and `use rule ... with` pattern. Move substantial analysis logic into focused scripts rather than expanding shell blocks or parse-time Python.

### Make Execution Reliable

- Quote individual paths with `{input.bam:q}` / `{output.bam:q}` and lists with `{input:q}`; escape literal shell/awk braces as `{{` and `}}`.
- Preserve Bash strict-mode failure propagation. Use unique per-job logs and scratch paths; never mask a failed tool with `|| true`.
- Prefer explicit output files. Use `directory()` only for a directory exclusively owned by one job; reruns can delete it. Use `temp()` only for disposable intermediates.
- Use `ensure(..., non_empty=True)` only when empty output is invalid; an empty peak set may be legitimate.
- Avoid adding parse-time filesystem writes or expensive computations. Perform processing inside rules so dry-runs remain useful; existing directory creation is documented above.
- Keep inputs immutable and jobs rerunnable. Never let concurrent samples append to one report or overwrite shared scratch files; aggregate their separate outputs in a dedicated rule.

### Resources, Configuration, And Reproducibility

- Pass effective `{threads}` to tools: Snakemake may reduce it. Budget concurrent subprocesses together. Declare per-job `mem_mb`, `disk_mb`, and `runtime` where useful; resource declarations guide scheduling, not tool memory enforcement.
- When changing piped commands or Bismark workers, verify that subprocess allocations respect the job's effective CPU budget, including reduced-core runs. Do not copy fixed config thread counts into commands independently of that budget.
- Validate new config options for type, range, and allowed values in the owning workflow, with actionable errors. Preserve the current lazy reference-file validation contract. Use a schema if configuration complexity warrants one; do not add a schema dependency for a trivial option.
- Keep scientific parameters in documented config and machine-specific scheduling settings in invocation options or profiles when needed. Continue using the committed Pixi environment; do not introduce a second environment manager incidentally.
- Expose random seeds for stochastic steps when supported, preserve sample ordering in aggregate outputs, and record meaningful tool options in `params` so rerun decisions can track changes.
- For behavior changes, choose checks that exercise the affected branch: PE/SE support, optional outputs, one versus multiple samples, or alternative peak callers as applicable. A small execution test with suitable data complements a dry-run when execution is in scope; do not launch production data merely to check syntax.

Official references: [best practices](https://snakemake.readthedocs.io/en/stable/snakefiles/best_practices.html), [rules](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html), [modularization](https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html), and [configuration](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html). Follow repository conventions where the general documentation offers alternative layouts or environment strategies.

## Samples And References

- A nonempty `<pipeline>.samples_info` replaces discovery; it is not merged with discovered files. Otherwise, `src/utils.py` scans only the immediate `raw_fastqs_dir` for compressed FASTQs.
- Discovery supports `_R1/_R2`, `_1/_2`, `.R1/.R2`, optional `_001` suffixes, and single-end `.fastq.gz`/`.fq.gz` files. It does not merge lanes; avoid ambiguous duplicate files for a sample/read.
- `rnaseq`, `wgbs`, and `atacseq` reject non-paired samples; `chip_cr` accepts PE and SE and infers `type` from presence of `R2` when omitted.
- ChIP/CUT&RUN input subtraction is name-driven: signal `<base>_repN` pairs with input `<base>_input_repN`; `<base>` may contain underscores.
- Workflows use `validate_required_path` to check that required reference values are nonempty, not that files exist. Missing files can fail during DAG input resolution; directory permissions can fail during module setup. The stricter helpers in `src/utils.py` are not the current workflow validation contract.
- The pipeline does not download references. RNA-seq can build `kallisto_index` and HISAT2 indexes, WGBS creates `Bisulfite_Genome/` next to `ref_genome`, and ChIP/ATAC create configured Bowtie2 index directories.

## Pipeline Gotchas

- RNA-seq gene body coverage is optional via `rnaseq.gene_body_coverage.enabled`; `refgene_bed` is required only when enabled, and its heatmap output exists only with at least 3 samples.
- RNA-seq HISAT2 alignments include per-sample read groups needed by Picard. Duplicate marking retains reads; do not silently change it to duplicate removal.
- RNA-seq, ChIP/CUT&RUN, and ATAC-seq correlation heatmaps are generated only when at least 2 BAMs are available; do not add unconditional MultiQC/report dependencies for those files.
- ChIP/CUT&RUN defaults to `results/chipseq_cutrun`, despite the config key and rule prefix being `chip_cr`. Fingerprint/correlation includes controls; heatmap selection currently excludes any sample name containing `input` and prefers a subtracted bigWig when available.
- WGBS uses `wgbs.log_dir` when set; other pipelines default logs under `logs/<pipeline>`.
- WGBS Bismark thread allocation derives Bowtie2 threads from `wgbs.bismark.threads / wgbs.bismark.parallel`; keep `threads >= parallel`.
- ATAC-seq supports only `peak_caller: "macs3"` or `"genrich"`; Genrich uses a queryname-sorted intermediate but both paths expose the shared `.narrowPeak` output.
- ATAC-seq marks duplicates for metrics, then removes duplicate-flagged reads during filtering. Downstream QC, bigWigs, and peaks use the filtered BAMs.
