# Agent Guidance for `snakemake-hydra`

This repository is a modular NGS pipeline built with Snakemake and Pixi.
Treat the workflow files as the source of truth and keep edits small, local,
and consistent with existing pipeline behavior.

## Repository Snapshot
- Entrypoint: `Snakefile`
- Environment manager: `pixi.toml`, `pixi.lock`
- Main config: `config/config.yaml`
- Shared Python helpers: `src/utils.py`
- Workflow modules: `workflows/rnaseq.smk`, `workflows/wgbs.smk`, `workflows/chip_cr.smk`, `workflows/atacseq.smk`
- Reusable rules live under `workflows/rules/**`
- Operational notes are tracked in `backlog.md`
- `Snakefile` dispatches to one subworkflow based on `config["pipeline"]` and re-exports that module's `all` rule as the default target.
- Supported pipelines: `rnaseq`, `wgbs`, `chip_cr`, `atacseq`

## Cursor / Copilot Rules
- No Cursor rules were found in `.cursor/rules/`.
- No `.cursorrules` file was found.
- No Copilot instructions were found in `.github/copilot-instructions.md`.
- If any of those files are added later, follow them in addition to this file.

## Build / Run / Test Commands
- Install dependencies: `pixi install`
- Run any tool inside the managed environment: `pixi run <command>`
- Dry-run the currently selected pipeline: `pixi run snakemake -n --cores 1`
- Run the selected pipeline: `pixi run snakemake --cores 8 --printshellcmds --latency-wait 60`
- Re-run incomplete jobs: `pixi run snakemake --cores 8 --rerun-incomplete`
- There is no dedicated formatter, linter, or unit-test suite configured.
- Do not invent repo-wide `pytest`, `ruff`, `black`, or similar workflows.
- In practice, validation is done with Snakemake dry-runs and targeted rule or target execution.

## Best Single-Test Equivalents
- Dry-run one rule path: `pixi run snakemake -n --cores 1 --until <rule_name>`
- Execute up to one rule: `pixi run snakemake --cores 4 --until <rule_name>`
- Dry-run one concrete output: `pixi run snakemake -n --cores 1 <path/to/output>`
- Build one concrete output: `pixi run snakemake --cores 4 <path/to/output>`
- Confirm that `config/config.yaml` selects the pipeline that owns the rule or output before running targeted commands.

## Useful Single-Test Examples
- RNA-seq rule path: `pixi run snakemake --cores 4 --until rnaseq_hisat2_align_pe`
- RNA-seq output target: `pixi run snakemake --cores 4 results/rnaseq/kallisto/sample/abundance.tsv`
- ChIP/CUT&RUN rule path: `pixi run snakemake --cores 4 --until chip_cr_fastqc_raw_pe`
- WGBS output target: `pixi run snakemake --cores 4 results/wgbs/methyldackel/sample_CpG.methylKit`
- ATAC-seq output target: `pixi run snakemake --cores 4 results/atacseq/peaks/sample_peaks.narrowPeak`
- Print the DAG: `pixi run snakemake --dag --cores 1`
- Render the DAG if Graphviz is installed: `pixi run snakemake --dag --cores 1 | dot -Tpng > dag.png`

## Validation Tips
- For config-only edits, a full dry-run is usually enough.
- For rule changes, prefer `--until <rule_name>` first.
- For path or naming changes, dry-run a concrete output target.
- When changing MultiQC gating, dry-run the final report target.
- Keep validations scoped to the active pipeline instead of probing unrelated workflows.

## Style Guidelines

### General Editing
- Prefer the smallest correct change.
- Preserve existing structure and naming unless there is a concrete reason to refactor.
- Do not reformat large sections of workflow code just for style.
- Keep path-sensitive changes localized; many rules assume stable directory layouts and filename suffixes.
- When implementing a notable fix or feature, add a short note to `backlog.md`.

### Imports
- Follow the existing pattern: standard library imports first.
- In workflow modules, `sys.path.insert(0, os.path.abspath("src"))` is the established way to import `utils.py`.
- Keep local imports from `utils` directly after the `sys.path.insert(...)` line.
- Do not add third-party Python imports unless they are truly required and fit the repo's current approach.

### Python Conventions
- Functions use `snake_case`.
- Module-level config-derived constants use uppercase names like `RAW_DIR` and `ALIGNMENT_DIR`.
- Keep helper functions small and focused.
- Use docstrings for helper functions in `src/utils.py`.
- The repository does not currently use type annotations broadly; keep new code consistent unless there is a strong reason to introduce types in a whole file.
- Prefer explicit returns and straightforward control flow over clever abstractions.

### Snakemake Workflow Structure
- Each top-level workflow pulls values from `config`, validates what it needs, and writes normalized values back into `config[...]` for included rule files.
- Preserve that pattern when adding new config keys needed by included rules.
- Use `include:` to compose workflows from rule modules.
- Use `use rule ... as ... with:` when specializing generic rules.
- Keep `rule all` as the workflow's final contract for expected outputs.

### Rule Authoring
- Rule names are `snake_case` and usually reflect scope plus action.
- Keep `input`, `output`, `params`, `threads`, and `log` explicit.
- Prefer `expand(...)` or small local helper functions for assembling output lists.
- Use `wildcard_constraints` where PE/SE naming would otherwise collide.
- When a tool builds an index, follow the existing sentinel-file approach where appropriate, such as `bismark_index_created.OK`.

### Shell Blocks
- Prefer `pixi run <tool>` inside `shell:` blocks for reproducibility.
- Use triple-quoted shell blocks for multi-line commands.
- Redirect command output to `{log}.out` and `{log}.err` in the existing style.
- Surface tunable values via `params` or nested config keys instead of hard-coding new magic values.
- Match the surrounding quoting style and avoid unnecessary shell complexity.

### Paths and Naming
- Build paths with `os.path.join(...)` rather than hard-coded separators.
- Keep output directory patterns consistent with the active pipeline, typically under `results/<pipeline>/...` and `logs/<pipeline>/...`.
- Preserve established suffixes such as `_pe`, `_se`, `_R1`, `_R2`, `.filtered.sorted.bam`, and `.markdup.metrics.txt`.
- MultiQC outputs are expected at `<multiqc_results_dir>/multiqc_report.html`.

### Config Conventions
- `config/config.yaml` is the main control surface.
- `pipeline` must be one of `rnaseq`, `wgbs`, `chip_cr`, or `atacseq`.
- Sample definitions live under `<pipeline>.samples_info` and are keyed by sample name.
- If `samples_info` is missing or empty, sample discovery may fall back to `raw_fastqs_dir` via `src/utils.py`.
- Tool-specific options generally live under nested keys like `fastp.extra_args`, `deeptools.plotCorrelation.extra_args`, or `picard.markduplicates.extra_args`.
- Reference FASTA and BED inputs are user-managed; the pipeline does not download them.

### Sample Handling Rules
- `rnaseq`, `wgbs`, and `atacseq` currently support paired-end samples only.
- `chip_cr` supports both paired-end and single-end samples.
- ChIP/CUT&RUN input subtraction expects names like `<base>_rep1` and `<base>_input_rep1`.
- Raw FASTQ discovery in `src/utils.py` supports common `_R1/_R2`, `_1/_2`, and `.R1/.R2` naming patterns.

### Error Handling
- Fail fast on invalid configuration or unsupported sample layouts.
- Prefer `ValueError` with short, actionable messages.
- Reuse helper validators from `src/utils.py` when possible, especially for required paths.
- Validate assumptions near the top of the workflow file so failures happen before Snakemake expands large portions of the DAG.

### Formatting Expectations
- Match the surrounding file's formatting style instead of applying a new one.
- Keep comments brief and only where they clarify non-obvious workflow logic.
- Prefer readable multi-line assignments for long config lookups or output lists.
- Avoid introducing new helper layers when one local helper function is enough.

## Workflow-Specific Notes
- `rnaseq` validates `transcriptome_fasta`, `kallisto_index`, `ref_genome`, `hisat2_index_dir`, and optionally `gene_body_coverage.refgene_bed`.
- `wgbs` expects `ref_genome` to exist and lets Bismark create `Bisulfite_Genome/` next to that reference.
- `chip_cr` validates `ref_genome`, `gene_bed`, and `bowtie2_index_dir`, then normalizes PE/SE sample types before rule expansion.
- `atacseq` validates `peak_caller` and currently supports only `macs3` or `genrich`.
- Each top-level workflow assembles `multiqc_input_files` and `multiqc_analysis_dirs` before including `rules/multiqc.smk`.

## Agent-Specific Guidance
- Read the target workflow and related rule modules before editing behavior.
- If you add a new output that should appear in MultiQC gating, update `multiqc_input_files` or the workflow's output aggregation accordingly.
- If you add a new config key, thread it through the owning workflow and document it in `README.md` when user-facing.
- Do not assume there is a conventional Python test suite hiding elsewhere.
- Prefer validating with a dry-run first, then a targeted rule or output if the change is execution-related.

## High-Value Files
- `Snakefile`, `config/config.yaml`, `src/utils.py`
- `workflows/rnaseq.smk`, `workflows/wgbs.smk`, `workflows/chip_cr.smk`, `workflows/atacseq.smk`
- `workflows/rules/qc.smk`, `workflows/rules/trimming_and_qc.smk`, `workflows/rules/multiqc.smk`
