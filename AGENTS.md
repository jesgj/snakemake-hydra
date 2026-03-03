# Agent Guidance for Snakemake_pipeline_TFM

This repository is a Snakemake + Pixi bioinformatics pipeline with multiple
workflows (RNA-seq, WGBS, ChIP-seq/CUT&RUN, ATAC-seq). Use the guidance below
when operating as an agentic coding assistant in this repo.

## Repository Facts
- Primary entrypoint: `Snakefile`
- Workflow modules: `workflows/*.smk` and `workflows/rules/**/*.smk`
- Shared Python utilities: `src/utils.py`
- Configuration: `config/config.yaml`
- Environment manager: Pixi (`pixi.toml`, `pixi.lock`)

## Build / Run / Lint / Test Commands
There is no explicit lint or test suite configured. Treat Snakemake dry-runs
or targeted rules as the practical “tests.”

### Environment
- Install dependencies: `pixi install`
- Run a command inside the env: `pixi run <command>`

### Snakemake runs
- Dry run: `pixi run snakemake -n --cores 1`
- Standard run: `pixi run snakemake --cores 8 --printshellcmds --latency-wait 60`
- Re-run incomplete: `pixi run snakemake --cores 8 --rerun-incomplete`

### “Single test” equivalents
Use one of the following to exercise a specific rule or output target:
- Run until a rule: `pixi run snakemake --cores 4 --until <rule_name>`
  - Example: `pixi run snakemake --cores 16 --until chip_cr_fastqc_raw_pe`
- Run a single output target: `pixi run snakemake --cores 4 <path/to/output>`
  - Example: `pixi run snakemake --cores 4 results/rnaseq/kallisto/sample/abundance.tsv`
- Run a subworkflow by changing `pipeline:` in `config/config.yaml` first.

### DAG / graph (optional)
- `pixi run snakemake --dag --cores 1 | dot -Tpng > dag.png`

### Lint / tests
- None configured. Do not invent a lint/test runner.
- Prefer a dry run or a targeted rule as a sanity check.

## Code Style and Conventions
Follow existing patterns. The workflows are the source of truth.

### Snakemake style
- Use `include:` for rule modules and keep workflows composable.
- Put config lookups at the top of a workflow file and mirror them back into
  `config[...]` so included rules can access them.
- Build paths with `os.path.join` and avoid hard-coded separators.
- Prefer `expand(...)` for enumerating outputs; avoid manual loops unless needed.
- Rule names are snake_case and indicate scope (e.g., `fastqc_raw_pe`).
- Rule outputs use consistent suffixes:
  - Paired-end: `_pe` or `_R1/_R2` in filenames
  - Single-end: `_se` or `_SE` in filenames
- Use sentinel files to mark index builds (e.g., `index_built.OK`).
- Keep `threads`, `params`, and `log` explicit on each rule.
- Logs are routed to `logs/<pipeline>/...` with separate `.out` and `.err`.
- Prefer `pixi run <tool>` in shells for reproducible tool execution.

### Python in workflows and `src/`
- Standard library imports first; no third-party imports are used in this repo.
- Use `sys.path.insert(0, os.path.abspath("src"))` in workflows when importing
  `src/utils.py` (existing pattern).
- Functions are `snake_case`, small, and focused.
- Use docstrings for helper functions in `src/utils.py`.
- Types are not annotated today; keep new code consistent unless there is a
  compelling reason to introduce typing across the file/module.
- Prefer explicit `ValueError` with helpful messages for invalid config/data.

### YAML config conventions (`config/config.yaml`)
- `pipeline` must be one of: `rnaseq`, `wgbs`, `chip_cr`, `atacseq`.
- Each pipeline config contains directories for raw inputs, outputs, and logs.
- Tool parameters live under `<tool>.extra_args` or similar nested keys.
- When adding new config keys, keep naming consistent and document in README.

### Naming conventions and paths
- Samples are identified by keys in `samples_info`.
- Raw FASTQ naming conventions are handled in `src/utils.py`.
- Output paths are consistent across workflows (e.g., `results/<pipeline>/...`).
- Logs are placed under `logs/<pipeline>/<step>/`.

### Error handling
- Fail fast when required inputs are missing (e.g., raise `ValueError`).
- Validate assumptions early (example: ATAC-seq only supports paired-end).
- Keep error messages short and actionable.

### Shell blocks in rules
- Use triple-quoted `shell` blocks for multi-line commands.
- Redirect stdout/stderr to `{log}.out` / `{log}.err` where used.
- Keep command parameters explicit; avoid magic defaults in shell strings.

## Cursor / Copilot Rules
- No Cursor rules found in `.cursor/rules/` or `.cursorrules`.
- No Copilot instructions found in `.github/copilot-instructions.md`.

## Editing Guidance
- Do not reformat large sections unless necessary; preserve existing style.
- Keep changes localized and minimal; pipeline code is sensitive to paths.
- If adding a rule, wire it into the appropriate workflow and update
  `multiqc_input_files` if the output should be reported.
- Log notable changes in `backlog.md` when implementing features or fixes.

## Useful Files
- `Snakefile`
- `config/config.yaml`
- `workflows/rnaseq.smk`
- `workflows/wgbs.smk`
- `workflows/chip_cr.smk`
- `workflows/atacseq.smk`
- `workflows/rules/*.smk`
- `src/utils.py`
