# Agent Guidance for `snakemake-hydra`

## Source Of Truth

- This is a Snakemake + Pixi NGS pipeline, not a Python package; workflow behavior lives in `Snakefile`, `workflows/*.smk`, `workflows/rules/**/*.smk`, and focused helpers under `src/` and `scripts/`.
- `Snakefile` reads the selected configuration, builds config objects for the original four modules and the optional amplicon module, then imports only the selected pipeline with prefixed rule names like `rnaseq_*`, `wgbs_*`, `chip_cr_*`, `atacseq_*`, or `amplicon_*`.
- Keep the original four top-level pipeline config mappings present: `Snakefile` accesses each before selecting rules. Inside a module, `config` is that pipeline's subsection plus `pipeline`; use `config["samples_info"]`, not `config["rnaseq"]["samples_info"]`.
- Supported `config["pipeline"]` values are `rnaseq`, `wgbs`, `chip_cr`, `atacseq`, and `amplicon`; focused rule/target runs must match the selected pipeline or the rule will not be imported.
- `pixi.toml` targets `linux-64`, has dependencies but no tasks, and is paired with committed `pixi.lock`; use explicit `pixi run ...` commands. Keep manifest and lockfile consistent when changing dependencies.
- GitHub Actions CI runs locked Pixi installation, Python syntax checks, the existing 16S smoke helper, and all eight real-data fixture variants on PRs and pushes to main; a manual input enables the public ENA fixture. Lint recommendations are retained in nonblocking diagnostic steps; the MACS3 job also runs discovery/empty-peak regressions. No pre-commit, formatter, typechecker, or general unit-test config is present; do not invent repo-wide `pytest`, `ruff`, `black`, etc. Snakemake's built-in `--lint` is available without a repository linter config.
- Inspect `git status --short` before editing and preserve existing user changes. Treat sample, reference, and output paths in `config/config.yaml` as user configuration; do not rewrite them to make validation pass.
- `config/config.local.yaml` is ignored by Git and automatically preferred by `Snakefile` when present; `HYDRA_CONFIG` overrides the selected config file.

## Where To Make Changes

- Pipeline setup, defaults, sample validation, output lists, and rule instantiation: `workflows/<pipeline>.smk`.
- Tool commands and their declared inputs/outputs: `workflows/rules/<pipeline>/*.smk`.
- Shared FASTQ QC/trimming: `workflows/rules/qc.smk` and `workflows/rules/trimming_and_qc.smk`.
- Shared BAM/deepTools QC and reporting: `workflows/rules/bam_qc.smk`, `workflows/rules/deeptools_qc.smk`, and `workflows/rules/multiqc.smk`. Check the including workflows when changing these shared interfaces.
- Sample discovery and path validation helpers: `src/utils.py`. The shell scripts under `src/` are Slurm launch examples with fixed resources/targets, not the workflow entrypoint.
- ATAC run normalization: `workflows/atacseq.smk`; FASTQ validation: `scripts/validate_fastq_pairs.py`; BAM read-group validation: `src/validate_bam_read_group.py` and `src/validate_merged_bam_read_groups.py`.
- Amplicon config validation: `src/amplicon_config.py`; read preparation/QC and R analysis: `scripts/amplicon_reads.py`, `scripts/amplicon_dada2.R`, and `scripts/amplicon_taxonomy.R`.
- CI: `.github/workflows/ci.yml`. Real-data fixture sources and preparation: `test_data/datasets.json` and `test_data/prepare.py`; execution verification: `test_data/finish_run.py`. See `test_data/README.md` for fixture limits and `docs/amplicon.md` for 16S configuration.

## Commands

Run from the repository root; config paths and the `src` import depend on the working directory.

- Install the committed environment: `pixi install --locked`. Use `pixi run --locked ...` for reproducible validation; deliberately regenerate the lockfile when changing dependencies.
- Dry-run the selected pipeline: `pixi run --locked snakemake -n --cores 1`.
- Run the selected pipeline: `pixi run --locked snakemake --cores 8 --printshellcmds --latency-wait 60`.
- Re-run incomplete jobs: `pixi run --locked snakemake --cores 8 --rerun-incomplete`.
- Dry-run one selected-pipeline rule path: `pixi run --locked snakemake -n --cores 1 --until <prefixed_rule_name>`.
- Dry-run one concrete output: `pixi run --locked snakemake -n --cores 1 <path/to/output>`.
- DAG only: `pixi run --locked snakemake --dag --cores 1`; render with Graphviz using `pixi run --locked snakemake --dag --cores 1 | dot -Tpng > dag.png`. Graphviz `dot` is optional and is not declared in `pixi.toml`.

## Validation Strategy

- Documentation-only changes: verify claims against source and run `git diff --check`; a pipeline run is unnecessary.
- For config-only edits, a selected-pipeline dry-run is usually the best verification.
- For rule edits, prefer `--until <prefixed_rule_name>` first, then a concrete output if paths or filenames changed.
- For workflow/rule edits, also use `pixi run --locked snakemake --lint` when the selected workflow can be parsed. Address findings introduced by the change; distinguish existing warnings and recommendations that do not fit the repository's Pixi environment strategy.
- Useful output templates: `<kallisto_output_dir>/<sample>/abundance.tsv`, `<methyldackel_dir>/<sample>_CpG.methylKit`, `<peaks_dir>/<sample>_peaks.narrowPeak`, `<featurecounts_dir>/consensus_peak_counts.txt`. Substitute actual selected-pipeline config values and sample IDs; do not copy placeholders literally.
- When changing MultiQC dependencies, dry-run the final report at `<multiqc_results_dir>/multiqc_report.html`, or `<amplicon.output_dir>/multiqc/multiqc_report.html` for 16S.
- Keep production-config verification scoped to the active pipeline. Exercise other branches through isolated fixture configs rather than editing the user's pipeline selection.
- Dry-runs and DAG generation still execute module-level Python, including sample discovery and `os.makedirs` for output, log, and index directories. Inspect configured paths first.
- A successful dry-run verifies DAG construction, not tool execution or biological results. If missing inputs or permissions block validation, report the command and blocker; do not create fake production inputs. Run real analysis only when it is part of the task.

### Isolated Real-Data Checks

- The eight original-pipeline variants are `rnaseq`, `wgbs`, `wgbs_full`, `chipseq`, `chipseq_se`, `cutrun`, `atacseq_macs3`, and `atacseq_genrich`. These use seeded subsets of real experimental reads; `chipseq_se` uses only R1 of paired ChIP libraries. They check execution and file interfaces, not biological accuracy.
- `atacseq_macs3` partitions each sample's real pairs into two disjoint software run fixtures with a shared library and distinct read groups; verification checks merged RG metadata and record counts against each run BAM. `atacseq_genrich` keeps legacy sample-level R1/R2 coverage. `--configs-only` regenerates config but does not create missing split FASTQs.
- Prepare only the needed variant: `pixi run --locked python test_data/prepare.py --variants rnaseq --provenance test_data/work/ci/rnaseq/provenance.json`. This downloads inputs/references and generates ignored configs; it does not launch the workflow. The preparer refuses to overwrite existing sampled FASTQs. Use `--configs-only --variants rnaseq` to regenerate config without downloading data.
- For fixture invocations, set `HYDRA_CONFIG` and pass the same file with `--replace-workflow-config --configfile`. This keeps subprocesses on the fixture config and prevents recursive merging of production samples. Example: `HYDRA_CONFIG=test_data/configs/rnaseq.yaml pixi run --locked snakemake -n --cores 2 --jobs 1 --replace-workflow-config --configfile test_data/configs/rnaseq.yaml`.
- When execution is needed, replace `-n` with `--printshellcmds --latency-wait 60` in that command. Keep local fixture runs sequential with `--jobs 1`; CI also uses `--resources ci_slots=1 --default-resources ci_slots=1`. CI runs `wgbs`, `chipseq_se`, and `atacseq_genrich` with one core, and the others with two. Scheduler limits do not enforce hard OS CPU/RAM caps, especially for Bismark subprocesses.
- After successful execution, verify and retain outputs with `pixi run --locked python test_data/finish_run.py rnaseq --cores 2 --keep-outputs --evidence test_data/work/ci/rnaseq/execution.json`. Without `--keep-outputs`, this helper deletes verified variant results and generated WGBS reference indexes. Use the actual variant/core count; do not run cleanup concurrently with jobs sharing a WGBS index.
- Write new evidence under ignored `test_data/work/` instead of overwriting committed historical `provenance.json` or `execution.json`. Historical records are not evidence that the current checkout passed. Report which checks actually ran.
- If Snakemake's user cache is unwritable, set `XDG_CACHE_HOME=/tmp/hydra-cache` for validation rather than changing pipeline paths.

## Workflow Conventions

- The original four top-level workflows read config, validate required keys, normalize values back into `config[...]`, build output/MultiQC lists, then include rule files. Amplicons call `src/amplicon_config.normalize()` and define their own targets and reporting inputs.
- Included rules expect parent workflows to set keys such as `samples_info`, output directories, `multiqc_input_files`, and `multiqc_analysis_dirs`.
- `multiqc_input_files` establishes dependencies that must finish before reporting; `multiqc_analysis_dirs` specifies where MultiQC searches. Update both as needed when adding reports or moving outputs; the report output directory alone is not a substitute for the analysis directories.
- Preserve explicit BAM-index dependencies for tools that need indexed BAMs. Match each producer's filename: RNA-seq Picard emits `.markdup.bai`, while samtools-indexed BAMs here use `.bam.bai`.
- Workflow modules import helpers via `sys.path.insert(0, os.path.abspath("src"))` followed by direct imports from `utils` or `amplicon_config`.
- Add user-facing config keys to the owning workflow/helper, `config/config.yaml`, and relevant documentation (`README.md` and, for 16S, `docs/amplicon.md`); add notable fixes/features to `backlog.md`.
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

- For the original four pipelines, a nonempty `<pipeline>.samples_info` replaces discovery; it is not merged with discovered files. Otherwise, `src/utils.py` scans only the immediate `raw_fastqs_dir` for compressed FASTQs. Amplicons use their separate manifest/discovery validator described below.
- Shared discovery supports `_R1/_R2`, `_1/_2`, `.R1/.R2`, optional `_001` suffixes, and single-end `.fastq.gz`/`.fq.gz` files. Duplicate sample/mate matches and orphan R2 files fail. Multiple `_L001`/`_L002`-style lanes remain separate samples with a warning; use explicit ATAC runs to group them.
- `rnaseq`, `wgbs`, and `atacseq` reject non-paired samples; `chip_cr` accepts PE and SE and infers `type` from presence of `R2` when omitted.
- ChIP/CUT&RUN input subtraction is name-driven: signal `<base>_repN` pairs with input `<base>_input_repN`; `<base>` may contain underscores.
- Reference values are required config keys but are generally only presence-checked at workflow parse time; selected ATAC prebuilt indexes are additionally checked for a complete `.bt2` or `.bt2l` family.
- The pipeline does not download references. RNA-seq can build `kallisto_index` and HISAT2 indexes, WGBS creates `Bisulfite_Genome/` next to `ref_genome`, and ChIP/ATAC can create configured Bowtie2 index directories.

## Pipeline Gotchas

- RNA-seq gene body coverage is optional via `rnaseq.gene_body_coverage.enabled`; `refgene_bed` is required only when enabled, and its heatmap output exists only with at least 3 samples.
- RNA-seq HISAT2 alignments include per-sample read groups needed by Picard. Duplicate marking retains reads; do not silently change it to duplicate removal.
- RNA-seq, ChIP/CUT&RUN, and ATAC-seq correlation heatmaps are generated only when at least 2 BAMs are available; do not add unconditional MultiQC/report dependencies for those files.
- ChIP/CUT&RUN defaults to `results/chipseq_cutrun`, despite the config key and rule prefix being `chip_cr`. Fingerprint/correlation includes controls; heatmaps exclude names matching `<base>_input_repN` and prefer a subtracted bigWig when available.
- WGBS and amplicons honor their section's `log_dir`, defaulting to `logs/<pipeline>`. RNA-seq, ChIP/CUT&RUN, and ATAC use fixed `logs/<pipeline>` paths; fixture output isolation does not isolate these logs.
- WGBS requires configured Bismark `threads >= parallel`. At execution, worker count is capped by effective Snakemake threads and per-worker Bowtie2 threads derive from the remaining budget. Omit `-p` for one thread: pinned Bismark rejects explicit `-p 1`. Bismark still has internal subprocess overhead; scheduler resources are not hard OS caps.
- ATAC-seq supports only `peak_caller: "macs3"` or `"genrich"`; Genrich uses a queryname-sorted intermediate but both paths expose the shared `.narrowPeak` output.
- ATAC-seq marks duplicates for metrics, then removes duplicate-flagged reads during filtering. Downstream QC, bigWigs, and peaks use the filtered BAMs.
- ATAC-seq `bowtie2_index_prefix` selects a read-only prebuilt index; otherwise Hydra builds six tracked files in `bowtie2_index_dir`. Generated `.bt2l` indexes require `bowtie2.large_index: true`.
- ATAC-seq index validation writes a prefix-specific marker under `reference_dir`; alignments depend directly on the six index files and that marker.
- ATAC-seq consensus peaks use sorted BED3 inputs and configurable `consensus_peaks.minimum_support`; featureCounts and FRiP depend on the cohort-wide consensus peak set.
- Empty per-sample peaks and empty consensus sets are valid. With no consensus intervals, emit a header-only SAF and zero-row count matrix/zero processed-fragment summary, skip featureCounts, and report FRiP as `NA`; keep the report dependencies intact.
- ATAC-seq MACS3 runs in true BAMPE fragment mode, while bigWigs are unshifted full-fragment CPM tracks; Tn5 insertion-site shifting is downstream analysis.

### ATAC Runs And Read Groups

- Preferred manifest: `samples_info.<sample>.runs.<run>` with original `R1`/`R2` and optional `read_group` fields (`id`, `sample`, `library`, `platform`, `platform_unit`). Legacy sample-level R1/R2 remains supported and normalizes to `run1`.
- Sample/run IDs allow letters, digits, dots, underscores, and hyphens, but forbid `__`, which separates them in filenames. Per-run files use `<sample>__<run>`; merged sample BAMs retain `<sample>_pe.sorted.bam`.
- Validate each original FASTQ pair before trimming; `fastq_validation.all_records` defaults to true. Sampled validation is explicitly weaker and cannot establish equal full-file counts.
- Align runs separately, validate each BAM's header and every alignment RG tag, then merge validated run BAMs per biological sample before duplicate marking. Preserve the validation dependencies and merged read-group report when changing this DAG.
- RG IDs must be unique across runs; RG sample must equal the biological sample. Defaults use `<sample>.<run>` for ID/platform unit, the sample (or sample-level `library`) for library, and `ILLUMINA` for platform. Preserve library identity across runs from the same library so Picard duplicate marking remains meaningful.

## Amplicons

- `workflows/amplicon.smk`, `workflows/rules/amplicon/*.smk`, `src/amplicon_config.py`, and `scripts/amplicon_*.R`/`scripts/amplicon_reads.py` own the configurable 16S branch. An amplicon section is optional for the original pipelines.
- Use either nonempty `samples_info` or a TSV `sample_sheet`; both together are rejected. Otherwise discover compressed or uncompressed `.fastq`/`.fq` in the immediate `raw_fastqs_dir`. Input paths, including TSV entries, resolve from the repository working directory. Normalization checks FASTQ existence without writing output directories.
- Set `primers.mode` explicitly to `remove` or `already_removed`. Removal uses configured IUPAC primers in read orientation, requires both 5' primers for paired inputs, and clips optional 3' read-through. Preparation audits mate synchronization before and after trimming; do not bypass this audit when changing positional pairing in DADA2 (`matchIDs = FALSE`).
- `mode: auto` chooses `merge` for PE and `forward` for SE. Explicit `forward`/`reverse` can select one mate from PE inputs; reverse-mode ASVs are reverse-complemented. Never silently fall back after merge failure or mix merge/forward/reverse modes in one cohort. Filter vectors always mean `[R1, R2]`.
- DADA2 owns quality filtering and learns errors separately by batch/read; pooling stays within those groups. Keep the seed, resolved config, retention table, and tool/reference metadata. Zero-filtered samples remain as zero count rows; entirely empty analyses fail.
- Keep `resolved_config.json`, scripts, prepared reads/stats, and optional taxonomy reference in Snakemake dependencies. `filtered/` and `errors/` are directories owned by the denoising rule and are final targets: deleting either must trigger reconstruction.
- Taxonomy requires an explicit DADA2-formatted training set and reference version when enabled; disabled taxonomy still emits a table and metadata with `not_requested`. ASV IDs follow lexical sequence order and can change when the cohort changes; compare sequences across runs.
- Amplicon MultiQC has its own explicit inputs/scans rather than the shared `multiqc_input_files`/`multiqc_analysis_dirs` interface. Preserve distinct raw/prepared FastQC names and the custom retention table. Report path: `<output_dir>/multiqc/multiqc_report.html`.
- Run `pixi run --locked python test_data/amplicon_smoke.py` for PE/SE, primers, TSV, reverse mode, taxonomy, batches, archive headers, expected overlap failure, repeatability, and deleted-directory rebuilds. It uses two real 1,500-pair samples and a small training reference bundled with pinned DADA2; primer tests add engineered sequences. No additional fixture downloads are needed. Outputs/evidence are ignored under `test_data/work/amplicon_smoke/`; `--workdir` can isolate another run.
- The separate public test uses 20,000 real 16S V4 pairs sampled from ENA ERR777695, experimental primers, and downloaded RDP 19 taxonomy. Prepare with `pixi run --locked python test_data/amplicon_public.py`, run with `HYDRA_CONFIG=test_data/configs/amplicon_public.yaml pixi run --locked snakemake --cores 2 --jobs 1`, then verify with `pixi run --locked python test_data/amplicon_public.py --verify`. Its primer/truncation choices are accession-specific, not defaults for other datasets. CI runs this only when manually enabled; the bundled smoke suite runs automatically.
