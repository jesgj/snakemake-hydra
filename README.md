# Modular NGS Pipeline

![Pipeline logo](logo.png)

Snakemake + Pixi pipeline for four NGS workflows: RNA-seq, WGBS, ChIP-seq/CUT&RUN, and ATAC-seq. A single value in `config/config.yaml` selects which workflow is imported and run.

## Repository Map
- `Snakefile`: top-level entrypoint; dispatches to the selected workflow.
- `config/config.yaml`: main configuration and sample/reference paths.
- `workflows/*.smk`: top-level workflow modules for `rnaseq`, `wgbs`, `chip_cr`, and `atacseq`.
- `workflows/rules/**/*.smk`: reusable and pipeline-specific Snakemake rules.
- `src/utils.py`: sample discovery and config helper functions.
- `backlog.md`: notable implemented changes and planned work.

## Requirements
- Linux.
- [Pixi](https://pixi.sh/) installed.
- Disk space for FASTQ, BAM, reference indexes, and results.
- Optional Graphviz `dot` to render DAG images.

## Setup
```bash
git clone <repo-url>
cd snakemake-hydra
pixi install
```

`pixi.toml` defines dependencies but no tasks, so run tools explicitly with `pixi run ...`.

## Configure
Edit `config/config.yaml` before running.

```yaml
pipeline: "wgbs"  # rnaseq, wgbs, chip_cr, or atacseq
```

At minimum, configure the selected pipeline block:
- `raw_fastqs_dir` and/or `samples_info`.
- Reference FASTA/BED/index paths required by that pipeline.
- Output directories if the defaults are not suitable.
- `<pipeline>.multiqc_results_dir` for the MultiQC report location.

The checked-in config contains safe placeholder paths. Replace them for your environment, but avoid changing unrelated pipeline blocks unless you intend to switch or update them.

Reference inputs are user-managed. The pipeline does not download genomes or annotations. Config values are checked for presence during workflow setup; missing files, permissions, or unwritable index directories fail later during Snakemake input resolution or tool execution.

## Samples
Manual sample declarations live under the selected pipeline's `samples_info`:

```yaml
chip_cr:
  raw_fastqs_dir: "data/chip_fastqs"
  samples_info:
    liver_rep1:
      R1: "data/chip_fastqs/liver_rep1_R1.fastq.gz"
      R2: "data/chip_fastqs/liver_rep1_R2.fastq.gz"
      type: "PE"
```

If `samples_info` is missing or empty, `src/utils.py` discovers FASTQs in `raw_fastqs_dir` using these patterns:
- Paired-end: `_R1/_R2`, `_1/_2`, `.R1/.R2`.
- Single-end: `<sample>.fastq.gz` or `<sample>.fq.gz`.

Sample support differs by pipeline:
- `rnaseq`, `wgbs`, and `atacseq` require paired-end samples with `R1` and `R2`.
- `chip_cr` accepts paired-end and single-end samples; `type` is inferred from `R2` when omitted.
- ChIP/CUT&RUN subtraction pairs signal `<base>_repN` with input `<base>_input_repN`; `<base>` may contain underscores.

## Run Commands
Dry-run the selected pipeline:

```bash
pixi run snakemake -n --cores 1
```

Run the selected pipeline:

```bash
pixi run snakemake --cores 8 --printshellcmds --latency-wait 60
```

Re-run incomplete jobs:

```bash
pixi run snakemake --cores 8 --rerun-incomplete
```

Dry-run one imported rule path:

```bash
pixi run snakemake -n --cores 1 --until <prefixed_rule_name>
```

Dry-run one concrete output:

```bash
pixi run snakemake -n --cores 1 <path/to/output>
```

Generate a DAG:

```bash
pixi run snakemake --dag --cores 1
pixi run snakemake --dag --cores 1 | dot -Tpng > dag.png
```

Focused rules must belong to the current `pipeline` value because `Snakefile` imports only the selected workflow. Imported rule names are prefixed, for example `rnaseq_hisat2_align_pe`, `wgbs_bismark_alignment`, `chip_cr_fastqc_raw_pe`, or `atacseq_bowtie2_align_pe`.

Useful concrete dry-run targets:
- RNA-seq: `results/rnaseq/kallisto/sample/abundance.tsv`.
- WGBS: `results/wgbs/methyldackel/sample_CpG.methylKit`.
- ATAC-seq: `results/atacseq/peaks/sample_peaks.narrowPeak`.
- MultiQC: `<multiqc_results_dir>/multiqc_report.html`.

## Validation
Small real-data fixtures for RNA-seq, WGBS, ChIP-seq/CUT&RUN, and both ATAC-seq
peak callers are documented in [test_data/README.md](test_data/README.md).
The preparer retains seeded random subsets, removes temporary source FASTQs,
and generates separate configs without editing production paths.

There is no configured CI, formatter, linter, typechecker, or unit-test suite in this repository. Use Snakemake dry-runs as the main validation path:
- Config-only changes: run `pixi run snakemake -n --cores 1` for the selected pipeline.
- Rule changes: dry-run through the affected prefixed rule with `--until`.
- Filename/path changes: dry-run a concrete output target.
- MultiQC dependency changes: dry-run the final `multiqc_report.html` path.

## Pipeline Notes

### RNA-seq
- Requires paired-end samples.
- Required config values include `transcriptome_fasta`, `kallisto_index`, `ref_genome`, and `hisat2_index_dir`.
- The pipeline can build the Kallisto index and HISAT2 indexes at the configured paths.
- Gene body coverage is controlled by `rnaseq.gene_body_coverage.enabled`; `refgene_bed` is required only when enabled.
- Gene body coverage heatmap output is tracked only when at least 3 samples are analyzed.
- The deepTools correlation heatmap is generated only when at least 2 aligned BAMs are available.

### WGBS
- Requires paired-end samples.
- `ref_genome` must point to the reference used by Bismark; Bismark creates `Bisulfite_Genome/` next to it.
- `wgbs.log_dir` overrides the default `logs/wgbs` location.
- Bismark Bowtie2 threads are derived from `wgbs.bismark.threads / wgbs.bismark.parallel`; keep `threads >= parallel`.
- MethylDackel outputs include `<sample>_CpG.methylKit` and merged-context `<sample>_CpG.bedGraph`.

### ChIP-seq/CUT&RUN
- Supports paired-end and single-end samples.
- Requires `ref_genome`, `gene_bed`, and `bowtie2_index_dir` config values.
- Bowtie2 indexes are created in `bowtie2_index_dir` when needed.
- `plotFingerprint` and correlation include input controls; heatmap `computeMatrix`/`plotHeatmap` excludes `_input` samples.
- Subtracted bigWigs are produced only for samples with matching `<base>_repN` and `<base>_input_repN` names.

### ATAC-seq
- Requires paired-end samples.
- Requires `ref_genome` and `bowtie2_index_dir` config values.
- `peak_caller` must be `macs3` or `genrich`.
- Both peak-caller paths expose the shared `peaks/<sample>_peaks.narrowPeak` output; Genrich uses an intermediate queryname-sorted BAM.
- Duplicate reads are marked for metrics, then removed during filtered BAM generation before downstream QC, bigWig generation, and peak calling.
- The deepTools correlation heatmap is generated only when at least 2 filtered BAMs are available.

## Outputs
Output paths are config-driven. Common defaults and examples include:
- MultiQC report: `<pipeline>.multiqc_results_dir/multiqc_report.html`.
- RNA-seq Kallisto: `results/rnaseq/kallisto/<sample>/abundance.tsv`.
- RNA-seq aligned BAM: `results/rnaseq/aligned_bams/<sample>_pe.sorted.bam`.
- WGBS MethylDackel: `results/wgbs/methyldackel/<sample>_CpG.methylKit`.
- ChIP/CUT&RUN bigWig: `results/chipseq_cutrun/bigwigs/<sample>_pe.bw` or `<sample>_se.bw`.
- ATAC-seq peaks: `results/atacseq/peaks/<sample>_peaks.narrowPeak`.

## Troubleshooting
- No samples found: check that `samples_info` is populated or `raw_fastqs_dir` exists with supported FASTQ names.
- Rule not found in a focused run: check `config["pipeline"]`; only the selected workflow's prefixed rules are imported.
- RNA-seq/WGBS/ATAC reject a sample: these workflows require paired-end `R1` and `R2` inputs.
- Missing correlation heatmap: the workflow skips correlation when fewer than 2 BAMs are available.
- Missing RNA-seq gene body coverage heatmap: it is emitted only with at least 3 samples and when gene body coverage is enabled.
- Missing ChIP/CUT&RUN subtraction: verify names follow `<base>_repN` and `<base>_input_repN`.
- Reference or index failure: configured reference files must exist, and index directories must be writable.

## License
See `LICENSE`.

![Pipeline](pipeline_v2.svg)
