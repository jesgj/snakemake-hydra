# CTA_TFM_UOC - Modular NGS Pipeline (Snakemake + Pixi)

![Pipeline logo](logo.png)

Language: [English](#english) | [Español](#espanol)

<a id="english"></a>
## English

Modular and reproducible NGS analysis pipeline built with Snakemake and Pixi. A single configuration file selects one of four supported workflows: RNA-seq, WGBS, ChIP-seq/CUT&RUN, or ATAC-seq.

### Main features
- Workflow selection from `config/config.yaml`.
- Reproducible software environment with Pixi and `pixi.lock`.
- Automatic sample discovery or manual sample declaration through `samples_info`.
- Shared QC and trimming with FastQC and `fastp`.
- Automatic reference indexing when required by the selected workflow.
- Unified reporting through MultiQC and structured per-step logs.

### Repository layout
- `Snakefile`: top-level entrypoint that dispatches to the selected subworkflow.
- `config/config.yaml`: central configuration file.
- `workflows/`: workflow modules.
- `workflows/rules/`: reusable rule files.
- `src/utils.py`: sample discovery and input handling helpers.
- `src/scripts/`: helper scripts used by specific workflows.
- `backlog.md`: tracked notable implementation and maintenance changes.

### Requirements
- Linux.
- [Pixi](https://pixi.sh/) installed.
- Enough disk space for FASTQ, BAM, and result files.
- Optional: Graphviz if you want to render DAGs with `dot`.

### Quick installation
```bash
git clone <repo-url>
cd Snakemake_pipeline_TFM
pixi install
```

### Basic configuration
The main configuration lives in `config/config.yaml`. At minimum you should:

1. Set the subworkflow in `pipeline`.
2. Define `raw_fastqs_dir` and/or `samples_info`.
3. Configure reference files and index directories required by the selected workflow.

#### Pipeline selection
```yaml
pipeline: "chip_cr"  # options: rnaseq, wgbs, chip_cr, atacseq
```

#### Manual sample input
```yaml
chip_cr:
  raw_fastqs_dir: "data/chip_fastqs"
  samples_info:
    sampleA:
      R1: "data/chip_fastqs/sampleA_R1.fastq.gz"
      R2: "data/chip_fastqs/sampleA_R2.fastq.gz"
      type: "PE"
```

#### Automatic sample discovery
If `samples_info` is omitted or empty, samples are discovered from `raw_fastqs_dir` using common naming patterns:
- Paired-end: `_R1/_R2`, `_1/_2`, `.R1/.R2`.
- Single-end: `<sample>.fastq.gz` or `<sample>.fq.gz`.

#### Workflow support notes
- `rnaseq`, `wgbs`, and `atacseq` currently support paired-end samples only.
- `chip_cr` supports both paired-end and single-end samples.
- Automatic ChIP/CUT&RUN bigWig subtraction expects names such as `liver_rep1` and `liver_input_rep1`.

### Running the pipeline
Dry run:
```bash
pixi run snakemake -n --cores 1
```

Standard execution:
```bash
pixi run snakemake --cores 8 --printshellcmds --latency-wait 60
```

Run until a specific rule:
```bash
pixi run snakemake --cores 16 --until chip_cr_fastqc_raw_pe
```

Run a single output target:
```bash
pixi run snakemake --cores 4 results/rnaseq/kallisto/sample/abundance.tsv
```

Re-run incomplete jobs:
```bash
pixi run snakemake --cores 8 --rerun-incomplete
```

Generate a DAG (requires Graphviz):
```bash
pixi run snakemake --dag --cores 1 | dot -Tpng > dag.png
```

### Step-by-step tutorials

#### 1) ChIP-seq/CUT&RUN (paired-end with input control)
1. Place paired-end FASTQ files in a directory such as `data/chip_fastqs/`.
2. Name the ChIP sample as `base_repX` and the matching input as `base_input_repX` to enable automatic bigWig subtraction.
3. Edit `config/config.yaml`:

```yaml
pipeline: "chip_cr"
chip_cr:
  raw_fastqs_dir: "data/chip_fastqs"
  samples_info:
    liver_rep1:
      R1: "data/chip_fastqs/liver_rep1_R1.fastq.gz"
      R2: "data/chip_fastqs/liver_rep1_R2.fastq.gz"
      type: "PE"
    liver_input_rep1:
      R1: "data/chip_fastqs/liver_input_rep1_R1.fastq.gz"
      R2: "data/chip_fastqs/liver_input_rep1_R2.fastq.gz"
      type: "PE"
  ref_genome: "ref/genome.fa"
  bowtie2_index_dir: "ref/bowtie2"
  gene_bed: "ref/genes.bed"
```

4. Run:
```bash
pixi run snakemake --cores 16 --printshellcmds --latency-wait 60
```

5. Key outputs:
- `results/chipseq_cutrun/multiqc_report.html`
- `results/chipseq_cutrun/bigwigs/*.bw`
- `results/chipseq_cutrun/subtracted_bigwigs/*.subtracted.bw`
- `results/chipseq_cutrun/deeptools/heatmap.png`

#### 2) ChIP-seq/CUT&RUN (single-end)
For single-end data, define only `R1` and set `type: "SE"`:

```yaml
pipeline: "chip_cr"
chip_cr:
  raw_fastqs_dir: "data/chip_fastqs"
  samples_info:
    sampleSE:
      R1: "data/chip_fastqs/sampleSE.fastq.gz"
      type: "SE"
```

Outputs use the `_se` suffix, for example `sampleSE_se.sorted.bam` and `sampleSE_se.bw`.

#### 3) RNA-seq (paired-end)
1. Place paired-end FASTQ files in `data/rnaseq_raw_fastqs/`.
2. Configure `transcriptome_fasta`, `kallisto_index`, `ref_genome`, and `hisat2_index_dir`.
3. Example:

```yaml
pipeline: "rnaseq"
rnaseq:
  raw_fastqs_dir: "data/rnaseq_raw_fastqs"
  samples_info:
    Ctl_1:
      R1: "data/rnaseq_raw_fastqs/Ctl_1_R1.fq.gz"
      R2: "data/rnaseq_raw_fastqs/Ctl_1_R2.fq.gz"
  transcriptome_fasta: "ref/transcriptome.fa.gz"
  kallisto_index: "ref/kallisto.idx"
  ref_genome: "ref/genome.fa"
  hisat2_index_dir: "ref/hisat2"
  marked_bam_dir: "results/rnaseq/marked_bams"
  duplication_qc_dir: "results/rnaseq/duplication_qc"
  gene_body_coverage_dir: "results/rnaseq/gene_body_coverage"
  picard:
    java_opts: "-Xmx4g"
    markduplicates:
      extra_args: ""
  gene_body_coverage:
    refgene_bed: "ref/duumy.bed"
    minimum_length: 100
    format: "pdf"
```

4. Run and review:
- `results/rnaseq/kallisto/<sample>/abundance.tsv`
- `results/rnaseq/aligned_bams/<sample>_pe.sorted.bam`
- `results/rnaseq/marked_bams/<sample>_pe.markdup.bam`
- `results/rnaseq/duplication_qc/<sample>_pe.markdup.metrics.txt`
- `results/rnaseq/gene_body_coverage/all_samples.geneBodyCoverage.txt`
- `results/rnaseq/gene_body_coverage/all_samples.geneBodyCoverage.curves.pdf`
- `results/rnaseq/gene_body_coverage/all_samples.geneBodyCoverage.heatMap.pdf` when 3 or more BAMs are analyzed
- `results/rnaseq/bigwigs/<sample>_pe.bw`
- `results/rnaseq/multiqc_report.html`

Duplicate marking is tracked in separate BAMs for QC and does not remove reads from the original RNA-seq alignment outputs.

Use a real BED12 gene model matched to the RNA-seq genome assembly for production runs. `ref/duumy.bed` is only a placeholder for dry-run testing.

#### 4) WGBS (paired-end)
1. Place paired-end FASTQ files in `data/wgbs_raw_fastqs/`.
2. Configure `ref_genome`. Bismark will create `Bisulfite_Genome/` next to that reference.
3. Example:

```yaml
pipeline: "wgbs"
wgbs:
  raw_fastqs_dir: "data/wgbs_raw_fastqs"
  samples_info:
    Ctr_1:
      R1: "data/wgbs_raw_fastqs/Ctr_1_R1.fastq.gz"
      R2: "data/wgbs_raw_fastqs/Ctr_1_R2.fastq.gz"
  ref_genome: "ref/danrer11_lambda.fa"
  bismark:
    threads: 66
    parallel: 8
```

4. Key outputs:
- `results/wgbs/methyldackel/<sample>_CpG.methylKit`
- `results/wgbs/methyldackel_mergecontext/<sample>_CpG.bedGraph`
- `results/wgbs/multiqc_report.html`

5. Optional resource tuning:
- `wgbs.bismark.threads` and `wgbs.bismark.parallel`
- `wgbs.sambamba.filter_threads` and `wgbs.sambamba.sort_threads`
- `wgbs.methyldackel_mbias.threads` and `wgbs.methyldackel_extract.threads`

#### 5) ATAC-seq (paired-end)
1. Place paired-end FASTQ files in `data/atacseq_raw_fastqs/`.
2. Configure `ref_genome` and `bowtie2_index_dir`.
3. Example:

```yaml
pipeline: "atacseq"
atacseq:
  raw_fastqs_dir: "data/atacseq_raw_fastqs"
  samples_info:
    atac_rep1:
      R1: "data/atacseq_raw_fastqs/atac_rep1_R1.fastq.gz"
      R2: "data/atacseq_raw_fastqs/atac_rep1_R2.fastq.gz"
      type: "PE"
  ref_genome: "ref/Danio_rerio.GRCz11.dna.primary_assembly.fa.gz"
  bowtie2_index_dir: "ref/bowtie2_atacseq"
  macs3:
    genome_size: "1.4e9"
```

4. Key outputs:
- `results/atacseq/peaks/<sample>_peaks.narrowPeak`
- `results/atacseq/bigwigs/<sample>_pe.bw`
- `results/atacseq/fragment_qc/<sample>_pe.insert_size_metrics.txt`
- `results/atacseq/multiqc_report.html`

#### 6) Automatic sample discovery
If you want to avoid `samples_info`, keep only `raw_fastqs_dir`. The pipeline will detect samples from common FASTQ naming conventions. Keep in mind that `rnaseq`, `wgbs`, and `atacseq` require paired-end inputs, while `chip_cr` accepts paired-end and single-end data.

### Main outputs by pipeline

#### RNA-seq
- QC: `results/rnaseq/qc_raw`, `results/rnaseq/qc_trimmed`
- Trimming: `results/rnaseq/trimmed_fastqs`
- Kallisto: `results/rnaseq/kallisto/<sample>/abundance.tsv`
- HISAT2: `results/rnaseq/aligned_bams/<sample>_pe.sorted.bam`
- Duplicate marking QC: `results/rnaseq/marked_bams/<sample>_pe.markdup.bam` and `results/rnaseq/duplication_qc/<sample>_pe.markdup.metrics.txt`
- Gene body coverage: `results/rnaseq/gene_body_coverage/all_samples.geneBodyCoverage.txt` plus curve and optional heatmap plots
- BigWig: `results/rnaseq/bigwigs/<sample>_pe.bw`
- MultiQC: `results/rnaseq/multiqc_report.html`

#### WGBS
- QC: `results/wgbs/qc`, `results/wgbs/qc_trimmed`
- Trimming: `results/wgbs/trimmed_fastqs`
- BAMs: `data/wgbs/raw_bams`, `data/wgbs/dedup_bams`, `data/wgbs/filtered_bams`, `data/wgbs/sorted_filtered_bams`
- M-bias: `results/wgbs/mbias`
- MethylDackel: `results/wgbs/methyldackel`, `results/wgbs/methyldackel_mergecontext`
- MultiQC: `results/wgbs/multiqc_report.html`

#### ChIP-seq/CUT&RUN
- QC: `results/chipseq_cutrun/qc_raw`, `results/chipseq_cutrun/qc_trimmed`
- Trimming: `results/chipseq_cutrun/trimmed_fastqs`
- BAMs: `results/chipseq_cutrun/aligned_bams`, `results/chipseq_cutrun/filtered_bams`
- BigWig: `results/chipseq_cutrun/bigwigs`, `results/chipseq_cutrun/subtracted_bigwigs`
- deepTools: `results/chipseq_cutrun/deeptools`
- MultiQC: `results/chipseq_cutrun/multiqc_report.html`

#### ATAC-seq
- QC: `results/atacseq/qc_raw`, `results/atacseq/qc_trimmed`
- Trimming: `results/atacseq/trimmed_fastqs`
- BAMs: `results/atacseq/aligned_bams`, `results/atacseq/filtered_bams`
- Fragment QC: `results/atacseq/fragment_qc`
- BigWig: `results/atacseq/bigwigs`
- Peaks: `results/atacseq/peaks`
- MultiQC: `results/atacseq/multiqc_report.html`

### Parameter customization
Tool parameters are configured in `config/config.yaml`:
- `fastp.extra_args`
- `hisat2.extra_args`
- `kallisto.extra_args`
- `gene_body_coverage.*`
- `picard.markduplicates.extra_args`
- `bismark.extra_args`
- `bismark.threads`, `bismark.parallel`
- `bowtie2.extra_args`
- `sambamba.extra_args`, `sambamba.filter_threads`, `sambamba.sort_threads`
- `methyldackel_mbias.*`
- `methyldackel_extract.*`
- `deeptools.*`
- `macs3.*`
- `picard.java_opts`

Adjust these fields to tune trimming, mapping, filtering, normalization, and peak-calling behavior.

### FAQ and troubleshooting
- **Single-end in ChIP/CUT&RUN:** use `type: "SE"` and define only `R1`.
- **RNA-seq/WGBS/ATAC-seq single-end:** not implemented; use paired-end data.
- **WGBS validation fails early:** every WGBS sample must include `R1` and `R2`; the workflow aborts if it finds empty or single-end samples.
- **No samples are detected:** check `raw_fastqs_dir` and FASTQ naming patterns.
- **No bigWig subtraction appears for ChIP/CUT&RUN:** make sure input samples are named as `base_input_repX`.
- **Logs and result paths:** logs are written under `logs/<pipeline>`; for `chip_cr`, result files are stored under `results/chipseq_cutrun`.
- **Reference index permissions:** Bismark creates `Bisulfite_Genome/` next to `ref_genome`, so that location must be writable.

### License
See `LICENSE`.

<a id="espanol"></a>
## Español

Pipeline bioinformática modular y reproducible para análisis de datos NGS con Snakemake y gestión de dependencias con Pixi. Un único archivo de configuración selecciona uno de cuatro workflows soportados: RNA-seq, WGBS, ChIP-seq/CUT&RUN o ATAC-seq.

### Características principales
- Selección del workflow desde `config/config.yaml`.
- Entorno de software reproducible con Pixi y `pixi.lock`.
- Detección automática de muestras o declaración manual mediante `samples_info`.
- QC y trimming compartidos con FastQC y `fastp`.
- Indexación automática de referencias cuando el workflow la necesita.
- Reporte unificado con MultiQC y logs estructurados por paso.

### Estructura del repositorio
- `Snakefile`: punto de entrada principal que deriva al subworkflow seleccionado.
- `config/config.yaml`: archivo de configuración central.
- `workflows/`: módulos de workflow.
- `workflows/rules/`: reglas reutilizables.
- `src/utils.py`: utilidades de detección de muestras y manejo de entradas.
- `src/scripts/`: scripts auxiliares usados por workflows concretos.
- `backlog.md`: registro de cambios de implementación y mantenimiento.

### Requisitos
- Linux.
- [Pixi](https://pixi.sh/) instalado.
- Espacio suficiente en disco para FASTQ, BAM y resultados.
- Opcional: Graphviz si quieres renderizar DAGs con `dot`.

### Instalación rápida
```bash
git clone <url-del-repo>
cd Snakemake_pipeline_TFM
pixi install
```

### Configuración básica
La configuración principal vive en `config/config.yaml`. Como mínimo debes:

1. Definir el subworkflow en `pipeline`.
2. Especificar `raw_fastqs_dir` y/o `samples_info`.
3. Configurar los archivos de referencia y directorios de índices necesarios para el workflow elegido.

#### Selección de pipeline
```yaml
pipeline: "chip_cr"  # opciones: rnaseq, wgbs, chip_cr, atacseq
```

#### Entrada manual de muestras
```yaml
chip_cr:
  raw_fastqs_dir: "data/chip_fastqs"
  samples_info:
    sampleA:
      R1: "data/chip_fastqs/sampleA_R1.fastq.gz"
      R2: "data/chip_fastqs/sampleA_R2.fastq.gz"
      type: "PE"
```

#### Autodetección de muestras
Si `samples_info` está vacío o no se define, las muestras se detectan desde `raw_fastqs_dir` usando convenciones comunes de nombres:
- Paired-end: `_R1/_R2`, `_1/_2`, `.R1/.R2`.
- Single-end: `<sample>.fastq.gz` o `<sample>.fq.gz`.

#### Notas de compatibilidad por workflow
- `rnaseq`, `wgbs` y `atacseq` soportan actualmente solo muestras paired-end.
- `chip_cr` soporta muestras paired-end y single-end.
- La sustracción automática de bigWig en ChIP/CUT&RUN espera nombres como `liver_rep1` y `liver_input_rep1`.

### Ejecución del pipeline
Dry run:
```bash
pixi run snakemake -n --cores 1
```

Ejecución estándar:
```bash
pixi run snakemake --cores 8 --printshellcmds --latency-wait 60
```

Ejecutar hasta una regla concreta:
```bash
pixi run snakemake --cores 16 --until chip_cr_fastqc_raw_pe
```

Ejecutar una salida concreta:
```bash
pixi run snakemake --cores 4 results/rnaseq/kallisto/sample/abundance.tsv
```

Reanudar trabajos incompletos:
```bash
pixi run snakemake --cores 8 --rerun-incomplete
```

Generar un DAG (requiere Graphviz):
```bash
pixi run snakemake --dag --cores 1 | dot -Tpng > dag.png
```

### Tutoriales paso a paso

#### 1) ChIP-seq/CUT&RUN (paired-end con control input)
1. Coloca los FASTQ paired-end en un directorio como `data/chip_fastqs/`.
2. Nombra la muestra ChIP como `base_repX` y el input correspondiente como `base_input_repX` para habilitar la sustracción automática de bigWig.
3. Edita `config/config.yaml`:

```yaml
pipeline: "chip_cr"
chip_cr:
  raw_fastqs_dir: "data/chip_fastqs"
  samples_info:
    liver_rep1:
      R1: "data/chip_fastqs/liver_rep1_R1.fastq.gz"
      R2: "data/chip_fastqs/liver_rep1_R2.fastq.gz"
      type: "PE"
    liver_input_rep1:
      R1: "data/chip_fastqs/liver_input_rep1_R1.fastq.gz"
      R2: "data/chip_fastqs/liver_input_rep1_R2.fastq.gz"
      type: "PE"
  ref_genome: "ref/genome.fa"
  bowtie2_index_dir: "ref/bowtie2"
  gene_bed: "ref/genes.bed"
```

4. Ejecuta:
```bash
pixi run snakemake --cores 16 --printshellcmds --latency-wait 60
```

5. Resultados clave:
- `results/chipseq_cutrun/multiqc_report.html`
- `results/chipseq_cutrun/bigwigs/*.bw`
- `results/chipseq_cutrun/subtracted_bigwigs/*.subtracted.bw`
- `results/chipseq_cutrun/deeptools/heatmap.png`

#### 2) ChIP-seq/CUT&RUN (single-end)
Para datos single-end, define solo `R1` y usa `type: "SE"`:

```yaml
pipeline: "chip_cr"
chip_cr:
  raw_fastqs_dir: "data/chip_fastqs"
  samples_info:
    sampleSE:
      R1: "data/chip_fastqs/sampleSE.fastq.gz"
      type: "SE"
```

Las salidas usan el sufijo `_se`, por ejemplo `sampleSE_se.sorted.bam` y `sampleSE_se.bw`.

#### 3) RNA-seq (paired-end)
1. Coloca los FASTQ paired-end en `data/rnaseq_raw_fastqs/`.
2. Configura `transcriptome_fasta`, `kallisto_index`, `ref_genome` y `hisat2_index_dir`.
3. Ejemplo:

```yaml
pipeline: "rnaseq"
rnaseq:
  raw_fastqs_dir: "data/rnaseq_raw_fastqs"
  samples_info:
    Ctl_1:
      R1: "data/rnaseq_raw_fastqs/Ctl_1_R1.fq.gz"
      R2: "data/rnaseq_raw_fastqs/Ctl_1_R2.fq.gz"
  transcriptome_fasta: "ref/transcriptome.fa.gz"
  kallisto_index: "ref/kallisto.idx"
  ref_genome: "ref/genome.fa"
  hisat2_index_dir: "ref/hisat2"
  marked_bam_dir: "results/rnaseq/marked_bams"
  duplication_qc_dir: "results/rnaseq/duplication_qc"
  gene_body_coverage_dir: "results/rnaseq/gene_body_coverage"
  picard:
    java_opts: "-Xmx4g"
    markduplicates:
      extra_args: ""
  gene_body_coverage:
    refgene_bed: "ref/duumy.bed"
    minimum_length: 100
    format: "pdf"
```

4. Ejecuta y revisa:
- `results/rnaseq/kallisto/<sample>/abundance.tsv`
- `results/rnaseq/aligned_bams/<sample>_pe.sorted.bam`
- `results/rnaseq/marked_bams/<sample>_pe.markdup.bam`
- `results/rnaseq/duplication_qc/<sample>_pe.markdup.metrics.txt`
- `results/rnaseq/gene_body_coverage/all_samples.geneBodyCoverage.txt`
- `results/rnaseq/gene_body_coverage/all_samples.geneBodyCoverage.curves.pdf`
- `results/rnaseq/gene_body_coverage/all_samples.geneBodyCoverage.heatMap.pdf` cuando se analizan 3 o mas BAMs
- `results/rnaseq/bigwigs/<sample>_pe.bw`
- `results/rnaseq/multiqc_report.html`

El marcado de duplicados se guarda en BAMs separados solo para QC y no elimina lecturas de los BAMs de alineamiento originales.

Usa un modelo génico BED12 real y compatible con el ensamblado del genoma para ejecuciones reales. `ref/duumy.bed` es solo un archivo de prueba para el dry-run.

#### 4) WGBS (paired-end)
1. Coloca los FASTQ paired-end en `data/wgbs_raw_fastqs/`.
2. Configura `ref_genome`. Bismark creará `Bisulfite_Genome/` junto a esa referencia.
3. Ejemplo:

```yaml
pipeline: "wgbs"
wgbs:
  raw_fastqs_dir: "data/wgbs_raw_fastqs"
  samples_info:
    Ctr_1:
      R1: "data/wgbs_raw_fastqs/Ctr_1_R1.fastq.gz"
      R2: "data/wgbs_raw_fastqs/Ctr_1_R2.fastq.gz"
  ref_genome: "ref/danrer11_lambda.fa"
  bismark:
    threads: 66
    parallel: 8
```

4. Resultados clave:
- `results/wgbs/methyldackel/<sample>_CpG.methylKit`
- `results/wgbs/methyldackel_mergecontext/<sample>_CpG.bedGraph`
- `results/wgbs/multiqc_report.html`

5. Ajustes opcionales de recursos:
- `wgbs.bismark.threads` y `wgbs.bismark.parallel`
- `wgbs.sambamba.filter_threads` y `wgbs.sambamba.sort_threads`
- `wgbs.methyldackel_mbias.threads` y `wgbs.methyldackel_extract.threads`

#### 5) ATAC-seq (paired-end)
1. Coloca los FASTQ paired-end en `data/atacseq_raw_fastqs/`.
2. Configura `ref_genome` y `bowtie2_index_dir`.
3. Ejemplo:

```yaml
pipeline: "atacseq"
atacseq:
  raw_fastqs_dir: "data/atacseq_raw_fastqs"
  samples_info:
    atac_rep1:
      R1: "data/atacseq_raw_fastqs/atac_rep1_R1.fastq.gz"
      R2: "data/atacseq_raw_fastqs/atac_rep1_R2.fastq.gz"
      type: "PE"
  ref_genome: "ref/Danio_rerio.GRCz11.dna.primary_assembly.fa.gz"
  bowtie2_index_dir: "ref/bowtie2_atacseq"
  macs3:
    genome_size: "1.4e9"
```

4. Resultados clave:
- `results/atacseq/peaks/<sample>_peaks.narrowPeak`
- `results/atacseq/bigwigs/<sample>_pe.bw`
- `results/atacseq/fragment_qc/<sample>_pe.insert_size_metrics.txt`
- `results/atacseq/multiqc_report.html`

#### 6) Autodetección de muestras
Si quieres evitar `samples_info`, deja solo `raw_fastqs_dir`. El pipeline detectará muestras usando convenciones comunes de nombres de FASTQ. Ten en cuenta que `rnaseq`, `wgbs` y `atacseq` requieren entradas paired-end, mientras que `chip_cr` acepta datos paired-end y single-end.

### Salidas principales por pipeline

#### RNA-seq
- QC: `results/rnaseq/qc_raw`, `results/rnaseq/qc_trimmed`
- Trimming: `results/rnaseq/trimmed_fastqs`
- Kallisto: `results/rnaseq/kallisto/<sample>/abundance.tsv`
- HISAT2: `results/rnaseq/aligned_bams/<sample>_pe.sorted.bam`
- QC de duplicados: `results/rnaseq/marked_bams/<sample>_pe.markdup.bam` y `results/rnaseq/duplication_qc/<sample>_pe.markdup.metrics.txt`
- Cobertura del cuerpo genico: `results/rnaseq/gene_body_coverage/all_samples.geneBodyCoverage.txt` mas la curva y el mapa de calor opcional
- BigWig: `results/rnaseq/bigwigs/<sample>_pe.bw`
- MultiQC: `results/rnaseq/multiqc_report.html`

#### WGBS
- QC: `results/wgbs/qc`, `results/wgbs/qc_trimmed`
- Trimming: `results/wgbs/trimmed_fastqs`
- BAMs: `data/wgbs/raw_bams`, `data/wgbs/dedup_bams`, `data/wgbs/filtered_bams`, `data/wgbs/sorted_filtered_bams`
- M-bias: `results/wgbs/mbias`
- MethylDackel: `results/wgbs/methyldackel`, `results/wgbs/methyldackel_mergecontext`
- MultiQC: `results/wgbs/multiqc_report.html`

#### ChIP-seq/CUT&RUN
- QC: `results/chipseq_cutrun/qc_raw`, `results/chipseq_cutrun/qc_trimmed`
- Trimming: `results/chipseq_cutrun/trimmed_fastqs`
- BAMs: `results/chipseq_cutrun/aligned_bams`, `results/chipseq_cutrun/filtered_bams`
- BigWig: `results/chipseq_cutrun/bigwigs`, `results/chipseq_cutrun/subtracted_bigwigs`
- deepTools: `results/chipseq_cutrun/deeptools`
- MultiQC: `results/chipseq_cutrun/multiqc_report.html`

#### ATAC-seq
- QC: `results/atacseq/qc_raw`, `results/atacseq/qc_trimmed`
- Trimming: `results/atacseq/trimmed_fastqs`
- BAMs: `results/atacseq/aligned_bams`, `results/atacseq/filtered_bams`
- Fragment QC: `results/atacseq/fragment_qc`
- BigWig: `results/atacseq/bigwigs`
- Peaks: `results/atacseq/peaks`
- MultiQC: `results/atacseq/multiqc_report.html`

### Personalización de parámetros
Los parámetros de herramientas se configuran en `config/config.yaml`:
- `fastp.extra_args`
- `hisat2.extra_args`
- `kallisto.extra_args`
- `gene_body_coverage.*`
- `picard.markduplicates.extra_args`
- `bismark.extra_args`
- `bismark.threads`, `bismark.parallel`
- `bowtie2.extra_args`
- `sambamba.extra_args`, `sambamba.filter_threads`, `sambamba.sort_threads`
- `methyldackel_mbias.*`
- `methyldackel_extract.*`
- `deeptools.*`
- `macs3.*`
- `picard.java_opts`

Modifica estos campos para ajustar trimming, alineamiento, filtrado, normalización y peak calling.

### FAQ y resolución de problemas
- **Single-end en ChIP/CUT&RUN:** usa `type: "SE"` y define solo `R1`.
- **RNA-seq/WGBS/ATAC-seq single-end:** no está implementado; usa datos paired-end.
- **La validación de WGBS falla al inicio:** cada muestra WGBS debe incluir `R1` y `R2`; el workflow se detiene si detecta muestras vacías o single-end.
- **No se detectan muestras:** revisa `raw_fastqs_dir` y las convenciones de nombres de los FASTQ.
- **No aparece la sustracción de bigWig en ChIP/CUT&RUN:** asegúrate de nombrar los inputs como `base_input_repX`.
- **Rutas de logs y resultados:** los logs se escriben en `logs/<pipeline>`; para `chip_cr`, los resultados se guardan en `results/chipseq_cutrun`.
- **Permisos del índice de referencia:** Bismark crea `Bisulfite_Genome/` junto a `ref_genome`, así que esa ubicación debe tener permisos de escritura.

### Licencia
Consulta `LICENSE`.

![Pipeline](pipeline_wbg.svg)
