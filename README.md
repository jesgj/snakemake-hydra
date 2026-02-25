# CTA_TFM_UOC - Pipeline modular para NGS (Snakemake + Pixi)

Pipeline bioinformática modular y reproducible para análisis de datos NGS con Snakemake y gestión de dependencias con Pixi. Soporta tres subpipelines seleccionables desde un único archivo de configuración: RNA-seq, WGBS y ChIP-seq/CUT&RUN.

## Características principales
- Selección de pipeline por configuración (sin tocar el código).
- Reproducibilidad garantizada con Pixi y `pixi.lock`.
- Detección automática de muestras o entrada manual por `samples_info`.
- QC y trimming centralizados (FastQC + fastp).
- Alineamientos e índices construidos automáticamente cuando aplica.
- Reporte unificado con MultiQC y logs por paso.

## Estructura del repositorio
- `Snakefile`: controlador que carga el subworkflow seleccionado.
- `config/config.yaml`: configuración central del pipeline y rutas.
- `workflows/`: workflows por subpipeline.
- `workflows/rules/`: reglas modulares (QC, trimming, BAM QC, etc.).
- `src/utils.py`: autodetección de muestras y lógica de entrada.
- `src/scripts/`: scripts auxiliares para WGBS.
- `DOC.md`: documentación metodológica detallada.

## Requisitos
- Linux.
- [Pixi](https://pixi.sh/) instalado.
- Espacio suficiente en disco para FASTQ, BAM y resultados.
- Opcional: Graphviz si quieres generar DAGs con `dot`.

## Instalación rápida
```bash
git clone <url-del-repo>
cd CTA_TFM_UOC
pixi install
```

## Configuración básica
La configuración vive en `config/config.yaml`. Lo mínimo es:

1. Seleccionar el subpipeline en `pipeline`.
2. Definir `raw_fastqs_dir` y/o `samples_info`.
3. Ajustar rutas a referencias (FASTA, índices, BED).

### Selección de pipeline
```yaml
pipeline: "chip_cr"  # opciones: rnaseq, wgbs, chip_cr
```

### Entrada de muestras (manual)
```yaml
chip_cr:
  raw_fastqs_dir: "data/chip_fastqs"
  samples_info:
    sampleA:
      R1: "data/chip_fastqs/sampleA_R1.fastq.gz"
      R2: "data/chip_fastqs/sampleA_R2.fastq.gz"
      type: "PE"
```

### Autodetección de muestras
Si `samples_info` está vacío o no se define, se detectan muestras desde `raw_fastqs_dir` con estas convenciones de nombres:
- Paired-end: `_R1/_R2`, `_1/_2`, `.R1/.R2`.
- Single-end: `<sample>.fastq.gz` o `<sample>.fq.gz`.

## Ejecución
Dry run:
```bash
pixi run snakemake -- -n
```

Ejecución estándar (ajusta cores):
```bash
pixi run snakemake --cores 8 --printshellcmds --latency-wait 60
```

Ejecutar hasta una regla concreta (ojo con el prefijo del módulo):
```bash
pixi run snakemake --cores 16 --until chip_cr_fastqc_raw_pe
```

Reanudar incompletos:
```bash
pixi run snakemake --cores 8 --rerun-incomplete
```

Generar DAG (requiere Graphviz):
```bash
pixi run snakemake --dag | dot -Tpng > dag.png
```

## Tutoriales paso a paso

### 1) ChIP-seq/CUT&RUN (paired-end con control de input)
1. Coloca FASTQ en un directorio, por ejemplo `data/chip_fastqs/`.
2. Nombra las muestras con el patrón `base_repX` y el input `base_input_repX` para activar la sustracción automática.
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

### 2) ChIP-seq/CUT&RUN (single-end)
Para single-end define solo `R1` y `type: "SE"`:
```yaml
pipeline: "chip_cr"
chip_cr:
  raw_fastqs_dir: "data/chip_fastqs"
  samples_info:
    sampleSE:
      R1: "data/chip_fastqs/sampleSE.fastq.gz"
      type: "SE"
```
Los outputs usan el sufijo `_se` (p. ej. `sampleSE_se.sorted.bam`, `sampleSE_se.bw`).

### 3) RNA-seq (paired-end)
1. Coloca FASTQ PE en `data/rnaseq_raw_fastqs/`.
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
```
4. Ejecuta y revisa:
- `results/rnaseq/kallisto/<sample>/abundance.tsv`
- `results/rnaseq/aligned_bams/<sample>_pe.sorted.bam`
- `results/rnaseq/multiqc_report.html`

### 4) WGBS (paired-end)
1. Coloca FASTQ PE en `data/wgbs_raw_fastqs/`.
2. Configura `ref_genome` (Bismark creará `Bisulfite_Genome` en el mismo directorio).
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
```
4. Resultados clave:
- `results/wgbs/methyldackel/<sample>_CpG.methylKit`
- `results/wgbs/methyldackel_mergecontext/<sample>_CpG.bedGraph`
- `results/wgbs/multiqc_report.html`

### 5) Autodetección de muestras
Si quieres evitar `samples_info`, deja solo `raw_fastqs_dir`. El pipeline detectará muestras con convenciones comunes. Ojo: `rnaseq` y `wgbs` requieren paired-end; `chip_cr` acepta PE o SE.

## Salidas principales por pipeline

### RNA-seq
- QC: `results/rnaseq/qc_raw`, `results/rnaseq/qc_trimmed`
- Trimming: `results/rnaseq/trimmed_fastqs`
- Kallisto: `results/rnaseq/kallisto/<sample>/abundance.tsv`
- HISAT2: `results/rnaseq/aligned_bams/<sample>_pe.sorted.bam`
- MultiQC: `results/rnaseq/multiqc_report.html`

### WGBS
- QC: `results/wgbs/qc`, `results/wgbs/qc_trimmed`
- Trimming: `results/wgbs/trimmed_fastqs`
- BAMs: `data/wgbs/raw_bams`, `data/wgbs/dedup_bams`, `data/wgbs/filtered_bams`
- MethylDackel: `results/wgbs/methyldackel`, `results/wgbs/methyldackel_mergecontext`
- MultiQC: `results/wgbs/multiqc_report.html`

### ChIP-seq/CUT&RUN
- QC: `results/chipseq_cutrun/qc_raw`, `results/chipseq_cutrun/qc_trimmed`
- Trimming: `results/chipseq_cutrun/trimmed_fastqs`
- BAMs: `results/chipseq_cutrun/aligned_bams`, `results/chipseq_cutrun/filtered_bams`
- BigWig: `results/chipseq_cutrun/bigwigs`, `results/chipseq_cutrun/subtracted_bigwigs`
- deepTools: `results/chipseq_cutrun/deeptools`
- MultiQC: `results/chipseq_cutrun/multiqc_report.html`

## Personalización de parámetros
Los parámetros de herramientas están en `config/config.yaml`:
- `fastp.extra_args`
- `hisat2.extra_args`
- `kallisto.extra_args`
- `bismark.extra_args`
- `sambamba.*_extra_args`
- `deeptools.*`
- `picard.java_opts`

Modifica esos campos para ajustar calidad, recorte, filtros, normalización de bigWig, etc.

## FAQ y resolución de problemas
- **Single-end en ChIP/CUT&RUN:** usa `type: "SE"` y solo `R1`.
- **RNA-seq/WGBS single-end:** no están implementados; usa paired-end.
- **No se detectan muestras:** revisa `raw_fastqs_dir` y el patrón de nombres.
- **No hay sustracción de bigWig:** asegúrate de nombrar input como `base_input_repX`.
- **Rutas de logs/resultados:** logs en `logs/<pipeline>`; para `chip_cr` los resultados están en `results/chipseq_cutrun`.
- **Permisos del índice:** Bismark crea `Bisulfite_Genome/` junto a `ref_genome`; necesitas permisos de escritura.

## Licencia
Ver `LICENSE`.

![Pipeline](pipeline_wbg.svg)