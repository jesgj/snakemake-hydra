# Small real-data fixtures

These fixtures cover all four pipeline values and the distinct ChIP-seq and
CUT&RUN assays supported by `chip_cr`. They test software behavior, not biological
accuracy or production-scale performance. No simulated reads are used.

| Assay | Real source | Retained size | Reference |
| --- | --- | --- | --- |
| RNA-seq | [GSE110004](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110004), SRR6357070/71/72 | 10,000 pairs each, 3 biological replicates | Yeast R64-1-1 chromosome I and matching transcriptome/GTF |
| ATAC-seq | [GSE66386](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE66386), SRR1822153/54 | 20,000 pairs each, 2 biological replicates | Yeast R64-1-1 genome (~12 Mb) |
| ChIP-seq | [SRR5204807](https://www.ncbi.nlm.nih.gov/sra/?term=SRR5204807) and [SRR5204809](https://www.ncbi.nlm.nih.gov/sra/?term=SRR5204809), Spt5 IP1 and matched Input1 | 20,000 pairs each | Yeast R64-1-1 genome and BED |
| CUT&RUN | [GSM2413650](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2413650), SRR5071545, Cse4 20s | 20,000 pairs | Yeast R64-1-1 genome and BED |
| WGBS | [GSE65685 / PRJNA274767](https://www.ebi.ac.uk/ena/browser/view/PRJNA274767), SRR1792821, methylC-seq F2__89 | 20,000 pairs | Arabidopsis TAIR10 chromosome 5 (~27 Mb), Ensembl Plants release 62 |

The complete TAIR10 genome (~120 Mb, including organelles) is also downloaded to
`references/tair10/genome.fa`. Use `wgbs_full.yaml` to align the same subset across
the full genome. Full and chromosome-only WGBS FASTAs are kept in separate
directories because Bismark indexes all FASTAs in the reference directory.
The complete yeast genome is shared by ATAC-seq, ChIP-seq, and CUT&RUN.

## Sampling and provenance

Run from the repository root:

```bash
pixi run python test_data/prepare.py
```

`datasets.json` records source URLs, pinned nf-core repository revisions, the
retained pair counts, and seed `20261002`. RNA-seq, ATAC-seq, and ChIP-seq use
small real-data subsets already hosted by nf-core; the RNA-seq source is enriched
for chromosome I. See the source [RNA-seq sampling method](https://github.com/nf-core/test-datasets/blob/626c8fab639062eade4b10747e919341cbf9b41a/README.md)
and [ATAC-seq sampling method](https://github.com/nf-core/test-datasets/blob/cd022b097372b078a68d8afadb172ad7342fd91f/README.md).

WGBS and CUT&RUN use their entire original libraries as the sampling population.
The preparer downloads one source pair into a temporary directory, validates mate
identifiers and FASTQ records, and uniformly samples pairs without replacement
using reservoir sampling. R1/R2 remain synchronized and retain their original
sequences and qualities. Both temporary source files are deleted when that sample
finishes, including on failure. Full source FASTQs are never pipeline inputs.
The script rejects downloads larger than 100 MiB per file.

Only the sampled FASTQs and small references remain under this directory.
`provenance.json` records actual source/retained counts, source and output SHA-256
checksums, compressed sizes, and sampling scope. WGBS source MD5s are also checked
against ENA metadata. Stored gzip outputs have deterministic headers. RAM for
sampling scales with the retained count, not library size.

This is uniform sampling within the specified source population. The nf-core
subsets and chromosome-restricted references introduce deliberate limitations:
WGBS reads from other chromosomes will not map; CUT&RUN spike-in reads will not
map to yeast; RNA-seq covers chromosome I. Small libraries can yield empty peak
sets or sparse QC plots. Do not compare their metrics with production benchmarks.

## Pipeline modes and bounded execution

The preparer generates independent configs in `test_data/configs/`, keeping all
four required module mappings. Outputs and indexes use `test_data/work/<variant>`.
WGBS logs use its configured test directory; other workflows have fixed
`logs/<pipeline>` paths. Use `--replace-workflow-config` on every fixture invocation: without it,
Snakemake recursively merges production sample mappings into the test config.
This option is available in the committed Snakemake environment. Existing
production configuration is not edited.

| Config name | Branch exercised |
| --- | --- |
| `rnaseq.yaml` | PE RNA-seq, three samples, correlation and enabled gene body coverage/heatmap |
| `wgbs.yaml` | PE WGBS, Bismark alignment and methylation extraction |
| `wgbs_full.yaml` | Same WGBS reads against complete TAIR10, with separate results |
| `chipseq.yaml` | PE ChIP-seq, matched input subtraction and correlation |
| `chipseq_se.yaml` | SE code path using just R1 of the same real ChIP libraries |
| `cutrun.yaml` | Actual PE CUT&RUN, with a 20 bp minimum to retain its 25 bp reads |
| `atacseq_macs3.yaml` | PE ATAC-seq, two samples, correlation and MACS3 |
| `atacseq_genrich.yaml` | Same ATAC inputs, with Genrich |

The SE variant exercises software handling of single reads; it is not a separate
biological SE experiment. RNA-seq BED12 is derived from the matching real GTF.
RNA-seq bootstraps are disabled for speed; Picard heaps are limited to 1 GB, and
ATAC MACS3 uses a yeast genome size. Configs inherit other current tool options
from `config/config.yaml`; regenerate them after changing defaults:

```bash
pixi run python test_data/prepare.py --configs-only
```

Set `HYDRA_CONFIG` to the same fixture config passed with `--configfile` so
Snakemake subprocesses also use the fixture instead of a production local config.
Inspect the DAG first; choose one variant per invocation:

```bash
HYDRA_CONFIG=test_data/configs/wgbs.yaml XDG_CACHE_HOME=/tmp/hydra-cache pixi run snakemake -n --cores 2 --replace-workflow-config --configfile test_data/configs/wgbs.yaml
HYDRA_CONFIG=test_data/configs/atacseq_genrich.yaml XDG_CACHE_HOME=/tmp/hydra-cache pixi run snakemake -n --cores 2 --replace-workflow-config --configfile test_data/configs/atacseq_genrich.yaml
```

For a real smoke run, use a single job at a time and two scheduler cores:

```bash
HYDRA_CONFIG=test_data/configs/wgbs.yaml XDG_CACHE_HOME=/tmp/hydra-cache pixi run snakemake --cores 2 --jobs 1 --replace-workflow-config --configfile test_data/configs/wgbs.yaml --printshellcmds --latency-wait 60
```

These options limit scheduled concurrency, not every tool's physical memory or
subprocess count. In particular, Bismark starts multiple aligner processes; do
not infer a hard OS CPU/RAM cap from `--cores`. No execution-time promise is made.
Keep runs sequential. `validation.json` records preparation and DAG checks;
`execution.json` records actual execution checks, artifact hashes, tool versions,
and cleanup. Dry-run success alone does not confirm tool execution.

After Snakemake exits successfully, verify the final report, BAM integrity, mapped
read counts, and peak record structure, record hashes, then delete that run:

```bash
pixi run python test_data/finish_run.py wgbs --cores 2
```

Use the variant name and actual scheduler core count. This helper refuses cleanup
if required checks fail. It deletes the variant work directory and generated
WGBS reference indexes, while retaining original reference FASTAs and FASTQs.

Generated FASTQs, references, configs, and results are locally ignored by Git;
the catalog, preparer, and provenance records remain reviewable. The preparer
refuses to overwrite existing sampled FASTQs; move an old fixture set aside if
you want to regenerate it.

## Checks performed

All 150,000 retained pairs passed FASTQ structure, synchronization, count, and
checksum checks. A sampling check on real reads confirmed that identical seeds
give identical outputs and a changed seed changes the selected subset.

All eight variant DAGs construct successfully. All eight complete execution runs
passed: RNA-seq (including three-sample gene body coverage), WGBS chromosome 5
and full TAIR10, ChIP PE and SE, actual CUT&RUN, and ATAC MACS3 and Genrich.
ChIP SE and ATAC Genrich ran with one core; the others used two cores. Every
run used one scheduled job at a time. Execution records are in
`execution.json`; each successful run was checked before its outputs were deleted.
The runs exposed three defects, now fixed: shared SE wildcard constraints,
ChIP/CUT&RUN duplicate-marking streaming, and missing ATAC read groups for Picard.
ChIP/CUT&RUN and ATAC filtering pipe directly into Sambamba sorting with an
allocation derived from effective job threads. One-core runs process the stages
sequentially using temporary files. Snakemake lint still reports existing layout,
path, and per-rule environment recommendations; the repository uses Pixi.

These small subsets test execution and file interfaces. Their coverage and peak
counts are unsuitable for biological conclusions.

After merging remote commit `4459588`, all eight fixture DAGs passed again and
the run-aware ATAC MACS3 workflow completed, including consensus featureCounts,
FRiP, read-group checks, and MultiQC. See `merge_validation.json` and the last
entry in `execution.json`. Its outputs were verified and deleted.
