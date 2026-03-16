# Backlog

Track notable changes and implementation work here.

## Implemented
- 2026-03-03: Add RNA-seq bigWig generation via deeptools bamCoverage (RPKM, binSize 10).
- 2026-03-03: Add RNA-seq config keys for bigWig output paths and deeptools settings.
- 2026-03-03: Wire bigWig outputs into the RNA-seq workflow targets and MultiQC inputs.
- 2026-03-05: Add explicit RNA-seq BAM indexing before bamCoverage and fail fast for empty or non-paired RNA-seq sample sets.
- 2026-03-05: Stabilize the WGBS workflow with PE-only validation, indexed filtered BAM dependencies, tracked Bismark/MethylDackel artifacts, and configurable WGBS runtime settings.
- 2026-03-05: Rewrite the README as a full bilingual English/Spanish guide and remove stale documentation references.
- 2026-03-16: Add RNA-seq RSeQC gene body coverage QC with configurable BED12 input, tracked outputs, and a dry-run placeholder BED file.
- 2026-03-16: Add RNA-seq Picard MarkDuplicates QC to mark duplicate reads without removal and report duplication metrics in MultiQC.

## Planned
- (Add upcoming work items here.)
