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
- 2026-03-17: Make RNA-seq gene body coverage optional via `gene_body_coverage.enabled` and skip its validation and tracked outputs when disabled.
- 2026-03-17: Add ATAC-seq Picard duplicate marking metrics and remove duplicate-flagged reads during filtered BAM generation.
- 2026-03-17: Refactor shared deepTools QC rules for fingerprint and correlation plots, and add ATAC-seq fingerprint/correlation QC on explicitly indexed filtered BAMs.
- 2026-03-17: Add RNA-seq deepTools multiBamSummary/plotCorrelation QC on explicitly indexed aligned BAMs, skipping the correlation step when fewer than 2 samples are available.
- 2026-03-17: Add configurable ATAC-seq peak calling with default MACS3 support and an optional Genrich path that uses a queryname-sorted filtered BAM intermediate.
- 2026-03-17: Harden RNA-seq HISAT2 index tracking, write per-sample HISAT2 summary files, and make MultiQC wait on config-aware RNA-seq report artifacts.

## Planned
- (Add upcoming work items here.)
