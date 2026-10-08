#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(jsonlite))
args <- commandArgs(trailingOnly = TRUE)
cfg <- fromJSON(args[[1]], simplifyVector = FALSE)
threads <- as.integer(args[[2]])
seqs <- readRDS(file.path(cfg$output_dir, "dada2.rds"))$sequences
t <- cfg$taxonomy
set.seed(cfg$dada2$seed)
if (t$enabled) {
    taxonomy <- assignTaxonomy(seqs$sequence, t$training_set, minBoot = t$min_boot,
                               tryRC = t$try_rc, multithread = if (threads > 1) threads else FALSE)
    table <- data.frame(asv_id = seqs$asv_id, taxonomy, check.names = FALSE)
    metadata <- list(status = "completed", reference = t$training_set, reference_version = t$reference_version,
                     reference_md5 = unname(tools::md5sum(t$training_set)), min_boot = t$min_boot, try_rc = t$try_rc)
} else {
    table <- data.frame(asv_id = seqs$asv_id, taxonomy_status = "not_requested")
    metadata <- list(status = "not_requested")
}
write.table(table, file.path(cfg$output_dir, "taxonomy.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write_json(metadata, file.path(cfg$output_dir, "taxonomy_metadata.json"), pretty = TRUE, auto_unbox = TRUE)
