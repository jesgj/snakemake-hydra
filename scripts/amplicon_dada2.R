#!/usr/bin/env Rscript
# One cohort per comparable amplicon region; learn errors separately by batch/read.
suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(jsonlite))
args <- commandArgs(trailingOnly = TRUE)
cfg <- fromJSON(args[[1]], simplifyVector = FALSE)
threads <- as.integer(args[[2]])
mt <- if (threads > 1) threads else FALSE
out <- cfg$output_dir
dir.create(file.path(out, "filtered"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "errors"), recursive = TRUE, showWarnings = FALSE)
samples <- names(cfg$samples_info)
settings <- cfg$dada2
set.seed(settings$seed)
prepared <- function(s, mate) file.path(out, "primers", paste0(s, "_R", mate, ".fastq.gz"))
filtered <- function(s, mate) file.path(out, "filtered", paste0(s, "_R", mate, ".fastq.gz"))
track <- data.frame(sample = samples, mode = NA_character_, batch = NA_character_,
                    input = 0, primer_retained = 0, filtered = 0,
                    denoised_forward = NA_real_, denoised_reverse = NA_real_,
                    merged = NA_real_, nonchimeric = 0)
rownames(track) <- samples
tasks <- list()
for (s in samples) {
    info <- cfg$samples_info[[s]]
    f <- info$filter
    stats <- fromJSON(file.path(out, "primers", paste0(s, ".reads.json")))
    track[s, c("mode", "batch")] <- c(info$mode, info$batch)
    track[s, c("input", "primer_retained")] <- c(stats$input, stats$primer_retained)
    mates <- if (info$mode == "merge") c(1, 2) else if (info$mode == "reverse") 2 else 1
    tl <- unlist(f$trunc_len)[mates]
    left <- unlist(f$trim_left)[mates]
    ee <- unlist(f$max_ee)[mates]
    if (info$mode == "merge") {
        result <- filterAndTrim(prepared(s, 1), filtered(s, 1), prepared(s, 2), filtered(s, 2),
                                truncLen = tl, trimLeft = left, maxEE = ee, maxN = 0,
                                truncQ = f$trunc_q, minLen = f$min_len, rm.phix = TRUE,
                                # Preparation audits every pair before and after trimming.
                                # Match by position here: DADA2's automatic ID detection
                                # rejects valid archive headers such as @SRR123.1/1.
                                compress = TRUE, multithread = mt, matchIDs = FALSE)
    } else {
        result <- filterAndTrim(prepared(s, mates), filtered(s, mates),
                                truncLen = tl, trimLeft = left, maxEE = ee, maxN = 0,
                                truncQ = f$trunc_q, minLen = f$min_len, rm.phix = TRUE,
                                compress = TRUE, multithread = mt)
    }
    track[s, "filtered"] <- result[1, "reads.out"]
    if (result[1, "reads.out"] > 0) {
        for (mate in mates) {
            tasks[[length(tasks) + 1]] <- list(sample = s, mate = mate, batch = info$batch)
        }
    }
}
if (!length(tasks)) stop("No reads passed DADA2 filtering; inspect read lengths/quality and filter settings.")
groups <- unique(vapply(tasks, function(t) paste(t$batch, t$mate, sep = "\r"), character(1)))
calls <- list()
for (i in seq_along(groups)) {
    group <- Filter(function(t) paste(t$batch, t$mate, sep = "\r") == groups[[i]], tasks)
    files <- vapply(group, function(t) filtered(t$sample, t$mate), character(1))
    err <- learnErrors(files, nbases = settings$nbases, randomize = TRUE, multithread = mt)
    saveRDS(err, file.path(out, "errors", paste0("batch_", i, ".rds")))
    derep <- lapply(files, derepFastq)
    names(derep) <- vapply(group, function(t) t$sample, character(1))
    inferred <- dada(derep, err = err, pool = settings$pool, multithread = mt)
    # DADA2 unwraps a one-sample input list to a single dada object.
    if (inherits(inferred, "dada")) inferred <- list(inferred)
    for (j in seq_along(group)) {
        task <- group[[j]]
        if (is.null(calls[[task$sample]])) calls[[task$sample]] <- list()
        calls[[task$sample]][[as.character(task$mate)]] <- list(dada = inferred[[j]], derep = derep[[j]])
        col <- if (task$mate == 1) "denoised_forward" else "denoised_reverse"
        track[task$sample, col] <- sum(getUniques(inferred[[j]]))
    }
}
abundances <- list()
for (s in names(calls)) {
    info <- cfg$samples_info[[s]]
    call <- calls[[s]]
    if (info$mode == "merge") {
        merged <- mergePairs(call[["1"]]$dada, call[["1"]]$derep,
                             call[["2"]]$dada, call[["2"]]$derep,
                             minOverlap = settings$min_overlap, maxMismatch = settings$max_mismatch)
        track[s, "merged"] <- sum(merged$abundance)
        if (!sum(merged$abundance)) {
            stop(paste("No merged reads for", s, "— check post-trimming overlap. Explicitly select mode: forward or reverse to analyse one mate; no automatic fallback is performed."))
        }
        abundances[[s]] <- getUniques(merged)
    } else {
        mate <- if (info$mode == "reverse") "2" else "1"
        abundance <- getUniques(call[[mate]]$dada)
        if (info$mode == "reverse") names(abundance) <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(names(abundance))))
        abundances[[s]] <- abundance
    }
}
seqtab <- makeSequenceTable(abundances)
seqtab <- removeBimeraDenovo(seqtab, method = settings$chimera_method, multithread = mt)
if (!ncol(seqtab) || sum(seqtab) == 0) stop("No nonchimeric ASVs remain; inspect filtering, error learning, and amplicon comparability.")
# Restore zero-retention samples and deterministic lexical sequence ordering.
full <- matrix(0, nrow = length(samples), ncol = ncol(seqtab), dimnames = list(samples, sort(colnames(seqtab))))
full[rownames(seqtab), colnames(seqtab)] <- seqtab
track$nonchimeric <- rowSums(full)
ids <- paste0("ASV", seq_len(ncol(full)))
seqs <- data.frame(asv_id = ids, sequence = colnames(full))
colnames(full) <- ids
write.table(data.frame(sample = samples, full, check.names = FALSE), file.path(out, "asv_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(as.vector(rbind(paste0(">", ids), seqs$sequence)), file.path(out, "asv_sequences.fasta"))
write.table(track, file.path(out, "read_tracking.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
saveRDS(list(counts = full, sequences = seqs), file.path(out, "dada2.rds"))
write_json(list(dada2_version = as.character(packageVersion("dada2")), R_version = R.version.string,
                seed = settings$seed, error_learning_groups = lapply(seq_along(groups), function(i) list(file = paste0("errors/batch_", i, ".rds"), batch = strsplit(groups[[i]], "\r", fixed = TRUE)[[1]][1], mate = strsplit(groups[[i]], "\r", fixed = TRUE)[[1]][2])),
                resolved_config = cfg), file.path(out, "run_metadata.json"), pretty = TRUE, auto_unbox = TRUE)
mqc <- file.path(out, "amplicon_read_retention_mqc.tsv")
connection <- file(mqc, "w")
writeLines(c("# id: amplicon_read_retention", "# section_name: Amplicon read retention", "# plot_type: table", "# description: Counts are reads for single-mate modes and read pairs for merge mode; NA means not applicable."), connection)
write.table(track, connection, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA", col.names = TRUE)
close(connection)
