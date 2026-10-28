#!/usr/bin/env Rscript
# 13a_complete_pipeline.R — resume from filtered fastq, complete through taxonomy.
suppressPackageStartupMessages({
  library(dada2); library(Biostrings)
})

ROOT        <- file.path(Sys.getenv("HOME"), "adesola_1")
WORK        <- file.path(ROOT, "reanalysis")
FILT_DIR    <- file.path(WORK, "intermediate", "filtered")
OUT_CSV     <- file.path(WORK, "outputs", "csv")
OUT_LOG     <- file.path(WORK, "outputs", "logs")
SILVA_TRAIN <- file.path(ROOT, "silva_nr99_v138.1_train_set.fa.gz")
SILVA_SP    <- file.path(ROOT, "silva_species_assignment_v138.1.fa.gz")

log_path <- file.path(OUT_LOG, sprintf("13a_complete_%s.log", format(Sys.time(), "%Y%m%d_%H%M%S")))
con <- file(log_path, "wt"); sink(con, split = TRUE)
on.exit({ sink(); close(con) }, add = TRUE)

cat("============================================================\n")
cat("13a_complete_pipeline.R |", format(Sys.time()), "\n")
cat("R", R.version.string, "| DADA2", as.character(packageVersion("dada2")), "\n")
cat("============================================================\n\n")

# --- 0. Locate cached filtered fastq ---
cat("[0] Loading filtered fastq\n")
filtFs <- sort(list.files(FILT_DIR, pattern = "_filtered\\.fastq\\.gz$", full.names = TRUE))
if (length(filtFs) != 7) stop("expected 7 filtered files, found ", length(filtFs))
sample.names <- sub("_filtered\\.fastq\\.gz$", "", basename(filtFs))
names(filtFs) <- sample.names
print(data.frame(sample = sample.names, file = basename(filtFs)))
stopifnot(all(file.exists(filtFs)))

# --- 1. Learn errors (all 7 samples, not just first 2) ---
cat("\n[1] learnErrors(): nbases = 1e10 (forces all samples)\n")
set.seed(42)
err <- learnErrors(filtFs, nbases = 1e10, multithread = TRUE, verbose = TRUE)

# --- 2. Denoise ---
cat("\n[2] dada(): pool = 'pseudo'\n")
dadaFs <- dada(filtFs, err = err, pool = "pseudo", multithread = TRUE, verbose = TRUE)
seqtab <- makeSequenceTable(dadaFs)
cat("ASVs before chimera removal:", ncol(seqtab), "\n")

# --- 3. Chimera removal ---
cat("\n[3] removeBimeraDenovo()\n")
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                    multithread = TRUE, verbose = TRUE)
cat("ASVs after chimera removal:", ncol(seqtab.nochim), "\n")

# --- 4. CHECKPOINT: save seqtab NOW ---
cat("\n[4] CHECKPOINT: writing seqtab_final.rds (pre-taxonomy)\n")
saveRDS(seqtab.nochim, file.path(ROOT, "seqtab_final.rds"))
writeXStringSet(DNAStringSet(getSequences(seqtab.nochim)),
                file.path(ROOT, "asv_sequences_final.fasta"))
cat("  seqtab_final.rds written\n")
cat("  asv_sequences_final.fasta written\n")

# --- 5. Read tracking (updated with new chimera pass) ---
track <- data.frame(
  sample   = sample.names,
  denoised = as.integer(rowSums(seqtab)),
  nonchim  = as.integer(rowSums(seqtab.nochim)),
  pct_of_denoised = round(100 * rowSums(seqtab.nochim) / rowSums(seqtab), 2),
  row.names = NULL
)
write.csv(track, file.path(OUT_CSV, "read_tracking_postchim.csv"), row.names = FALSE)
cat("\n[5] Post-chimera retention:\n"); print(track)

# --- 6. Taxonomy (LONG step, ~5-15 min) ---
cat("\n[6] assignTaxonomy() + addSpecies() — this takes 5-15 minutes, do not interrupt\n")
tax <- assignTaxonomy(seqtab.nochim, SILVA_TRAIN, multithread = TRUE, verbose = TRUE)
cat("  assignTaxonomy done\n")
tax <- addSpecies(tax, SILVA_SP, verbose = TRUE)
cat("  addSpecies done\n")

# --- 7. Drop non-Bacteria ---
keep   <- !is.na(tax[, "Kingdom"]) & tax[, "Kingdom"] == "Bacteria"
n_drop <- sum(!keep)
cat("\n[7] Dropping non-Bacteria ASVs:", n_drop, "\n")
if (n_drop > 0) {
  seqtab.nochim <- seqtab.nochim[, keep, drop = FALSE]
  tax           <- tax[keep, , drop = FALSE]
}

# --- 8. Save everything ---
cat("\n[8] Writing final objects\n")
saveRDS(seqtab.nochim, file.path(ROOT, "seqtab_final.rds"))
saveRDS(tax,           file.path(ROOT, "taxonomy_final.rds"))
writeXStringSet(DNAStringSet(getSequences(seqtab.nochim)),
                file.path(ROOT, "asv_sequences_final.fasta"))

cat("\nDONE.\n")
cat("  seqtab_final.rds:",    nrow(seqtab.nochim), "samples x", ncol(seqtab.nochim), "ASVs\n")
cat("  taxonomy_final.rds:",  nrow(tax), "ASVs\n")
cat("  Log:", log_path, "\n")
sessionInfo()
