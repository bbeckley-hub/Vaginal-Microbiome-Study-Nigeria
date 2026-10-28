#!/usr/bin/env Rscript
# 12_dada2_reprocess.R  (cutadapt, revcomp-aware, corrected REVCOMP)
suppressPackageStartupMessages({
  library(dada2); library(Biostrings)
})

ROOT         <- file.path(Sys.getenv("HOME"), "adesola_1")
WORK         <- file.path(ROOT, "reanalysis")
FASTQ_STAGE  <- file.path(WORK, "intermediate", "staged_fastq")
PRIMER_DIR   <- file.path(WORK, "intermediate", "primer_removed")
FILT_DIR     <- file.path(WORK, "intermediate", "filtered")
OUT_CSV      <- file.path(WORK, "outputs", "csv")
OUT_LOG      <- file.path(WORK, "outputs", "logs")
SILVA_TRAIN  <- file.path(ROOT, "silva_nr99_v138.1_train_set.fa.gz")
SILVA_SP     <- file.path(ROOT, "silva_species_assignment_v138.1.fa.gz")

for (d in c(FASTQ_STAGE, PRIMER_DIR, FILT_DIR, OUT_CSV, OUT_LOG))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

log_path <- file.path(OUT_LOG, sprintf("12_dada2_%s.log", format(Sys.time(), "%Y%m%d_%H%M%S")))
con <- file(log_path, "wt"); sink(con, split = TRUE)
on.exit({ sink(); close(con) }, add = TRUE)

cat("==============================================================\n")
cat("DADA2 reprocess | started", format(Sys.time()), "\n")
cat("R", R.version.string, "| DADA2", as.character(packageVersion("dada2")), "\n")
cat("==============================================================\n\n")

# --- Stage raw reads ---
cat("[0] Staging raw reads\n")
raw_gz <- sort(list.files(ROOT, pattern = "hifi_reads\\.fastq\\.gz$", full.names = TRUE))
stopifnot(length(raw_gz) == 7)
sample.names <- sub("-16s.*$", "", basename(raw_gz))
fnFs <- file.path(FASTQ_STAGE, paste0(sample.names, ".fastq"))
for (i in seq_along(raw_gz)) {
  if (!file.exists(fnFs[i]))
    system(sprintf("gunzip -c %s > %s", shQuote(raw_gz[i]), shQuote(fnFs[i])))
}
names(fnFs) <- sample.names
print(data.frame(sample = sample.names, file = basename(fnFs)))

count_reads <- function(f) {
  if (!file.exists(f)) return(0L)
  as.integer(system2("bash", c("-c", shQuote(sprintf("wc -l < %s", shQuote(f)))), stdout = TRUE)) %/% 4L
}
count_reads_gz <- function(f) {
  if (!file.exists(f)) return(0L)
  as.integer(system2("bash", c("-c", shQuote(sprintf("zcat %s | wc -l", shQuote(f)))), stdout = TRUE)) %/% 4L
}
n_raw <- sapply(fnFs, count_reads)

# --- Primer removal with cutadapt (revcomp-aware, corrected REVCOMP) ---
cat("\n[1] Primer removal with cutadapt: 27F / 1492R, --revcomp\n")
FWD     <- "AGAGTTTGATTCTGGCTCAG"
REVCOMP <- "AAGTCGTAACAAGGTAACC"   # correct 19-bp RC of GGTTACCTTGTTACGACTT

fn_noprimer <- file.path(PRIMER_DIR, paste0(sample.names, "_noprimer.fastq.gz"))
for (i in seq_along(fnFs)) {
  message("  -> ", sample.names[i])
  cmd <- sprintf(
    "cutadapt --revcomp -g %s -a %s -e 0.1 --minimum-length 500 --discard-untrimmed -j 0 -o %s %s > %s 2>&1",
    FWD, REVCOMP,
    shQuote(fn_noprimer[i]),
    shQuote(fnFs[i]),
    shQuote(file.path(OUT_LOG, paste0("cutadapt_", sample.names[i], ".log")))
  )
  status <- system(cmd)
  if (status != 0) stop(sprintf("cutadapt failed for %s", sample.names[i]))
}
n_noprimer <- sapply(fn_noprimer, count_reads_gz)

# --- filterAndTrim ---
cat("\n[2] filterAndTrim(): truncLen = 0, maxEE = 2\n")
filtFs <- file.path(FILT_DIR, paste0(sample.names, "_filtered.fastq.gz"))
names(filtFs) <- sample.names
filterAndTrim(
  fwd = fn_noprimer, filt = filtFs,
  truncLen = 0, maxEE = 2, truncQ = 2, maxN = 0,
  rm.phix = TRUE, compress = TRUE,
  multithread = TRUE, verbose = TRUE
)
n_filt <- sapply(filtFs, count_reads_gz)

# --- learnErrors + dada ---
cat("\n[3] learnErrors() + dada()\n")
err    <- learnErrors(filtFs, multithread = TRUE, verbose = TRUE)
dadaFs <- dada(filtFs, err = err, pool = "pseudo", multithread = TRUE, verbose = TRUE)
seqtab <- makeSequenceTable(dadaFs)
cat("ASVs before chimera removal:", ncol(seqtab), "\n")

# --- chimera removal ---
cat("\n[4] removeBimeraDenovo()\n")
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                    multithread = TRUE, verbose = TRUE)
cat("ASVs after chimera removal:", ncol(seqtab.nochim), "\n")

# --- tracking ---
track <- data.frame(
  sample        = sample.names,
  raw_reads     = as.integer(n_raw),
  after_primers = as.integer(n_noprimer),
  after_filter  = as.integer(n_filt),
  denoised      = as.integer(rowSums(seqtab)),
  nonchim       = as.integer(rowSums(seqtab.nochim)),
  pct_final     = round(100 * rowSums(seqtab.nochim) / as.integer(n_raw), 2),
  row.names     = NULL
)
write.csv(track, file.path(OUT_CSV, "read_tracking.csv"), row.names = FALSE)
cat("\n[5] Read tracking table:\n"); print(track)

# --- taxonomy ---
cat("\n[6] assignTaxonomy() + addSpecies()\n")
tax <- assignTaxonomy(seqtab.nochim, SILVA_TRAIN, multithread = TRUE, verbose = TRUE)
tax <- addSpecies(tax, SILVA_SP, verbose = TRUE)

keep <- !is.na(tax[, "Kingdom"]) & tax[, "Kingdom"] == "Bacteria"
n_drop <- sum(!keep)
cat("\n[7] Dropping non-Bacteria ASVs:", n_drop, "\n")
if (n_drop > 0) {
  seqtab.nochim <- seqtab.nochim[, keep, drop = FALSE]
  tax           <- tax[keep, , drop = FALSE]
}

cat("\n[8] Writing final objects\n")
saveRDS(seqtab.nochim, file.path(ROOT, "seqtab_final.rds"))
saveRDS(tax,           file.path(ROOT, "taxonomy_final.rds"))
writeXStringSet(DNAStringSet(getSequences(seqtab.nochim)),
                file.path(ROOT, "asv_sequences_final.fasta"))

cat("\nDONE. Log:", log_path, "\n")
sessionInfo()
