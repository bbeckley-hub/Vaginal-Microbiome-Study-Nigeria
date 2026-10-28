#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(Biostrings) })

ROOT      <- file.path(Sys.getenv("HOME"), "adesola_1")
BLAST_DIR <- file.path(ROOT, "reanalysis/blast")
OUT_CSV   <- file.path(ROOT, "reanalysis/outputs/csv")

blast <- read.csv(file.path(BLAST_DIR, "lacto_nt_final.csv"),
                  stringsAsFactors = FALSE, check.names = FALSE)
blast$qseqid <- sub("^>", "", trimws(blast$qseqid))
blast$sseqid <- trimws(blast$sseqid)

sp <- read.delim(file.path(BLAST_DIR, "accession_species_master.tsv"),
                 header = FALSE, col.names = c("sseqid", "species"),
                 stringsAsFactors = FALSE, quote = "")
sp$sseqid  <- trimws(sp$sseqid)
sp$species <- trimws(sp$species)

blast <- merge(blast, sp, by = "sseqid", all.x = TRUE)
blast$species[is.na(blast$species) | blast$species %in% c("NA", "")] <- "unclassified"

seqtab <- readRDS(file.path(ROOT, "seqtab_final.rds"))
tax    <- readRDS(file.path(ROOT, "taxonomy_final.rds"))
genus  <- setNames(tax[, "Genus"], rownames(tax))[colnames(seqtab)]

lacto_idx  <- which(!is.na(genus) & genus == "Lactobacillus")
lacto_seqs <- colnames(seqtab)[lacto_idx]
lacto_ids  <- paste0("LactoASV_", seq_along(lacto_seqs), "_len", nchar(lacto_seqs))
query_len  <- setNames(nchar(lacto_seqs), lacto_ids)

blast_asvs   <- unique(blast$qseqid)
asvs_present <- intersect(lacto_ids, blast_asvs)
asvs_missing <- setdiff(lacto_ids, blast_asvs)

blast$qcov <- round(100 * (blast$qend - blast$qstart + 1) / query_len[blast$qseqid], 2)

# Named species pattern: genus + species, both alphabetic, not "uncultured"/"bacterium"
is_named <- function(x) {
  x <- trimws(as.character(x))
  !is.na(x) & x != "" &
  !grepl("uncultured|unclassified|bacterium|metagenome|environmental|synthetic|clone",
         x, ignore.case = TRUE) &
  grepl("^[A-Z][a-z]+ [a-z]+", x)
}

calls <- do.call(rbind, lapply(asvs_present, function(id) {
  hits <- blast[blast$qseqid == id, , drop = FALSE]
  if (nrow(hits) == 0) return(NULL)
  hits <- hits[order(-hits$bitscore), , drop = FALSE]
  hits$species_clean <- trimws(as.character(hits$species))

  # Identify absolute top hit (for transparency) and best named hit (for call)
  abs_top     <- hits[1, ]
  abs_top_unc <- !is_named(abs_top$species_clean)

  named_hits <- hits[is_named(hits$species_clean), , drop = FALSE]
  # Deduplicate named hits by species, keep best per species
  if (nrow(named_hits) > 0) {
    named_hits <- named_hits[!duplicated(named_hits$species_clean), , drop = FALSE]
  }

  call_val   <- "Lactobacillus sp. (unassigned)"
  reason_val <- "no named hit"

  if (nrow(named_hits) == 0) {
    if (abs_top_unc) reason_val <- "only uncultured hits"
    else             reason_val <- "no named species hit"
  } else {
    top    <- named_hits[1, ]
    second <- if (nrow(named_hits) >= 2) named_hits[2, ] else NULL
    prefix <- if (abs_top_unc) "[top hit uncultured; best named] " else ""

    if (!is.na(top$pident) && top$pident >= 99 && !is.na(top$qcov) && top$qcov >= 95) {
      second_ambig <- (!is.null(second) &&
                       !is.na(second$pident) && second$pident >= 99 &&
                       !is.na(second$qcov) && second$qcov >= 95 &&
                       !identical(second$species_clean, top$species_clean))
      if (second_ambig) {
        call_val   <- sprintf("Lactobacillus sp. (aff. %s / %s)",
                              top$species_clean, second$species_clean)
        reason_val <- sprintf("%stop two named both >=99%% (%.2f vs %.2f)",
                              prefix, top$pident, second$pident)
      } else {
        call_val   <- as.character(top$species_clean)
        reason_val <- sprintf("%s%.2f%% identity, qcov %.1f%%",
                              prefix, top$pident, top$qcov)
      }
    } else if (!is.na(top$pident) && top$pident >= 97) {
      call_val   <- sprintf("Lactobacillus sp. (nearest: %s)", top$species_clean)
      reason_val <- sprintf("%s%.2f%% identity - genus only", prefix, top$pident)
    } else {
      reason_val <- sprintf("%s%.2f%% identity - too low", prefix,
                            ifelse(is.na(top$pident), 0, top$pident))
    }
  }

  data.frame(
    asv_id         = as.character(id),
    seq_length     = as.integer(query_len[id]),
    abs_top_species = as.character(abs_top$species_clean),
    abs_top_pident  = as.numeric(abs_top$pident),
    best_named_species = if (nrow(named_hits) > 0) as.character(named_hits[1, "species_clean"]) else NA_character_,
    best_named_pident  = if (nrow(named_hits) > 0) as.numeric(named_hits[1, "pident"]) else NA_real_,
    best_named_qcov    = if (nrow(named_hits) > 0) as.numeric(named_hits[1, "qcov"]) else NA_real_,
    second_named_species = if (nrow(named_hits) >= 2) as.character(named_hits[2, "species_clean"]) else NA_character_,
    second_named_pident  = if (nrow(named_hits) >= 2) as.numeric(named_hits[2, "pident"]) else NA_real_,
    call           = call_val,
    reason         = reason_val,
    stringsAsFactors = FALSE
  )
}))

if (length(asvs_missing) > 0) {
  missing_df <- data.frame(
    asv_id = asvs_missing, seq_length = as.integer(query_len[asvs_missing]),
    abs_top_species = NA_character_, abs_top_pident = NA_real_,
    best_named_species = NA_character_, best_named_pident = NA_real_,
    best_named_qcov = NA_real_, second_named_species = NA_character_,
    second_named_pident = NA_real_,
    call = "Lactobacillus sp. (no nt hit)",
    reason = "no BLAST hit at 1e-10",
    stringsAsFactors = FALSE
  )
  calls <- rbind(calls, missing_df)
}
rownames(calls) <- NULL
write.csv(calls, file.path(OUT_CSV, "table_lactobacillus_species_calls.csv"),
          row.names = FALSE)
cat("=== Species call table written ===\n"); print(table(calls$call)); cat("\n")

# Species % per sample
rel <- sweep(seqtab, 1, rowSums(seqtab), "/") * 100
species_of_asv <- setNames(calls$call, calls$asv_id)
lacto_rel <- rel[, lacto_seqs, drop = FALSE]
colnames(lacto_rel) <- lacto_ids
all_calls <- sort(unique(calls$call))
species_mat <- sapply(all_calls, function(s) {
  idx <- intersect(names(species_of_asv)[species_of_asv == s], lacto_ids)
  if (!length(idx)) return(rep(0, nrow(lacto_rel)))
  rowSums(lacto_rel[, idx, drop = FALSE])
})
species_mat <- round(species_mat, 2)
species_df  <- data.frame(sample_name = rownames(rel), species_mat, check.names = FALSE)
write.csv(species_df, file.path(OUT_CSV, "audit_lactobacillus_species_per_sample.csv"),
          row.names = FALSE)
cat("=== Species % per sample ===\n"); print(species_df); cat("\n")

summary_df <- data.frame(
  call      = all_calls,
  n_asvs    = as.integer(table(factor(calls$call, levels = all_calls))),
  total_pct = round(colSums(species_mat), 2),
  row.names = NULL
)
write.csv(summary_df, file.path(OUT_CSV, "species_call_summary.csv"), row.names = FALSE)
cat("=== Species call summary ===\n"); print(summary_df); cat("\n")
cat("=== DONE ===\n")
