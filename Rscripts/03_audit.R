#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(Biostrings) })

ROOT     <- file.path(Sys.getenv("HOME"), "adesola_1")
WORK     <- file.path(ROOT, "reanalysis")
OUT_CSV  <- file.path(WORK, "outputs", "csv")
OUT_LOG  <- file.path(WORK, "outputs", "logs")
dir.create(OUT_CSV, showWarnings = FALSE, recursive = TRUE)

log_path <- file.path(OUT_LOG, sprintf("14_audit_%s.log", format(Sys.time(), "%Y%m%d_%H%M%S")))
con <- file(log_path, "wt"); sink(con, split = TRUE)
on.exit({ sink(); close(con) }, add = TRUE)

cat("=========================================================\n")
cat("14_audit_clean.R |", format(Sys.time()), "\n")
cat("=========================================================\n\n")

seqtab <- readRDS(file.path(ROOT, "seqtab_final.rds"))
tax    <- readRDS(file.path(ROOT, "taxonomy_final.rds"))
meta   <- read.csv(file.path(ROOT, "metadata_corrected.csv"),
                   stringsAsFactors = FALSE, check.names = FALSE)
rownames(meta) <- meta$sample_name

meta_order <- c("01PND","02PND","005PND","003PND","06DO","007NC","008NC")
missing <- setdiff(meta_order, rownames(seqtab))
if (length(missing)) stop("samples missing from seqtab: ", paste(missing, collapse = ", "))
seqtab <- seqtab[meta_order, , drop = FALSE]
meta   <- meta[meta_order, , drop = FALSE]

tot   <- rowSums(seqtab)
rel   <- sweep(seqtab, 1, tot, "/") * 100
genus <- setNames(tax[, "Genus"], rownames(tax))[colnames(seqtab)]

shannon  <- apply(rel, 1, function(p) { p <- p[p > 0]; -sum((p/100)*log(p/100)) })
observed <- rowSums(seqtab > 0)

lacto_idx <- which(!is.na(genus) & genus == "Lactobacillus")
lacto_pct <- if (length(lacto_idx)) rowSums(rel[, lacto_idx, drop=FALSE]) else rep(0, nrow(rel))

top_genus <- sapply(seq_len(nrow(rel)), function(i) {
  g <- tapply(rel[i, ], genus, sum); names(sort(g, decreasing=TRUE))[1]
})

per_sample <- data.frame(
  sample_name = rownames(rel), participant_id = meta$participant_id,
  time_point = meta$time_point, hiv_status = meta$hiv_status,
  age = meta$age, total_reads = as.integer(tot),
  observed_asvs = as.integer(observed), shannon_index = round(shannon, 3),
  lactobacillus_pct = round(lacto_pct, 2), top_genus = top_genus, row.names = NULL
)
write.csv(per_sample, file.path(OUT_CSV, "audit_per_sample_summary.csv"), row.names = FALSE)
cat("--- Per-sample summary ---\n"); print(per_sample); cat("\n")

genus_mat <- sapply(unique(genus[!is.na(genus)]), function(g) {
  idx <- which(!is.na(genus) & genus == g)
  if (!length(idx)) return(rep(0, nrow(rel)))
  rowSums(rel[, idx, drop=FALSE])
})
genus_mat <- round(genus_mat, 2)
write.csv(data.frame(sample_name = rownames(rel), genus_mat, check.names=FALSE),
          file.path(OUT_CSV, "audit_genus_per_sample.csv"), row.names = FALSE)
cat("--- Top 6 genera overall ---\n")
print(round(sort(colMeans(genus_mat), decreasing = TRUE)[1:6], 2)); cat("\n")

pathogens <- c("Gardnerella","Prevotella","Atopobium","Sneathia",
               "Mobiluncus","Mycoplasma","Ureaplasma","Megasphaera")
path_df <- data.frame(sample_name = rownames(rel))
for (p in pathogens) {
  idx <- which(!is.na(genus) & genus == p)
  path_df[[p]] <- if (length(idx)) round(rowSums(rel[, idx, drop=FALSE]), 2) else 0
}
write.csv(path_df, file.path(OUT_CSV, "audit_pathogen_per_sample.csv"), row.names = FALSE)
cat("--- Pathogen % per sample ---\n"); print(path_df); cat("\n")

grp <- with(meta, ifelse(hiv_status == "Neg", "HIV_neg_control",
                         ifelse(art_months == 0, "HIV_pos_baseline", "HIV_pos_1mo")))
desc <- do.call(rbind, lapply(split(seq_len(nrow(rel)), grp), function(ix) {
  data.frame(group = grp[ix][1], n_samples = length(ix),
    n_participants = length(unique(meta$participant_id[ix])),
    lactobacillus_mean = round(mean(lacto_pct[ix]), 2),
    lactobacillus_min  = round(min(lacto_pct[ix]), 2),
    lactobacillus_max  = round(max(lacto_pct[ix]), 2),
    shannon_mean       = round(mean(shannon[ix]), 3))
}))
rownames(desc) <- NULL
write.csv(desc, file.path(OUT_CSV, "audit_group_means.csv"), row.names = FALSE)
cat("--- Descriptive group means ---\n"); print(desc); cat("\n")

traj <- per_sample[per_sample$participant_id %in% c("PND02","PND03"),
                   c("participant_id","time_point","lactobacillus_pct","shannon_index")]
traj <- traj[order(traj$participant_id, match(traj$time_point, c("baseline","1_month_post"))), ]
write.csv(traj, file.path(OUT_CSV, "audit_trajectories.csv"), row.names = FALSE)
cat("--- Paired trajectories ---\n"); print(traj); cat("\n")

if (length(lacto_idx)) {
  lacto_seqs <- colnames(seqtab)[lacto_idx]
  writeXStringSet(DNAStringSet(lacto_seqs, names = lacto_seqs),
                  file.path(OUT_CSV, "lactobacillus_asvs_for_blast.fasta"))
  cat("--- Lactobacillus ASVs exported for BLASTn:", length(lacto_seqs), "sequences ---\n")
}
cat("\n=== DONE ===\n")
