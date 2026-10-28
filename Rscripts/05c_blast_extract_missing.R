#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(Biostrings))

FASTA <- file.path(Sys.getenv("HOME"),
  "adesola_1/reanalysis/outputs/csv/lactobacillus_asvs_for_blast.fasta")
BLAST <- file.path(Sys.getenv("HOME"),
  "adesola_1/reanalysis/DHKUWNGA014-Alignment-HitTable.csv")
OUT   <- file.path(Sys.getenv("HOME"),
  "adesola_1/reanalysis/blast/missing_29.fasta")

fasta <- readDNAStringSet(FASTA)
blast <- read.csv(BLAST, header = FALSE, stringsAsFactors = FALSE)
have  <- unique(blast$V1)
missing_ids <- setdiff(names(fasta), have)
missing <- fasta[names(fasta) %in% missing_ids]

writeXStringSet(missing, OUT)
cat("Total ASVs in FASTA:  ", length(fasta), "\n")
cat("ASVs with nt hits:    ", length(have), "\n")
cat("Missing ASVs:         ", length(missing), "\n")
cat("Written to:           ", OUT, "\n\n")
cat("Missing ASV IDs:\n")
cat(names(missing), sep = "\n")
