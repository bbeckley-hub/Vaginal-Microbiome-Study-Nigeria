library(phyloseq)
library(readxl)
library(dplyr)

ps <- readRDS("vaginal_microbiome_phyloseq.rds")
lacto <- subset_taxa(ps, Genus == "Lactobacillus")
lacto_rel <- transform_sample_counts(lacto, function(x) x/sum(x)*100)
lacto_abund <- as.data.frame(otu_table(lacto_rel))
write.csv(lacto_abund, "lactobacillus_abundances.csv")

blast_lacto <- read_excel("NBLAST_asvs.fasta (104 ASVs).xlsx", sheet = "Sheet1")
colnames(blast_lacto)[1:6] <- c("SN", "ASV_ID", "Top_hit", "Pct_id", "Align_len", "Eval")
blast_lacto$Species <- gsub("^(Lactobacillus [a-z]+).*", "\\1", blast_lacto$Top_hit)
blast_lacto$Species[!grepl("Lactobacillus", blast_lacto$Species)] <- "unclassified"
blast_lacto$pident <- as.numeric(blast_lacto$Pct_id) * 100
blast_lacto <- blast_lacto[blast_lacto$pident >= 99, ]

species_table <- data.frame(
  ASV_Number = 1:nrow(blast_lacto),
  Species = blast_lacto$Species,
  Percent_identity = blast_lacto$pident,
  Top_BLAST_hit = blast_lacto$Top_hit
)
write.csv(species_table, "Supplementary_Table_S2_Lactobacillus_species.csv", row.names = FALSE)
