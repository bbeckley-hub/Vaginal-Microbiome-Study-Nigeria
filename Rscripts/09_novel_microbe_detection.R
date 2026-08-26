library(readxl)

blast_all <- read_excel("NBLAST_asvs.fasta_342 ASVs.xlsx", sheet = 1)
pct <- blast_all$`Percent identity`
pct <- ifelse(pct <= 1, pct * 100, pct)
novel <- pct < 97
cat("Total ASVs:", nrow(blast_all), "\n")
cat("Novel candidates (<97% identity):", sum(novel, na.rm = TRUE), "\n")

novel_df <- blast_all[novel, ]
write.csv(novel_df, "Supplementary_Table_S6_novel_candidates.csv", row.names = FALSE)
