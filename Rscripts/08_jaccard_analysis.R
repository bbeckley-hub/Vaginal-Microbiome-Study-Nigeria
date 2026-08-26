library(phyloseq)

ps <- readRDS("vaginal_microbiome_phyloseq.rds")
ps_pa <- transform_sample_counts(ps, function(x) ifelse(x > 0, 1, 0))
jacc_dist <- phyloseq::distance(ps_pa, method = "jaccard")
jacc_mat <- as.matrix(jacc_dist)
write.csv(jacc_mat, "Supplementary_Table_S5_Jaccard_distances.csv")

sample_data(ps_pa)$group <- interaction(sample_data(ps_pa)$hiv_status, sample_data(ps_pa)$age_group)
groups <- unique(sample_data(ps_pa)$group)
within_means <- sapply(groups, function(g) {
  idx <- which(sample_data(ps_pa)$group == g)
  if(length(idx) > 1) mean(jacc_mat[idx, idx][lower.tri(jacc_mat[idx, idx])]) else NA
})
write.csv(data.frame(Group = groups, Mean_Jaccard = within_means), "jaccard_within_means.csv", row.names = FALSE)
