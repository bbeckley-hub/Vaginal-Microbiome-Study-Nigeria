library(phyloseq)
library(ggplot2)
library(pheatmap)
library(dplyr)

ps <- readRDS("vaginal_microbiome_phyloseq.rds")
ps_rel <- transform_sample_counts(ps, function(x) x/sum(x)*100)
ps_genus <- tax_glom(ps_rel, taxrank = "Genus")
genus_melt <- psmelt(ps_genus)

p1a <- ggplot(genus_melt, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(~ hiv_status + time_point, scales = "free_x") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(y = "Relative abundance (%)", title = "Genus-level composition (all genera)")
ggsave("Figure1A.png", p1a, width = 14, height = 8, limitsize = FALSE)

genus_melt_grouped <- genus_melt %>%
  group_by(Sample) %>%
  mutate(Abundance2 = ifelse(Abundance < 1, "Other", Genus)) %>%
  group_by(Sample, Abundance2) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")
p1b <- ggplot(genus_melt_grouped, aes(x = Sample, y = Abundance, fill = Abundance2)) +
  geom_bar(stat = "identity", position = "stack") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(y = "Relative abundance (%)", title = "Genus composition (<1% grouped as 'Other')")
ggsave("Figure1B.png", p1b, width = 10, height = 6)

genus_avg <- tapply(genus_melt$Abundance, genus_melt$Genus, mean)
top20 <- names(sort(genus_avg, decreasing = TRUE)[1:20])
ps_top20 <- subset_taxa(ps_genus, Genus %in% top20)
otu_top20 <- as(otu_table(ps_top20), "matrix")
rownames(otu_top20) <- as(tax_table(ps_top20), "matrix")[, "Genus"]
pheatmap(otu_top20, scale = "none", main = "Top 20 genera (raw abundances)",
         filename = "Supplementary_Figure_S1.png", width = 8, height = 10)
