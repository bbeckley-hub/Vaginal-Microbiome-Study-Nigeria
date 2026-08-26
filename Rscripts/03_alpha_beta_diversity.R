library(phyloseq)
library(ggplot2)
library(vegan)

ps <- readRDS("vaginal_microbiome_phyloseq.rds")
ps_rel <- transform_sample_counts(ps, function(x) x/sum(x)*100)
sample_data(ps)$age_group <- ifelse(sample_data(ps)$age < 30, "<30", ">=30")

p_alpha <- plot_richness(ps, x = "age_group", measures = c("Observed", "Shannon")) +
  geom_boxplot() + theme_bw() + facet_wrap(~hiv_status)
ggsave("Figure2A.png", p_alpha, width = 8, height = 5)

ord <- ordinate(ps_rel, method = "PCoA", distance = "bray")
p_beta <- plot_ordination(ps_rel, ord, color = "hiv_status", shape = "age_group") +
  stat_ellipse() + theme_bw() +
  labs(title = "PCoA – Bray-Curtis (colored by HIV, shaped by age)")
ggsave("Figure2B.png", p_beta, width = 6, height = 5)

dist <- phyloseq::distance(ps_rel, method = "bray")
adonis_res <- adonis2(dist ~ hiv_status + age_group, data = data.frame(sample_data(ps)))
capture.output(print(adonis_res), file = "permanova_age_hiv.txt")
