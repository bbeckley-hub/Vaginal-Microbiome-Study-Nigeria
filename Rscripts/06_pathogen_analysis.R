library(phyloseq)
library(ggplot2)
library(reshape2)
library(dplyr)

ps <- readRDS("vaginal_microbiome_phyloseq.rds")

pathogens <- c("Gardnerella", "Prevotella", "Atopobium", "Sneathia",
               "Mobiluncus", "Mycoplasma", "Ureaplasma", "Megasphaera")
ps_path <- subset_taxa(ps, Genus %in% pathogens)
if(ntaxa(ps_path) > 0) {
  ps_path_rel <- transform_sample_counts(ps_path, function(x) x/sum(x)*100)
  path_melt <- psmelt(ps_path_rel)
  p_path <- ggplot(path_melt, aes(x = Sample, y = Abundance, fill = Genus)) +
    geom_bar(stat = "identity") + facet_grid(.~hiv_status, scales = "free_x") +
    theme_bw() + theme(axis.text.x = element_text(angle = 90)) +
    labs(y = "Relative abundance (%)", title = "Pathogenic genera")
  ggsave("Figure3.png", p_path, width = 10, height = 5)

  path_wide <- dcast(path_melt, Sample + hiv_status ~ Genus, value.var = "Abundance", fill = 0)
  path_long <- melt(path_wide, id.vars = c("Sample", "hiv_status"), variable.name = "Pathogen", value.name = "Abundance")
  path_long$Present <- ifelse(path_long$Abundance > 0, 1, 0)
  table1 <- path_long %>%
    group_by(hiv_status, Pathogen) %>%
    summarise(Count = sum(Present), .groups = "drop") %>%
    dcast(Pathogen ~ hiv_status, value.var = "Count", fill = 0)
  write.csv(table1, "Table1_pathogen_prevalence.csv", row.names = FALSE)
}
