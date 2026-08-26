library(phyloseq)
library(ggplot2)
library(dplyr)

ps <- readRDS("vaginal_microbiome_phyloseq.rds")
lacto <- subset_taxa(ps, Genus == "Lactobacillus")
lacto_rel <- transform_sample_counts(lacto, function(x) x/sum(x)*100)
lacto_df <- psmelt(lacto_rel)
lacto_df <- lacto_df[, c("Sample", "time_point", "Abundance")]
paired <- lacto_df[lacto_df$Sample %in% c("01PND", "02PND", "003PND", "06DO"), ]
paired$participant <- ifelse(paired$Sample %in% c("01PND","02PND"), "PND02", "PND03")
paired_sum <- paired %>%
  group_by(participant, time_point) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

p_art <- ggplot(paired_sum, aes(x = time_point, y = Abundance, fill = participant)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) + theme_bw() +
  labs(y = "Relative abundance of Lactobacillus (%)", x = "Time point",
       title = "Effect of 1 month ART on Lactobacillus abundance") +
  theme(legend.position = "bottom")
ggsave("Figure5.png", p_art, width = 5, height = 4)
write.csv(paired_sum, "Figure5_data.csv", row.names = FALSE)
