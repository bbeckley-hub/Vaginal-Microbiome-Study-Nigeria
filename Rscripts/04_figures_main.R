#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

ROOT     <- file.path(Sys.getenv("HOME"), "adesola_1")
OUT_CSV  <- file.path(ROOT, "reanalysis/outputs/csv")
OUT_FIG  <- file.path(ROOT, "reanalysis/outputs/figures")
dir.create(OUT_FIG, showWarnings = FALSE, recursive = TRUE)

sample_order <- c("01PND","02PND","005PND","003PND","06DO","007NC","008NC")

genus_df <- read.csv(file.path(OUT_CSV, "audit_genus_per_sample.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
genus_long <- pivot_longer(genus_df, cols = -sample_name,
                           names_to = "genus", values_to = "abundance")
genus_long$sample_name <- factor(genus_long$sample_name, levels = sample_order)

per_sample <- read.csv(file.path(OUT_CSV, "audit_per_sample_summary.csv"),
                       stringsAsFactors = FALSE)
per_sample$sample_name <- factor(per_sample$sample_name, levels = sample_order)
per_sample <- per_sample[order(per_sample$sample_name), ]

species_df <- read.csv(file.path(OUT_CSV, "audit_lactobacillus_species_per_sample.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
species_df$sample_name <- factor(species_df$sample_name, levels = sample_order)

# -----------------------------------------------------------------
# 16-colour palette: Okabe-Ito + extensions
# -----------------------------------------------------------------
pal16 <- c(
  "#E69F00","#56B4E9","#009E73","#F0E442",
  "#0072B2","#D55E00","#CC79A7","#000000",
  "#8B4513","#6A5ACD","#FF7F50","#2E8B57",
  "#B8860B","#8B008B","#20B2AA","#999999"
)

# =================================================================
# FIGURE 1A - all genera
# =================================================================
cat("[Fig 1A] Building...\n")

top_genera <- genus_long %>%
  group_by(genus) %>%
  summarise(total = sum(abundance, na.rm = TRUE)) %>%
  arrange(desc(total)) %>%
  head(15) %>%
  pull(genus)

genus_long$genus_grp <- ifelse(genus_long$genus %in% top_genera,
                               genus_long$genus, "Other")

legend_order <- genus_long %>%
  group_by(genus_grp) %>%
  summarise(total = sum(abundance, na.rm = TRUE)) %>%
  arrange(desc(total)) %>%
  pull(genus_grp)
genus_long$genus_grp <- factor(genus_long$genus_grp, levels = legend_order)
fill_map <- setNames(pal16[seq_along(legend_order)], legend_order)

p1a <- ggplot(genus_long, aes(x = sample_name, y = abundance, fill = genus_grp)) +
  geom_col(position = "stack", width = 0.85) +
  scale_fill_manual(values = fill_map) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100.5)) +
  labs(
    x = NULL,
    y = "Relative abundance (%)",
    title = "A. Genus composition, all samples",
    subtitle = paste("Adjacent pairs 02PND/005PND (PND02) and 003PND/06DO (PND03)",
                     "are baseline / post-ART from the same participant"),
    fill = "Genus"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    legend.text = element_text(size = 9, face = "italic"),
    legend.title = element_text(size = 10, face = "bold"),
    legend.key.size = unit(0.42, "cm"),
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 8.5, colour = "grey30"),
    plot.margin = margin(8, 18, 8, 8)
  ) +
  guides(fill = guide_legend(ncol = 1))

ggsave(file.path(OUT_FIG, "Figure1A_genus_composition.png"), p1a,
       width = 11.5, height = 7.5, dpi = 300)
ggsave(file.path(OUT_FIG, "Figure1A_genus_composition.pdf"), p1a,
       width = 11.5, height = 7.5)

# =================================================================
# FIGURE 1B - collapsed <1%
# =================================================================
cat("[Fig 1B] Building...\n")

genus_long_b <- genus_long %>%
  group_by(sample_name, genus) %>%
  summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(sample_name) %>%
  mutate(genus_grp = ifelse(abundance < 1, "Other (<1% per sample)", as.character(genus))) %>%
  group_by(sample_name, genus_grp) %>%
  summarise(abundance = sum(abundance), .groups = "drop")

legend_b <- genus_long_b %>%
  group_by(genus_grp) %>%
  summarise(total = sum(abundance)) %>%
  arrange(desc(total)) %>%
  pull(genus_grp)
genus_long_b$genus_grp <- factor(genus_long_b$genus_grp, levels = legend_b)
fill_map_b <- setNames(pal16[seq_along(legend_b)], legend_b)

p1b <- ggplot(genus_long_b, aes(x = sample_name, y = abundance, fill = genus_grp)) +
  geom_col(position = "stack", width = 0.85) +
  scale_fill_manual(values = fill_map_b) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100.5)) +
  labs(
    x = NULL,
    y = "Relative abundance (%)",
    title = "B. Genus composition (genera <1% grouped as 'Other')",
    fill = "Genus"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    legend.text = element_text(size = 9, face = "italic"),
    legend.title = element_text(size = 10, face = "bold"),
    legend.key.size = unit(0.5, "cm"),
    plot.title = element_text(size = 12, face = "bold"),
    plot.margin = margin(8, 18, 8, 8)
  ) +
  guides(fill = guide_legend(ncol = 1))

ggsave(file.path(OUT_FIG, "Figure1B_genus_grouped.png"), p1b,
       width = 9.5, height = 6.5, dpi = 300)
ggsave(file.path(OUT_FIG, "Figure1B_genus_grouped.pdf"), p1b,
       width = 9.5, height = 6.5)

# =================================================================
# FIGURE 1 composite
# =================================================================
if (requireNamespace("gridExtra", quietly = TRUE)) {
  library(gridExtra)
  png(file.path(OUT_FIG, "Figure1_composite.png"),
      width = 3450, height = 4200, res = 300)
  grid.arrange(p1a, p1b, ncol = 1, heights = c(1, 1))
  dev.off()
  cat("  Figure 1 composite written\n")
}

# =================================================================
# FIGURE 2 - trajectories with annotations
# =================================================================
cat("[Fig 2] Building...\n")

traj <- per_sample %>%
  filter(participant_id %in% c("PND02","PND03")) %>%
  mutate(time_point = factor(time_point,
                             levels = c("baseline","1_month_post"),
                             labels = c("Baseline", "1 month\npost-ART")))

p2 <- ggplot(traj, aes(x = time_point, y = lactobacillus_pct,
                       group = participant_id, colour = participant_id)) +
  geom_line(linewidth = 1.4) +
  geom_point(size = 5) +
  geom_text(aes(label = sprintf("%.1f%%", lactobacillus_pct)),
            vjust = -1.4, size = 4.2, show.legend = FALSE) +
  annotate("text", x = 2.4, y = 18,
           label = "PND02: lost Lactobacillus\n(Shuttleworthia-dominant)",
           hjust = 0, size = 3.3, colour = "#D55E00") +
  annotate("text", x = 2.4, y = 82,
           label = "PND03: recovered Lactobacillus\n(L. iners-dominant)",
           hjust = 0, size = 3.3, colour = "#0072B2") +
  scale_colour_manual(values = c("PND02" = "#D55E00", "PND03" = "#0072B2")) +
  scale_y_continuous(limits = c(-5, 105), breaks = seq(0, 100, 25)) +
  scale_x_discrete(expand = expansion(add = c(0.4, 1.4))) +
  labs(
    x = "Time point",
    y = expression(paste(italic("Lactobacillus"), " relative abundance (%)")),
    title = "Divergent 1-month ART trajectories in paired participants",
    colour = "Participant"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 12, face = "bold"))

ggsave(file.path(OUT_FIG, "Figure2_trajectories.png"), p2,
       width = 9.5, height = 6.5, dpi = 300)
ggsave(file.path(OUT_FIG, "Figure2_trajectories.pdf"), p2,
       width = 9.5, height = 6.5)

# =================================================================
# FIGURE 3 - species with shortened labels
# =================================================================
cat("[Fig 3] Building...\n")

species_long <- pivot_longer(species_df, cols = -sample_name,
                             names_to = "species", values_to = "abundance")

shorten_species <- function(x) {
  x <- gsub("Lactobacillus sp. \\(nearest: Lactobacillus iners\\)",
            "L. sp. (near L. iners)", x)
  x <- gsub("Lactobacillus sp. \\(nearest: Lactobacillus crispatus\\)",
            "L. sp. (near L. crispatus)", x)
  x <- gsub("Lactobacillus sp. \\(no nt hit\\)", "L. sp. (no BLAST hit)", x)
  x <- gsub("Lactobacillus sp. \\(unassigned\\)", "L. sp. (unassigned)", x)
  x <- gsub("Lactobacillus iners", "L. iners", x)
  x <- gsub("Lactobacillus crispatus", "L. crispatus", x)
  x
}
species_long$species_short <- shorten_species(species_long$species)
species_long <- species_long %>%
  group_by(sample_name, species_short) %>%
  summarise(abundance = sum(abundance), .groups = "drop")

species_order <- c("L. iners","L. crispatus",
                   "L. sp. (near L. iners)","L. sp. (near L. crispatus)",
                   "L. sp. (no BLAST hit)","L. sp. (unassigned)")
present <- species_order[species_order %in% unique(species_long$species_short)]
species_long$species_short <- factor(species_long$species_short, levels = present)
fill_map_species <- setNames(pal16[seq_along(present)], present)

p3 <- ggplot(species_long, aes(x = sample_name, y = abundance, fill = species_short)) +
  geom_col(position = "stack", width = 0.85) +
  scale_fill_manual(values = fill_map_species) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    x = NULL,
    y = expression(paste(italic("Lactobacillus"), " species abundance (%)")),
    title = expression(paste(italic("Lactobacillus"), " species composition per sample")),
    fill = "Species"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    legend.text = element_text(size = 9, face = "italic"),
    legend.title = element_text(size = 10, face = "bold"),
    legend.key.size = unit(0.5, "cm"),
    plot.title = element_text(size = 12, face = "bold"),
    plot.margin = margin(8, 15, 8, 8)
  ) +
  guides(fill = guide_legend(ncol = 1))

ggsave(file.path(OUT_FIG, "Figure3_species.png"), p3,
       width = 10, height = 6.5, dpi = 300)
ggsave(file.path(OUT_FIG, "Figure3_species.pdf"), p3,
       width = 10, height = 6.5)

cat("\n=== DONE ===\n")
cat("Figures in:", OUT_FIG, "\n")
print(list.files(OUT_FIG))
