#!/usr/bin/env Rscript
# 18b_build_heatmap.R — regenerate Supplementary Figure S1 (heatmap)
suppressPackageStartupMessages({
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("pheatmap not installed. Run: install.packages('pheatmap')")
  }
  library(pheatmap)
})

CSV  <- file.path(Sys.getenv("HOME"), "adesola_1/reanalysis/outputs/csv")
SUPP <- file.path(Sys.getenv("HOME"), "adesola_1/submission/Supplementary")
dir.create(SUPP, showWarnings = FALSE, recursive = TRUE)

gen <- read.csv(file.path(CSV, "audit_genus_per_sample.csv"),
                check.names = FALSE, row.names = 1)

totals <- sort(colSums(gen, na.rm = TRUE), decreasing = TRUE)
keep   <- names(totals)[1:min(20, length(totals))]
m      <- t(as.matrix(gen[, keep, drop = FALSE]))
m      <- m[order(-rowSums(m, na.rm = TRUE)), , drop = FALSE]
m[is.na(m)] <- 0

sample_order <- c("01PND","02PND","005PND","003PND","06DO","007NC","008NC")
m <- m[, sample_order]

ann_col <- data.frame(
  Group = factor(c("HIV+ baseline","HIV+ baseline","HIV+ 1-month",
                   "HIV+ baseline","HIV+ 1-month",
                   "Control","Control"),
                 levels = c("HIV+ baseline","HIV+ 1-month","Control")),
  row.names = sample_order
)
ann_colors <- list(
  Group = c("HIV+ baseline" = "#56B4E9",
            "HIV+ 1-month"  = "#D55E00",
            "Control"       = "#009E73")
)

palette <- colorRampPalette(c("#f7f7f7", "#fee8c8", "#fdbb84",
                              "#e34a33", "#b30000"))(100)

common_args <- list(
  mat               = m,
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  scale             = "none",
  color             = palette,
  border_color      = "white",
  fontsize          = 10,
  fontsize_row      = 9,
  fontsize_col      = 9,
  angle_col         = 45,
  main              = "Top 20 genera: relative abundance (%)",
  treeheight_row    = 40,
  treeheight_col    = 40,
  cellwidth         = 40,
  cellheight        = 22,
  annotation_col    = ann_col,
  annotation_colors = ann_colors
)

# PNG
png(file.path(SUPP, "Supplementary_Figure_S1_heatmap_top20_genera.png"),
    width = 2200, height = 2600, res = 220)
do.call(pheatmap, common_args)
dev.off()

# PDF
pdf(file.path(SUPP, "Supplementary_Figure_S1_heatmap_top20_genera.pdf"),
    width = 10, height = 12)
do.call(pheatmap, common_args)
dev.off()

cat("Wrote Supplementary Figure S1 (PNG + PDF)\n")
cat("Genera plotted (", nrow(m), "):\n")
cat(paste0("  ", rownames(m), collapse = "\n"), "\n")
