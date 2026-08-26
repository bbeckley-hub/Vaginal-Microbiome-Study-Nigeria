library(phyloseq)
library(DESeq2)
library(ggplot2)

ps <- readRDS("vaginal_microbiome_phyloseq.rds")
sample_data(ps)$hiv_status <- factor(sample_data(ps)$hiv_status, levels = c("Neg", "Pos"))
diagdds <- phyloseq_to_deseq2(ps, ~ hiv_status)
gm <- apply(counts(diagdds), 1, function(x) exp(mean(log(x[x>0]))))
diagdds <- estimateSizeFactors(diagdds, geoMeans = gm)
diagdds <- DESeq(diagdds, test = "Wald", fitType = "parametric")
res <- results(diagdds, cooksCutoff = FALSE)
res_df <- as.data.frame(res)
res_df$ASV <- rownames(res_df)
tax_df <- as.data.frame(tax_table(ps))
res_df$Genus <- tax_df[rownames(res_df), "Genus"]

res_df$diffexpressed <- "Not significant"
res_df$diffexpressed[res_df$padj < 0.1 & res_df$log2FoldChange > 0] <- "Enriched in HIV+"
res_df$diffexpressed[res_df$padj < 0.1 & res_df$log2FoldChange < 0] <- "Enriched in HIV-"

p_volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), col = diffexpressed)) +
  geom_point(alpha = 0.6) + theme_bw() +
  scale_color_manual(values = c("Enriched in HIV+" = "red", "Enriched in HIV-" = "blue", "Not significant" = "grey")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = "Volcano plot: HIV+ vs HIV-", x = "log2 Fold Change", y = "-log10 p-value")
ggsave("Figure4.png", p_volcano, width = 7, height = 6)
write.csv(res_df, "differential_abundance_full.csv", row.names = FALSE)
