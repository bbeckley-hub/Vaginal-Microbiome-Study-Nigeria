library(phyloseq)

seqtab <- readRDS("seqtab_nochim.rds")
tax <- readRDS("taxonomy_with_species.rds")
metadata <- read.csv("metadata.csv", stringsAsFactors = FALSE)

rownames(seqtab) <- gsub("_filtered\\.fastq\\.gz$", "", rownames(seqtab))
rownames(metadata) <- metadata$sample_name
metadata$sample_name <- NULL
rownames(metadata)[rownames(metadata) == "006DO"] <- "06DO"

common <- intersect(rownames(seqtab), rownames(metadata))
seqtab <- seqtab[common, ]
metadata <- metadata[common, ]

ps <- phyloseq(otu_table(seqtab, taxa_are_rows = FALSE),
               sample_data(metadata),
               tax_table(tax))
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
saveRDS(ps, "vaginal_microbiome_phyloseq.rds")
