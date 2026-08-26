library(dada2)
library(Biostrings)

path <- "~/adesola_1"
setwd(path)

fnFs <- sort(list.files(path, pattern = ".fastq.gz$", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "-"), `[`, 1)
sample.names <- gsub(".fastq.gz$", "", sample.names)

filtFs <- file.path(path, "filtered", paste0(sample.names, "_filtered.fastq.gz"))
filterAndTrim(fnFs, filtFs, truncLen = 1450, maxEE = 2, truncQ = 2,
              rm.phix = TRUE, compress = TRUE, multithread = TRUE)

err <- learnErrors(filtFs, multithread = TRUE)
dadaFs <- dada(filtFs, err = err, multithread = TRUE)
seqtab <- makeSequenceTable(dadaFs)
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE)
saveRDS(seqtab.nochim, "seqtab_nochim.rds")

tax <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138.1_train_set.fa.gz", multithread = TRUE)
tax <- addSpecies(tax, "silva_species_assignment_v138.1.fa.gz")
saveRDS(tax, "taxonomy_with_species.rds")

asv_seqs <- DNAStringSet(getSequences(seqtab.nochim))
writeXStringSet(asv_seqs, "asv_sequences.fasta")
