library(dada2); packageVersion("dada2")

path <- "E:/Anderson-BackUp/TP_SA/analysis/seqs/unpacked"
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_1.fq and SAMPLENAME_2.fq

fnFs <- sort(list.files(path, pattern="_1.fq", full.names = TRUE))
#fnFS
fnRs <- sort(list.files(path, pattern="_2.fq", full.names = TRUE))
#fnRs

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft = c(20,20),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out)

#Errors

errF <- learnErrors(filtFs, multithread=FALSE)

errR <- learnErrors(filtRs, multithread=FALSE)

#Dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)
gc()
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err=errF, multithread=FALSE)
gc()
dadaRs <- dada(derepRs, err=errR, multithread=FALSE)

dadaFs[[1]]

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)

dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

taxa <- assignTaxonomy(seqtab.nochim, "E:/Anderson-BackUp/Rice/silva_nr99_v138.1_train_set.fa.gz", multithread=FALSE)

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

saveRDS(seqtab.nochim, "seqtab_final.rds")
saveRDS(taxa, "taxonomy_genera.rds")


library(phyloseq); packageVersion("phyloseq")

library(ggplot2); packageVersion("ggplot2")

seqtab = readRDS("seqtab_final.rds")
taxa = readRDS("taxonomy_genera.rds")
mapfile <- "E:/Anderson-BackUp/TP_SA/map_tpsa.csv"

map = import_qiime_sample_data(mapfile)

# Building a phyloseq object

ps <- phyloseq(tax_table(taxa), otu_table(seqtab, taxa_are_rows=FALSE))

input = merge_phyloseq(ps,map) 

input
