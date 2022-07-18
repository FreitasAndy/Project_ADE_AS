## After dada2 pipeline ##

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

microbiome::summarize_phyloseq(input)

inputR = rarefy_even_depth(input, sample.size = 82270, replace = FALSE)

plot_richness(inputR, x="Treatment", measures=c("InvSimpson", "Observed"), color="Substrate")

input.clr = microbiome::transform(input, "clr")

input_ord = ordinate(input.clr, "NMDS" , "euclidean") 
p3 = plot_ordination(input.clr, input_ord, color = "Substrate")
p1 = p3 + geom_point(aes(shape = Period), size = 6, alpha = 0.8) +
  theme(legend.position = "bottom", legend.title = element_blank())

p1

#permanova
library(vegan)

df = as(sample_data(input.clr), "data.frame")
ds = phyloseq::distance(input.clr, method = "euclidean")
permanova = vegan::adonis2(ds ~ Substrate*Period, data = df, permutations = 999)

