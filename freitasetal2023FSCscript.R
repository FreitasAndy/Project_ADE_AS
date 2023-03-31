---
title: "ADEs and Agricolae Soils"
author: "Anderson Freitas"
date: "Jul/14/2022"
output: html_notebook
---
  
#How we made our analysis and figures:
  
## Dirty work
  
  #Libraries
  
library(hues)
library(knitr)
library(readr)
library(vegan)
library(dplyr)
library(ggpubr)
library(ggplot2)
library(magrittr)
library(microeco)
library(corrplot)
library(phyloseq)
library(agricolae)
library(tidyverse)
library(microbiome)
library(factoextra)
library(dunn.test)

#############################################

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- plyr::rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

#############################################

theme_set(theme_bw())
pal = "Set1"
scale_colour_discrete <-  function(palname=pal, ...){
  scale_colour_brewer(palette=palname, ...)
}
scale_fill_discrete <-  function(palname=pal, ...){
  scale_fill_brewer(palette=palname, ...)
}

ze=c("#79b643", "#bb59bd", "#52ac7a", "#cf426b", "#4dadcf","#cf5234",
     "#7475cb", "#cea13d", "#be6c91", "#748037", "#be774e")

map_tpsa <- read_delim("map_tpsa.csv", delim = "\t",
                       escape_double = FALSE, na = "NA", trim_ws = TRUE)
map = data.frame(map_tpsa)

seqtab = readRDS("dada2/seqtab_final.rds")
taxa = readRDS("dada2/taxonomy_genera.rds")
mapfile <- "map_tpsa.csv"

map = import_qiime_sample_data(mapfile)

# Building a phyloseq object

ps <- phyloseq(tax_table(taxa), otu_table(seqtab, taxa_are_rows=FALSE))

input = merge_phyloseq(ps,map) 

input

microbiome::summarize_phyloseq(input)

#Rarefaction
inputR = rarefy_even_depth(input, sample.size = 82270, replace = FALSE)

#Subsetting
Cedrela  = subset_samples(inputR, Tree == "C. fissilis")
Pelto    = subset_samples(inputR, Tree == "P. dubium")
Cecropia = subset_samples(inputR, Tree == "C. pachystachya")
Bulk     = subset_samples(inputR, Tree == "Bulk")

#Creating clr objects
inputR.clr = microbiome::transform(inputR, "clr")

## for Cedrela ##

otu_ced = as.data.frame(t(otu_table(Cedrela)))
tax_ced = as.data.frame(tax_table(Cedrela))
map_ced = data.frame(sample_data(Cedrela))

dataset.ced <- microtable$new(sample_table = map_ced,
                              otu_table    = otu_ced,
                              tax_table    = tax_ced)
dataset.ced$tidy_dataset()
dataset.ced

## for Peltophorum ##

otu_pel = as.data.frame(t(otu_table(Pelto)))
tax_pel = as.data.frame(tax_table(Pelto))
map_pel = data.frame(sample_data(Pelto))

dataset.pel <- microtable$new(sample_table = map_pel,
                              otu_table    = otu_pel,
                              tax_table    = tax_pel)
dataset.pel$tidy_dataset()
dataset.pel

## for Cecropia ##

otu_cec = as.data.frame(t(otu_table(Cecropia)))
tax_cec = as.data.frame(tax_table(Cecropia))
map_cec = data.frame(sample_data(Cecropia))

dataset.cec <- microtable$new(sample_table = map_cec,
                              otu_table    = otu_cec,
                              tax_table    = tax_cec)
dataset.cec$tidy_dataset()
dataset.cec

## for bulk ##

otu_bk = as.data.frame(t(otu_table(Bulk)))
tax_bk = as.data.frame(tax_table(Bulk))
map_bk = data.frame(sample_data(Bulk))

dataset.bk <- microtable$new(sample_table  = map_bk,
                             otu_table    = otu_bk,
                             tax_table    = tax_bk)
dataset.bk$tidy_dataset()
dataset.bk

## for the whole dataset ##

otu_tp = as.data.frame(t(otu_table(inputR)))
tax_tp = as.data.frame(tax_table(inputR))
map_tp = data.frame(sample_data(inputR))

# Let's create a microtable object with more information
dataset <- microtable$new(sample_table = map_tp,
                          otu_table    = otu_tp,
                          tax_table    = tax_tp)
dataset

dataset$tax_table %<>% tidy_taxonomy

#Calculating Alpha Diversity
observed=microbiome::alpha(inputR, index = "all")
meta=microbiome::meta(inputR)

#Creating a file to plot a graph
alpha= cbind(observed,meta)

#Testing differences
dunn.test::dunn.test(alpha$observed, alpha$Treatment,
                     method="bh", kw=TRUE, label=TRUE,
                     wrap=FALSE, table=F, list=T, rmc=FALSE, 
                     alpha=0.05, altp=FALSE)

alpha$Substrate <- factor(alpha$Substrate, levels = c("Control", "20%ADE", "100%ADE"))
alpha$Tree <- factor(alpha$Tree, levels = c("Bulk", "C. pachystachya",
                                            "P. dubium", "C. fissilis"))

alpha.plot = 
  ggplot(data = alpha, aes(x=Tree, y=observed, fill = Substrate)) +
  geom_boxplot() +
  labs(x = "Tree", y= "Number of Taxa") +
  theme(axis.text.x = element_text(angle = -15, vjust = 0.5, hjust=0.5)) +
  labs(color='Treatment') + guides(color = "none")+
  scale_fill_manual(values = c("#cb6751", "#9e6ebd", "#7aa457"))

alpha.plot

#Don't forget setwd
#dev.print(tiff, "./Figures/alpha_div.tiff", compression = "lzw", res=600, height=5, width=8, units="in")

#Firstly, a perMANOVA

df        = as(sample_data(inputR.clr), "data.frame")
ds        = phyloseq::distance(inputR.clr, method = "euclidean")
permanova = adonis2(ds ~ Tree*Substrate, data = df, permutations = 999)

permanova

#And the plot

input_ord = ordinate(inputR.clr, "PCoA" , "euclidean") 
#p3 = plot_ordination(inputR.clr, input_ord, color = "Tree", shape = "Substrate")
#p1.subs_shape = p3 + geom_point(aes(shape = Substrate), size = 6, alpha = 0.8) +
#  theme(legend.position = "right") +
#  scale_fill_manual(values = ze)

p4 = plot_ordination(inputR.clr, input_ord, color = "Substrate", shape = "Tree")
p1.tree_shape = p4 + geom_point(aes(shape = Tree), size = 6, alpha = 0.8) +
  theme(legend.position = "right") +
  annotate("text", x = -60, y = -48, hjust = 0.2 , 
           label = bquote('Substrate:'~R^2~'= 0.11  |  p = 0.001'), size = 3)+
  annotate("text", x = -60, y = -52, hjust = 0.2 , 
           label = bquote('Tree:'~R^2~'= 0.09  |  p = 0.004'), size = 3)+
  scale_fill_manual(values = ze)


#p1.subs_shape
p1.tree_shape

#dev.print(tiff, "./Figures/Beta_div.tiff", compression = "lzw", res=600, height=5, width=8, units="in")

# create trans_abund object
# using 10 Phyla with the highest abundance in the dataset.
dataset$cal_abund()
t1 <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 10)

# two facets example
# require package ggh4x, first run install.packages("ggh4x") if not installed
t1$plot_bar(others_color = "grey70", facet = "Tree",
            facet2 = "Substrate", xtext_keep = FALSE, legend_text_italic = FALSE, barwidth = 1)

#dev.print(tiff, "./Figures/Overall_phyla.tiff", compression = "lzw", res=600, height=5, width=8, units="in")

#Excluding initial samples
#DM is the object with only T1 samples

DM <- filter(map_tpsa, U_brizantha_weight != "NA")
DM = DM %>%
  arrange(U_brizantha_weight) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
  mutate(Substrate=factor(Substrate, levels=c("Control", "20%ADE", "100%ADE")))

comparison = list(c("Control", "100%ADE"), c("Control", "20%ADE"), c("20%ADE", "100%ADE"))

## Urochloa brizantha ##

DM.plot = ggplot(data = DM, mapping = aes(x = Substrate, y = U_brizantha_weight, fill = Substrate)) +
  geom_boxplot() +
  theme(legend.position = "none") +
  stat_compare_means(comparisons = comparison) +
  labs(title = "U. brizantha production", y = "Dry mass (g)") +
  scale_fill_manual(values = c("#7aa457", "#9e6ebd", "#cb6751"))
DM.plot

#dev.print(tiff, "./Figures/Urochloa.tiff", compression = "lzw", res=600, height=4, width=6, units="in")

## Plants' dry matter ##
TM.plot = ggplot(data = DM, mapping = aes(x = Substrate, y = Dry_Mass, fill = Treatment)) +
  geom_boxplot() +
  facet_grid(.~Tree, scales="free", space="free_x") +
  theme(legend.position="none") +
  labs(title="Dry matter production",
       x ="", y = "DM production (g)")

## Plant height ##
H.plot = ggplot(data = DM, mapping = aes(x = Substrate, y = Height, fill = Treatment)) +
  geom_boxplot() +
  facet_grid(.~Tree, scales="free", space="free_x") + 
  theme(legend.position="none") +
  labs(title="Plant height",
       x ="", y = "Height (cm)")

## Root size ##
RS.plot = ggplot(data = DM, mapping = aes(x = Substrate, y = Root_size, fill = Treatment)) +
  geom_boxplot() +
  facet_grid(.~Tree, scales="free", space="free_x") + 
  theme(legend.position="none") +
  labs(title="Root size",
       x ="Substrate", y = "Root size (cm)")

#Binding graphics
plants.g = ggarrange(TM.plot, H.plot, RS.plot, ncol=1, nrow = 3, align = "hv", common.legend = F)
plants.g

#dev.print(tiff, "./Figures/Growing.tiff", compression = "lzw", res=600, height=8, width=6, units="in")

cec.chem = filter(map_tpsa) 

cec.chem[is.na(cec.chem)] <- 0

my.variables <- cec.chem[c(4,7:9, 11:38)]

for(i in 1:length(my.variables)) {
  aaa = kruskal.test(x = as.matrix(my.variables[,i]), g = cec.chem$Treatment)
  if(aaa$p.value <= 0.05) {
    print(colnames(my.variables[i]))
    print(aaa)  
  } else {
  }
}

#RDA
#Analysing
#colnames(map_tp)

t1 <- trans_env$new(dataset = dataset, add_data = map_tp[, c(11:16, 24, 27, 35:38)])

t1$cal_diff(method = "KW_dunn", group = "Substrate")

t1$cal_ordination(method = "RDA", taxa_level = "Family")
t1$trans_ordination(show_taxa = 10, adjust_arrow_length = TRUE,
                    max_perc_env = 3, max_perc_tax = 0.5,
                    min_perc_env = 3, min_perc_tax = 0.5)

#Plotting
t1$plot_ordination(plot_color = "Substrate",
                   plot_shape = "Tree", plot_type = c("point", "ellipse"))

#dev.print(tiff, "RDA.tiff", compression = "lzw", res=600, height=6, width=8, units="in")

#FAPROTAX
dataset <- microtable$new(sample_table = map_tp,
                          otu_table    = otu_tp,
                          tax_table    = tax_tp)
dataset

dataset$tidy_dataset()

t1 <- trans_func$new(dataset)
t1$cal_spe_func(prok_database = "FAPROTAX")
t1$cal_spe_func_perc(abundance_weighted = FALSE)

#Saving the output
#genes = cbind(as.data.frame(t1$res_spe_func_perc), map_tp)
#write.csv(genes, file = "genes.csv")

genes <- read_csv("genes.csv")

genes[is.na(genes)] <- 0

my.variables <- genes[c(2:17)]

for(i in 1:length(my.variables)) {
  aaa = kruskal.test(x = as.matrix(my.variables[,i]), g = genes$Treatment)
  if(aaa$p.value <= 0.05) {
    print(colnames(my.variables[i]))
    print(aaa)
  } else {
  }
}