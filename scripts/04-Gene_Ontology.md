Gene Ontology Analysis for DMGs
================
EL Strand

# Gene Ontology for Differentially Methylated Genes

I’m following Hollie Putnam and Kevin Wong’s pipelines. -
<https://github.com/kevinhwong1/Thermal_Transplant_Molecular/blob/main/scripts/WGBS_GM.Rmd>

## Load libraries

<https://github.com/hputnam/Geoduck_Meth/blob/master/RAnalysis/Scripts/GM.Rmd>

``` r
# BiocManager::install("goseq")

library(plyr)
library(dplyr)
library(purrr)
library(ggplot2)
library(tidyverse)
library(readxl)
library(plotrix) 
library(gridExtra)
library(seacarb) 
library(pheatmap)
library(tidyr)
library(goseq)
library(genefilter)
library(gplots)
library(cowplot)
library(lsmeans)
library(data.table)
library(RColorBrewer)
library(randomcoloR)
#library(GSEABase)
library(ggpubr)
library(hrbrthemes)
library(viridis)
library(factoextra)
#library(ropls)
#library(mixOmics)
library(vegan)
library(Rmisc)
```

## Input data from other script

### Methylation data

``` r
load("data/WGBS/meth_table5x_filtered2.RData") 
#meth_table5x_filtered ## OR run this prior to GOseq large function 

load("data/WGBS/meth_table5x_filtered2_sigDMG.RData")
#meth_table5x_filtered2_sigDMG <- read.csv("data/WGBS/output/meth_table5x_filtered2_sigDMG.csv") %>% dplyr::select(-X)
## OR run other script (DMG-Pacuta.Rmd) to create this file and then run this script

hypermeth_trip2 <- read.csv("data/WGBS/hyper-hypo methylation/DMG_hypermeth_triploidy2.csv") 
hypermeth_trip1 <- read.csv("data/WGBS/hyper-hypo methylation/DMG_hypermeth_triploidy1.csv") 
hypermeth_dip <- read.csv("data/WGBS/hyper-hypo methylation/DMG_hypermeth_diploidy.csv") 
hypometh_dip <- read.csv("data/WGBS/hyper-hypo methylation/DMG_hypometh_diploidy.csv")
hypometh_trip1 <- read.csv("data/WGBS/hyper-hypo methylation/DMG_hypometh_triploidy1.csv")
hypometh_trip2 <- read.csv("data/WGBS/hyper-hypo methylation/DMG_hypometh_triploidy2.csv")

hypermeth_trip2_GOlist <- meth_table5x_filtered2_sigDMG[meth_table5x_filtered2_sigDMG$gene %in% hypermeth_trip2$gene,]
hypermeth_trip1_GOlist <- meth_table5x_filtered2_sigDMG[meth_table5x_filtered2_sigDMG$gene %in% hypermeth_trip1$gene,]
hypermeth_dip_GOlist <- meth_table5x_filtered2_sigDMG[meth_table5x_filtered2_sigDMG$gene %in% hypermeth_dip$gene,]
hypometh_dip_GOlist <- meth_table5x_filtered2_sigDMG[meth_table5x_filtered2_sigDMG$gene %in% hypometh_dip$gene,]
hypometh_trip1_GOlist <- meth_table5x_filtered2_sigDMG[meth_table5x_filtered2_sigDMG$gene %in% hypometh_trip1$gene,]
hypometh_trip2_GOlist <- meth_table5x_filtered2_sigDMG[meth_table5x_filtered2_sigDMG$gene %in% hypometh_trip2$gene,]
```

### Read in gff file to get lengths and merge with gene_meta created file

``` r
eggNOGG <- read.delim(file = "data/gene expression/Pocillopora_acuta_HIv2.genes.EggNog_results.txt",
                      sep = "\t", header=TRUE) %>% dplyr::rename(gene = X.query)

# OE <- read.csv(file = "data/WGBS/Pacuta_CpGOE_full.csv", sep = ",",
#                  header = TRUE) 

# there is a space after each gene name.. get rid of. This was some error in CpG_OE.Rmd script 
# OE$Gene <- gsub(" ","",OE$Gene)

#merge with annotation information
# gene_meta <- dplyr::full_join(eggNOGG, OE, by = "Gene") %>% dplyr::rename(gene = Gene)

gff <- read.delim("data/gene expression/Pocillopora_acuta_HIv2.genes.gff3", 
                  sep = "\t", header = FALSE) %>% dplyr::rename(gene = V9) %>% 
  dplyr::rename(start = V4) %>% 
  dplyr::rename(stop = V5) %>% 
  dplyr::select(gene, start, stop) %>%
  mutate(Length = stop-start) %>%
  dplyr::filter(grepl("ID",gene)) %>%
  mutate(across(everything(), gsub, pattern = "ID=", replacement = "")) %>%
  dplyr::select(gene, Length)
```

    ## Warning: There was 1 warning in `mutate()`.
    ## ℹ In argument: `across(everything(), gsub, pattern = "ID=", replacement = "")`.
    ## Caused by warning:
    ## ! The `...` argument of `across()` is deprecated as of dplyr 1.1.0.
    ## Supply arguments directly to `.fns` through an anonymous function instead.
    ## 
    ##   # Previously
    ##   across(a:b, mean, na.rm = TRUE)
    ## 
    ##   # Now
    ##   across(a:b, \(x) mean(x, na.rm = TRUE))

``` r
gff$Length <- as.numeric(gff$Length)

# gene_meta <- left_join(gene_meta, gff, by = "gene")
gene_meta <- left_join(eggNOGG, gff, by = "gene")
```

### Setting up files for go enrichment

Generate list of all genes identified and background list for GO
enrichment analysis

``` r
# 5x
all5x.genes <- as.data.frame(unique(meth_table5x_filtered2$gene)) #1,796 genes total
colnames(all5x.genes) <- "gene"
all5x.genes_annot <- left_join(all5x.genes, gene_meta, by = "gene")
```

### Vector with all genes after filtering

``` r
# 5x 
ALL.5x.vector <-c(t(all5x.genes_annot$gene))
ID.5x.vector <- all5x.genes_annot$gene
LENGTH.5x.vector <-all5x.genes_annot$Length

#goslim <- getOBOCollection("http://current.geneontology.org/ontology/subsets/goslim_generic.obo")
#goslimdf <- data.frame(goSlim(BPGO_collection, goslim))
  
goslim <- read.csv("https://raw.githubusercontent.com/kevinhwong1/Thermal_Transplant_Molecular/main/data/WGBS/GO-GOslim.csv")
goslim <- goslim %>% dplyr::select(-term)
```

### Creating GO ID file

``` r
meta_ID <- gene_meta %>% dplyr::select(gene, GOs) %>% 
  dplyr::rename(GO_IDs = GOs) %>%
  dplyr::mutate(across(everything(), gsub, pattern = ",", replacement = "; "),
         across(everything(), gsub, pattern = "-", replacement = NA))

splitted <- strsplit(as.character(meta_ID$GO_IDs), "; ") #slit into multiple GO ids to get them all in a single row
GO.terms <- data.frame(v1 = rep.int(meta_ID$gene, sapply(splitted, length)), 
                       v2 = unlist(splitted)) #list all genes with each of their GO terms in a single row 
colnames(GO.terms) <- c("gene", "GO_IDs")
```

### Writing goenrich function

<https://bioconductor.org/packages/devel/bioc/vignettes/goseq/inst/doc/goseq.pdf>

PWF = Probability Weighting Function

``` r
goenrich5x <- function(filename, 
                     identifier){
  #filename=methgroup_GO_sig
                DMG <- as.character(filename$gene) #set the enrichment test list
                DMG.vector <-c(t(DMG)) #change to vectors
                
                gene.vector=as.integer(ALL.5x.vector%in%DMG.vector) #Construct new vector with 1 for DEG and 0 for others
                names(gene.vector)=ALL.5x.vector #set names
                DEG.pwf<-nullp(gene.vector, ID.5x.vector, bias.data=LENGTH.5x.vector) #weight vector by length of gene
                ## col names of DEG.pwf (3 total): gene as row name, DEgenes, bias.data, pwf 
                
                #Find enriched GO terms
                GO.wall<-goseq(DEG.pwf, ID.5x.vector, gene2cat=GO.terms, 
                               test.cats=c("GO:CC", "GO:BP", "GO:MF"), 
                               method="Wallenius", use_genes_without_cat=TRUE)
                ### col names of GO.wall (5 total): category, over_represented_pvalue, under_represented_pvalue, 
                ### numDEInCat, numInCat
                
                GO <- GO.wall[order(GO.wall$over_represented_pvalue),]
                colnames(GO)[1] <- "GO_ID"
                
                #GOslim 
                GO.slim <- merge(GO, goslim, by = "GO_ID")
                GO.slim <- GO.slim[!duplicated(GO.slim$GO_ID), ]
                
                #Filtering for p > 0.05
                filename_sig.GO <- GO.slim %>%
                  dplyr::filter(over_represented_pvalue <0.05) %>%
                  arrange(., ontology, term, over_represented_pvalue)
                
                # Formatting GOs column 
                filename <- filename %>%
                  dplyr::rename(GO_IDs = GOs) %>%
                  mutate(across(everything(), gsub, pattern = ",", replacement = "; "),
                         across(everything(), gsub, pattern = "-", replacement = NA))
                
                #Formatting sig gene file with a goterm per row
                split <- strsplit(as.character(filename$GO_IDs), "; ") 
                split2 <- data.frame(v1 = rep.int(filename$gene, sapply(split, length)), 
                                     v2 = unlist(split)) 
                colnames(split2) <- c("gene", "GO_IDs")
                
                filename2 <- filename %>% dplyr::select(-GO_IDs, -Length)
                filename_GO <- merge(split2, filename2, by = "gene")
                
                colnames(filename_GO)[2] <- "GO_ID"
                
                # Merge sig meth genes with goslim
                filename_GOslim <- dplyr::left_join(filename_GO, filename_sig.GO, by = "GO_ID")
                
                return(filename_GOslim)
}
```

### Merge files of interest with annotation

``` r
# ## want to run the above function on: meth_table5x_filtered_sigDMG
# ## this set is too big to run this dataframe. 
# # methgroup_GO_sig <- left_join(meth_table5x_filtered2_sigDMG, gene_meta, by = "gene")
# # GOenrich_sig5x <- goenrich5x(methgroup_GO_sig, meth_exp_group)
# # write.csv(GOenrich_sig5x , file = "data/WGBS/output/GOresults_methgroup_sig5x.csv")
# 
# GOenrich_sig5x_annotated <- GOenrich_sig5x %>%
#   select(-GO_ID) %>%
#   filter(over_represented_pvalue != "NA") %>%
#   filter(!is.na(term)) %>%
#   distinct(gene, .keep_all = TRUE) %>% #keep_all = keep all columns; 1,539 genes annotated with term 
#   mutate(GeneRatio = numDEInCat/numInCat)
# 
# GOplot <- GOenrich_sig5x_annotated %>%
#   ggplot(aes(x = (GeneRatio*100), y = reorder(GOSlim_bin, GeneRatio), 
#              fill = over_represented_pvalue, size=numDEInCat)) +
#   xlab("% of genes involved") +
#   scale_y_discrete(position = "right") +
#   geom_hline(aes(yintercept = GOSlim_bin), linetype = "dotted", color = "grey") +
#   geom_point(shape = 21, color = "black") +
#   theme_bw() + 
#   theme(panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(), 
#     axis.line = element_line(colour = "black"),
#     legend.text=element_text(size=10), legend.title=element_text(size=10),
#     strip.text.y = element_text(size = 15),
#     strip.text.x = element_text(size = 15),
#     axis.text = element_text(size = 10, color = "black"),
#     axis.title.x = element_text(size = 15),
#     axis.title.y = element_blank(), 
#     legend.position = "right") 
# 
# ggsave(filename="data/WGBS/output/figures/GOenrich_ploidy_bin_plot.png", 
#        plot=GOplot, dpi=300, width=8, height=8, units="in")
# 
# GOplot2 <- GOenrich_sig5x_annotated %>%
#   #ggplot(aes(x = numDEInCat, y = reorder(term, numDEInCat), fill = over_represented_pvalue)) +
#   ggplot(aes(x = (GeneRatio*100), y = reorder(term, GeneRatio), 
#              fill = over_represented_pvalue, size=numDEInCat)) +
#   xlab("% of genes involved") +
#   scale_y_discrete(position = "right") +
#   geom_hline(aes(yintercept = term), linetype = "dotted", color = "grey") +
#   geom_point(shape = 21, color = "black") +
#   theme_bw() + theme(panel.grid.major = element_blank(),
#   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
#   strip.text.y = element_text(size = 12),
#   strip.text.x = element_text(size = 12),
#   axis.text = element_text(size = 12, color = "black"),
#   axis.title.x = element_text(size = 12),
#   axis.title.y = element_blank(), 
#   legend.position = "right") #+
#   #facet_grid(GOSlim_bin ~ ., scale = "free_y",
#             # labeller = label_wrap_gen(width = 2, multi_line = TRUE))
# 
# ggsave(filename="data/WGBS/output/figures/GOenrich_ploidy_term_plot.png", 
#        plot=GOplot2, dpi=300, width=12, height=12, units="in")
```

### Subset to genes of interest

``` r
# unique(GOenrich_sig5x_annotated$GOSlim_bin)
# ERgolgi <- GOenrich_sig5x_annotated %>% subset(GOSlim_bin=="ER/Golgi")
# cellORG <- GOenrich_sig5x_annotated %>% subset(GOSlim_bin=="cell organization and biogenesis")
# RNAmetab <- GOenrich_sig5x_annotated %>% subset(GOSlim_bin=="RNA metabolism")
# cellpro <- GOenrich_sig5x_annotated %>% subset(GOSlim_bin=="cell cycle and proliferation")
# cyto <- GOenrich_sig5x_annotated %>% subset(GOSlim_bin=="cytoskeleton")
# nuc <- GOenrich_sig5x_annotated %>% subset(GOSlim_bin=="nucleus")
# othercyto <- GOenrich_sig5x_annotated %>% subset(GOSlim_bin=="other cytoplasmic organelle")
# dev <- GOenrich_sig5x_annotated %>% subset(GOSlim_bin=="developmental processes")
# othcell <- GOenrich_sig5x_annotated %>% subset(GOSlim_bin=="other cellular component")
# othmetab <- GOenrich_sig5x_annotated %>% subset(GOSlim_bin=="other metabolic processes")
# protmetab <- GOenrich_sig5x_annotated %>% subset(GOSlim_bin=="protein metabolism")
# othmolec <- GOenrich_sig5x_annotated %>% subset(GOSlim_bin=="other molecular function")
# othbio <- GOenrich_sig5x_annotated %>% subset(GOSlim_bin=="other biological processes")
# transport <- GOenrich_sig5x_annotated %>% subset(GOSlim_bin=="transport")
# sig <- GOenrich_sig5x_annotated %>% subset(GOSlim_bin=="signal transduction")
```

## Graph methylation levels

``` r
# ## 699 genes in developmental processes 
# GOenrich_sig5x_dev <- GOenrich_sig5x %>%
#   subset(GOSlim_bin == "developmental processes") %>%
#   dplyr::select(gene) %>% distinct() %>%
#   left_join(., meth_table5x_filtered, by = "gene") %>%
#   group_by(Sample.ID, gene) %>%
#   mutate(mean = mean(per.meth)) %>% ungroup() %>%
#   group_by(Sample.ID) %>%
#   mutate(mean.sample = mean(mean)) %>% ungroup() %>%
#   group_by(treesplit, gene) %>%
#   mutate(mean.ploidy.gene = mean(mean)) %>% ungroup() %>%
#   filter(!is.na(treesplit)) %>%
#   dplyr::select(treesplit, mean.sample) %>% distinct()
# GOenrich_sig5x_dev$GOSlim_bin <- "developmental processes"
# 
# ## 430 genes in stress response 
# GOenrich_sig5x_stress <- GOenrich_sig5x %>%
#   subset(GOSlim_bin == "stress response") %>%
#   dplyr::select(gene) %>% distinct() %>%
#   left_join(., meth_table5x_filtered, by = "gene") %>%
#   group_by(Sample.ID, gene) %>%
#   mutate(mean = mean(per.meth)) %>% ungroup() %>%
#   group_by(Sample.ID) %>%
#   mutate(mean.sample = mean(mean)) %>% ungroup() %>%
#   group_by(treesplit, gene) %>%
#   mutate(mean.ploidy.gene = mean(mean)) %>% ungroup() %>%
#   filter(!is.na(treesplit)) %>%
#   dplyr::select(treesplit, mean.sample) %>% distinct()
# GOenrich_sig5x_stress$GOSlim_bin <- "stress response"
# 
# GOenrich_subset <- bind_rows(GOenrich_sig5x_dev, GOenrich_sig5x_stress)
# 
# GOenrich_subset %>% 
#   ggplot(., aes(x=treesplit, y=mean.sample, color=treesplit)) + 
#   geom_boxplot() + 
#   geom_jitter(alpha=0.1, width = 0.25) + 
#   facet_wrap(~GOSlim_bin, scales="free_y") +
#   theme_classic() + ylab("Methylation % per dinucleotide") + xlab("Ploidy") +
#   scale_colour_manual(values = c("skyblue3", "olivedrab4", "darkgreen")) +
#   theme(strip.text.x = element_text(size = 12, color = "black", face = "bold"),
#         strip.text.y = element_text(size = 12, color = "black", face = "bold")) +
#   theme(strip.background = element_rect(color="black", fill="white", linewidth=0.5, linetype="solid"))
# 
# aov <-aov(mean.sample ~ treesplit*GOSlim_bin, data=GOenrich_subset)
# summary(aov)
# TukeyHSD(aov)
```

``` r
# GOenrich_sig5x_subset <- GOenrich_sig5x %>%
#   subset(GOSlim_bin == "stress response" | GOSlim_bin == "developmental processes")
# 
# DMG_ploidysig ### gene counts data from Targeted-GE script (Vst transformed)
# meth_table5x_filtered ### meth counts from DMG-Pacuta script 
# 
# DMG_ploidysig_GE <- meth_table5x_filtered[meth_table5x_filtered$gene %in% GOenrich_sig5x_subset$gene,]
# DMG_ploidysig_GE <- DMG_ploidysig_GE %>%
#   dplyr::select(gene, Sample.ID, per.meth, treesplit, ) %>%
#   group_by(gene, Sample.ID) %>%
#   mutate(mean.per.meth = mean(per.meth)) %>% ungroup() %>% dplyr::select(-per.meth) %>%
#   distinct() 
# DMG_ploidysig_GE$Sample.ID <- as.character(DMG_ploidysig_GE$Sample.ID)
# 
# DMG_ploidysig2 <- DMG_ploidysig[DMG_ploidysig$gene %in% DMG_ploidysig_GE$gene,]
# DMG_ploidysig2 <- DMG_ploidysig2 %>% dplyr::select(gene, Sample.ID, count)
# 
# DMG_ploidysig_GE <- DMG_ploidysig_GE %>% 
#   left_join(., DMG_ploidysig2, by = c("gene", "Sample.ID"))
# 
# DMG_ploidysig_GE_summary <- summarySE(data=DMG_ploidysig_GE, measurevar = "count", 
#                                       groupvars = c("Sample.ID", "treesplit"), na.rm = TRUE)
#   
# DMG_ploidysig_GE_summary <- DMG_ploidysig_GE_summary %>%
#   mutate(CV = (sd/count)*100,
#          inv_CV = (1/CV),
#          log_expression = log10(count+1))
# 
# DMG_ploidysig_GE_summary %>%
#   filter(!is.na(treesplit)) %>%
#   gather("measurement", "value", 8:10) %>%
#   subset(!measurement=="inv_CV") %>%
#   ggplot(., aes(x=treesplit, y=value, color=treesplit)) + 
#   facet_wrap(~measurement, scales="free_y") +
#   geom_boxplot() + 
#   geom_jitter(alpha=0.1, width = 0.25) +
#   theme_classic() + 
#   scale_colour_manual(values = c("skyblue3", "olivedrab4", "darkgreen")) 
# 
# aov <-aov(CV ~ treesplit, data=DMG_ploidysig_GE_summary)
# summary(aov)
# TukeyHSD(aov)
```

## Run GO term analysis on hypermeth subset from DMG-Pacuta.Rmd

``` r
## want to run the above function on: hypermeth_trip2
hypermeth_trip2_GOlist <- hypermeth_trip2_GOlist %>% dplyr::select(gene) %>% distinct()
hypermeth_trip2_GO <- left_join(hypermeth_trip2_GOlist, gene_meta, by = "gene")
GOenrich_hypermeth_trip2 <- goenrich5x(hypermeth_trip2_GO, triploidy2)
```

    ## Warning in pcls(G): initial point very close to some inequality constraints

![](04-Gene_Ontology_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

    ## Using manually entered categories.

    ## Calculating the p-values...

    ## 'select()' returned 1:1 mapping between keys and columns

``` r
write.csv(GOenrich_hypermeth_trip2 , file = "data/WGBS/GOenrich/GOenrich_hypermeth_trip2.csv", row.names = FALSE)

hypometh_trip2_GOlist <- hypometh_trip2_GOlist %>% dplyr::select(gene) %>% distinct()
hypometh_trip2_GO <- left_join(hypometh_trip2_GOlist, gene_meta, by = "gene")
GOenrich_hypometh_trip2 <- goenrich5x(hypometh_trip2_GO, triploidy2)
```

![](04-Gene_Ontology_files/figure-gfm/unnamed-chunk-12-2.png)<!-- -->

    ## Using manually entered categories.

    ## Calculating the p-values...

    ## 'select()' returned 1:1 mapping between keys and columns

``` r
write.csv(GOenrich_hypometh_trip2 , file = "data/WGBS/GOenrich/GOenrich_hypometh_trip2.csv", row.names = FALSE)

hypermeth_trip1_GOlist <- hypermeth_trip1_GOlist %>% dplyr::select(gene) %>% distinct()
hypermeth_trip1_GO <- left_join(hypermeth_trip1, gene_meta, by = "gene")
GOenrich_hypermeth_trip1 <- goenrich5x(hypermeth_trip1_GO, triploidy1)
```

![](04-Gene_Ontology_files/figure-gfm/unnamed-chunk-12-3.png)<!-- -->

    ## Using manually entered categories.

    ## Calculating the p-values...

    ## 'select()' returned 1:1 mapping between keys and columns

``` r
write.csv(GOenrich_hypermeth_trip1 , file = "data/WGBS/GOenrich/GOenrich_hypermeth_trip1.csv", row.names = FALSE)

hypometh_trip1_GOlist <- hypometh_trip1_GOlist %>% dplyr::select(gene) %>% distinct()
hypometh_trip1_GO <- left_join(hypometh_trip1_GOlist, gene_meta, by = "gene")
GOenrich_hypometh_trip1 <- goenrich5x(hypometh_trip1_GO, triploidy1)
```

    ## Warning in pcls(G): initial point very close to some inequality constraints

![](04-Gene_Ontology_files/figure-gfm/unnamed-chunk-12-4.png)<!-- -->

    ## Using manually entered categories.

    ## Calculating the p-values...

    ## 'select()' returned 1:1 mapping between keys and columns

``` r
write.csv(GOenrich_hypometh_trip1 , file = "data/WGBS/GOenrich/GOenrich_hypometh_trip1.csv", row.names = FALSE)

hypermeth_dip_GOlist <- hypermeth_dip_GOlist %>% dplyr::select(gene) %>% distinct()
hypermeth_dip_GO <- left_join(hypermeth_dip, gene_meta, by = "gene")
GOenrich_hypermeth_dip <- goenrich5x(hypermeth_dip_GO, diploidy)
```

![](04-Gene_Ontology_files/figure-gfm/unnamed-chunk-12-5.png)<!-- -->

    ## Using manually entered categories.

    ## Calculating the p-values...

    ## 'select()' returned 1:1 mapping between keys and columns

``` r
write.csv(GOenrich_hypermeth_dip , file = "data/WGBS/GOenrich/GOenrich_hypermeth_dip.csv", row.names = FALSE)

hypometh_dip_GOlist <- hypometh_dip_GOlist %>% dplyr::select(gene) %>% distinct()
hypometh_dip_GO <- left_join(hypometh_dip, gene_meta, by = "gene")
GOenrich_hypometh_dip <- goenrich5x(hypometh_dip_GO, diploidy)
```

![](04-Gene_Ontology_files/figure-gfm/unnamed-chunk-12-6.png)<!-- -->

    ## Using manually entered categories.

    ## Calculating the p-values...

    ## 'select()' returned 1:1 mapping between keys and columns

``` r
write.csv(GOenrich_hypometh_dip , file = "data/WGBS/GOenrich/GOenrich_hypometh_dip.csv", row.names = FALSE)

# GOenrich_hypermeth_dip <- read.csv("data/WGBS/output/GOenrich_hypermeth_dip.csv") 
# GOenrich_hypermeth_trip1 <- read.csv("data/WGBS/output/GOenrich_hypermeth_trip1.csv") 
# GOenrich_hypermeth_trip2 <- read.csv("data/WGBS/output/GOenrich_hypermeth_trip2.csv") 

### no. annotated / total no. genes 
## HYPER trip2 = 53 / 101  
## HYPER trip1 = 35 / 70 
## HYPER dip = 100 / 415 
## HYPO trip2 = 56 / 272
## HYPO trip1 =  62 / 140
## HYPO dip = 56 / 143 

#length(unique(GOenrich_hypometh_dip_annotated$gene))

## TRIP 2
GOenrich_hypermeth_trip2_annotated <- GOenrich_hypermeth_trip2 %>%
  dplyr::select(-GO_ID) %>%
  filter(over_represented_pvalue != "NA") %>%
  filter(!is.na(term)) %>%
  distinct(gene, .keep_all = TRUE) %>% #keep_all = keep all columns
  mutate(GeneRatio = numDEInCat/numInCat) %>%
  mutate(category = "T2")

GOenrich_hypometh_trip2_annotated <- GOenrich_hypometh_trip2 %>%
  dplyr::select(-GO_ID) %>%
  filter(over_represented_pvalue != "NA") %>%
  filter(!is.na(term)) %>%
  distinct(gene, .keep_all = TRUE) %>% #keep_all = keep all columns
  mutate(GeneRatio = numDEInCat/numInCat) %>%
  mutate(category = "T2")

## TRIP 1
GOenrich_hypermeth_trip1_annotated <- GOenrich_hypermeth_trip1 %>%
  dplyr::select(-GO_ID) %>%
  filter(over_represented_pvalue != "NA") %>%
  filter(!is.na(term)) %>%
  distinct(gene, .keep_all = TRUE) %>% #keep_all = keep all columns; 139/399 genes annotated with term 
  mutate(GeneRatio = numDEInCat/numInCat) %>%
  mutate(category = "T1")

GOenrich_hypometh_trip1_annotated <- GOenrich_hypometh_trip1 %>%
  dplyr::select(-GO_ID) %>%
  filter(over_represented_pvalue != "NA") %>%
  filter(!is.na(term)) %>%
  distinct(gene, .keep_all = TRUE) %>% #keep_all = keep all columns; 139/399 genes annotated with term 
  mutate(GeneRatio = numDEInCat/numInCat) %>%
  mutate(category = "T1")

## DIPLOIDY 
GOenrich_hypermeth_dip_annotated <- GOenrich_hypermeth_dip %>%
  dplyr::select(-GO_ID) %>%
  filter(over_represented_pvalue != "NA") %>%
  filter(!is.na(term)) %>%
  distinct(gene, .keep_all = TRUE) %>% #keep_all = keep all columns; 470/1,056 genes annotated with term 
  mutate(GeneRatio = numDEInCat/numInCat) %>%
  mutate(category = "D")

GOenrich_hypometh_dip_annotated <- GOenrich_hypometh_dip %>%
  dplyr::select(-GO_ID) %>%
  filter(over_represented_pvalue != "NA") %>%
  filter(!is.na(term)) %>%
  distinct(gene, .keep_all = TRUE) %>% #keep_all = keep all columns; 470/1,056 genes annotated with term 
  mutate(GeneRatio = numDEInCat/numInCat) %>%
  mutate(category = "D")

# making large dataframe of them altogether 
GOenrich_hypermeth_annotated_total <- 
  bind_rows(GOenrich_hypermeth_trip2_annotated, 
            GOenrich_hypermeth_trip1_annotated, GOenrich_hypermeth_dip_annotated)

GOenrich_hypometh_annotated_total <- 
  bind_rows(GOenrich_hypometh_trip2_annotated, 
            GOenrich_hypometh_trip1_annotated, GOenrich_hypometh_dip_annotated)
```

Hypermeth plot

``` r
GOenrich_hypermeth_annotated_total %>%
  filter(ontology == "BP") %>% 
  dplyr::select(term, GOSlim_bin, category, GeneRatio, over_represented_pvalue) %>%
  distinct() %>%
  ggplot(aes(x = category, y = reorder(term, GOSlim_bin), 
             fill = over_represented_pvalue, size=GeneRatio)) +
  scale_y_discrete(position = "right") + xlab("") +
  geom_hline(aes(yintercept = term), linetype = "dotted", color = "grey") +
  geom_point(position = position_jitter(w = 0, h = 0), shape=21, color="black") +
  geom_vline(aes(xintercept = 1.5), color="grey", lty = "dotted") +
  geom_vline(aes(xintercept = 2.5), color="grey", lty = "dotted") +
  facet_grid(GOSlim_bin ~ ., scales = "free", space = "free",
             labeller = label_wrap_gen(width = 8, multi_line = TRUE)) +
  labs(
    fill = "Over represented p-value",
    size = "Gene ratio"
  ) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black"),
    legend.text=element_text(size=10), legend.title=element_text(size=10),
    strip.text.y = element_text(size = 10, angle=0, face="bold"),
    strip.clip = 'off',
    strip.placement = "outside",
    axis.text = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_blank(), 
    legend.position = "right") 
```

    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA

![](04-Gene_Ontology_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
ggsave("data/figures/Fig3_GOenrich_term_hypermeth.png", dpi=300, width=9, height=10.5, units="in")
```

    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA

Hypometh plot

``` r
GOenrich_hypometh_annotated_total %>%
  filter(ontology == "BP") %>% 
  dplyr::select(term, GOSlim_bin, category, GeneRatio, over_represented_pvalue) %>%
  distinct() %>%
  ggplot(aes(x = category, y = reorder(term, GOSlim_bin), 
             fill = over_represented_pvalue, size=GeneRatio)) +
  scale_y_discrete(position = "right") + xlab("") +
  geom_hline(aes(yintercept = term), linetype = "dotted", color = "grey") +
  geom_point(position = position_jitter(w = 0, h = 0), shape=21, color="black") +
  geom_vline(aes(xintercept = 1.5), color="grey", lty = "dotted") +
  geom_vline(aes(xintercept = 2.5), color="grey", lty = "dotted") +
  facet_grid(GOSlim_bin ~ ., scales = "free", space = "free",
             labeller = label_wrap_gen(width = 8, multi_line = TRUE)) +
  labs(
    fill = "Over represented p-value",
    size = "Gene ratio"
  ) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black"),
    legend.text=element_text(size=10), legend.title=element_text(size=10),
    strip.text.y = element_text(size = 10, angle=0, face="bold"),
    strip.clip = 'off',
    strip.placement = "outside",
    axis.text = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_blank(), 
    legend.position = "right") 
```

    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA

![](04-Gene_Ontology_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

``` r
ggsave("data/figures/Fig3_GOenrich_term_hypometh.png", dpi=300, width=9.67, height=10.5, units="in")
```

    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
    ## Warning in mean.default(X[[i]], ...): argument is not numeric or logical:
    ## returning NA
