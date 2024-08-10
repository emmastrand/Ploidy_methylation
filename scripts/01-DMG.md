Identifying differentially methylated genes
================
EL Strand

# Differentially Methylated Genes

Input from this script:
<https://github.com/emmastrand/EmmaStrand_Notebook/blob/master/_posts/2021-10-21-KBay-Bleaching-Pairs-WGBS-Analysis-Pipeline.md#kbay-wgbs-methylation-analysis-pipeline>

I’m following Hollie Putnam and Kevin Wong’s pipelines. -
<https://github.com/kevinhwong1/Thermal_Transplant_Molecular/blob/main/scripts/WGBS_GM.Rmd>

## Load libraries

<https://github.com/hputnam/Geoduck_Meth/blob/master/RAnalysis/Scripts/GM.Rmd>

``` r
# BiocManager::install("goseq")

library(plyr)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(readxl)
library(plotrix) 
library(gridExtra)
library(seacarb) 
library(pheatmap)
library(tidyr)
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

## Load DNA methylation files

Below commands were only run use to produce data file and then not run
again for computational purposes.

``` r
# meth5x <- list.files(path = "data/WGBS/output/meth_counts_5x", 
#                      pattern = ".bed$", full.names=TRUE) %>% set_names(.) %>%
#   map_dfr(read.csv,.id="Sample.ID", header=FALSE, sep="\t", na.string="NA", stringsAsFactors = FALSE) %>%
#   dplyr::select(-c(V3,V7:V14)) 
# 
# colnames(meth5x) <- c("Sample.ID", "scaffold", "position","per.meth","meth","unmeth","gene")
# meth5x$gene <- gsub("ID=","",meth5x$gene) #remove extra characters
# meth5x$gene <- gsub("Parent=","",meth5x$gene) #remove extra characters
# meth5x$Sample.ID <- gsub("data/WGBS/output/meth_counts_5x/","",meth5x$Sample.ID) #remove extra characters
# meth5x$Sample.ID <- gsub("_5x_sorted.tab_gene_CpG_5x_enrichment.bed","",meth5x$Sample.ID) #remove extra characters 
# meth5x$Sample.ID <- as.character(meth5x$Sample.ID)
```

Load df created in code above

``` r
load("data/WGBS/meth5x.RData")
```

### Metadata files

``` r
meta <- read.csv("data/metadata/Molecular_metadata.csv") %>%
  mutate(Timepoint = if_else(Plug_ID == "2153", "12 hour", Timepoint)) %>%
  dplyr::rename(Sample.ID = Plug_ID) %>% 
  dplyr::select(Sample.ID, Species, Tank, Treatment, Temperature, CO2, Timepoint, Sample.Date) 
  
meta$Sample.ID <- as.character(meta$Sample.ID)

meta$Timepoint <- factor(meta$Timepoint, 
                           levels = c("0 hour", "6 hour", "12 hour", "30 hour",
                                      "1 week", "2 week", "4 week", "6 week",
                                      "8 week", "12 week", "16 week"))

pacuta_clade <- read.delim2("data/metadata/clade_annotations_Pacuta.txt", header = TRUE, sep="\t") %>%
  separate(., Sample, c("Species", "Treatment", "Timepoint", "Sample.ID")) %>% 
  dplyr::select(-Treatment, -Timepoint) %>% dplyr::select(Sample.ID, Clade)

ploidyinfo <- read.delim2("data/metadata/samples_Pacuta.annotations.txt", 
                          header=TRUE, sep="\t", na.strings=c("","NA")) %>% 
  dplyr::rename(., Sample.ID = plugid) %>% 
  mutate_if(is.integer, as.character) %>% 
  dplyr::select(Sample.ID, ploidy, group, treesplit) %>%
  #mutate(treesplit = coalesce(treesplit, "D")) %>% 
  #changing NA values to "D" in treesplit -- this is just so triploidy is split in 2 groups and diploidy is altogether 
  left_join(., pacuta_clade, by="Sample.ID")

meta <- dplyr::left_join(meta, ploidyinfo, by="Sample.ID")
```

### Merge dataframes together

``` r
df5x <- full_join(meth5x, meta, by = "Sample.ID") %>%
  mutate(meth_exp_group = case_when(
    group == "Group1" ~ "T1",
    group == "Group2" ~ "T1",
    group == "Group3" ~ "T2",
    group == "Group4" ~ "T2",
    group == "Group5" ~ "D",
    group == "Group6" ~ "D",
    group == "Group7" ~ "ungroup",
    group == "Ungroup" ~ "ungroup")) %>% 
  subset(!group == "Group7" & !group == "Ungroup") %>%
  filter(!Sample.ID == "1225") %>% 
  filter(!Sample.ID == "2197") 

df5x %>% 
  dplyr::select(Sample.ID, meth_exp_group, per.meth) %>% 
  filter(!is.na(per.meth)) %>% dplyr::select(-per.meth) %>%
  filter(!Sample.ID == "1225") %>% 
  filter(!Sample.ID == "2197") %>%
  distinct() %>% write.csv("data/metadata/meth_pattern_groups.csv")
```

## Binomial GLM to test for differentially methylated genes

### Filter for genes with \>5 methylated positions

Come back to this value of a cut off of 5 methylated positions per gene.

``` r
meth_table5x <- df5x
meth_table5x$sample_gene <- paste0(meth_table5x$Sample.ID, meth_table5x$gene)
```

``` r
meth_table5x_position_counts <- dplyr::count(meth_table5x, vars = c(sample_gene))
meth_table5x_position_counts <- meth_table5x_position_counts[which(meth_table5x_position_counts$n > 5), ]

meth_table5x_filtered <- meth_table5x[meth_table5x$sample_gene %in% meth_table5x_position_counts$vars,]
```

### Run loop for GLM

Create data frame to store results

``` r
results5x <- data.frame()
```

First subset the unique dataframes and second run the GLMs

``` r
gs5x <- unique(meth_table5x_filtered$gene)
```

#### 5x loop

``` r
# for(i in 1:length(meth_table5x_filtered$gene)){
#     
#     #subset the dataframe gene by gene
#     meth_table_filt1 <- subset(meth_table5x_filtered, gene ==gs5x[1])
#     
#     # fit glm position model
#     fit2 <- glm(matrix(c(meth, unmeth), ncol=2) ~ meth_exp_group, 
#                data=meth_table_filt1, family=binomial)
#     a2 <- anova(fit2, test="LRT")
# 
#     # capture summary stats to data frame
#     df2 <- data.frame(gene = meth_table_filt1[1,7], ###for treesplit [1,8] or at least double check this
#                      pval.methgroup = a2$`Pr(>Chi)`[2],
#                      stringsAsFactors = F)
#     # bind rows of temporary data frame to the results data frame
#     results5x <- rbind(results5x, df2)
# }
```

Export df created. Run the read.csv() function below the export code to
read back in those df. When running the script multiple times, you don’t
need to run the binomial function above. Skip that code and just read in
results5x or 10x df.

``` r
#results5x %>% write.csv("data/WGBS/results5x_methgroup.csv") #9,647 genes - this should match length of gs5x

results5x <- read.csv("data/WGBS/results5x_methgroup.csv", header = TRUE) %>% 
  dplyr::select(-X) %>% distinct() ## distinct() necessary in case the model above went through multiple times
```

Some genes return the error:
`"glm.fit: fitted probabilities numerically 0 or 1 occurredglm.fit: fitted probabilities numerically 0 or 1". These will be filtered out of analysis later on`
and
`Error in`contrasts\<-`(`*tmp*`, value = contr.funs[1 + isOF[nn]]) : contrasts can be applied only to factors with 2 or more levels`.

## Calculating adjusted p-values for glm

``` r
results5x[is.na(results5x)] <- 0
results5x$adj.pval.methgroup <- p.adjust(results5x$pval.methgroup, method='BH')
```

### subsetting new df to only those adjusted p value columns

methgroup: 3,870 genes significantly DMG out of 9,467 genes that had at
least 5 methylated positions and appeared in all samples.

``` r
DMG_5x_sig <- results5x %>%
  dplyr::select(gene, adj.pval.methgroup) %>% 
  mutate(across(2, round, 8)) %>%
  filter(adj.pval.methgroup < 0.050) 
```

    ## Warning: There was 1 warning in `mutate()`.
    ## ℹ In argument: `across(2, round, 8)`.
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
nrow(DMG_5x_sig) # 3,870 for methgroup
```

    ## [1] 3870

``` r
DMG_5x_sig %>% write.csv("data/WGBS/DMG_5x_sig_methgroup.csv")
```

#### Dataframes at this point

- `eggNOGG` = original gene annotation file
- `OE` = output from CpG OE script with weakly and heavily methylated
  gene sets
- `gene_meta` = merged metadata of genes
- `df5x / df10x` = .bed files and metadata together
- `meth_table5x_filtered / meth_table10x_filtered` = output from
  filtering all genes with 5+ methylated positions; input for the glm
  model
- `gs5x / gs10x` = list of unique genes from filtered methylation
  table  
- `results5x / results10x` = output from glm model; contains gene names,
  pvalues and adjusted pvalues for differentially methylated  
- `DMG_5x_sig / DMG_5x_sig` = only adjusted pvalues from the results 5
  and 10x dfs

### Adding “loci status” to each loci as a category and filtering to only significant DMGs

``` r
meth_table5x_filtered <- meth_table5x_filtered %>% 
  mutate(loci_status = case_when(
          per.meth <= 10 ~ "unmethylated",
          per.meth > 10 & per.meth < 50 ~ "sparsely methylated",
          per.meth >= 50 ~ "methylated"))

## hist for per.meth distribution for DMGs
meth_table5x_filtered_sigDMG <- meth_table5x_filtered[meth_table5x_filtered$gene %in% DMG_5x_sig$gene,]
```

## Methylation level calculations

### GENERAL HISTOGRAM OF CPG LOCI METHYLATION LEVEL

``` r
### number of loci per sample = 103,629 unique loci for all; 44,366 for sig DMG
## sample size of D = 16; T1=20; T2=17 

#### GENERAL HISTOGRAM OF CPG LOCI METHYLATION LEVEL 
## hist for per.meth distribution for all genes
histall <- meth_table5x_filtered %>%
  ggplot(aes(x=per.meth, fill=meth_exp_group)) + theme_bw() +
  facet_grid(factor(loci_status, levels=c('unmethylated','sparsely methylated','methylated'))~meth_exp_group, scales="free") +
  geom_histogram(alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("skyblue3", "olivedrab4", "darkgreen")) +
  labs(fill="") +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold"),
        strip.text.y = element_text(size = 12, color = "black", face = "bold")) +
  theme(strip.background = element_rect(color="black", fill="white", linewidth=0.5, linetype="solid"))

ggsave(filename="data/figures/Loci_Percent_methylation_all.jpeg", width=8, height=7, units="in")
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

``` r
## hist for per.meth distribution for significant DMGs 
histDMG <- meth_table5x_filtered_sigDMG %>%
  ggplot(aes(x=per.meth, fill=meth_exp_group)) + theme_bw() +
  facet_grid(factor(loci_status, levels=c('unmethylated','sparsely methylated','methylated'))~meth_exp_group, scales="free") +
  geom_histogram(alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("skyblue3", "olivedrab4", "darkgreen")) +
  labs(fill="") +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold"),
        strip.text.y = element_text(size = 12, color = "black", face = "bold")) +
  theme(strip.background = element_rect(color="black", fill="white", linewidth=0.5, linetype="solid"))

ggsave(filename="data/figures/Loci_Percent_methylation_sigDMGs.jpeg", width=8, height=7, units="in")
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

``` r
##### without unmethylated 
histun <- meth_table5x_filtered %>%
  filter(!loci_status == "unmethylated") %>%
  ggplot(aes(x=per.meth, fill=meth_exp_group)) + theme_bw() +
  facet_grid(~meth_exp_group, scales="free") +
  geom_histogram(alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("skyblue3", "olivedrab4", "darkgreen")) +
  labs(fill="") +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold"),
        strip.text.y = element_text(size = 12, color = "black", face = "bold")) +
  theme(strip.background = element_rect(color="black", fill="white", linewidth=0.5, linetype="solid"))

ggsave(filename="data/figures/Loci_Percent_methylation_filtered.jpeg", width=6, height=3, units="in")
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

### Filtering out loci with zero’s

Filter out loci that have median \< 10% across all samples.

meth_table5x_filtered = 5,346,733 loci meth_table5x_filtered_sigDMG =
2,278,227 loci

meth_table5x_filtered2 = 229,661 loci in 1,976 genes (4,503 loci per
sample) meth_table5x_filtered2_sigDMG = 198,647 loci in 1,447 genes
(3,895 loci per sample)

``` r
meth_table5x_filtered2 <- meth_table5x_filtered %>% group_by(position) %>%
  mutate(position.median = median(per.meth)) %>%
  filter(!position.median < 10) %>% ungroup()

## viewing one loci to see if the above function worked 
meth_table5x_filtered2 %>% filter(position == "3396")
```

    ## # A tibble: 0 × 22
    ## # ℹ 22 variables: Sample.ID <chr>, scaffold <chr>, position <int>,
    ## #   per.meth <dbl>, meth <int>, unmeth <int>, gene <chr>, Species <chr>,
    ## #   Tank <int>, Treatment <chr>, Temperature <chr>, CO2 <chr>, Timepoint <fct>,
    ## #   Sample.Date <int>, ploidy <chr>, group <chr>, treesplit <chr>, Clade <chr>,
    ## #   meth_exp_group <chr>, sample_gene <chr>, loci_status <chr>,
    ## #   position.median <dbl>

``` r
meth_table5x_filtered %>% filter(position == "3396")
```

    ##    Sample.ID                           scaffold position per.meth meth unmeth
    ## 1       1047 Pocillopora_acuta_HIv2___Sc0000000     3396 0.000000    0      9
    ## 2       1051 Pocillopora_acuta_HIv2___Sc0000000     3396 0.000000    0     20
    ## 3       1090 Pocillopora_acuta_HIv2___Sc0000000     3396 0.000000    0     26
    ## 4       1103 Pocillopora_acuta_HIv2___Sc0000000     3396 0.000000    0     11
    ## 5       1147 Pocillopora_acuta_HIv2___Sc0000000     3396 0.000000    0      7
    ## 6       1159 Pocillopora_acuta_HIv2___Sc0000000     3396 5.882353    1     16
    ## 7       1168 Pocillopora_acuta_HIv2___Sc0000000     3396 0.000000    0     23
    ## 8       1184 Pocillopora_acuta_HIv2___Sc0000000     3396 0.000000    0     19
    ## 9       1238 Pocillopora_acuta_HIv2___Sc0000000     3396 0.000000    0     12
    ## 10      1281 Pocillopora_acuta_HIv2___Sc0000000     3396 0.000000    0     30
    ## 11      1296 Pocillopora_acuta_HIv2___Sc0000000     3396 0.000000    0     24
    ## 12      1303 Pocillopora_acuta_HIv2___Sc0000000     3396 0.000000    0     24
    ## 13      1329 Pocillopora_acuta_HIv2___Sc0000000     3396 0.000000    0     26
    ## 14      1416 Pocillopora_acuta_HIv2___Sc0000000     3396 3.333333    1     29
    ## 15      1427 Pocillopora_acuta_HIv2___Sc0000000     3396 0.000000    0     21
    ## 16      1459 Pocillopora_acuta_HIv2___Sc0000000     3396 0.000000    0     13
    ## 17      1487 Pocillopora_acuta_HIv2___Sc0000000     3396 0.000000    0     18
    ## 18      1536 Pocillopora_acuta_HIv2___Sc0000000     3396 0.000000    0     23
    ## 19      1559 Pocillopora_acuta_HIv2___Sc0000000     3396 4.761905    1     20
    ## 20      1563 Pocillopora_acuta_HIv2___Sc0000000     3396 0.000000    0     28
    ## 21      1571 Pocillopora_acuta_HIv2___Sc0000000     3396 0.000000    0     20
    ## 22      1582 Pocillopora_acuta_HIv2___Sc0000000     3396 0.000000    0     18
    ## 23      1596 Pocillopora_acuta_HIv2___Sc0000000     3396 0.000000    0     30
    ## 24      1641 Pocillopora_acuta_HIv2___Sc0000000     3396 0.000000    0     25
    ## 25      1647 Pocillopora_acuta_HIv2___Sc0000000     3396 0.000000    0     11
    ## 26      1707 Pocillopora_acuta_HIv2___Sc0000000     3396 0.000000    0      5
    ## 27      1709 Pocillopora_acuta_HIv2___Sc0000000     3396 0.000000    0     77
    ## 28      1728 Pocillopora_acuta_HIv2___Sc0000000     3396 0.000000    0     23
    ## 29      1732 Pocillopora_acuta_HIv2___Sc0000000     3396 0.000000    0     21
    ## 30      1755 Pocillopora_acuta_HIv2___Sc0000000     3396 0.000000    0     29
    ## 31      1757 Pocillopora_acuta_HIv2___Sc0000000     3396 0.000000    0     21
    ## 32      1765 Pocillopora_acuta_HIv2___Sc0000000     3396 0.000000    0     14
    ## 33      1777 Pocillopora_acuta_HIv2___Sc0000000     3396 0.000000    0      8
    ## 34      1820 Pocillopora_acuta_HIv2___Sc0000000     3396 0.000000    0     23
    ## 35      2012 Pocillopora_acuta_HIv2___Sc0000000     3396 0.000000    0     24
    ## 36      2064 Pocillopora_acuta_HIv2___Sc0000000     3396 0.000000    0     17
    ## 37      2072 Pocillopora_acuta_HIv2___Sc0000000     3396 0.000000    0     12
    ## 38      2087 Pocillopora_acuta_HIv2___Sc0000000     3396 0.000000    0     25
    ## 39      2212 Pocillopora_acuta_HIv2___Sc0000000     3396 0.000000    0      8
    ## 40      2300 Pocillopora_acuta_HIv2___Sc0000000     3396 0.000000    0     25
    ## 41      2304 Pocillopora_acuta_HIv2___Sc0000000     3396 0.000000    0     22
    ## 42      2306 Pocillopora_acuta_HIv2___Sc0000000     3396 0.000000    0     17
    ## 43      2409 Pocillopora_acuta_HIv2___Sc0000000     3396 0.000000    0     18
    ## 44      2413 Pocillopora_acuta_HIv2___Sc0000000     3396 2.702703    1     36
    ## 45      2513 Pocillopora_acuta_HIv2___Sc0000000     3396 0.000000    0     27
    ## 46      2564 Pocillopora_acuta_HIv2___Sc0000000     3396 0.000000    0     17
    ## 47      2668 Pocillopora_acuta_HIv2___Sc0000000     3396 0.000000    0     17
    ## 48      2861 Pocillopora_acuta_HIv2___Sc0000000     3396 0.000000    0     20
    ## 49      2877 Pocillopora_acuta_HIv2___Sc0000000     3396 0.000000    0     13
    ## 50      2878 Pocillopora_acuta_HIv2___Sc0000000     3396 0.000000    0     35
    ## 51      2879 Pocillopora_acuta_HIv2___Sc0000000     3396 0.000000    0     12
    ##                                     gene Species Tank Treatment Temperature
    ## 1  Pocillopora_acuta_HIv2___TS.g10153.t1  Pacuta    1      ATAC     Ambient
    ## 2  Pocillopora_acuta_HIv2___TS.g10153.t1  Pacuta    1      ATAC     Ambient
    ## 3  Pocillopora_acuta_HIv2___TS.g10153.t1  Pacuta    2      HTHC        High
    ## 4  Pocillopora_acuta_HIv2___TS.g10153.t1  Pacuta    1      ATAC     Ambient
    ## 5  Pocillopora_acuta_HIv2___TS.g10153.t1  Pacuta    3      ATHC     Ambient
    ## 6  Pocillopora_acuta_HIv2___TS.g10153.t1  Pacuta    1      ATAC     Ambient
    ## 7  Pocillopora_acuta_HIv2___TS.g10153.t1  Pacuta    2      HTHC        High
    ## 8  Pocillopora_acuta_HIv2___TS.g10153.t1  Pacuta    2      HTHC        High
    ## 9  Pocillopora_acuta_HIv2___TS.g10153.t1  Pacuta    2      HTHC        High
    ## 10 Pocillopora_acuta_HIv2___TS.g10153.t1  Pacuta    3      ATHC     Ambient
    ## 11 Pocillopora_acuta_HIv2___TS.g10153.t1  Pacuta    3      ATHC     Ambient
    ## 12 Pocillopora_acuta_HIv2___TS.g10153.t1  Pacuta    4      HTAC        High
    ## 13 Pocillopora_acuta_HIv2___TS.g10153.t1  Pacuta    4      HTAC        High
    ## 14 Pocillopora_acuta_HIv2___TS.g10153.t1  Pacuta    5      HTHC        High
    ## 15 Pocillopora_acuta_HIv2___TS.g10153.t1  Pacuta    5      HTHC        High
    ## 16 Pocillopora_acuta_HIv2___TS.g10153.t1  Pacuta    3      ATHC     Ambient
    ## 17 Pocillopora_acuta_HIv2___TS.g10153.t1  Pacuta    7      HTAC        High
    ## 18 Pocillopora_acuta_HIv2___TS.g10153.t1  Pacuta    7      HTAC        High
    ## 19 Pocillopora_acuta_HIv2___TS.g10153.t1  Pacuta    8      ATAC     Ambient
    ## 20 Pocillopora_acuta_HIv2___TS.g10153.t1  Pacuta    8      ATAC     Ambient
    ## 21 Pocillopora_acuta_HIv2___TS.g10153.t1  Pacuta    7      HTAC        High
    ## 22 Pocillopora_acuta_HIv2___TS.g10153.t1  Pacuta    9      HTAC        High
    ## 23 Pocillopora_acuta_HIv2___TS.g10153.t1  Pacuta    4      HTAC        High
    ## 24 Pocillopora_acuta_HIv2___TS.g10153.t1  Pacuta    6      ATAC     Ambient
    ## 25 Pocillopora_acuta_HIv2___TS.g10153.t1  Pacuta    7      HTAC        High
    ## 26 Pocillopora_acuta_HIv2___TS.g10153.t1  Pacuta    9      HTAC        High
    ## 27 Pocillopora_acuta_HIv2___TS.g10153.t1  Pacuta   10      HTHC        High
    ## 28 Pocillopora_acuta_HIv2___TS.g10153.t1  Pacuta    4      HTAC        High
    ## 29 Pocillopora_acuta_HIv2___TS.g10153.t1  Pacuta   10      HTHC        High
    ## 30 Pocillopora_acuta_HIv2___TS.g10153.t1  Pacuta    6      ATAC     Ambient
    ## 31 Pocillopora_acuta_HIv2___TS.g10153.t1  Pacuta    6      ATAC     Ambient
    ## 32 Pocillopora_acuta_HIv2___TS.g10153.t1  Pacuta    7      HTAC        High
    ## 33 Pocillopora_acuta_HIv2___TS.g10153.t1  Pacuta    8      ATAC     Ambient
    ## 34 Pocillopora_acuta_HIv2___TS.g10153.t1  Pacuta   10      HTHC        High
    ## 35 Pocillopora_acuta_HIv2___TS.g10153.t1  Pacuta    8      ATAC     Ambient
    ## 36 Pocillopora_acuta_HIv2___TS.g10153.t1  Pacuta    9      HTAC        High
    ## 37 Pocillopora_acuta_HIv2___TS.g10153.t1  Pacuta    9      HTAC        High
    ## 38 Pocillopora_acuta_HIv2___TS.g10153.t1  Pacuta   10      HTHC        High
    ## 39 Pocillopora_acuta_HIv2___TS.g10153.t1  Pacuta   11      ATHC     Ambient
    ## 40 Pocillopora_acuta_HIv2___TS.g10153.t1  Pacuta    5      HTHC        High
    ## 41 Pocillopora_acuta_HIv2___TS.g10153.t1  Pacuta    5      HTHC        High
    ## 42 Pocillopora_acuta_HIv2___TS.g10153.t1  Pacuta    6      ATAC     Ambient
    ## 43 Pocillopora_acuta_HIv2___TS.g10153.t1  Pacuta   11      ATHC     Ambient
    ## 44 Pocillopora_acuta_HIv2___TS.g10153.t1  Pacuta    8      ATAC     Ambient
    ## 45 Pocillopora_acuta_HIv2___TS.g10153.t1  Pacuta    9      HTAC        High
    ## 46 Pocillopora_acuta_HIv2___TS.g10153.t1  Pacuta   12      ATHC     Ambient
    ## 47 Pocillopora_acuta_HIv2___TS.g10153.t1  Pacuta   11      ATHC     Ambient
    ## 48 Pocillopora_acuta_HIv2___TS.g10153.t1  Pacuta   11      ATHC     Ambient
    ## 49 Pocillopora_acuta_HIv2___TS.g10153.t1  Pacuta   12      ATHC     Ambient
    ## 50 Pocillopora_acuta_HIv2___TS.g10153.t1  Pacuta   12      ATHC     Ambient
    ## 51 Pocillopora_acuta_HIv2___TS.g10153.t1  Pacuta   12      ATHC     Ambient
    ##        CO2 Timepoint Sample.Date   ploidy  group treesplit    Clade
    ## 1  Ambient    2 week    20181006  Diploid Group6         D CladePA5
    ## 2  Ambient    4 week    20181020 Triploid Group1        T1 CladePA2
    ## 3     High    2 week    20181006 Triploid Group3        T2 CladePA3
    ## 4  Ambient   12 week    20181215 Triploid Group2        T1 CladePA1
    ## 5     High   12 week    20181215 Triploid Group2        T1 CladePA1
    ## 6  Ambient    8 week    20181117  Diploid Group5      <NA> CladePA6
    ## 7     High   30 hour    20180923  Diploid Group6         D CladePA5
    ## 8     High    4 week    20181020 Triploid Group4        T2 CladePA4
    ## 9     High    8 week    20181117 Triploid Group3        T2 CladePA3
    ## 10    High    2 week    20181006  Diploid Group6         D CladePA5
    ## 11    High   30 hour    20180923  Diploid Group6         D CladePA5
    ## 12 Ambient   30 hour    20180923  Diploid Group6         D CladePA5
    ## 13 Ambient    4 week    20181020  Diploid Group6         D CladePA5
    ## 14    High   12 week    20181215  Diploid Group6         D CladePA5
    ## 15    High    2 week    20181006 Triploid Group4        T2 CladePA4
    ## 16    High    4 week    20181020 Triploid Group1        T1 CladePA2
    ## 17 Ambient    2 week    20181006  Diploid Group6         D CladePA5
    ## 18 Ambient    8 week    20181117 Triploid Group3        T2 CladePA3
    ## 19 Ambient    8 week    20181117 Triploid Group3        T2 CladePA3
    ## 20 Ambient   30 hour    20180923 Triploid Group3        T2 CladePA3
    ## 21 Ambient   30 hour    20180923  Diploid Group6         D CladePA5
    ## 22 Ambient   12 week    20181215 Triploid Group3        T2 CladePA3
    ## 23 Ambient   12 week    20181215 Triploid Group2        T1 CladePA1
    ## 24 Ambient    8 week    20181117 Triploid Group3        T2 CladePA3
    ## 25 Ambient   12 week    20181215 Triploid Group3        T2 CladePA3
    ## 26 Ambient   30 hour    20180923 Triploid Group2        T1 CladePA1
    ## 27    High    4 week    20181020 Triploid Group2        T1 CladePA1
    ## 28 Ambient    2 week    20181006 Triploid Group2        T1 CladePA1
    ## 29    High    8 week    20181117 Triploid Group3        T2 CladePA3
    ## 30 Ambient    4 week    20181020  Diploid Group5      <NA> CladePA6
    ## 31 Ambient   30 hour    20180923 Triploid Group3        T2 CladePA3
    ## 32 Ambient    4 week    20181020 Triploid Group3        T2 CladePA3
    ## 33 Ambient   12 week    20181215 Triploid Group1        T1 CladePA2
    ## 34    High    2 week    20181006 Triploid Group3        T2 CladePA3
    ## 35 Ambient    4 week    20181020 Triploid Group2        T1 CladePA1
    ## 36 Ambient    8 week    20181117 Triploid Group2        T1 CladePA1
    ## 37 Ambient    2 week    20181006 Triploid Group1        T1 CladePA2
    ## 38    High   30 hour    20180923 Triploid Group2        T1 CladePA1
    ## 39    High   30 hour    20180923  Diploid Group5      <NA> CladePA6
    ## 40    High    8 week    20181117 Triploid Group2        T1 CladePA1
    ## 41    High    4 week    20181020 Triploid Group2        T1 CladePA1
    ## 42 Ambient   12 week    20181215 Triploid Group2        T1 CladePA1
    ## 43    High    2 week    20181006 Triploid Group3        T2 CladePA3
    ## 44 Ambient    2 week    20181006 Triploid Group2        T1 CladePA1
    ## 45 Ambient    4 week    20181020 Triploid Group2        T1 CladePA1
    ## 46    High    4 week    20181020 Triploid Group3        T2 CladePA3
    ## 47    High   12 week    20181215  Diploid Group6         D CladePA5
    ## 48    High    4 week    20181020  Diploid Group6         D CladePA5
    ## 49    High   30 hour    20180923 Triploid Group2        T1 CladePA1
    ## 50    High    2 week    20181006 Triploid Group2        T1 CladePA1
    ## 51    High   12 week    20181215  Diploid Group6         D CladePA5
    ##    meth_exp_group                               sample_gene  loci_status
    ## 1               D 1047Pocillopora_acuta_HIv2___TS.g10153.t1 unmethylated
    ## 2              T1 1051Pocillopora_acuta_HIv2___TS.g10153.t1 unmethylated
    ## 3              T2 1090Pocillopora_acuta_HIv2___TS.g10153.t1 unmethylated
    ## 4              T1 1103Pocillopora_acuta_HIv2___TS.g10153.t1 unmethylated
    ## 5              T1 1147Pocillopora_acuta_HIv2___TS.g10153.t1 unmethylated
    ## 6               D 1159Pocillopora_acuta_HIv2___TS.g10153.t1 unmethylated
    ## 7               D 1168Pocillopora_acuta_HIv2___TS.g10153.t1 unmethylated
    ## 8              T2 1184Pocillopora_acuta_HIv2___TS.g10153.t1 unmethylated
    ## 9              T2 1238Pocillopora_acuta_HIv2___TS.g10153.t1 unmethylated
    ## 10              D 1281Pocillopora_acuta_HIv2___TS.g10153.t1 unmethylated
    ## 11              D 1296Pocillopora_acuta_HIv2___TS.g10153.t1 unmethylated
    ## 12              D 1303Pocillopora_acuta_HIv2___TS.g10153.t1 unmethylated
    ## 13              D 1329Pocillopora_acuta_HIv2___TS.g10153.t1 unmethylated
    ## 14              D 1416Pocillopora_acuta_HIv2___TS.g10153.t1 unmethylated
    ## 15             T2 1427Pocillopora_acuta_HIv2___TS.g10153.t1 unmethylated
    ## 16             T1 1459Pocillopora_acuta_HIv2___TS.g10153.t1 unmethylated
    ## 17              D 1487Pocillopora_acuta_HIv2___TS.g10153.t1 unmethylated
    ## 18             T2 1536Pocillopora_acuta_HIv2___TS.g10153.t1 unmethylated
    ## 19             T2 1559Pocillopora_acuta_HIv2___TS.g10153.t1 unmethylated
    ## 20             T2 1563Pocillopora_acuta_HIv2___TS.g10153.t1 unmethylated
    ## 21              D 1571Pocillopora_acuta_HIv2___TS.g10153.t1 unmethylated
    ## 22             T2 1582Pocillopora_acuta_HIv2___TS.g10153.t1 unmethylated
    ## 23             T1 1596Pocillopora_acuta_HIv2___TS.g10153.t1 unmethylated
    ## 24             T2 1641Pocillopora_acuta_HIv2___TS.g10153.t1 unmethylated
    ## 25             T2 1647Pocillopora_acuta_HIv2___TS.g10153.t1 unmethylated
    ## 26             T1 1707Pocillopora_acuta_HIv2___TS.g10153.t1 unmethylated
    ## 27             T1 1709Pocillopora_acuta_HIv2___TS.g10153.t1 unmethylated
    ## 28             T1 1728Pocillopora_acuta_HIv2___TS.g10153.t1 unmethylated
    ## 29             T2 1732Pocillopora_acuta_HIv2___TS.g10153.t1 unmethylated
    ## 30              D 1755Pocillopora_acuta_HIv2___TS.g10153.t1 unmethylated
    ## 31             T2 1757Pocillopora_acuta_HIv2___TS.g10153.t1 unmethylated
    ## 32             T2 1765Pocillopora_acuta_HIv2___TS.g10153.t1 unmethylated
    ## 33             T1 1777Pocillopora_acuta_HIv2___TS.g10153.t1 unmethylated
    ## 34             T2 1820Pocillopora_acuta_HIv2___TS.g10153.t1 unmethylated
    ## 35             T1 2012Pocillopora_acuta_HIv2___TS.g10153.t1 unmethylated
    ## 36             T1 2064Pocillopora_acuta_HIv2___TS.g10153.t1 unmethylated
    ## 37             T1 2072Pocillopora_acuta_HIv2___TS.g10153.t1 unmethylated
    ## 38             T1 2087Pocillopora_acuta_HIv2___TS.g10153.t1 unmethylated
    ## 39              D 2212Pocillopora_acuta_HIv2___TS.g10153.t1 unmethylated
    ## 40             T1 2300Pocillopora_acuta_HIv2___TS.g10153.t1 unmethylated
    ## 41             T1 2304Pocillopora_acuta_HIv2___TS.g10153.t1 unmethylated
    ## 42             T1 2306Pocillopora_acuta_HIv2___TS.g10153.t1 unmethylated
    ## 43             T2 2409Pocillopora_acuta_HIv2___TS.g10153.t1 unmethylated
    ## 44             T1 2413Pocillopora_acuta_HIv2___TS.g10153.t1 unmethylated
    ## 45             T1 2513Pocillopora_acuta_HIv2___TS.g10153.t1 unmethylated
    ## 46             T2 2564Pocillopora_acuta_HIv2___TS.g10153.t1 unmethylated
    ## 47              D 2668Pocillopora_acuta_HIv2___TS.g10153.t1 unmethylated
    ## 48              D 2861Pocillopora_acuta_HIv2___TS.g10153.t1 unmethylated
    ## 49             T1 2877Pocillopora_acuta_HIv2___TS.g10153.t1 unmethylated
    ## 50             T1 2878Pocillopora_acuta_HIv2___TS.g10153.t1 unmethylated
    ## 51              D 2879Pocillopora_acuta_HIv2___TS.g10153.t1 unmethylated

``` r
meth_table5x_filtered2_sigDMG <- meth_table5x_filtered2[meth_table5x_filtered2$gene %in% DMG_5x_sig$gene,]
length(unique(meth_table5x_filtered2_sigDMG$gene)) ##1,447 genes left 
```

    ## [1] 1447

``` r
save(meth_table5x_filtered2, file = "data/WGBS/meth_table5x_filtered2.RData")
save(meth_table5x_filtered2_sigDMG, file = "data/WGBS/meth_table5x_filtered2_sigDMG.RData")
```
