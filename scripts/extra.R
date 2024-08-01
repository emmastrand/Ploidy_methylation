### P. acuta annotation from Rutgers 

http://cyanophora.rutgers.edu/Pocillopora_acuta/
  
  GO Slim filed created from GOSlim.Rmd script in same repo. 

```{r}
eggNOGG <- read.delim(file = "data/Pocillopora_acuta_HIv2.genes.EggNog_results.txt",
                      sep = "\t", header=TRUE) %>% dplyr::rename(gene = X.query)
```



### CpG OE calculated from prior scripts 

Generated in this script: https://github.com/emmastrand/EmmaStrand_Notebook/blob/master/_posts/2022-11-19-CpG-OE-Analysis-for-DNA-Methylation.md

```{r}
OE <- read.csv(file = "data/WGBS/Pacuta_CpGOE_full.csv", sep = ",",
               header = TRUE) %>% dplyr::rename(gene = Gene)

# there is a space after each gene name.. get rid of. This was some error in CpG_OE.Rmd script 
OE$gene <- gsub(" ","",OE$gene)

#merge with annotation information
gene_meta <- dplyr::full_join(eggNOGG, OE, by = "gene")  
OE_NA <- gene_meta %>% subset(.,is.na(status)) ##26 genes found with no CpG O/E calculated (this gene was not in genome but in gff)
```