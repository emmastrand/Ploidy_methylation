Methylation statistics
================
EL Strand

Based on Kevin Wong’s methylseq analysis scripts.

## Resources

Ignore names with “\_val”; <https://github.com/ewels/MultiQC/issues/626>

``` r
library(tidyverse)
library(ggplot2)
library(readxl)
library(ggtext)
```

## Load data

``` r
# load data
bismark_data <- read_excel("data/WGBS/methylseqstats_total.xlsx", sheet = "Qualimap_Bismark") %>%
  mutate_at(., vars(contains('Percent')), ~ (. * 100)) %>%
  dplyr::rename(
    `Q: GC (%)` = `Percent GC`,
    `Q: Insert Size (bp)` = `Ins. size`,
    `Q: >= 5X (%)` = `Percent ≥ 5X`,
    `Q: >= 10X (%)` = `Percent ≥ 10X`,
    `Q: Median Coverage (X)` = `Median cov`,
    `Q: Mean Coverage (X)` = `Mean cov`,
    `Q: Aligned (%)` = `Percent Aligned...8`,
    `Q: Aligned (M)` = `M Aligned...9`,
    `Q: Total Reads (M)` = `M Total reads`,
    `Q: Error rate (%)`= `Percent Error rate`,
    `B: mCpG (%)` = `Percent mCpG`,
    `B: mCHG (%)` = `Percent mCHG`,
    `B: mCHH (%)` = `Percent mCHH`,
    `B: Cytosines (M)` = `M C's`,
    `B: Duplicate Reads (%)` = `Percent Dups`,
    `B: Unique Reads (M)` = `M Unique`,
    `B: Aligned (M)` = `M Aligned...18`,
    `B: Aligned (%)` = `Percent Aligned...19`
  )
```

    ## New names:
    ## • `Percent Aligned` -> `Percent Aligned...8`
    ## • `M Aligned` -> `M Aligned...9`
    ## • `M Aligned` -> `M Aligned...18`
    ## • `Percent Aligned` -> `Percent Aligned...19`

``` r
bismark_data$Plug_ID = substr(bismark_data$`Sample Name`, 1, 4)

fastqc_data <- read_excel("data/WGBS/methylseqstats_total.xlsx", sheet = "cutadapt_fastqc") %>%
  mutate_at(., vars(contains('Percent')), ~ (. * 100)) %>%
  dplyr::rename(
    `FQ: Trimmed (% bp)` = `Percent BP Trimmed`,
    `FQ: Duplicates (%)` = `Percent Dups`,
    `FQ: GC (%)` = `Percent GC`,
    `FQ: Length (bp)` = `Length`,
    `FQ: Failed (%)` = `Percent Failed`,
    `FQ: Total Reads (M)` = `M Seqs`
  )

fastqc_data$Plug_ID = substr(fastqc_data$`Sample Name`, 1, 4)

loci_numbers <- read_excel("data/WGBS/methylseqstats_total.xlsx", sheet = "Loci_pre_allfilt")
loci_numbers$Plug_ID <- as.character(loci_numbers$Plug_ID)
```

## Metadata

``` r
load("data/metadata/meta.RData") ## created in script 01-DMG]
meta <- meta %>% dplyr::rename(Plug_ID = Sample.ID)

data <- left_join(bismark_data, meta, by = "Plug_ID") %>% 
  left_join(., loci_numbers, by = "Plug_ID")

fastqc_data <- left_join(fastqc_data, meta, by = "Plug_ID")
```

## Plotting

``` r
data %>% 
  gather("metric", "value", 2:19) %>%
  ggplot(., aes(x=ploidy, y=value)) + 
  geom_boxplot(aes(color = ploidy), outlier.shape=NA, fill=NA) + 
  geom_jitter(aes(fill = ploidy), width=0.15, color='black', shape=21, shape=2, alpha=0.4) +
  labs(
    fill = "Ploidy Group",
    color = "Ploidy Group",
    y="",
    x=""
  ) +
  
    facet_wrap(~metric, scales = "free_y") + 
  
  theme_bw() +
    scale_colour_manual(values = c("skyblue3", "olivedrab4")) +
    scale_fill_manual(values = c("skyblue3", "olivedrab4")) +
  
  theme(strip.text.x = element_text(size = 8, color = "black", face = "bold"),
        strip.text.y = element_text(size = 8, color = "black", face = "bold"),
        legend.position = "none",
        strip.background = element_rect(color="black", fill="white", linewidth=0.5, linetype="solid"))
```

    ## Warning: Duplicated aesthetics after name standardisation: shape

![](00-QC-statistics_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
ggsave("data/figures/methylseq_stats.png", width=10, height=6.5)
```

``` r
fastqc_data %>% 
  gather("metric", "value", 3:8) %>%
  filter(!is.na(value)) %>%
  ggplot(., aes(x=ploidy, y=value, color=ploidy)) + 
    geom_boxplot(aes(color = ploidy), outlier.shape=NA, fill=NA) + 
    geom_jitter(aes(fill = ploidy), width=0.15, color='black', shape=21, alpha=0.4) +
    facet_wrap(Trimming~metric, scales = "free_y",
               labeller = labeller(
                 Trimming = c(
                   "pre" = "Pre-Filtering",
                   "post" = "Post-Filtering"
                 ))) +
    labs(
      y="",
      x=""
    ) +
    theme_bw() +
    scale_colour_manual(values = c("skyblue3", "olivedrab4")) +
    scale_fill_manual(values = c("skyblue3", "olivedrab4")) +
  
    theme(
      strip.text.x = element_markdown(size = 8, color = "black", face = "bold"),
      strip.text.y = element_markdown(size = 8, color = "black", face = "bold"),
      legend.position = "none",
      strip.background = element_rect(color="black", fill="white", linewidth=0.5, linetype="solid")
    )
```

![](00-QC-statistics_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
ggsave("data/figures/FQ_stats.png", width=7, height=5)
```
