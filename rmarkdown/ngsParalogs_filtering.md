fst ngsParalog filtering
================

## Load relevant libraries and data

``` r
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.0     ✔ readr     2.1.4
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.0
    ## ✔ ggplot2   3.4.1     ✔ tibble    3.1.8
    ## ✔ lubridate 1.9.2     ✔ tidyr     1.3.0
    ## ✔ purrr     1.0.1     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the ]8;;http://conflicted.r-lib.org/conflicted package]8;; to force all conflicts to become errors

``` r
# table for chromosome name and number
# two columns: chr (numbered) and chrName (what your chromosomes are called in your angsd output file)
chrom_df <- read.table("/fs/cbsubscb16/storage/rkc/sample_lists/chrom_meta_data.txt", header = TRUE, col.names = c("chrName","chrNewName")) %>% mutate(chr = seq(1,104,1)) %>% dplyr::select(-chrNewName)

# read in position files
AI_wgph_pos <- read_delim("/fs/cbsubscb16/storage/rkc/angsd/counts/PCAM-PPLA_AI_wgph.pos") %>% dplyr::rename(chrName = chr) %>% left_join(chrom_df, by = "chrName", ) %>%
    subset(select = -c(chrName))
```

    ## Rows: 2954767 Columns: 3
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (1): chr
    ## dbl (2): pos, totDepth
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
EBS_wgph_pos <- read_delim("/fs/cbsubscb16/storage/rkc/angsd/counts/PCAM-PPLA_EBS_wgph.pos") %>% dplyr::rename(chrName = chr) %>% left_join(chrom_df, by = "chrName", ) %>%
    subset(select = -c(chrName))
```

    ## Rows: 4619509 Columns: 3
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (1): chr
    ## dbl (2): pos, totDepth
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
GOA_wgph_pos <- read_delim("/fs/cbsubscb16/storage/rkc/angsd/counts/PCAM-PPLA_GOA_wgph.pos") %>% dplyr::rename(chrName = chr) %>% left_join(chrom_df, by = "chrName", ) %>%
    subset(select = -c(chrName))
```

    ## Rows: 5252065 Columns: 3
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (1): chr
    ## dbl (2): pos, totDepth
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
NBS_wgph_pos <- read_delim("/fs/cbsubscb16/storage/rkc/angsd/counts/PCAM-PPLA_NBS_wgph.pos") %>% dplyr::rename(chrName = chr) %>% left_join(chrom_df, by = "chrName", ) %>%
    subset(select = -c(chrName))
```

    ## Rows: 3029137 Columns: 3
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (1): chr
    ## dbl (2): pos, totDepth
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
SEAK_wgph_pos <- read_delim("/fs/cbsubscb16/storage/rkc/angsd/counts/PCAM-PPLA_SEAK_wgph.pos") %>% dplyr::rename(chrName = chr) %>% left_join(chrom_df, by = "chrName", ) %>%
    subset(select = -c(chrName))
```

    ## Rows: 4168293 Columns: 3
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (1): chr
    ## dbl (2): pos, totDepth
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
AI_EBS_fst <- read_tsv("/fs/cbsubscb16/storage/rkc/angsd/fst/AI-EastBering.fst.txt", col_names = c("region", "chrName", "pos", "Nsites","fst"), skip = 1) %>%  mutate(correct_fst = ifelse(fst < 0, 0, fst)) %>%  
  left_join(chrom_df, by = "chrName") %>%
    subset(select = -c(chrName)) %>%
    dplyr::select(chr, pos, correct_fst)
```

    ## Rows: 2980157 Columns: 5
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (2): region, chrName
    ## dbl (3): pos, Nsites, fst
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
AI_GOA_fst <- read_tsv("/fs/cbsubscb16/storage/rkc/angsd/fst/AI-GOA.fst.txt", col_names = c("region", "chrName", "pos", "Nsites","fst"), skip = 1) %>%  mutate(correct_fst = ifelse(fst < 0, 0, fst)) %>%  
  left_join(chrom_df, by = "chrName") %>%
    subset(select = -c(chrName)) %>%
    dplyr::select(chr, pos, correct_fst)
```

    ## Rows: 3024251 Columns: 5
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (2): region, chrName
    ## dbl (3): pos, Nsites, fst
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
NBS_AI_fst <- read_tsv("/fs/cbsubscb16/storage/rkc/angsd/fst/AI-NorthBering.fst.txt", col_names = c("region", "chrName", "pos", "Nsites","fst"), skip = 1) %>%  mutate(correct_fst = ifelse(fst < 0, 0, fst)) %>%  
  left_join(chrom_df, by = "chrName") %>%
    subset(select = -c(chrName)) %>%
    dplyr::select(chr, pos, correct_fst)
```

    ## Rows: 2331023 Columns: 5
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (2): region, chrName
    ## dbl (3): pos, Nsites, fst
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
AI_SEAK_fst <- read_tsv("/fs/cbsubscb16/storage/rkc/angsd/fst/AI-SEAK.fst.txt", col_names = c("region", "chrName", "pos", "Nsites","fst"), skip = 1) %>%  mutate(correct_fst = ifelse(fst < 0, 0, fst)) %>%  
  left_join(chrom_df, by = "chrName") %>%
    subset(select = -c(chrName)) %>%
    dplyr::select(chr, pos, correct_fst)
```

    ## Rows: 2813926 Columns: 5
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (2): region, chrName
    ## dbl (3): pos, Nsites, fst
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
EBS_GOA_fst <- read_tsv("/fs/cbsubscb16/storage/rkc/angsd/fst/EastBering-GOA.fst.txt", col_names = c("region", "chrName", "pos", "Nsites","fst"), skip = 1) %>%  mutate(correct_fst = ifelse(fst < 0, 0, fst)) %>%  
  left_join(chrom_df, by = "chrName") %>%
    subset(select = -c(chrName)) %>%
    dplyr::select(chr, pos, correct_fst)
```

    ## Rows: 4352166 Columns: 5
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (2): region, chrName
    ## dbl (3): pos, Nsites, fst
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
NBS_EBS_fst <- read_tsv("/fs/cbsubscb16/storage/rkc/angsd/fst/EastBering-NorthBering.fst.txt", col_names = c("region", "chrName", "pos", "Nsites","fst"), skip = 1) %>%  mutate(correct_fst = ifelse(fst < 0, 0, fst)) %>%  
  left_join(chrom_df, by = "chrName") %>%
    subset(select = -c(chrName)) %>%
    dplyr::select(chr, pos, correct_fst)
```

    ## Rows: 2933753 Columns: 5
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (2): region, chrName
    ## dbl (3): pos, Nsites, fst
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
EBS_SEAK_fst <- read_tsv("/fs/cbsubscb16/storage/rkc/angsd/fst/EastBering-SEAK.fst.txt", col_names = c("region", "chrName", "pos", "Nsites","fst"), skip = 1) %>%  mutate(correct_fst = ifelse(fst < 0, 0, fst)) %>%  
  left_join(chrom_df, by = "chrName") %>%
    subset(select = -c(chrName)) %>%
    dplyr::select(chr, pos, correct_fst)
```

    ## Rows: 3727073 Columns: 5
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (2): region, chrName
    ## dbl (3): pos, Nsites, fst
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
NBS_GOA_fst <- read_tsv("/fs/cbsubscb16/storage/rkc/angsd/fst/GOA-NorthBering.fst.txt", col_names = c("region", "chrName", "pos", "Nsites","fst"), skip = 1) %>%  mutate(correct_fst = ifelse(fst < 0, 0, fst)) %>%  
  left_join(chrom_df, by = "chrName") %>%
    subset(select = -c(chrName)) %>%
    dplyr::select(chr, pos, correct_fst)
```

    ## Rows: 3033644 Columns: 5
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (2): region, chrName
    ## dbl (3): pos, Nsites, fst
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
GOA_SEAK_fst <- read_tsv("/fs/cbsubscb16/storage/rkc/angsd/fst/GOA-SEAK.fst.txt", col_names = c("region", "chrName", "pos", "Nsites","fst"), skip = 1) %>%  mutate(correct_fst = ifelse(fst < 0, 0, fst)) %>%  
  left_join(chrom_df, by = "chrName") %>%
    subset(select = -c(chrName)) %>%
    dplyr::select(chr, pos, correct_fst)
```

    ## Rows: 4001653 Columns: 5
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (2): region, chrName
    ## dbl (3): pos, Nsites, fst
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
NBS_SEAK_fst <- read_tsv("/fs/cbsubscb16/storage/rkc/angsd/fst/NorthBering-SEAK.fst.txt", col_names = c("region", "chrName", "pos", "Nsites","fst"), skip = 1) %>%  mutate(correct_fst = ifelse(fst < 0, 0, fst)) %>%  
  left_join(chrom_df, by = "chrName") %>%
    subset(select = -c(chrName)) %>%
    dplyr::select(chr, pos, correct_fst)
```

    ## Rows: 2707622 Columns: 5
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (2): region, chrName
    ## dbl (3): pos, Nsites, fst
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

#### Test script for filtering fst files

``` r
AI_EBS_wgph_fst <- AI_EBS_fst %>% inner_join(AI_wgph_pos, by = c("chr" = "chr", "pos" = "pos")) %>% inner_join(EBS_wgph_pos, by = c("chr" = "chr", "pos" = "pos")) %>% dplyr::select(c(-totDepth.x,-totDepth.y))
#write_tsv(AI_EBS_wgph_fst, "/fs/cbsubscb16/storage/rkc/angsd/fst/AI-EBS_wgph.fst.txt")

AI_GOA_wgph_fst <- AI_GOA_fst %>% inner_join(AI_wgph_pos, by = c("chr" = "chr", "pos" = "pos")) %>% inner_join(GOA_wgph_pos, by = c("chr" = "chr", "pos" = "pos")) %>% dplyr::select(c(-totDepth.x,-totDepth.y))
#write_tsv(AI_GOA_wgph_fst, "/fs/cbsubscb16/storage/rkc/angsd/fst/AI-GOA_wgph.fst.txt")

NBS_AI_wgph_fst <- NBS_AI_fst %>% inner_join(AI_wgph_pos, by = c("chr" = "chr", "pos" = "pos")) %>% inner_join(NBS_wgph_pos, by = c("chr" = "chr", "pos" = "pos")) %>% dplyr::select(c(-totDepth.x,-totDepth.y))
#write_tsv(NBS_AI_wgph_fst, "/fs/cbsubscb16/storage/rkc/angsd/fst/AI-NBS_wgph.fst.txt")

AI_SEAK_wgph_fst <- AI_SEAK_fst %>% inner_join(AI_wgph_pos, by = c("chr" = "chr", "pos" = "pos")) %>% inner_join(SEAK_wgph_pos, by = c("chr" = "chr", "pos" = "pos")) %>% dplyr::select(c(-totDepth.x,-totDepth.y))
#write_tsv(AI_SEAK_wgph_fst, "/fs/cbsubscb16/storage/rkc/angsd/fst/AI-SEAK_wgph.fst.txt")

EBS_GOA_wgph_fst <- EBS_GOA_fst %>% inner_join(EBS_wgph_pos, by = c("chr" = "chr", "pos" = "pos")) %>% inner_join(GOA_wgph_pos, by = c("chr" = "chr", "pos" = "pos")) %>% dplyr::select(c(-totDepth.x,-totDepth.y))
#write_tsv(EBS_GOA_wgph_fst, "/fs/cbsubscb16/storage/rkc/angsd/fst/EBS-GOA_wgph.fst.txt")

NBS_EBS_wgph_fst <- NBS_EBS_fst %>% inner_join(EBS_wgph_pos, by = c("chr" = "chr", "pos" = "pos")) %>% inner_join(NBS_wgph_pos, by = c("chr" = "chr", "pos" = "pos")) %>% dplyr::select(c(-totDepth.x,-totDepth.y))
#write_tsv(NBS_EBS_wgph_fst, "/fs/cbsubscb16/storage/rkc/angsd/fst/EBS-NBS_wgph.fst.txt")

EBS_SEAK_wgph_fst <- EBS_SEAK_fst %>% inner_join(EBS_wgph_pos, by = c("chr" = "chr", "pos" = "pos")) %>% inner_join(SEAK_wgph_pos, by = c("chr" = "chr", "pos" = "pos")) %>% dplyr::select(c(-totDepth.x,-totDepth.y))
#write_tsv(EBS_SEAK_wgph_fst, "/fs/cbsubscb16/storage/rkc/angsd/fst/EBS-SEAK_wgph.fst.txt")

NBS_GOA_wgph_fst <- NBS_GOA_fst %>% inner_join(GOA_wgph_pos, by = c("chr" = "chr", "pos" = "pos")) %>% inner_join(NBS_wgph_pos, by = c("chr" = "chr", "pos" = "pos")) %>% dplyr::select(c(-totDepth.x,-totDepth.y))
#write_tsv(NBS_GOA_wgph_fst, "/fs/cbsubscb16/storage/rkc/angsd/fst/GOA-NBS_wgph.fst.txt")

GOA_SEAK_wgph_fst <- GOA_SEAK_fst %>% inner_join(GOA_wgph_pos, by = c("chr" = "chr", "pos" = "pos")) %>% inner_join(SEAK_wgph_pos, by = c("chr" = "chr", "pos" = "pos")) %>% dplyr::select(c(-totDepth.x,-totDepth.y))
#write_tsv(GOA_SEAK_wgph_fst, "/fs/cbsubscb16/storage/rkc/angsd/fst/GOA-SEAK_wgph.fst.txt")

NBS_SEAK_wgph_fst <- NBS_SEAK_fst %>% inner_join(NBS_wgph_pos, by = c("chr" = "chr", "pos" = "pos")) %>% inner_join(SEAK_wgph_pos, by = c("chr" = "chr", "pos" = "pos")) %>% dplyr::select(c(-totDepth.x,-totDepth.y))
#write_tsv(NBS_SEAK_wgph_fst, "/fs/cbsubscb16/storage/rkc/angsd/fst/NBS-SEAK_wgph.fst.txt")
```
