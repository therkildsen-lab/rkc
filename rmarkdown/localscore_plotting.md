Localscore plotting
================

## Load required packages

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

## Load data

``` r
pop_list <- c("NBS", "AI", "EBS", "GOA", "SEAK")
pairs_list <- c()
for(a in 1:4){
  for(b in 1:5){ # sorry for the nested for loops...
    tmp <- paste0("/fs/cbsubscb16/storage/rkc/localscore/output/",pop_list[a],"-",pop_list[b],"_xi2.ls") # local score values
    tmp2 <- paste0("/fs/cbsubscb16/storage/rkc/localscore/output/",pop_list[a],"-",pop_list[b],"_xi2_01sig.txt") # list of significant regions
    tmp3 <- paste0("/fs/cbsubscb16/storage/rkc/localscore/output/",pop_list[a],"-",pop_list[b],"_xi2_chr_sig.txt") # per chromosome significance values
    if(file.exists(tmp) == T){
      tmp3_df <- read_delim(tmp3, col_names = T, delim = " ") %>% dplyr::rename(chr = Chr)
      assign(paste0(pop_list[a],"_",pop_list[b],"_out"), read_delim(tmp, col_names = c("chr", "pos", "ls"), delim = " ") %>% left_join(tmp3_df, by = "chr")) # join significance values with local score values
      assign(paste0(pop_list[a],"_",pop_list[b],"_sig"), read_delim(tmp2, col_names = T, delim = "\t")) # read in significant regions
      
    }
  }
}

## chromosome names and numbers
chrom_df <- read.table("/fs/cbsubscb16/storage/rkc/sample_lists/chrom_meta_data.txt", header = TRUE, col.names = c("chrName","chrNewName")) %>% mutate(chr = seq(1,104,1)) %>% dplyr::select(-chrNewName)
```

## Filter out data for chrs with significant peaks

``` r
sig_list <- list(NBS_AI_sig, NBS_EBS_sig, NBS_GOA_sig, NBS_SEAK_sig, AI_EBS_sig, AI_GOA_sig, AI_SEAK_sig, EBS_GOA_sig, EBS_SEAK_sig, GOA_SEAK_sig)
out_list <- list(NBS_AI_out, NBS_EBS_out, NBS_GOA_out, NBS_SEAK_out, AI_EBS_out, AI_GOA_out, AI_SEAK_out, EBS_GOA_out, EBS_SEAK_out, GOA_SEAK_out)
name_list <- c("NBS_AI", "NBS_EBS", "NBS_GOA", "NBS_SEAK", "AI_EBS", "AI_GOA", "AI_SEAK", "EBS_GOA", "EBS_SEAK", "GOA_SEAK")

for(s in 1:10){
  tmp <- unique(sig_list[[s]]$chr)
  assign(paste0(name_list[s],"_outSig"), out_list[[s]] %>% filter(chr %in% tmp))
}

all_sig_regions <- rbind(sig_list[[1]], sig_list[[2]], sig_list[[3]], sig_list[[4]], sig_list[[5]], sig_list[[6]], sig_list[[7]], sig_list[[8]], sig_list[[9]], sig_list[[10]])

# write_tsv(all_sig_regions, "/fs/cbsubscb16/storage/rkc/localscore/all_sig_regions_wgph.tsv")
```

## Merge local score with Fst file

``` r
## Read in Fst files
## correct negative values of fst to 0
## rename chromosomes to chromosome numbers
## join fst to localscore dataset. 

AI_EBS_fst <- read_tsv("/fs/cbsubscb16/storage/rkc/angsd/fst/AI-EBS_wgph.fst.txt")
AI_EBS_outSig <- AI_EBS_outSig %>% left_join(AI_EBS_fst, by = c("chr" = "chr", "pos" = "pos"))

AI_GOA_fst <- read_tsv("/fs/cbsubscb16/storage/rkc/angsd/fst/AI-GOA_wgph.fst.txt")
AI_GOA_outSig <- AI_GOA_outSig %>% left_join(AI_GOA_fst, by = c("chr" = "chr", "pos" = "pos"))

NBS_AI_fst <- read_tsv("/fs/cbsubscb16/storage/rkc/angsd/fst/AI-NBS_wgph.fst.txt")
NBS_AI_outSig <- NBS_AI_outSig %>% left_join(NBS_AI_fst, by = c("chr" = "chr", "pos" = "pos"))

AI_SEAK_fst <- read_tsv("/fs/cbsubscb16/storage/rkc/angsd/fst/AI-SEAK_wgph.fst.txt")
AI_SEAK_outSig <- AI_SEAK_outSig %>% left_join(AI_SEAK_fst, by = c("chr" = "chr", "pos" = "pos"))

EBS_GOA_fst <- read_tsv("/fs/cbsubscb16/storage/rkc/angsd/fst/EBS-GOA_wgph.fst.txt")
EBS_GOA_outSig <- EBS_GOA_outSig %>% left_join(EBS_GOA_fst, by = c("chr" = "chr", "pos" = "pos"))

NBS_EBS_fst <- read_tsv("/fs/cbsubscb16/storage/rkc/angsd/fst/EBS-NBS_wgph.fst.txt")
NBS_EBS_outSig <- NBS_EBS_outSig %>% left_join(NBS_EBS_fst, by = c("chr" = "chr", "pos" = "pos"))

EBS_SEAK_fst <- read_tsv("/fs/cbsubscb16/storage/rkc/angsd/fst/EBS-SEAK_wgph.fst.txt")
EBS_SEAK_outSig <- EBS_SEAK_outSig %>% left_join(EBS_SEAK_fst, by = c("chr" = "chr", "pos" = "pos"))

NBS_GOA_fst <- read_tsv("/fs/cbsubscb16/storage/rkc/angsd/fst/GOA-NBS_wgph.fst.txt")
NBS_GOA_outSig <- NBS_GOA_outSig %>% left_join(NBS_GOA_fst, by = c("chr" = "chr", "pos" = "pos"))

GOA_SEAK_fst <- read_tsv("/fs/cbsubscb16/storage/rkc/angsd/fst/GOA-SEAK_wgph.fst.txt")
GOA_SEAK_outSig <- GOA_SEAK_outSig %>% left_join(GOA_SEAK_fst, by = c("chr" = "chr", "pos" = "pos"))

NBS_SEAK_fst <- read_tsv("/fs/cbsubscb16/storage/rkc/angsd/fst/NBS-SEAK_wgph.fst.txt")
NBS_SEAK_outSig <- NBS_SEAK_outSig %>% left_join(NBS_SEAK_fst, by = c("chr" = "chr", "pos" = "pos"))
```

## Facet plot each pair

#### NBS vs others

``` r
p_NBS_AI <- ggplot(NBS_AI_outSig, aes(x = pos/10^6, y = ls)) + facet_wrap(chr ~.) +
  coord_cartesian(ylim=c(0, max(NBS_AI_outSig$sig.01*1.1)), xlim = c(0,92), expand = F) +
  xlab("position (Mbp)") +
  ylab("local score") +
  cowplot::theme_cowplot() +
  geom_line(linewidth=0.4, alpha=0.5) +
  geom_point(aes(y = correct_fst*max(NBS_AI_outSig$sig.01*1.1)), alpha = 0.3, color = "dodgerblue") +
  scale_y_continuous(sec.axis = sec_axis(~./max(NBS_AI_outSig$sig.01*1.1), name = "Fst")) +
  geom_hline(aes(yintercept = NBS_AI_outSig$sig.01), color = "blue", lty = 2)
p_NBS_AI
```

    ## Warning: Removed 8 rows containing missing values (`geom_point()`).

![](localscore_plotting_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
p_NBS_EBS <- ggplot(NBS_EBS_outSig, aes(x = pos/10^6, y = ls)) + facet_wrap(chr ~.) +
  coord_cartesian(ylim=c(0, max(NBS_EBS_outSig$sig.01*1.1)), xlim = c(0,92), expand = F) +
  xlab("position (Mbp)") +
  ylab("local score") +
  cowplot::theme_cowplot() +
  geom_line(linewidth=0.4, alpha=0.5) +
  geom_point(aes(y = correct_fst*max(NBS_EBS_outSig$sig.01*1.1)), alpha = 0.3, color = "dodgerblue") +
  scale_y_continuous(sec.axis = sec_axis(~./max(NBS_EBS_outSig$sig.01*1.1), name = "Fst")) +
  geom_hline(aes(yintercept = NBS_EBS_outSig$sig.01), color = "blue", lty = 2)
p_NBS_EBS
```

    ## Warning: Removed 14 rows containing missing values (`geom_point()`).

![](localscore_plotting_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->

``` r
p_NBS_GOA <- ggplot(NBS_GOA_outSig, aes(x = pos/10^6, y = ls)) + facet_wrap(chr ~.) +
  coord_cartesian(ylim=c(0, max(NBS_GOA_outSig$sig.01*1.1)), xlim = c(0,92), expand = F) +
  xlab("position (Mbp)") +
  ylab("local score") +
  cowplot::theme_cowplot() +
  geom_line(linewidth=0.4, alpha=0.5) +
  geom_point(aes(y = correct_fst*max(NBS_GOA_outSig$sig.01*1.1)), alpha = 0.3, color = "dodgerblue") +
  scale_y_continuous(sec.axis = sec_axis(~./max(NBS_GOA_outSig$sig.01*1.1), name = "Fst")) +
  geom_hline(aes(yintercept = NBS_GOA_outSig$sig.01), color = "blue", lty = 2)
p_NBS_GOA
```

    ## Warning: Removed 24 rows containing missing values (`geom_point()`).

![](localscore_plotting_files/figure-gfm/unnamed-chunk-5-3.png)<!-- -->

``` r
p_NBS_SEAK <- ggplot(NBS_SEAK_outSig, aes(x = pos/10^6, y = ls)) + facet_wrap(chr ~.) +
  coord_cartesian(ylim=c(0, max(NBS_SEAK_outSig$sig.01*1.1)), xlim = c(0,92), expand = F) +
  xlab("position (Mbp)") +
  ylab("local score") +
  cowplot::theme_cowplot() +
  geom_line(linewidth=0.4, alpha=0.5) +
  geom_point(aes(y = correct_fst*max(NBS_SEAK_outSig$sig.01*1.1)), alpha = 0.3, color = "dodgerblue") +
  scale_y_continuous(sec.axis = sec_axis(~./max(NBS_SEAK_outSig$sig.01*1.1), name = "Fst")) +
  geom_hline(aes(yintercept = NBS_SEAK_outSig$sig.01), color = "blue", lty = 2)
p_NBS_SEAK
```

    ## Warning: Removed 9 rows containing missing values (`geom_point()`).

![](localscore_plotting_files/figure-gfm/unnamed-chunk-5-4.png)<!-- -->

#### AI vs others

``` r
p_AI_EBS <- ggplot(AI_EBS_outSig, aes(x = pos/10^6, y = ls)) + facet_wrap(chr ~.) +
  coord_cartesian(ylim=c(0, max(AI_EBS_outSig$sig.01*1.1)), xlim = c(0,92), expand = F) +
  xlab("position (Mbp)") +
  ylab("local score") +
  cowplot::theme_cowplot() +
  geom_point(aes(y = correct_fst*max(AI_EBS_outSig$sig.01*1.1)), alpha = 0.3, color = "dodgerblue") +
  geom_line(linewidth=0.4, alpha=0.5) +
  scale_y_continuous(sec.axis = sec_axis(~./max(AI_EBS_outSig$sig.01*1.1), name = "Fst")) +
  geom_hline(aes(yintercept = AI_EBS_outSig$sig.01), color = "blue", lty = 2)
p_AI_EBS
```

    ## Warning: Removed 12 rows containing missing values (`geom_point()`).

![](localscore_plotting_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
p_AI_GOA <- ggplot(AI_GOA_outSig, aes(x = pos/10^6, y = ls)) + facet_wrap(chr ~.) +
  coord_cartesian(ylim=c(0, max(AI_GOA_outSig$sig.01*1.1)), xlim = c(0,92), expand = F) +
  xlab("position (Mbp)") +
  ylab("local score") +
  cowplot::theme_cowplot() +
  geom_point(aes(y = correct_fst*max(AI_GOA_outSig$sig.01*1.1)), alpha = 0.3, color = "dodgerblue") +
  geom_line(linewidth=0.4, alpha=0.5) +
  scale_y_continuous(sec.axis = sec_axis(~./max(AI_GOA_outSig$sig.01*1.1), name = "Fst")) +
  geom_hline(aes(yintercept = AI_GOA_outSig$sig.01), color = "blue", lty = 2)
p_AI_GOA
```

    ## Warning: Removed 10 rows containing missing values (`geom_point()`).

![](localscore_plotting_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->

``` r
p_AI_SEAK <- ggplot(AI_SEAK_outSig, aes(x = pos/10^6, y = ls)) + facet_wrap(chr ~.) +
  coord_cartesian(ylim=c(0, max(AI_SEAK_outSig$sig.01*1.1)), xlim = c(0,92), expand = F) +
  xlab("position (Mbp)") +
  ylab("local score") +
  cowplot::theme_cowplot() +
  geom_point(aes(y = correct_fst*max(AI_SEAK_outSig$sig.01*1.1)), alpha = 0.3, color = "dodgerblue") +
  geom_line(linewidth=0.4, alpha=0.5) +
  scale_y_continuous(sec.axis = sec_axis(~./max(AI_SEAK_outSig$sig.01*1.1), name = "Fst")) +
  geom_hline(aes(yintercept = AI_SEAK_outSig$sig.01), color = "blue", lty = 2)
p_AI_SEAK
```

    ## Warning: Removed 9 rows containing missing values (`geom_point()`).

![](localscore_plotting_files/figure-gfm/unnamed-chunk-6-3.png)<!-- -->

#### EBS vs others

``` r
p_EBS_GOA <- ggplot(EBS_GOA_outSig, aes(x = pos/10^6, y = ls)) + facet_wrap(chr ~.) +
  coord_cartesian(ylim=c(0, max(EBS_GOA_outSig$sig.01*1.1)), xlim = c(0,92), expand = F) +
  xlab("position (Mbp)") +
  ylab("local score") +
  cowplot::theme_cowplot() +
  geom_point(aes(y = correct_fst*max(EBS_GOA_outSig$sig.01*1.1)), alpha = 0.3, color = "dodgerblue") +
  geom_line(linewidth=0.4, alpha=0.5) +
  scale_y_continuous(sec.axis = sec_axis(~./max(EBS_GOA_outSig$sig.01*1.1), name = "Fst")) +
  geom_hline(aes(yintercept = EBS_GOA_outSig$sig.01), color = "blue", lty = 2)
p_EBS_GOA
```

    ## Warning: Removed 8 rows containing missing values (`geom_point()`).

![](localscore_plotting_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
p_EBS_SEAK <- ggplot(EBS_SEAK_outSig, aes(x = pos/10^6, y = ls)) + facet_wrap(chr ~.) +
  coord_cartesian(ylim=c(0, max(EBS_SEAK_outSig$sig.01*1.1)), xlim = c(0,92), expand = F) +
  xlab("position (Mbp)") +
  ylab("local score") +
  cowplot::theme_cowplot() +
  geom_point(aes(y = correct_fst*max(EBS_SEAK_outSig$sig.01*1.1)), alpha = 0.3, color = "dodgerblue") +
  geom_line(linewidth=0.4, alpha=0.5) +
  scale_y_continuous(sec.axis = sec_axis(~./max(EBS_SEAK_outSig$sig.01*1.1), name = "Fst")) +
  geom_hline(aes(yintercept = EBS_SEAK_outSig$sig.01), color = "blue", lty = 2)
p_EBS_SEAK
```

    ## Warning: Removed 14 rows containing missing values (`geom_point()`).

![](localscore_plotting_files/figure-gfm/unnamed-chunk-7-2.png)<!-- -->

## GOA vs others

``` r
p_GOA_SEAK <- ggplot(GOA_SEAK_outSig, aes(x = pos/10^6, y = ls)) + facet_wrap(chr ~.) +
  coord_cartesian(ylim=c(0, max(GOA_SEAK_outSig$sig.01*1.1)), xlim = c(0,92), expand = F) +
  xlab("position (Mbp)") +
  ylab("local score") +
  cowplot::theme_cowplot() +
  geom_point(aes(y = correct_fst*max(GOA_SEAK_outSig$sig.01*1.1)), alpha = 0.3, color = "dodgerblue") +
  geom_line(linewidth=0.4, alpha=0.5) +
  scale_y_continuous(sec.axis = sec_axis(~./max(GOA_SEAK_outSig$sig.01*1.1), name = "Fst")) +
  geom_hline(aes(yintercept = GOA_SEAK_outSig$sig.01), color = "blue", lty = 2)
p_GOA_SEAK
```

    ## Warning: Removed 9 rows containing missing values (`geom_point()`).

![](localscore_plotting_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

## Export all plots as png

``` r
ggsave("/fs/cbsubscb16/storage/rkc/localscore/plots/NBS_AI_localscore_sig.png", plot = p_NBS_AI, device = "png", height = 2.5, width = 7, units = "in")
```

    ## Warning: Removed 8 rows containing missing values (`geom_point()`).

``` r
ggsave("/fs/cbsubscb16/storage/rkc/localscore/plots/NBS_EBS_localscore_sig.png", plot = p_NBS_EBS, device = "png", height = 5, width = 7, units = "in")
```

    ## Warning: Removed 14 rows containing missing values (`geom_point()`).

``` r
ggsave("/fs/cbsubscb16/storage/rkc/localscore/plots/NBS_GOA_localscore_sig.png", plot = p_NBS_GOA, device = "png", height = 5, width = 7, units = "in")
```

    ## Warning: Removed 24 rows containing missing values (`geom_point()`).

``` r
ggsave("/fs/cbsubscb16/storage/rkc/localscore/plots/NBS_SEAK_localscore_sig.png", plot = p_NBS_SEAK, device = "png", height = 2.5, width = 7, units = "in")
```

    ## Warning: Removed 9 rows containing missing values (`geom_point()`).

``` r
ggsave("/fs/cbsubscb16/storage/rkc/localscore/plots/AI_EBS_localscore_sig.png", plot = p_AI_EBS, device = "png", height = 5, width = 7, units = "in")
```

    ## Warning: Removed 12 rows containing missing values (`geom_point()`).

``` r
ggsave("/fs/cbsubscb16/storage/rkc/localscore/plots/AI_GOA_localscore_sig.png", plot = p_AI_GOA, device = "png", height = 5, width = 7, units = "in")
```

    ## Warning: Removed 10 rows containing missing values (`geom_point()`).

``` r
ggsave("/fs/cbsubscb16/storage/rkc/localscore/plots/AI_SEAK_localscore_sig.png", plot = p_AI_SEAK, device = "png", height = 5, width = 5.5, units = "in")
```

    ## Warning: Removed 9 rows containing missing values (`geom_point()`).

``` r
ggsave("/fs/cbsubscb16/storage/rkc/localscore/plots/EBS_GOA_localscore_sig.png", plot = p_EBS_GOA, device = "png", height = 2.5, width = 5.5, units = "in")
```

    ## Warning: Removed 8 rows containing missing values (`geom_point()`).

``` r
ggsave("/fs/cbsubscb16/storage/rkc/localscore/plots/EBS_SEAK_localscore_sig.png", plot = p_EBS_SEAK, device = "png", height = 2.5, width = 5.5, units = "in")
```

    ## Warning: Removed 14 rows containing missing values (`geom_point()`).

``` r
ggsave("/fs/cbsubscb16/storage/rkc/localscore/plots/GOA_SEAK_localscore_sig.png", plot = p_GOA_SEAK, device = "png", height = 2.5, width = 7, units = "in")
```

    ## Warning: Removed 9 rows containing missing values (`geom_point()`).

## Plot all interesting chromosomes

#### Prepare data

``` r
tot_sig <- unique(c(sig_list[[1]]$chr, sig_list[[2]]$chr, sig_list[[3]]$chr, sig_list[[4]]$chr, sig_list[[5]]$chr, sig_list[[6]]$chr, sig_list[[7]]$chr, sig_list[[8]]$chr, sig_list[[9]]$chr,sig_list[[10]]$chr)) %>% sort()

## filter all out files for all significant chr
for(s in 1:10){
  assign(paste0(name_list[s],"_outSigAll"), out_list[[s]] %>% filter(chr %in% tot_sig) %>% mutate(outlier = F))
}

## create outlier column to id outliers in each population
for(s in 1:length(NBS_AI_sig$chr)){
  NBS_AI_outSigAll$outlier <- ifelse(NBS_AI_outSigAll$chr == NBS_AI_sig$chr[s] & NBS_AI_outSigAll$pos > NBS_AI_sig$beg[s] & NBS_AI_outSigAll$pos < NBS_AI_sig$end[s], T, NBS_AI_outSigAll$outlier)
}

for(s in 1:length(NBS_EBS_sig$chr)){
  NBS_EBS_outSigAll$outlier <- ifelse(NBS_EBS_outSigAll$chr == NBS_EBS_sig$chr[s] & NBS_EBS_outSigAll$pos > NBS_EBS_sig$beg[s] & NBS_EBS_outSigAll$pos < NBS_EBS_sig$end[s], T, NBS_EBS_outSigAll$outlier)
}

for(s in 1:length(NBS_GOA_sig$chr)){
  NBS_GOA_outSigAll$outlier <- ifelse(NBS_GOA_outSigAll$chr == NBS_GOA_sig$chr[s] & NBS_GOA_outSigAll$pos > NBS_GOA_sig$beg[s] & NBS_GOA_outSigAll$pos < NBS_GOA_sig$end[s], T, NBS_GOA_outSigAll$outlier)
}

for(s in 1:length(NBS_SEAK_sig$chr)){
  NBS_SEAK_outSigAll$outlier <- ifelse(NBS_SEAK_outSigAll$chr == NBS_SEAK_sig$chr[s] & NBS_SEAK_outSigAll$pos > NBS_SEAK_sig$beg[s] & NBS_SEAK_outSigAll$pos < NBS_SEAK_sig$end[s], T, NBS_SEAK_outSigAll$outlier)
}

for(s in 1:length(AI_EBS_sig$chr)){
  AI_EBS_outSigAll$outlier <- ifelse(AI_EBS_outSigAll$chr == AI_EBS_sig$chr[s] & AI_EBS_outSigAll$pos > AI_EBS_sig$beg[s] & AI_EBS_outSigAll$pos < AI_EBS_sig$end[s], T, AI_EBS_outSigAll$outlier)
}

for(s in 1:length(AI_GOA_sig$chr)){
  AI_GOA_outSigAll$outlier <- ifelse(AI_GOA_outSigAll$chr == AI_GOA_sig$chr[s] & AI_GOA_outSigAll$pos > AI_GOA_sig$beg[s] & AI_GOA_outSigAll$pos < AI_GOA_sig$end[s], T, AI_GOA_outSigAll$outlier)
}

for(s in 1:length(AI_SEAK_sig$chr)){
  AI_SEAK_outSigAll$outlier <- ifelse(AI_SEAK_outSigAll$chr == AI_SEAK_sig$chr[s] & AI_SEAK_outSigAll$pos > AI_SEAK_sig$beg[s] & AI_SEAK_outSigAll$pos < AI_SEAK_sig$end[s], T, AI_SEAK_outSigAll$outlier)
}

for(s in 1:length(EBS_GOA_sig$chr)){
  EBS_GOA_outSigAll$outlier <- ifelse(EBS_GOA_outSigAll$chr == EBS_GOA_sig$chr[s] & EBS_GOA_outSigAll$pos > EBS_GOA_sig$beg[s] & EBS_GOA_outSigAll$pos < EBS_GOA_sig$end[s], T, EBS_GOA_outSigAll$outlier)
}

for(s in 1:length(EBS_SEAK_sig$chr)){
  EBS_SEAK_outSigAll$outlier <- ifelse(EBS_SEAK_outSigAll$chr == EBS_SEAK_sig$chr[s] & EBS_SEAK_outSigAll$pos > EBS_SEAK_sig$beg[s] & EBS_SEAK_outSigAll$pos < EBS_SEAK_sig$end[s], T, EBS_SEAK_outSigAll$outlier)
}

for(s in 1:length(GOA_SEAK_sig$chr)){
  GOA_SEAK_outSigAll$outlier <- ifelse(GOA_SEAK_outSigAll$chr == GOA_SEAK_sig$chr[s] & GOA_SEAK_outSigAll$pos > GOA_SEAK_sig$beg[s] & GOA_SEAK_outSigAll$pos < GOA_SEAK_sig$end[s], T, GOA_SEAK_outSigAll$outlier)
}

## Merge with Fst
AI_EBS_outSigAll <- AI_EBS_outSigAll %>% left_join(AI_EBS_fst, by = c("chr" = "chr", "pos" = "pos"))
AI_GOA_outSigAll <- AI_GOA_outSigAll %>% left_join(AI_GOA_fst, by = c("chr" = "chr", "pos" = "pos"))
NBS_AI_outSigAll <- NBS_AI_outSigAll %>% left_join(NBS_AI_fst, by = c("chr" = "chr", "pos" = "pos"))
AI_SEAK_outSigAll <- AI_SEAK_outSigAll %>% left_join(AI_SEAK_fst, by = c("chr" = "chr", "pos" = "pos"))
EBS_GOA_outSigAll <- EBS_GOA_outSigAll %>% left_join(EBS_GOA_fst, by = c("chr" = "chr", "pos" = "pos"))
NBS_EBS_outSigAll <- NBS_EBS_outSigAll %>% left_join(NBS_EBS_fst, by = c("chr" = "chr", "pos" = "pos"))
EBS_SEAK_outSigAll <- EBS_SEAK_outSigAll %>% left_join(EBS_SEAK_fst, by = c("chr" = "chr", "pos" = "pos"))
NBS_GOA_outSigAll <- NBS_GOA_outSigAll %>% left_join(NBS_GOA_fst, by = c("chr" = "chr", "pos" = "pos"))
GOA_SEAK_outSigAll <- GOA_SEAK_outSigAll %>% left_join(GOA_SEAK_fst, by = c("chr" = "chr", "pos" = "pos"))
NBS_SEAK_outSigAll <- NBS_SEAK_outSigAll %>% left_join(NBS_SEAK_fst, by = c("chr" = "chr", "pos" = "pos"))

## How many outliers for each population 
AI_sigs <- rbind(AI_EBS_sig, AI_GOA_sig, NBS_AI_sig, AI_SEAK_sig) %>% arrange(chr)
EBS_sigs <- rbind(AI_EBS_sig, EBS_GOA_sig, NBS_EBS_sig, EBS_SEAK_sig) %>% arrange(chr)
NBS_sigs <- rbind(NBS_AI_sig, NBS_EBS_sig, NBS_GOA_sig, NBS_SEAK_sig) %>% arrange(chr)
GOA_sigs <- rbind(GOA_SEAK_sig, EBS_GOA_sig, AI_GOA_sig, NBS_GOA_sig) %>% arrange(chr)
SEAK_sigs <- rbind(AI_SEAK_sig, EBS_SEAK_sig, GOA_SEAK_sig, NBS_SEAK_sig) %>% arrange(chr)
```

## Plot Fst manhattan plots

``` r
get_cum_pos <- function(x){
  x_new <- group_by(x, chr) %>%
    # Compute lg size
    dplyr::summarise(lg_length=max(pos)) %>%
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(lg_length)-lg_length) %>%
    # Add this info to the initial dataset
    left_join(x, ., by=c("chr")) %>%
    # Add a cumulative position of each SNP
    dplyr::arrange(chr, pos) %>%
    group_by(chr) %>%
    mutate(pos_cum=pos+tot)
  return(x_new)
}

get_cum_pos_test <- function(x){
  x_new <- group_by(x, chr) %>%
    # Compute lg size
    dplyr::summarise(lg_length=max(pos)) %>%
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(lg_length)-lg_length)#%>%
    # Add this info to the initial dataset
    #left_join(x, ., by=c("chr")) %>%
    # Add a cumulative position of each SNP
    #arrange(chr, pos) %>%
    #group_by(chr) %>%
    #mutate(pos_cum=pos+tot)
  return(x_new)
}

# add option to define y  limit #    
my_color = c("#7fcdbb", "#2c7fb8")
plot_gwas_per_snp <- function(x){
  ggplot(x, aes(x=pos_cum/10^6, y=correct_fst, color=as.factor(chr))) +
    geom_point(size=0.2, alpha=0.5) +
    scale_color_manual(values=rep(my_color, 100000)) +
    geom_point(data = x %>% filter(outlier == T), aes(x = pos_cum/10^6, y = correct_fst), color = "red", size = 0.4) +
    coord_cartesian(ylim=c(0, 0.8), expand = F) +
    scale_y_continuous(breaks=c(0,0.3,0.6)) +
    cowplot::theme_cowplot() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line = element_line(linewidth = 0.2),
          legend.position="none",
          text = element_text(size=20),
          axis.text = element_text(size=8),
          strip.text.x = element_blank())
}



get_cum_pos(NBS_AI_outSigAll) %>% plot_gwas_per_snp()
```

    ## Warning: Removed 62 rows containing missing values (`geom_point()`).

![](localscore_plotting_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
ggsave("/fs/cbsubscb16/storage/rkc/figures/NBS_AI_ls_wgph_fst_outlier_scan.png", width = 7.5, height = 1, units = "in")
```

    ## Warning: Removed 62 rows containing missing values (`geom_point()`).

``` r
get_cum_pos(NBS_EBS_outSigAll) %>% plot_gwas_per_snp()
```

    ## Warning: Removed 64 rows containing missing values (`geom_point()`).

![](localscore_plotting_files/figure-gfm/unnamed-chunk-11-2.png)<!-- -->

``` r
ggsave("/fs/cbsubscb16/storage/rkc/figures/NBS_EBS_ls_wgph_fst_outlier_scan.png", width = 7.5, height = 1, units = "in")
```

    ## Warning: Removed 64 rows containing missing values (`geom_point()`).

``` r
get_cum_pos(NBS_GOA_outSigAll) %>% plot_gwas_per_snp()
```

    ## Warning: Removed 72 rows containing missing values (`geom_point()`).

![](localscore_plotting_files/figure-gfm/unnamed-chunk-11-3.png)<!-- -->

``` r
ggsave("/fs/cbsubscb16/storage/rkc/figures/NBS_GOA_ls_wgph_fst_outlier_scan.png", width = 7.5, height = 1, units = "in")
```

    ## Warning: Removed 72 rows containing missing values (`geom_point()`).

``` r
get_cum_pos(NBS_SEAK_outSigAll) %>% plot_gwas_per_snp()
```

    ## Warning: Removed 61 rows containing missing values (`geom_point()`).

![](localscore_plotting_files/figure-gfm/unnamed-chunk-11-4.png)<!-- -->

``` r
ggsave("/fs/cbsubscb16/storage/rkc/figures/NBS_SEAK_ls_wgph_fst_outlier_scan.png", width = 7.5, height = 1, units = "in")
```

    ## Warning: Removed 61 rows containing missing values (`geom_point()`).

``` r
get_cum_pos(AI_EBS_outSigAll) %>% plot_gwas_per_snp()
```

    ## Warning: Removed 61 rows containing missing values (`geom_point()`).

![](localscore_plotting_files/figure-gfm/unnamed-chunk-11-5.png)<!-- -->

``` r
ggsave("/fs/cbsubscb16/storage/rkc/figures/AI_EBS_ls_wgph_fst_outlier_scan.png", width = 7.5, height = 1, units = "in")
```

    ## Warning: Removed 61 rows containing missing values (`geom_point()`).

``` r
get_cum_pos(AI_GOA_outSigAll) %>% plot_gwas_per_snp()
```

    ## Warning: Removed 69 rows containing missing values (`geom_point()`).

![](localscore_plotting_files/figure-gfm/unnamed-chunk-11-6.png)<!-- -->

``` r
ggsave("/fs/cbsubscb16/storage/rkc/figures/AI_GOA_ls_wgph_fst_outlier_scan.png", width = 7.5, height = 1, units = "in")
```

    ## Warning: Removed 69 rows containing missing values (`geom_point()`).

``` r
get_cum_pos(AI_SEAK_outSigAll) %>% plot_gwas_per_snp()
```

    ## Warning: Removed 61 rows containing missing values (`geom_point()`).

![](localscore_plotting_files/figure-gfm/unnamed-chunk-11-7.png)<!-- -->

``` r
ggsave("/fs/cbsubscb16/storage/rkc/figures/AI_SEAK_ls_wgph_fst_outlier_scan.png", width = 7.5, height = 1, units = "in")
```

    ## Warning: Removed 61 rows containing missing values (`geom_point()`).

``` r
get_cum_pos(EBS_GOA_outSigAll) %>% plot_gwas_per_snp()
```

    ## Warning: Removed 72 rows containing missing values (`geom_point()`).

![](localscore_plotting_files/figure-gfm/unnamed-chunk-11-8.png)<!-- -->

``` r
ggsave("/fs/cbsubscb16/storage/rkc/figures/EBS_GOA_ls_wgph_fst_outlier_scan.png", width = 7.5, height = 1, units = "in")
```

    ## Warning: Removed 72 rows containing missing values (`geom_point()`).

``` r
get_cum_pos(EBS_SEAK_outSigAll) %>% plot_gwas_per_snp()
```

    ## Warning: Removed 61 rows containing missing values (`geom_point()`).

![](localscore_plotting_files/figure-gfm/unnamed-chunk-11-9.png)<!-- -->

``` r
ggsave("/fs/cbsubscb16/storage/rkc/figures/EBS_SEAK_ls_wgph_fst_outlier_scan.png", width = 7.5, height = 1, units = "in")
```

    ## Warning: Removed 61 rows containing missing values (`geom_point()`).

``` r
get_cum_pos(GOA_SEAK_outSigAll) %>% plot_gwas_per_snp()
```

    ## Warning: Removed 67 rows containing missing values (`geom_point()`).

![](localscore_plotting_files/figure-gfm/unnamed-chunk-11-10.png)<!-- -->

``` r
ggsave("/fs/cbsubscb16/storage/rkc/figures/GOA_SEAK_ls_wgph_fst_outlier_scan.png", width = 7.5, height = 1, units = "in")
```

    ## Warning: Removed 67 rows containing missing values (`geom_point()`).

## plot chr 100 for all pops

``` r
p100_NBS_GOA <- ggplot(NBS_GOA_outSigAll %>% filter(chr == 100), aes(x = pos/10^6, y = ls)) + 
  coord_cartesian(ylim=c(0, max(NBS_GOA_outSig$sig.01*1.1)), expand = F, xlim = c(0,67)) +
  xlab("position (Mbp)") +
  ylab("local score") +
  cowplot::theme_cowplot() +
  geom_point(aes(y = correct_fst*max(NBS_GOA_outSig$sig.01*1.1)), alpha = 0.3, color = "#2c7fb8", size = 3) + 
  geom_line(linewidth=0.4) +
  geom_point(data = NBS_GOA_outSigAll %>% filter(chr == 100 & outlier == T), aes(x = pos/10^6, y = correct_fst*max(NBS_GOA_outSig$sig.01*1.1)), color = "red", size = 3) +
  geom_line(linewidth=0.4) +
  scale_y_continuous(sec.axis = sec_axis(~./max(NBS_GOA_outSig$sig.01*1.1), name = "Fst")) +
  geom_hline(aes(yintercept = NBS_GOA_outSigAll$sig.01[1]), color = "blue", lty = 2) +
  theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position="none",
          text = element_text(size=20),
          axis.text = element_text(size=20),
          strip.text.x = element_blank())
p100_NBS_GOA
```

    ## Warning: Removed 3 rows containing missing values (`geom_point()`).

![](localscore_plotting_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

``` r
ggsave("/fs/cbsubscb16/storage/rkc/figures/NBS_GOA_ls_wgph_fst_outlier_scan_100.png", width = 9, height = 3, units = "in")
```

    ## Warning: Removed 3 rows containing missing values (`geom_point()`).

``` r
p100_AI_GOA <- ggplot(AI_GOA_outSigAll %>% filter(chr == 100), aes(x = pos/10^6, y = ls)) + 
  coord_cartesian(ylim=c(0, max(AI_GOA_outSig$sig.01*1.1)), expand = F, xlim = c(0,67)) +
  xlab("position (Mbp)") +
  ylab("local score") +
  cowplot::theme_cowplot() +
  geom_point(aes(y = correct_fst*max(AI_GOA_outSig$sig.01*1.1)), alpha = 0.3, color = "#2c7fb8", size = 3) + 
  geom_line(linewidth=0.4) +
  geom_point(data = AI_GOA_outSigAll %>% filter(chr == 100 & outlier == T), aes(x = pos/10^6, y = correct_fst*max(AI_GOA_outSig$sig.01*1.1)), color = "red", size = 3) +
  geom_line(linewidth=0.4) +
  scale_y_continuous(sec.axis = sec_axis(~./max(AI_GOA_outSig$sig.01*1.1), name = "Fst")) +
  geom_hline(aes(yintercept = AI_GOA_outSigAll$sig.01[1]), color = "blue", lty = 2) +
  theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position="none",
          text = element_text(size=20),
          axis.text = element_text(size=20),
          strip.text.x = element_blank())
p100_AI_GOA
```

    ## Warning: Removed 3 rows containing missing values (`geom_point()`).

![](localscore_plotting_files/figure-gfm/unnamed-chunk-12-2.png)<!-- -->

``` r
ggsave("/fs/cbsubscb16/storage/rkc/figures/AI_GOA_ls_wgph_fst_outlier_scan_100.png", width = 9, height = 3, units = "in")
```

    ## Warning: Removed 3 rows containing missing values (`geom_point()`).

``` r
p100_EBS_GOA <- ggplot(EBS_GOA_outSigAll %>% filter(chr == 100), aes(x = pos/10^6, y = ls)) + 
  coord_cartesian(ylim=c(0, max(EBS_GOA_outSig$sig.01*1.1)), expand = F, xlim = c(0,67)) +
  xlab("position (Mbp)") +
  ylab("local score") +
  cowplot::theme_cowplot() +
  geom_point(aes(y = correct_fst*max(EBS_GOA_outSig$sig.01*1.1)), alpha = 0.3, color = "#2c7fb8", size = 3) + 
  geom_line(linewidth=0.4) +
  geom_point(data = EBS_GOA_outSigAll %>% filter(chr == 100 & outlier == T), aes(x = pos/10^6, y = correct_fst*max(EBS_GOA_outSig$sig.01*1.1)), color = "red", size = 3) +
  geom_line(linewidth=0.3) +
  scale_y_continuous(sec.axis = sec_axis(~./max(EBS_GOA_outSig$sig.01*1.1), name = "Fst")) +
  geom_hline(aes(yintercept = EBS_GOA_outSigAll$sig.01[1]), color = "blue", lty = 2) +
  theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position="none",
          text = element_text(size=20),
          axis.text = element_text(size=20),
          strip.text.x = element_blank())
p100_EBS_GOA
```

    ## Warning: Removed 4 rows containing missing values (`geom_point()`).

![](localscore_plotting_files/figure-gfm/unnamed-chunk-12-3.png)<!-- -->

``` r
ggsave("/fs/cbsubscb16/storage/rkc/figures/EBS_GOA_ls_wgph_fst_outlier_scan_100.png", width = 9, height = 3, units = "in")
```

    ## Warning: Removed 4 rows containing missing values (`geom_point()`).

``` r
p100_GOA_SEAK <- ggplot(GOA_SEAK_outSigAll %>% filter(chr == 100), aes(x = pos/10^6, y = ls)) + 
  coord_cartesian(ylim=c(0, max(GOA_SEAK_outSig$sig.01*1.1)), expand = F, xlim = c(0,67)) +
  xlab("position (Mbp)") +
  ylab("local score") +
  cowplot::theme_cowplot() +
  geom_point(aes(y = correct_fst*max(GOA_SEAK_outSig$sig.01*1.1)), alpha = 0.3, color = "#2c7fb8", size = 3) + 
  geom_line(linewidth=0.4) +
  geom_point(data = GOA_SEAK_outSigAll %>% filter(chr == 100 & outlier == T), aes(x = pos/10^6, y = correct_fst*max(GOA_SEAK_outSig$sig.01*1.1)), color = "red", size = 3) +
  geom_line(linewidth=0.4) +
  scale_y_continuous(sec.axis = sec_axis(~./max(GOA_SEAK_outSig$sig.01*1.1), name = "Fst")) +
  geom_hline(aes(yintercept = GOA_SEAK_outSigAll$sig.01[1]), color = "blue", lty = 2) +
  theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position="none",
          text = element_text(size=20),
          axis.text = element_text(size=20),
          strip.text.x = element_blank())
p100_GOA_SEAK
```

    ## Warning: Removed 3 rows containing missing values (`geom_point()`).

![](localscore_plotting_files/figure-gfm/unnamed-chunk-12-4.png)<!-- -->

``` r
ggsave("/fs/cbsubscb16/storage/rkc/figures/GOA_SEAK_ls_wgph_fst_outlier_scan_100.png", width =9, height = 3, units = "in")
```

    ## Warning: Removed 3 rows containing missing values (`geom_point()`).
