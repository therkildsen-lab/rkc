Candidate region annotation
================

## Annotation hits for candiate regions

- We will search only within the candidate regions

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
library(knitr)
library(Rgff)
library(ape)
```

    ## 
    ## Attaching package: 'ape'
    ## 
    ## The following object is masked from 'package:dplyr':
    ## 
    ##     where

``` r
gff <- read.gff("/fs/cbsubscb16/storage/rkc/genome/king_crab_annotation.gff3") %>% tibble()
gff <- gff %>% mutate(chr = as.numeric(str_sub(seqid,14,18), .keep = "unused", .before = source)) %>% filter(chr <= 104) %>% arrange(chr)
outliers_list <- read_tsv("/fs/cbsubscb16/storage/rkc/localscore/output/sig01.txt")
```

    ## Warning: One or more parsing issues, call `problems()` on your data frame for details,
    ## e.g.:
    ##   dat <- vroom(...)
    ##   problems(dat)

    ## Rows: 41 Columns: 6
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (2): pop1, pop2
    ## dbl (4): chr, beg, end, peak
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
gff_outliers <- data.frame(seqid=character(),
                           beg=double(),
                           end=double(),
                           hit=character())

for(i in 1:length(outliers_list$chr)){
  tmp <- gff %>% filter(chr == outliers_list$chr[i] & start >= outliers_list$beg[i] & end <= outliers_list$end[i])
      gff_outliers <- rbind(gff_outliers,tmp)
      #print(lambda_tab[i,])
}

# gff_outliers %>% select(-chr) %>% write_tsv("/fs/cbsubscb16/storage/rkc/genome/RKC_annotation_outliers.gff3")
```

#### Get only parent feature

``` r
gff_outliers_p <- gff_outliers %>% filter(!grepl("Parent",attributes))
kable(gff_outliers_p)
```

| seqid            | source | type |    start |      end | score | strand | phase | attributes                                                                  | chr |
|:-----------------|:-------|:-----|---------:|---------:|------:|:-------|:------|:----------------------------------------------------------------------------|----:|
| HiC_scaffold_100 | EVM    | gene | 10808013 | 10842399 |    NA | \-     | NA    | ID=evm.TU.HiC_scaffold_100.43;Name=EVM%20prediction%20HiC_scaffold_100.43   | 100 |
| HiC_scaffold_100 | EVM    | gene | 10852258 | 10852487 |    NA | \-     | NA    | ID=evm.TU.HiC_scaffold_100.44;Name=EVM%20prediction%20HiC_scaffold_100.44   | 100 |
| HiC_scaffold_100 | EVM    | gene | 10898514 | 10907795 |    NA | \-     | NA    | ID=evm.TU.HiC_scaffold_100.335;Name=EVM%20prediction%20HiC_scaffold_100.335 | 100 |
| HiC_scaffold_100 | EVM    | gene | 10936501 | 10940638 |    NA | \+     | NA    | ID=evm.TU.HiC_scaffold_100.336;Name=EVM%20prediction%20HiC_scaffold_100.336 | 100 |
| HiC_scaffold_100 | EVM    | gene | 10929856 | 10931118 |    NA | \+     | NA    | ID=evm.TU.HiC_scaffold_100.337;Name=EVM%20prediction%20HiC_scaffold_100.337 | 100 |
| HiC_scaffold_29  | EVM    | gene | 29027351 | 29042655 |    NA | \+     | NA    | ID=evm.TU.HiC_scaffold_29.19;Name=EVM%20prediction%20HiC_scaffold_29.19     |  29 |
| HiC_scaffold_29  | EVM    | gene | 29066194 | 29068752 |    NA | \-     | NA    | ID=evm.TU.HiC_scaffold_29.20;Name=EVM%20prediction%20HiC_scaffold_29.20     |  29 |
| HiC_scaffold_29  | EVM    | gene | 29052519 | 29052714 |    NA | \+     | NA    | ID=evm.TU.HiC_scaffold_29.21;Name=EVM%20prediction%20HiC_scaffold_29.21     |  29 |
| HiC_scaffold_29  | EVM    | gene | 29125822 | 29127404 |    NA | \-     | NA    | ID=evm.TU.HiC_scaffold_29.22;Name=EVM%20prediction%20HiC_scaffold_29.22     |  29 |
| HiC_scaffold_29  | EVM    | gene | 28824117 | 28858842 |    NA | \-     | NA    | ID=evm.TU.HiC_scaffold_29.75;Name=EVM%20prediction%20HiC_scaffold_29.75     |  29 |
| HiC_scaffold_77  | EVM    | gene | 32263023 | 32263544 |    NA | \-     | NA    | ID=evm.TU.HiC_scaffold_77.75;Name=EVM%20prediction%20HiC_scaffold_77.75     |  77 |
| HiC_scaffold_77  | EVM    | gene | 32714830 | 32714996 |    NA | \+     | NA    | ID=evm.TU.HiC_scaffold_77.89;Name=EVM%20prediction%20HiC_scaffold_77.89     |  77 |
| HiC_scaffold_96  | EVM    | gene | 38143011 | 38144592 |    NA | \+     | NA    | ID=evm.TU.HiC_scaffold_96.135;Name=EVM%20prediction%20HiC_scaffold_96.135   |  96 |
| HiC_scaffold_96  | EVM    | gene | 38178885 | 38179434 |    NA | \+     | NA    | ID=evm.TU.HiC_scaffold_96.170;Name=EVM%20prediction%20HiC_scaffold_96.170   |  96 |
| HiC_scaffold_77  | EVM    | gene | 32263023 | 32263544 |    NA | \-     | NA    | ID=evm.TU.HiC_scaffold_77.75;Name=EVM%20prediction%20HiC_scaffold_77.75     |  77 |
| HiC_scaffold_96  | EVM    | gene | 38143011 | 38144592 |    NA | \+     | NA    | ID=evm.TU.HiC_scaffold_96.135;Name=EVM%20prediction%20HiC_scaffold_96.135   |  96 |
| HiC_scaffold_96  | EVM    | gene | 38178885 | 38179434 |    NA | \+     | NA    | ID=evm.TU.HiC_scaffold_96.170;Name=EVM%20prediction%20HiC_scaffold_96.170   |  96 |
| HiC_scaffold_59  | EVM    | gene | 42693085 | 42722969 |    NA | \-     | NA    | ID=evm.TU.HiC_scaffold_59.452;Name=EVM%20prediction%20HiC_scaffold_59.452   |  59 |
| HiC_scaffold_100 | EVM    | gene | 10808013 | 10842399 |    NA | \-     | NA    | ID=evm.TU.HiC_scaffold_100.43;Name=EVM%20prediction%20HiC_scaffold_100.43   | 100 |
| HiC_scaffold_100 | EVM    | gene | 10852258 | 10852487 |    NA | \-     | NA    | ID=evm.TU.HiC_scaffold_100.44;Name=EVM%20prediction%20HiC_scaffold_100.44   | 100 |
| HiC_scaffold_100 | EVM    | gene | 10898514 | 10907795 |    NA | \-     | NA    | ID=evm.TU.HiC_scaffold_100.335;Name=EVM%20prediction%20HiC_scaffold_100.335 | 100 |
| HiC_scaffold_100 | EVM    | gene | 10936501 | 10940638 |    NA | \+     | NA    | ID=evm.TU.HiC_scaffold_100.336;Name=EVM%20prediction%20HiC_scaffold_100.336 | 100 |
| HiC_scaffold_100 | EVM    | gene | 10929856 | 10931118 |    NA | \+     | NA    | ID=evm.TU.HiC_scaffold_100.337;Name=EVM%20prediction%20HiC_scaffold_100.337 | 100 |
| HiC_scaffold_77  | EVM    | gene | 32263023 | 32263544 |    NA | \-     | NA    | ID=evm.TU.HiC_scaffold_77.75;Name=EVM%20prediction%20HiC_scaffold_77.75     |  77 |
| HiC_scaffold_96  | EVM    | gene | 38143011 | 38144592 |    NA | \+     | NA    | ID=evm.TU.HiC_scaffold_96.135;Name=EVM%20prediction%20HiC_scaffold_96.135   |  96 |
| HiC_scaffold_96  | EVM    | gene | 38178885 | 38179434 |    NA | \+     | NA    | ID=evm.TU.HiC_scaffold_96.170;Name=EVM%20prediction%20HiC_scaffold_96.170   |  96 |
