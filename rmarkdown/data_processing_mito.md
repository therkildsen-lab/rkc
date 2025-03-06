Data Processing mito
================

## load relevent libraries and data

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
sample_table <- read_tsv("/fs/cbsubscb16/storage/rkc/sample_lists/sample_table.tsv")
```

    ## Rows: 183 Columns: 9
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (6): population, Loc, GeneralLoc, k3_inferred_pop, k4_inferred_pop, k5_i...
    ## dbl (3): ABLG, StartLatDD, StartLonDD
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

#### Sample lists and tables

``` r
bam_list <- read.delim("/fs/cbsubscb16/storage/rkc/sample_lists/bam_list_mito.txt", header = F) %>% mutate(ABLG = as.double(substr(V1,46,49))) %>% inner_join(sample_table, by = "ABLG") %>% dplyr::select(V1) %>% as.matrix()

write_lines(bam_list, "/fs/cbsubscb16/storage/rkc/sample_lists/bam_list_mito.txt")

sample_table_mito <- read.delim("/fs/cbsubscb16/storage/rkc/sample_lists/bam_list_mito.txt", header = F) %>% mutate(ABLG = as.double(substr(V1,46,49))) %>% inner_join(sample_table, by = "ABLG") %>% select(-1) 

write_tsv(sample_table_mito, "/fs/cbsubscb16/storage/rkc/sample_lists/sample_table_mito.tsv")
```

#### No need to deduplicate and overlap clip. Laura did this already. Data is ready for angsd
