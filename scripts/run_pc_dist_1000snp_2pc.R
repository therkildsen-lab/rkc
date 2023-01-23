library(tidyverse)
library(lostruct)
## Read the input
pca_summary <- read_tsv('/fs/cbsubscb16/storage/rkc/angsd/local_pca/window1000/pca_summary_1000snp_2pc.tsv', col_names = F) 
## Run pc_dist with pca_summary
pca_summary <- as.matrix(pca_summary)
attr(pca_summary, 'npc') <- 2
dist <- pc_dist(pca_summary)
write_tsv(as.data.frame(dist), '/fs/cbsubscb16/storage/rkc/angsd/local_pca/window1000/window_dist_1000snp_2pc.tsv', col_names = F)
