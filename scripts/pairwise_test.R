library(tidyverse)

paired_samples <- read_tsv("paired_by_distance.txt")
wilcox.test(log2(paired_samples$mtrD_azi_mic), log2(paired_samples$wt_azi_mic), paired = TRUE)
