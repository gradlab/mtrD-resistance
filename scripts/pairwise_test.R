library(tidyverse)

paired_samples <- read_tsv("paired_by_distance.txt")
wilcox.test(paired_samples$mtrD_azi_mic, paired_samples$wt_azi_mic, paired = TRUE)
