# code to analyze RNA-seq from CEEHRC
# by Kirti Mittal, PhD
# The input and output tables are available from GEO (GSE232939).

library(DESeq2)
library(tidyverse)
#setwd("current_dir/")

###### read in raw counts ######
counts <- read_delim("count_matrix_repeats.txt", delim = "\t",
col_names = c("GeneID", "Naive_HESC_RSeT_DT_1", "Naive_HESC_RSeT_DT_2", "Naive_HESC_RSeT_1", "Naive_HESC_RSeT_2", "Primed_HESC_1", "Primed_HESC_2"), comment = "#") %>%
dplyr::slice(-1) %>%
mutate_at(vars(Naive_HESC_RSeT_DT_1:Primed_HESC_2), as.numeric) %>%
filter(!grepl("-Mar", GeneID))
head(counts)
counts <- counts[,c(1,2,3,4,5,6,7)]
counts_df <- counts %>% select(-GeneID) %>% as.data.frame()
rownames(counts_df) <- counts$GeneID

## Normalization using Deseq##
samples <- c("RSeT_DT1", "RSeT_DT2", "RSeT_1", "RSeT_2", "Primed_1", "Primed_2")
condition <- factor(c(1,1,2,2,3,3),labels = c("RSeT_DT" , "RSeT" , "Primed" ))
info <- as.data.frame(cbind(samples, condition))
matrix <- DESeqDataSetFromMatrix(countData = counts_df, colData = info, design = ~condition)
matrix <- matrix[rowSums(counts(matrix)) >= 10,]
DEmatrix <- DESeq(matrix)
table_counts_normalised <- cpm(DEmatrix)
