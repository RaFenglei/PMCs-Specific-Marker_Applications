library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggrepel)
library("ggpubr")
library("cli")
library(edgeR)


setwd("/dir")

####################   READ IN RAW COUNTS and merge in one table  ####################
 ## The raw counts tables are available from GEO (GSE232939).

Primed_L1KD <- read.table("Primed_RAW COUNTS TABLE.txt", header=TRUE, sep="\t")
RSeT_L1KD <- read.table("RSeT_RAW COUNTS TABLE.txt", header=TRUE, sep="\t")
RSeTDT_L1KD <- read.table("RSeT+DT_RAW COUNTS TABLE.txt", header=TRUE, sep="\t")
PXGL_L1KD <- read.table("PXGL_RAW COUNTS TABLE.txt", header=TRUE, sep="\t")

Merge_Naive <- merge(merge(RSeT_L1KD ,PXGL_L1KD,by.x ='Gene', all.x = TRUE),RSeTDT_L1KD,by.x ='Gene', all.x = TRUE)
Merge_All <- merge(Merge_Naive,Primed_L1KD,by.x ='Gene', all.x = TRUE)
Merge_All <- na.omit(Merge_All)

write.table(Merge_All, file = "ALL-hESC_L1KD,Combin_raw_counsts.txt", row.names = T, sep = "\t", quote = F)

####################  generate MDS plot with edgeR  ####################
raw.data <- read.delim("ALL-hESC_L1KD,Combin_raw_counsts.txt", header=TRUE, sep="\t")
counts<-raw.data[,-1]
row.names(counts)<-raw.data[,1]


head(counts)
sapply(counts, class)

snames <- colnames(counts)
snames

group <- factor(substr(snames,1,9))
group

batch <- c(1,2,3,4,1,2,3,4,1,2,3,4,5,6,7,5,6,7,8,9,10,11,8,9,10,11,8,9,10,11)

d0 <- DGEList(counts)
d0 <- calcNormFactors(d0)
d0

cutoff <- 0.5   #21097/ 53713

drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d0)
dim(d) # number of genes left

plotMDS(d0, col = as.numeric(group))
plotMDS(d, col = as.numeric(group))
plotMDS(d, col = as.numeric(group),pch=1, cex = 1, dim.plot = c(1,2),plot = TRUE)

