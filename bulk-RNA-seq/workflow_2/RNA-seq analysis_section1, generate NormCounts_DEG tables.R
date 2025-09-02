# This code is to input the raw counts table, to out put
  ## 1. the generate differential gene expression (DEG) table (no batch correction).  
  ## 2. the normalized counts table (ComBat batch corrected). 

# The input and out put tables are available from GEO (GSE232939) for experiments 1. to 6. set below.

# This code applys to the experiments below: 
  ## 1.RNA-seq from ASO-L1 KD in Primed hESCs
  ## 2.RNA-seq from ASO-L1 KD in RSeT hESCs
  ## 3.RNA-seq from ASO-L1 KD in RSeT+DT hESCs
  ## 4.RNA-seq from ASO-L1 KD in PXGL hESCs
  ## 5.RNA-seq from EZH2i in RSeT+DT hESCs
  ## 6.RNA-seq from ASO-L1 and TPRX1-siRNA transfected co-KD in RSeT+DT hESCs.
# The raw counts are generated following the RNA-seq analysis pipeline created by M. Percharde (https://github.com/mpercharde/RNAseq), also see methods.
# The raw counts generated from feature counts (the pipeline above) are formatted (e.g. remove duplicates) as needed for input read. 

# This code also applys to The experiments:
## 5.Hendrickson, et al. (PGH)_EmbryoStages_data analysis for Genes expression (Fig.4f). 
### for 5., The input raw counts table "PGHendrickson_EmbryoStage_Genes,raw_counts.txt" and output "PGHendrickson_EmbryoStage_Genes,N_counts.txt" are included in this repository.


library(tidyverse)
library(dplyr)
library(edgeR)
library(ggplot2)
library(sva)

#setwd("/dir")

####################   READ IN RAW COUNTS  ####################
counts <- read.delim("RAW COUNTS TABLE.txt", header=TRUE, sep="\t")            # The raw counts table is available from GEO (.....)
head(counts)
sapply(counts, class)
snames <- colnames(counts)                      
snames
group <- factor(substr(snames,1,6))       # sub-string sample names from 1 to 6 (or 3, 5, as needed), so samples from different batches are grouped by name.  
group
batch <- c(1,1,1,2,2,2,3,3,3,4,4,4)       # set sample batches

d0 <- DGEList(counts)
d0 <- calcNormFactors(d0)
d0

####################   Filtering low count genes    ####################
cutoff <- 0.2                  # or cutoff <- 1, or cutoff <- 0.5, lower value cutoff more genes, but also more noisy background in normalization.  
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d0)  # number of counts input genes left
dim(d)   # number of genes left after cutoff filter


#########   Building the model for limma-voom _ Read depth normaliztion   #########
mm <- model.matrix(~0 + group + batch)         
y <- voom(d, mm, plot = T)
fit <- lmFit(y, mm)
head(coef(fit))
dim(coef(fit))


####################   Differential Gene Expression analysis  ####################

# Differential analysis for ASO-L1 to SO-L1 cells. 
# Or do other comparisons as need by changing the groupNAME. The groupNAME should be the same as the group names from above. 
contr <- makeContrasts(groupASO-L1 - groupSO-L1., levels = colnames(coef(fit)))           # ASO-L1 and SO-L1. is compared. 
contr
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC)>log2(1.5)))    # to inquire how many genes are: adj.P.Val < 0.05 & abs(logFC)>log2(1.5)
top.table$Gene <- rownames(top.table)
top.table <- top.table[,c("Gene", names(top.table)[1:6])]
write.table(top.table, file = "hESC_L1KD,ASO-L1KD to ASO-Ctr_DEG,cutoff=0.2.txt", row.names = F, sep = "\t", quote = F)  # exported the differential gene expression table 

# The exported DEG_table applies for MA, volcano and violin plot, and the "t" column in the table is used in GSEA analysis. 

View(y$E)
write.table(y$E, file = "hESC_L1KD,lo2Nom_counts,cutoff=0.2.txt", row.names = T, sep = "\t", quote = F)  # export the normalized counts table (no batch correction).  


####################   ComBat batch correction function  ####################
##https://github.com/zhangyuqing/ComBat-seq
#The normalized counts table after ComBat batch correction applies for MDS plot and heatmap plot.

adjusted <- ComBat_seq(as.matrix(counts), batch=batch, group=group)   
d1 <- DGEList(adjusted)                   # data matrix with ComBat batch adjusted counts table
d1 <- calcNormFactors(d1)
d1

cutoff <- 0.2                                 # cutoff low count genes
drop <- which(apply(cpm(d1), 1, max) < cutoff)
d <- d1[-drop,] 
dim(d1)
dim(d) # number of genes left

####################  generate MDS plot  ####################

plotMDS(d, col = as.numeric(group))    

mm <- model.matrix(~0 + group + batch)        
y <- voom(d, mm, plot = T)              # Normalization counts on ComBat batch adjusted table

################   Generate log2 normalized Counts table  ####################
write.table(y$E, file = "hESC_L1KD,lo2Nom_counts,cutoff=0.2_BatchCo.txt", row.names = T, sep = "\t", quote = F)

#The batch-corrected normalized counts table is utilized for heatmap plot.





############################################################################## 
#############     Generate RefGene-deep normalized and      ################## 
#############   batch-corrected Counts table for repeats    ##################
##############################################################################

# This code applys to the experiments: 
## 1. ASO-L1 KD in RSeT+DT hESCs, for plotting of repeats expression (Figure S2C). 
### The repeats raw counts table "RSeT+DT,hL1KD-raw_counts,s1_16_repeat.txt" is included in this repository. The RefGene_raw_counts for this experiment is available from GEO (GSE232939)
### The output "RSeT+DT hESC_L1KD,lo2Nom_counts_repeats,Ref-N_BathCo.txt" is included in this repository, and utilized for the heatmap plot in Figure S2C. 

## 2. Hendrickson, et al. (PGH)_EmbryoStages_data analysis for TEs/LINE1s expression (Figure S1A).  
### The repeats raw counts table "PGHendrickson_EmbryoStage_Repeats, raw_counts.txt" and "PGHendrickson_EmbryoStage_Genes,raw_counts.txt" is included in this repository 
### The output "PGHendrickson_EmbryoStage_Repeats,RefN_counts.txt" is included in this repository, and utilized for the heatmap plot in Figure S1A. 


# Input raw counts table of the repeats 
repeats <- read.delim("repeats RAW COUNTS TABLE.txt", header=TRUE, sep="\t")
repeats <-na.omit(repeats)     # remove NA rows
head(repeats)
sapply(repeats, class)

snames <- colnames(repeats)                      
snames
group <- factor(substr(snames,1,6))       # sub-string sample names from 1 to 6 (or 3, 5, as needed), so samples from different batches are grouped by name.  
group
batch <- c(1,1,1,2,2,2,3,3,3,4,4,4)       # set sample batches

# ComBat batch correcte the repeats table 
adjusted_repeats <- ComBat_seq(as.matrix(repeats), batch=batch, group=group)    #https://github.com/zhangyuqing/ComBat-seq

# generate data matrix with the ComBat batch-corrected repeats table, and cutoff 
d0 <- DGEList(adjusted_repeats)
d0 <- calcNormFactors(d0)
d0
  
cutoff <- 1                                       # cutoff low count repeats
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
d <- calcNormFactors(d)
d
dim(d0)
dim(d) # number of genes left


# Imput raw counts table of the ref genes, ComBat batch corrected, and then calculate the normalization factor. 
Refcounts  <- read.delim("RefGenes RAW COUNTS TABLE.txt", header=TRUE, sep="\t")            # The raw counts table is available from GEO (.....). This is count of the ref genes, used for calculate normalization factors.
head(Refcounts)

adjusted_Refcounts <- ComBat_seq(as.matrix(Refcounts), batch=batch, group=group)    #https://github.com/zhangyuqing/ComBat-seq

N <-colSums(adjusted_repeats)                          # count sum of reads from each colum -> each smaple       
nf <- calcNormFactors(adjusted_Refcounts,lib.size = N)      # count normalization factor based on the ref sample reads

plotMDS(d, col = as.numeric(group))         #MDS the ComBat batch-corrected repeats table matrix 

y <- voom(d, mm,lib.size = N * nf, plot = T)  # calculate the normalization factor on ComBat batch corrected Ref-gene matrix 
mm <- model.matrix(~0 + group + batch)    

View(y$E)
write.table(y$E, file = "RSeT+DT hESC_L1KD,lo2Nom_counts_repeats,Ref-N_BathCo.txt", row.names = T, sep = "\t", quote = F)

