# This code is to input the normalized counts table (with or without ComBat batch correction), to out put heatmap plot
# To generate normalized counts table (with or without ComBat batch correction), see the "RNA-seq analysis_section1,generate NormCounts_DEG tables.R" script 

# The input tables above is available from GEO (GSE232939) for each experimental set below.
# The DEG table (.txt) is coverted to (.CTP) for input, just by changing the file type.
# Gene sets used in this study are listed in the Supplementary Table 2. Save the gene sets/names in a .csv file with header "Gene"), see the "8cell_genes_Stowers.csv" file in the repository for example. 


# This code applys to the experiments below: 
  ## 1.RNA-seq from ASO-L1 KD in Primed hESCs
  ## 2.RNA-seq from ASO-L1 KD in RSeT hESCs
  ## 3.RNA-seq from ASO-L1 KD in RSeT+DT hESCs
  ## 4.RNA-seq from ASO-L1 KD in PXGL hESCs
  ## 5.RNA-seq from EZH2i in RSeT+DT hESCs
  ## 6.RNA-seq from ASO-L1 and TPRX1-siRNA transfected co-KD in RSeT+DT hESCs.


install.packages("gplots", dependencies = TRUE)
install.packages("RColorBrewer", dependencies = TRUE)
library(gplots)
library(RColorBrewer)
library(tidyverse)
library(dplyr)
library(edgeR)
library(ggplot2)
library(ggrepel)

library(metaMA)

#setwd("/dir")

####################   READ the normalized counts.TABLE and gene_sets.TABLE  ####################
topT <- read.table("THE NORMALIZED COUNTS TABLE.txt", header=T, sep="\t",comment.char="#")  

Morulae8cell <- read.table("8cell_Morulae_1687.csv",sep=",",header=TRUE)          #Yan paper, 1687 gene
Eightcell_Yan <- read.table("8cell_gene.csv",sep=",",header=TRUE)                 #Yan paper, 533 gene
Eightcell_Yu <- read.table("8cell_genes_Yu.csv",sep=",",header=TRUE)              #Yan paper, 743 gene
Eightcell_Stowers <- read.table("8cell_genes_Stowers.csv",sep=",",header=TRUE)         #Stowers paper, 45 gene

####################   Subset 8C genes from the gene_sets.TABLE expressed in the normalized counts.TABLE   ####################
topT<- topT[topT$GeneID %in% Eightcell_Stowers$gene, ]      #21
topT<- topT[topT$GeneID %in% Eightcell_Yan$Gene, ]       #449
topT<- topT[topT$GeneID %in% Eightcell_Yu$Gene, ]          #338
topT<- topT[topT$GeneID %in% Morulae8cell$Gene, ]         #1370

####################   Subset genes   #####################
#subtopT <- subset(topT, rowSums(topT[ ,2:16])>10)   # remove sum-each-row <10
#topT <- topT[order(abs(rowVars(topT[,2:13])), decreasing=TRUE)[seq_len(5000)],]  # subset top5000 most variable gene in the table


rnames <- topT[,1]
mat_topT <- data.matrix(topT[,2:ncol(topT)])

test<- mat_topT[1:50, ]
head(mat_topT,5)
rownames(mat_topT) <- rnames

# creates a own color palette from red to blue
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)

# (optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(seq(-1,0,length=100),               # for blue
               seq(0.01,0.8,length=100),           # for white
               seq(0.81,1,length=100))             # for red

## Row clustering (adjust here distance/linkage methods to what you need!)
#hr1 <- hclust(as.dist(1-cor(t(mat_topT), method="pearson")), method="complete")
hr2 <- hclust(as.dist(1-cor(t(mat_topT), method="spearman")), method="complete")
## Column clustering (adjust here distance/linkage methods to what you need!)
hc <- hclust(as.dist(1-cor(mat_topT, method="spearman")), method="complete")

heatmap.2(mat_topT,
          #cellnote = mat_data,         # same data set for cell labels
          main = "TITLE",               # heat map title
          #notecol="black",             # change font color of cell labels to black
          notecol="NA",
          density.info="none",      # turns off density plot inside color legend
          trace="none",             # turns off trace lines inside the heat map
          margins =c(12,9),         # widens margins around plot
          col=my_palette,           # use on color palette defined earlier
          scale="row",
          #breaks=col_breaks,       # enable color transition at specified limits
          dendrogram="col",         # only draw a col dendrogram
          Rowv = as.dendrogram(hr2),        # apply default clustering method
          Colv=as.dendrogram(hc))           # turn off column clustering




#####################  apply to plot genes/TEs expression during embryonic development  ############################3 

# the heatmap.R script above is also applied to plot genes or TEs of interest for their expression comparison in developing embryos. 
# The data resource from Hendrickson et al., 2017; Yan et al., 2013 and Yandım & Karakülah, 2019 are available from their Supplemental Table. Format the date table to the "THE NORMALIZED COUNTS TABLE.txt" for input (Detials are described in methods and/or each figure legend)

# The formated "THE NORMALIZED COUNTS TABLE.txt" table for input are available in the repository:
## 1.scRNAseq_Yan, human embryo stages.txt
## 2.scRNA_Yandim_Yan, EmbryoStages_repeats_avg.CPM.txt
## 3.PGHendrickson_EmbryoStage_Genes,N_counts.txt (To generate this table, see "RNA-seq analysis_section1, generate NormCounts_DEG tables.R")
## 4.PGHendrickson_EmbryoStage_Repeats,RefN_counts.txt (To generate this table, see "RNA-seq analysis_section1, generate NormCounts_DEG tables.R")

# The "genes/TEs.csv" to subset from the normalized counts table to plot in this study is shown in each figure. 




