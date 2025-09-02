# This code is to input the differential gene expression table table, to out put MA or volcano plot
# To generate differential gene expression (DEG) table, see the "RNA-seq analysis_section1,generate NormCounts_DEG tables.R" script 

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

library(tidyverse)
library(dplyr)
library(edgeR)
library(ggplot2)
library(ggrepel)

#setwd("/dir")

####################   READ DEG.TABLE  ####################
KD <- read.table("DEG TABLE.CTP", header=TRUE, sep="\t")

####################  set Categorical Variables  ####################
KD$cat <- 0                                                    # create a new "cat" column with 0 in all the rows
KD[KD$logFC > log2(1.5) & KD$adj.P.Val < 0.05, "cat"] <- 1     # genes with logFC > log2(1.5) & adj.P.Val < 0.05 are changed to to be 1 in the "cat" column
KD[KD$logFC < -log2(1.5) & KD$adj.P.Val < 0.05, "cat"] <- 2    # genes with logFC < -log2(1.5) & adj.P.Val < 0.05 are changed to to be 2 in the "cat" column
KD$cat <- as.factor(KD$cat)                                    # change the number in "cat" column into factor

################ read in gene sets of interest   ###################
Eightcell_Stowers <- read.table("8cell_genes_Stowers.csv",sep=",",header=TRUE) 


################ Generate MA plot   ###################

ggplot(KD, aes(x = AveExpr, y = logFC), label=Gene) + 
  geom_point(alpha = 0.2, color = "grey50") +                                           #alpha for transparency
  geom_point(data = subset(KD, cat == '1'), alpha=0.2, color = "red4") +                     #for colored dots (~cat), "red4" for cat=1 genes. 
  geom_point(data = subset(KD, cat == '2'), alpha=0.2, color = "blue4") +                     #for colored dots (~cat), "blue4" for cat=2 genes.
  geom_text_repel(data = subset(KD, KD$Gene %in% Eightcell_Stowers$Gene), aes(label=Gene),hjust=2,vjust=1, box.padding = 0.5, point.padding = 0.5,segment.size=0.2,size=4, color = 'black', max.overlaps = Inf) +   # Label genes in the "8cell_genes_Stowers.csv" list.
  #geom_text_repel(data = subset(KD, KD$Gene %in% c("Gene1","Gene2")), aes(label=Gene),hjust=0.1,vjust=0.5, box.padding = 0.5, point.padding = 0.5,segment.size=0.5,segment.color = 'grey4') +                      # Label "Gene1" "",Gene2" of your interest
  #geom_line (y=log2(1.5),color = 'grey4') +               # add line at y= log2(1.5)                      
  theme_bw() +                                             # A theme with white background and black gridlines. 
  theme(legend.position = "none") +
  ylab("Log2 fold change (ASO-L1 / SO-L1)") + xlab("log2 Counts per Million") +    #Change name of axis as needed. ylab is the compared log2 fold changes between samples in the input DEG table, here "ASO-L1 / SO-L1" for example. xlab is the "AveExpr" column values in the DEG table, the "log2 Counts per Million" value of the gene.  
  ggtitle("YOUR TITLE") +                                  #Change name of axis as needed
  theme(plot.title = element_text(lineheight=1.2, face="bold",hjust = 0.5, size = 12)) +      #change element text sizes as needed
  coord_cartesian(ylim=c(-3,5))        #change scale limits as needed


################ Generate Volcano plot   ###################

ggplot(KD, aes(x = logFC, y = -log10(adj.P.Val), col = cat)) +
  geom_point(alpha = 0.2, color = "grey50") +                                           #alpha for transparency
  geom_point(data = subset(KD, cat == '1'), alpha=0.2, color = "red4") +                     #for colored dots (~cat), "red4" for cat=1 genes. 
  geom_point(data = subset(KD, cat == '2'), alpha=0.2, color = "blue4") +                     #for colored dots (~cat), "blue4" for cat=2 genes.
  theme_bw() +                                                                   # A theme with white background and black gridlines. 
  geom_text_repel(data = subset(KD, KD$Gene %in% Eightcell_Stowers$Gene), aes(label=Gene),hjust=2,vjust=1, box.padding = 0.5, point.padding = 0.5,segment.size=0.2,size=4, color = 'black', max.overlaps = Inf) +   # Label genes in the "8cell_genes_Stowers.csv" list.
  #geom_text_repel(data = subset(KD, KD$Gene %in% c("Gene1","Gene2")), aes(label=Gene),hjust=0.1,vjust=0.5, box.padding = 0.5, point.padding = 0.5,segment.size=0.5,segment.color = 'grey4') +                      # Label "Gene1" "",Gene2" of your interest
  ggtitle("YOUR TITLE") +
  theme(plot.title = element_text(lineheight=1.2, face="bold",hjust = 0.5, size = 12)) +
  ylab("-log10 FDR") + xlab("log2 fold-change (ASO-L1 / SO-L1)")     #Change name of axis as needed. xlab is the adj.P.Val of each gene. ylab is the compared log2 fold changes between samples in the input DEG table, here "ASO-L1 / SO-L1" for example


################ Correlation plot   ###################
# compare gene expression change by EZH2i and L1KD, for correlation efficiency

install.packages("ggpubr")
library("ggpubr")

EZH2i <- read.table("RSeT+DT_EZH2i,DEG TABLE.CTP", header=TRUE, sep="\t")
L1KD <- read.table("RSeT+DT_L1KD,DEG TABLE.CTP", header=TRUE, sep="\t")

EZH2i <- EZH2i[,c('Gene', 'logFC')]
L1KD <- L1KD[,c('Gene', 'logFC')]
Merge <- merge(EZH2i,L1KD,by='Gene', all = TRUE)
Merge <- na.omit(Merge)

write.table(Merge, file = "L1KD_EZH2i-logFC_Merged.txt", row.names = T, sep = "\t", quote = F)

ggscatter(Merge, x = "logFC.x", y = "logFC.y", 
          add = "reg.line", conf.int = TRUE, 
          add.params = list(color = "red2", fill = "lightgray"),
          cor.coeff.args = list(method = "spearman", label.y = 4, label.sep = "\n"),
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "logFC EZH2i", ylab = "logFC L1KD")
cor.test(Merge$logFC.x, Merge$logFC.y, alternative = c("two.sided"), method = c("spearman"))

