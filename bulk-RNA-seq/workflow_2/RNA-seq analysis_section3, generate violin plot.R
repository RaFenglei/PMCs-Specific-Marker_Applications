# This code is to input the differential gene expression table, to out put violin plot
# To generate differential gene expression (DEG) table, see the "RNA-seq analysis_section1,generate NormCounts_DEG tables.R" script 


# The input tables above is available from GEO (GSE232939) for each experimental set below.
# The DEG table (.txt) is coverted to (.CTP) for input, just by changing the file type.
# Gene sets used in this study are listed in the Supplementary Table 2. Save the gene sets/names in a .csv file with header "Gene"), see the "8cell_genes_Stowers.csv" file in the repository for example. 



library(tidyverse)
library(dplyr)
library(edgeR)
library(ggplot2)
library(ggrepel)

#################################################################
########## This code applys to the experiments below: ###########
########## 1.RNA-seq from ASO-L1 KD in RSeT+DT hESCs  ###########
################## For Extended Data Fig.5f    ################## 
#################################################################

#setwd("/dir")
####################   READ DEG.TABLE and gene_sets.TABLE  ####################
RSeTwtDT_ASOint <- read.table("ASO-L1intORF_sO-L1intORF,cutoff=0.2.CTP", header=TRUE, sep="\t")       #"ASO-L1intORF_sO-L1intORF,cutoff=0.2.txt" is available from GEO (.....)

Morulae8cell <- read.table("8cell_Morulae_1687.csv",sep=",",header=TRUE)          #Yan paper, 1687 gene
Eightcell_Yan <- read.table("8cell_gene.csv",sep=",",header=TRUE)                 #Yan paper, 533 gene
Eightcell_Yu <- read.table("8cell_genes_Yu.csv",sep=",",header=TRUE)              #Yan paper, 743 gene
Eightcell_Stowers <- read.table("8cell_genes_Stowers.csv",sep=",",header=TRUE)         #Stowers paper, 45 gene

####################   Subset 8C genes from the gene_sets.TABLE expressed in the DEG.TABLE   ####################

RSeTwtDT_ASOint_8cell_Stowers<- RSeTwtDT_ASOint[RSeTwtDT_ASOint$Gene %in% Eightcell_Stowers$gene, ]      #19
RSeTwtDT_ASOint_8cell_Yan<- RSeTwtDT_ASOint[RSeTwtDT_ASOint$Gene %in% Eightcell_Yan$gene, ]       #443
RSeTwtDT_ASOint_8cell_Yu<- RSeTwtDT_ASOint[RSeTwtDT_ASOint$Gene %in% Eightcell_Yu$gene, ]             #341
RSeTwtDT_ASOint_Morulae8cell<- RSeTwtDT_ASOint[RSeTwtDT_ASOint$Gene %in% Morulae8cell$gene, ]         #1343

RSeTwtDT_random <- sample_n(RSeTwtDT_ASOint, 3000, replace = FALSE)    #take a RANDOM 3000 SUBSET of genes


####################    Combine the gene sets expression tables to one, set Categorical Variables   ####################
RSeTwtDT_ASOint <- bind_rows("RSeTwtDT_random3000" = RSeTwtDT_random,"RSeTwtDT_ASOint_Morulae8cell" = RSeTwtDT_ASOint_Morulae8cell, "RSeTwtDT_ASOint_8cell_Yan" = RSeTwtDT_ASOint_8cell_Yan,"RSeTwtDT_ASOint_8cell_Yu" = RSeTwtDT_ASOint_8cell_Yu, "RSeTwtDT_ASOint_8cell_Reik" = RSeTwtDT_ASOint_8cell_Reik, .id = "change")   # the combined table called "RSeTwtDT_ASOint" 
RSeTwtDT_ASOint$change <- factor(RSeTwtDT_ASOint$change, levels = c("RSeTwtDT_random3000", "RSeTwtDT_ASOint_Morulae8cell", "RSeTwtDT_ASOint_8cell_Yan","RSeTwtDT_ASOint_8cell_Yu", "RSeTwtDT_ASOint_8cell_Reik"))    # set Categorical Variables in the "change" column

####################    plot each categorical variables as violin plot for comparison  ####################
ggplot(RSeTwtDT_ASOint, aes(change, logFC, fill = change)) +
  geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim=c(-2,2)) +
  theme_bw() +                                             # A theme with white background and black gridlines. 
  geom_hline(yintercept=c(0)) +
  geom_hline(yintercept=c(-log2(1.5),log2(1.5)), linetype = "dashed") +
  scale_fill_manual(values = c("grey65", "indianred2", "red2","green2","blue2")) +          # set the color for each category  
  theme(axis.text.x = element_text(color = "black",size= 12), axis.text.y = element_text(color = "black", size = 12)) +
  ggtitle("RSeT+DT H9, L1KD") + 
  ylab("Log2 FC (ASO-L1/SO-L1)") + 
  theme(plot.title = element_text(lineheight=1.2, face="bold",hjust = 0.5, size = 12)) +
  theme(axis.text.x = element_text(color = "black",size= 5), axis.text.y = element_text(color = "black", size = 10)) +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width =10)) +
  theme(panel.background = element_rect(fill = "transparent", color = "black", size = 0.5), plot.background = element_rect(fill = "transparent"), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  coord_cartesian(ylim=c(-1,3))        #change scale limits as needed

####################    Statistic: wilcox.test for expression changes between the gene sets   ####################
wilcox.test(RSeTwtDT_ASOint_Morulae8cell$logFC,RSeTwtDT_random $logFC)   
wilcox.test(RSeTwtDT_ASOint_8cell_Yan$logFC,RSeTwtDT_random$logFC) 
wilcox.test(RSeTwtDT_ASOint_8cell_Yu$logFC,RSeTwtDT_random$logFC) 
wilcox.test(RSeTwtDT_ASOint_8cell_Stowers$logFC,RSeTwtDT_random$logFC) 








#################################################################################################
################################### This code applys to the experiments below: ##################
########## 1.RNA-seq from ASO-L1 and TPRX1-siRNA transfected co-KD in RSeT+DT hESCs.  ###########
############################################    For Fig.2j    ###################################
#################################################################################################

####################   READ DEG.TABLE and gene_sets.TABLE  ####################

L1KD_NTsi<- read.table("L1KD.NTsiRNA_CtrA.NTsi_cutoff=0.2.CTP", header=TRUE, sep="\t")       # Table is available from GEO (.....)
L1KD_TPRX1si<- read.table("L1KD.TPRX1si_CtrA.NTsi_cutoff=0.2.CTP", header=TRUE, sep="\t")    # Table is available from GEO (.....)

Eightcell_Stowers <- read.table("8cell_genes_Stowers.csv",sep=",",header=TRUE)         #Stowers paper, 45 gene

######## subset for protein coding gene sets #######
PcG <- read.table("human_NM_hg38_Genes.csv",sep=",",header=TRUE)   # 19378 protein coding genes,hg38_NM. ## The "human_NM_hg38_Genes.csv" file is in the repository  
L1KD_NTsi<- L1KD_NTsi[L1KD_NTsi$Gene %in% PcG$gene, ]              #15264
L1KD_TPRX1si<- L1KD_TPRX1si[L1KD_TPRX1si$Gene %in% PcG$gene, ]     #15264

####################   Subset 8cell_Stowers genes expressed in each DEG.TABLE   ####################
L1KD_NTsi_8cell_Stowers<- L1KD_NTsi[L1KD_NTsi$Gene %in% Eightcell_Stowers$gene, ]      #20
L1KD_TPRX1si_8cell_Stowers<- L1KD_TPRX1si[L1KD_TPRX1si$Gene %in% Eightcell_Stowers$gene, ]      #20

L1KD_NTsi_random <- sample_n(L1KD_NTsi, 45, replace = FALSE)     #take a RANDOM 45 SUBSET of genes 

####################    Combine 8cell_Stowers genes in all the DEG.Table to one table, set Categorical Variables   ####################
Eightcell_Stowers <-bind_rows("L1KD_NTsi_random"=L1KD_NTsi_random,"L1KD_NTsi_8cell_Stowers"= L1KD_NTsi_8cell_Stowers,"L1KD_TPRX1si_8cell_Stowers"= L1KD_TPRX1si_8cell_Stowers, .id = "change")    # the combined table called "Eightcell_Stowers"
Eightcell_Stowers$change <- factor(Eightcell_Stowers$change, levels = c("L1KD_NTsi_random", "L1KD_NTsi_8cell_Stowers", "L1KD_TPRX1si_8cell_Stowers"))      # set Categorical Variables in the "change" column

####################    plot each categorical variables as violin plot for comparison  ####################
ggplot(Eightcell_Stowers, aes(change, logFC, fill = change)) +
  geom_violin(outlier.shape = NA) +
  coord_cartesian(ylim=c(-2,2)) +
  theme_bw() +                                              # A theme with white background and black gridlines. 
  geom_hline(yintercept=c(0)) +
  geom_hline(yintercept=c(-log2(1.5),log2(1.5)), linetype = "dashed") +
  scale_fill_manual(values = c("grey80", "#595A5A", "#EC7D30", "#4971B8","#70AD45")) +                                # set the color for each category  
  theme(axis.text.x = element_text(color = "black",size= 12), axis.text.y = element_text(color = "black", size = 12)) +
  ggtitle("RSeT+DT H9,8cell_Stowers,cutoff=0.2") + 
  ylab("Log2 FC (~ /Ctr_NTsi)") + 
  theme(plot.title = element_text(lineheight=1.2, face="bold",hjust = 0.5, size = 12)) +
  theme(axis.text.x = element_text(color = "black",size= 5), axis.text.y = element_text(color = "black", size = 10)) +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width =10)) +
  theme(panel.background = element_rect(fill = "transparent", color = "black", size = 0.5), plot.background = element_rect(fill = "transparent"), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  coord_cartesian(ylim=c(-2,5))        #change scale limits as needed

####################    Statistic: wilcox.test for expression changes between the gene sets   ####################
wilcox.test(L1KD_NTsi_8cell_Stowersk$logFC,L1KD_TPRX1si_8cell_Stowers$logFC) 
wilcox.test(L1KD_NTsi_8cell_Stowers$logFC,L1KD_NTsi_random$logFC)
