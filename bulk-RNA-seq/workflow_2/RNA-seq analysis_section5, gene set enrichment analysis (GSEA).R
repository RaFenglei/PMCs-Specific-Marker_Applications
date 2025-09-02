# This code is to input the t-ranked gene list from the differential gene expression (DEG) table.
# To generate DEG table, see the "RNA-seq analysis_section1,generate NormCounts_DEG tables.R" script, 
# The DEG tables above is available from GEO (GSE232939) for each experimental set below.

# The "Gene" and "t" columns (ranked by descending t values) from the DEG table are saved to an new .RNK table without header for input, just by changing the file type.See the "RSeT+DT hESC_L1KD,ASO-L1intORF_sO-L1intORF,-t.RNK" file in the repository for example.
# Gene sets used in this study are listed in the Supplementary Table 2. Save the gene sets/names in a .gmt file without header, see the "8cell_genes_Stowers.gmt" file in the repository for example. 

# This code applys to the experiments below: 
  ## 1.RNA-seq from ASO-L1 KD in Primed hESCs
  ## 2.RNA-seq from ASO-L1 KD in RSeT hESCs
  ## 3.RNA-seq from ASO-L1 KD in RSeT+DT hESCs
  ## 4.RNA-seq from ASO-L1 KD in PXGL hESCs
  ## 5.RNA-seq from EZH2i in RSeT+DT hESCs
  ## 6.RNA-seq from ASO-L1 and TPRX1-siRNA transfected co-KD in RSeT+DT hESCs.
  ## 7.RNA-seq from Collinson et al 2016, Cell Reports_GSE76626_RNA-Seq ("GSE76626_ESC-EZH2KO_raw_counts_rmDup.txt" and "GSE76626_EZH2null_WT,cutoff=0.2,Norm_counts.txt" are analyzed following RNA-seq analysis methods in this paper and are included in the repository.) 


library(fgsea)
library(tidyverse)

#setwd("/dir")
####################   READ RNK.TABLE  ####################
rnk.file <- "THE RNK FILE.RNK"
ranks <- read.table(rnk.file,header=FALSE, colClasses = c("character", "numeric"))
colnames(ranks) <- c("ID", "logFC")
ranks <- setNames(ranks$logFC, ranks$ID)
str(ranks)

################ read in gene sets of interest and save as "pathway"   ###################
List1 <- read.table("Maternal.gmt", header = F, stringsAsFactors = F)
List1 <- as.vector(List1$V1)
List2 <- read.delim("Zygotic.gmt", header = F, stringsAsFactors = F)
List2 <- as.vector(List2$V1)
List3 <- read.table("8cell_genes.gmt", header = F, stringsAsFactors = F)
List3 <- as.vector(List3$V1)
## TO CONTINUE the lists (gene sets)....

pathways <- list("Maternal" = List1, "Zygotic" = List2,"8cell" = List3)  # Name each gene list, and list them together as pathway
str(head(pathways))

########################### run the fgsea function ###############################
fgseaRes <- fgsea(pathways, ranks, minSize=15, maxSize=3000, nperm=6000)
head(fgseaRes)          
  ## this retuns a table, in which the NES and padj indicates the enrichment score and adjusted-P values of each gene sets for the input file (RNK. file from the DEG table)
  ## plot the NES of each gene sets in xlsx. as bar plot. 

##########################  GSEA plot #################################
GSEA.down <- plotEnrichment(pathways[["Maternal"]],       # Gene name of list1, or list2, or...
                            ranks,
                            ticksSize = 0.01) + #width of ticks
  ylab("Enrichment score") + xlab("Rank") +
  ggtitle("ASO-L1intORF/sO,Maternal") + 
  theme(plot.title = element_text(lineheight=1.2, face="bold",hjust = 0.5, size = 12)) 

GSEA.down                # plot shows in RStudio.

tiff("TIFF FILE TITLE.tiff", units="in", width=3.5, height=2, res=300)
GSEA.down
dev.off()                # plot.tiff file saved to the folder.

