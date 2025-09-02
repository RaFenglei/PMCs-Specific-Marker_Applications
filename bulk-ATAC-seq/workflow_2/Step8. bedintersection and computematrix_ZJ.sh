#SBATCH --job-name=BED_computeMatrixPlot
#SBATCH -N 1
#SBATCH -n 6


#1.	bed intersction to generate BED files of Ref genes nearby H3K27me3-ChIP peaks, and H3K27me3-ChIP peaks in LADsNADs compartments. 
#2.	in R, subset Ref genes BED in each chromosome.
#3.     average BigWig of replicates for ploting.
#4.     computeMAtrix and Plot



#1. bed intersction to generate BED files of Ref genes nearby H3K27me3-ChIP peaks, and H3K27me3-ChIP peaks in LADsNADs compartments. 
# ctr-ChIP_merge_peaks_blacklisted.bed is used for computeMatrix below.
# ctr-ChIP_merge_peaks_blacklisted.bed is generated from step3.(macs3 callpeak on ctr-ChIP_merge.bam) and step4 (rm-blacklisting).

bedtools window -a ctr-ChIP_merge_peaks_blacklisted.bed -b hg38_genes_refseq.bed -w 3000 > K27me3peaks_hg38ref.bed      # bedtools window for hg38 ref genes within 3kb of K27me3 peaks

bedtools intersect -a ctr-ChIP_merge_peaks_blacklisted.bed -b 2312simple_naive_LAD.bed > naive_LAD_K27me3peaks.bed   
bedtools intersect -a ctr-ChIP_merge_peaks_blacklisted.bed -b 2312simple_naive_NAD.bed > naive_NAD_K27me3peaks.bed 
bedtools intersect -a ctr-ChIP_merge_peaks_blacklisted.bed -b 2312naive_LAD-NAD_Overlaps.bed > naive_coLADNAD_K27me3peaks.bed 
bedtools intersect -a ctr-ChIP_merge_peaks_blacklisted.bed -b 2312naive_nonLAD-NAD.bed > naive_nonLADNAD_K27me3peaks.bed 

bedtools intersect -a ctr-ChIP_merge_peaks_blacklisted.bed -b 2312naive_LAD-only.bed > naive_LAD-only_K27me3peaks.bed 
bedtools intersect -a ctr-ChIP_merge_peaks_blacklisted.bed -b 2312naive_NAD-only.bed > naive_NAD-only_K27me3peaks.bed 

#2. in R, subset Ref genes BED in each chromosome(subtract Chr.19 for example)

library(tidyverse)
library(dplyr)

hg38<- read.table("hg38.refGene_sort_rmDup.bed", sep="\t") 
hg38_chr19 <- subset(hg38,hg38$V1 == 'chr19')   
write.table(hg38_chr19, "hg38.refGene_chr19.bed", sep="\t", row.names = F, col.names = F, quote  = F)

#3. average BigWig of replicates (from step.6) for ploting.

bigwigAverage -b ctr-ChIP-rep1_CPM.bw ctr-ChIP-rep2_CPM.bw -o Average_ctr-ChIP_rep12.bw
bigwigAverage -b L1KD-ChIP-rep1_CPM.bw L1KD-ChIP-rep2_CPM.bw -o Average_L1KD-ChIP_rep12.bw


#4. computeMAtrix and Plot

##########################################
## plot at center point of K27me3 peaks ##
##########################################
computeMatrix reference-point --referencePoint center \
       -b 3000 -a 3000 \
       --regionsFileName ctr-ChIP_merge_peaks_blacklisted.bed \
       --scoreFileName Average_ctr-ChIP_rep12.bw Average_L1KD-ChIP_rep12.bw \
       --missingDataAsZero \
       --skipZeros \
       --outFileName Average_rep12-K27me3-ChIP_K27me3peaks.gz
plotProfile -m Average_rep12-K27me3-ChIP_K27me3peaks.gz --perGroup -out Average_rep12-K27me3-ChIP_K27me3peaks.pdf
plotHeatmap -m Average_rep12-K27me3-ChIP_K27me3peaks.gz -out Average_rep12-K27me3_K27me3peaks.pdf --colorMap Blues --missingDataColor "#DDEAF6"


############################################
## plot at center point of L1-ChIRP peaks ##
############################################

computeMatrix reference-point --referencePoint center \
       -b 3000 -a 3000 \
       --regionsFileName broad_l1-chirp_merge_peaks_blacklisted.bed \
       --scoreFileName Average_ctr-ChIP_rep12.bw Average_L1KD-ChIP_rep12.bw \
       --missingDataAsZero \
       --skipZeros \
       --outFileName Average_rep12-K27me3-ChIP_L1ChiRPpeaks.gz
plotProfile -m Average_rep12-K27me3-ChIP_L1ChiRPpeaks.gz --perGroup -out Average_rep12-K27me3-ChIP_L1ChiRPpeaks.pdf

#################################################################
## plot at TTS and TES of all genes within 3kb of K27me3 peaks ##
#################################################################

computeMatrix scale-regions  -S Average_ctr-ChIP_rep12.bw \
				Average_L1KD-ChIP_rep12.bw \
			     -R hg38ref_within3kb_K27me3peaks_rmDup.bed \
	--beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 \
        --missingDataAsZero \
        --skipZeros \
	-o Average_rep12-K27me3-ChIP_K27me3peaks_hg38ref.gz 
plotProfile -m Average_rep12-K27me3-ChIP_K27me3peaks_hg38ref.gz --perGroup -out Average_rep12-K27me3-ChIP_K27me3peaks_hg38ref.pdf

#################################################
##  plot at TTS and TES of genes by chromosomes##
#################################################
#(subtract Chr.19 for example)

computeMatrix scale-regions  -S Average_ctr-ChIP_rep12.bw \
				Average_L1KD-ChIP_rep12.bw \
			     -R hg38.refGene_chr19.bed  \
	--beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 \
        --missingDataAsZero \
        --skipZeros \
	-o Average_rep12-K27me3-ChIP_chr19.gz 
plotProfile -m Average_rep12-K27me3-ChIP_chr19.gz --perGroup -out Average_rep12-K27me3-ChIP_chr19.pdf


