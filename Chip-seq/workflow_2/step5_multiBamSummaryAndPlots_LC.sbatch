#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=5GB
#SBATCH --cpus-per-task=4
#SBATCH --job-name=multiBamSummary

###these tools are just part of bedtools still
export PATH=/home/wranalab/lcaldwel/programs/bedtools:$PATH


#Example
#multiBamSummary bins \
 #--bamfiles testFiles/*bam \ # using all BAM files in the folder
 #--minMappingQuality 30 \
 #--region 19 \ # limiting the binning of the genome to chromosome 19
 #--labels H3K27me3 H3K4me1 H3K4me3 HeK9me3 input \
 #-out readCounts.npz --outRawCounts readCounts.tab

#should be able to just call these, don't need to reference bedtools specifically?
#this set of commands doesn't work with number of processors command, it just gives a bin error ????

echo "starting multibam summary of ctr-ChIP..."
multiBamSummary bins \
--bamfiles filtered_s3.bam filtered_s5.bam \
--labels ctr-ChIP-rep1 ctr-ChIP-rep2 -p 6 \
-o ctr-ChIP.npz --outRawCounts ctr-ChIP-readCounts.tab  #-p 6#--numberOfProcessors 8

echo "starting multibam summary of L1KD-ChIP..."
multiBamSummary bins \
--bamfiles filtered_s4.bam filtered_s6.bam \
--labels L1KD-ChIP-rep1 L1KD-ChIP-rep2 \
-o L1KD-ChIP.npz --outRawCounts L1KD-ChIP-readCounts.tab -p 6  #-p 6#--numberOfProcessors 8

echo "starting multibam summary of K27me3-ChIP-all..."
multiBamSummary bins \
--bamfiles filtered_S1.bam filtered_s2.bam filtered_s3.bam filtered_s4.bam filtered_s5.bam \
filtered_s6.bam \
--labels ctr-input L1KD-input ctr-ChIP-rep1 L1KD-ChIP-rep1 ctr-ChIP-rep2 L1KD-ChIP-rep2 \
-o K27me3-ChIP-all.npz --outRawCounts K27me3-ChIP-all-readCounts.tab -p 6 #-p 6#--numberOfProcessors 8
#can use this .npz file later for the PCA


echo "plotting PCA..."
#this is one that I added
plotPCA -in K27me3-ChIP-all.npz \
-o K27me3-ChIP-all_PCA_readCounts.png \
-T "PCA of read counts, K27me3-ChIP-all"


###plotFingerprint
#This quality control will most likely be of interest for you if you are dealing with ChIP-seq
#samples as a pressing question in ChIP-seq experiments is “Did my ChIP work?”, i.e. did the antibody-treatment
#enrich sufficiently so that the ChIP signal can be separated from the background signal?
echo "plotting fingerprints..."
#(After all, around 90% of all DNA fragments in a ChIP experiment will represent the genomic background).

plotFingerprint -b filtered_S1.bam filtered_s2.bam filtered_s3.bam filtered_s4.bam \
filtered_s5.bam filtered_s6.bam \
--labels ctr-input L1KD-input ctr-ChIP-rep1 L1KD-ChIP-rep1 ctr-ChIP-rep2 L1KD-ChIP-rep2 \
-T "Fingerprints of K27me3-ChIP all samples"  \
--plotFile plotFingerprint_K27me3-ChIP-all.png \
--outRawCounts plotFingerprint_K27me3-ChIP-all.tab

echo "plotting correlation scatterplot..."
#scatterplot
plotCorrelation \
-in K27me3-ChIP-all.npz \
--corMethod pearson --skipZeros \
--plotTitle "Pearson Correlation of K27me3-ChIP-all" \
--whatToPlot scatterplot \
-o scatterplot_K27me3-ChIP-all_PearsonCorr_bigwigScores.png   \
--outFileCorMatrix K27me3-ChIP-all_PearsonCorr_bigwigScores.tab

echo "plotting correlation heatmap..."
#heatmap
plotCorrelation \
    -in K27me3-ChIP-all.npz \
    --corMethod spearman --skipZeros \
    --plotTitle "Spearman Correlation of K27me3-ChIP-all" \
    --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
    -o heatmap_K27me3-ChIP-all_SpearmanCorr_readCounts.png   \
    --outFileCorMatrix K27me3-ChIP-all_SpearmanCorr_readCounts.tab



echo "Finished chipseq step 5: multi-bam summary and plots with exit code $? at: `date`"
