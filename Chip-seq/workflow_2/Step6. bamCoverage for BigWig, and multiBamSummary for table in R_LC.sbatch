#!/bin/bash
#SBATCH --job-name=bamOpt1
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=20GB


#1.	BamCoverage on all individual files with CPM norm. Output as bigwig.
#2.	MultibigwigSummary, with all individual files (so rep1 and 2 separately), with bins = 50, outputting the raw table with --outRawCounts
#3.	Then I can put the resulting table from MBWS into R, and manually get mean values between reps and make the plots with either individual replicates or as mean values


export PATH=/home/wranalab/lcaldwel/programs/bedtools:$PATH
#uses deeptools which I installed with pip
blacklist=/home/wranalab/lcaldwel/ChipSeq2024/blacklist/hg38.Kundaje.GRCh38_unified_blacklist.bed


#Bam cov is the same for both so only do it once.

echo "Starting BamCoverage..."
bamCoverage -b filtered_S1.bam -o ctr-input_CPM.bw --outFileFormat bigwig --normalizeUsing CPM -bl $blacklist --effectiveGenomeSize 2862010578 --ignoreDuplicates -p 8
bamCoverage -b filtered_s2.bam -o L1KD-input_CPM.bw --outFileFormat bigwig --normalizeUsing CPM -bl $blacklist --effectiveGenomeSize 2862010578 --ignoreDuplicates -p 8
bamCoverage -b filtered_s3.bam -o ctr-ChIP-rep1_CPM.bw --outFileFormat bigwig --normalizeUsing CPM -bl $blacklist --effectiveGenomeSize 2862010578 --ignoreDuplicates -p 8
bamCoverage -b filtered_s4.bam -o L1KD-ChIP-rep1_CPM.bw --outFileFormat bigwig --normalizeUsing CPM -bl $blacklist --effectiveGenomeSize 2862010578 --ignoreDuplicates -p 8
bamCoverage -b filtered_s5.bam -o ctr-ChIP-rep2_CPM.bw --outFileFormat bigwig --normalizeUsing CPM -bl $blacklist --effectiveGenomeSize 2862010578 --ignoreDuplicates -p 8
bamCoverage -b filtered_s6.bam -o L1KD-ChIP-rep2_CPM.bw --outFileFormat bigwig --normalizeUsing CPM -bl $blacklist --effectiveGenomeSize 2862010578 --ignoreDuplicates -p 8

##option 1 - 	MultibigwigSummary, with all individual files (so rep1 and 2 separately), with bins = 50, outputting the raw table with --outRawCounts
echo "Starting MultiBigWigSummary..."
multiBigwigSummary bins -bs 50 -b ctr-input_CPM.bw ctr-ChIP-rep1_CPM.bw ctr-ChIP-rep2_CPM.bw L1KD-input_CPM.bw L1KD-ChIP-rep1_CPM.bw L1KD-ChIP-rep2_CPM.bw -o Average_CPM_bins50bp_MBWSummary_Option1.bw.npz -bl $blacklist -p 8 --outRawCounts Average_CPM_bins50bp_MBWSummary_Option1.txt


echo "Finished chipseq option 1 comparison: 50 bp bin size analysis, with exit code $? at: `date`"
