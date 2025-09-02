#!/bin/bash
#PBS -N ATAC_omasum_1
#PBS -q high1
#PBS -l nodes=1:ppn=32
#PBS -j oe
#PBS -o omasum_1.pbs.log

source /public/home/wangwen_lab/zhoubotong/.bashrc
echo "Starting analysis for omasum_1"
date

# Step 1: Quality control with fastp
echo "Running fastp"
bash 01.omasum_1.fastp.sh
date

# Step 2: Alignment with bowtie2
echo "Running bowtie2"
bash 02.omasum_1.bowtie2.sh
date

# Step 3: Remove duplicates with picard
echo "Running picard"
bash 03.omasum_1.picard.sh
date

# Step 4: Generate coverage with bamCoverage
echo "Running bamCoverage"
bash 04.omasum_1.bamCoverage.sh
date

# Step 5: Convert to BED format
echo "Running bedtools"
bash 05.omasum_1.bedtools.sh
date

# Step 6: Peak calling with macs3
echo "Running macs3"
bash 06.omasum_1.callPeak.sh
date

# Step 7: Fragment length analysis
echo "Running fragment length analysis"
bash 07.omasum_1.fragmentsLength.sh
date

# Step 8: Compute matrix around TSS (optional)
echo "Running computeMatrix (TSS)"
bash 08.omasum_1.computeMatrix.sh
date

# Step 9: Plot heatmap of TSS signal (optional)
echo "Running plotHeatmap (TSS)"
bash 09.omasum_1.plotHeatmap.sh
date

# Step 10: Calculate FRiP score
echo "Running FRiP calculation"
bash 10.omasum_1.FRiP.sh
date

echo "Analysis completed for omasum_1"
date
