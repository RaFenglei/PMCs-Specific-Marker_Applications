#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=10GB
#SBATCH --cpus-per-task=4
#SBATCH --job-name=blacklisting

export PATH=/home/wranalab/lcaldwel/programs/bedtools:$PATH

# Define path to blacklist
blacklist=/home/wranalab/lcaldwel/ChipSeq2024/blacklist/hg38.Kundaje.GRCh38_unified_blacklist.bed

# Remove peak file header
for file in *_peaks.xls; do
            if [[ "$file" == *"$broad"* ]];then #run normal peak code
               echo "checking peaks in broad files: $file"
               awk 'NR > 30 { print }' < $file > temp #
               awk -v OFS='\t' {'print $1,$2,$3,$9,$5,"+"'} temp > ${file%.xls}.bed

            else #run broad peak code
              echo "checking peaks in non-broad files: $file"
              awk 'NR > 23 { print }' < $file > temp #
              awk -v OFS='\t' {'print $1,$2,$3,$10,$6,"+"'} temp > ${file%.xls}.bed
            fi

done

rm temp

# Count number of peaks called in each sample
echo "counting peaks per sample, outputting to peak_counts..."
wc -l *peaks.bed > peak_counts

# Remove the blacklisted regions from peaks and quantify number of peaks lost
# Use bedtools v2.30.0 to intersect peaks with blacklisted regions and exclude
echo "running intersects with bedtools..."

for file in *peaks.bed; do
               echo "running intersect on $file"
               bedtools.static.binary intersect -v -wa -a $file -b $blacklist > ${file%.bed}_blacklisted.bed
done

echo "outputting blacklisted counts to blacklisted_counts..."
wc -l *blacklisted* > blacklisted_counts

echo "Finished chipseq step 4: blacklist removal with bedtools with exit code $? at: `date`"
