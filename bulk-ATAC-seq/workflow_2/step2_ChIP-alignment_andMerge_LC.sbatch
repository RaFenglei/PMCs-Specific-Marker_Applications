#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=20GB
#SBATCH --cpus-per-task=8
#SBATCH --job-name chipstep1

##bowtie2 cant run in parallel so don't request more nodes

export PATH=/home/wranalab/lcaldwel/programs/bowtie2-2.5.3:$PATH
export PATH=/home/wranalab/lcaldwel/programs/fastqc:$PATH
export PATH=/home/wranalab/lcaldwel/programs/samtools-1.18:$PATH

which bowtie2
which fastqc
which samtools

# Define path to genome indices
genome_index=/home/wranalab/lcaldwel/ChipSeq2024/genome_files/hg38

for file in *_L001_R1_001.fastq.gz; do
  sample=${file%_L001_R1_001.fastq.gz}

  echo $sample

  bowtie2 -p 6 --local -X 1000 -x $genome_index -1 $file -2 ${sample}_L001_R2_001.fastq.gz -S ${sample}.sam 2>${sample}.align.stats
  # -p is threads;
  # --local means Bowtie 2 does not require that the entire read align from one end to the other. instead it will soft clip from the ends in order to achieve the greatest possible alignment score.
  # -X is The maximum fragment length for valid paired-end alignments. E.g. if -X 100 is specified and a paired-end alignment consists of two 20-bp alignments in the proper orientation with a 60-bp gap between them, that alignment is considered valid (as long as -I is also satisfied).
  # -x is index base name
  # -S is output sam name

	 echo "wc -l aligned"
	 wc -l ${sample}.sam

# Count reads assigned to chrM, Un, random, then remove -- this info is outputted to slurm out file!!
	echo "wc -l chrM"
	grep chrM ${sample}.sam | wc -l
	echo "wc -l chrUn"
	grep chrUn ${sample}.sam | wc -l
	echo "wc -l random"
	grep random ${sample}.sam | wc -l
	sed '/chrM/d;/random/d;/chrUn/d' < ${sample}.sam > ${sample}removedchrs.sam

# Convert to BAM, sort
	samtools view -bu ${sample}removedchrs.sam | samtools sort -o sorted_${sample}.bam

# Remove duplicates - don’t need any of these files later, so name intermediate files redundantly
	samtools collate -o namecollate.bam sorted_${sample}.bam
	samtools fixmate -m namecollate.bam fixmate.bam
  samtools sort -o positionsort.bam fixmate.bam
	samtools markdup -rs positionsort.bam filtered_${sample}.bam
	#markdup flag the duplicates and -r argument remove them

# Index bam
	samtools index filtered_${sample}.bam

# For completeness sake, get a line count of the number of reads passing all filters, then remove the file
  samtools view filtered_${sample}.bam -o filtered_${sample}.sam
  echo "wc -l filtered"
  wc -l filtered_${sample}.sam
	rm filtered_${sample}.sam

# Remove intermediate sam files
	# Only uncomment this if you are sure this script is working for your samples; this is just here to save space
  rm ${sample}.sam
	rm ${sample}removedchrs.sam
  rm namecollate.bam
  rm fixmate.bam
  rm positionsort.bam

done

# Compile alignment statistics with multiqc
multiqc . --ignore "FastQC/"

echo "Finished chipseq step 1 run with exit code $? at: `date`"

### Here we merge and run normalize with CPM
#export PATH=/home/wranalab/lcaldwel/.local/bin/macs3:$PATH
export PATH=/home/wranalab/lcaldwel/.local/bin/deeptools:$PATH


#which macs3 #macs3 doesn't work on cluster because it has 3.8 and macs3 needs 3.9+
which deeptools

###merge with samtools
echo "merging with samtools..."
samtools merge ctr-ChIP_merge.bam filtered_s3.bam filtered_s5.bam 
samtools merge L1KD-ChIP_merge.bam filtered_s4.bam filtered_s6.bam 
samtools index ctr-ChIP_merge.bam
samtools index L1KD-ChIP_merge.bam


#get coverage
for file in *merge.bam; do
  bamCoverage -b ctr-ChIP_merge.bam -o ctr-ChIP_merge.bw -v --normalizeUsing CPM
  bamCoverage -b L1KD-ChIP_merge.bam -o L1KD-ChIP_merge.bw -v --normalizeUsing CPM
done

echo "Finished chipseq step 3: Merging and deeptools run with exit code $? at: `date`"
