#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=5GB
#SBATCH --cpus-per-task=8
#SBATCH --job-name chipQC

export PATH=/home/wranalab/lcaldwel/programs/bowtie2-2.5.3:$PATH
export PATH=/home/wranalab/lcaldwel/programs/fastqc:$PATH

fastqc *.gz

multiqc . --ignore "FastQC/"

echo "completed FastQC and MultiQC with exit code $? at: `date`"
