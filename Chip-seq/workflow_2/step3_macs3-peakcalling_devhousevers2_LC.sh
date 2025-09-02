#!/bin/bash


echo "calling peaks round 1 with macs..."
# -t is treatment.
# -n is name of outfile.

# for broad cutoff: here you use the control samples as a background.
echo "broad calling..."
macs3 callpeak -t filtered_s3.bam -c filtered_S1.bam --broad -g hs --broad-cutoff 0.1 -n broad_ctr-ChIP-rep1 -f BAMPE
macs3 callpeak -t filtered_s4.bam -c filtered_s2.bam --broad -g hs --broad-cutoff 0.1 -n broad_L1KD-ChIP-rep1 -f BAMPE
macs3 callpeak -t filtered_s5.bam -c filtered_S1.bam --broad -g hs --broad-cutoff 0.1 -n broad_ctr-ChIP-rep2 -f BAMPE
macs3 callpeak -t filtered_s6.bam -c filtered_s2.bam --broad -g hs --broad-cutoff 0.1 -n broad_L1KD-ChIP-rep2 -f BAMPE

#broad merge now.
echo "broad merged calling..."
macs3 callpeak -t ctr-ChIP_merge.bam -c filtered_S1.bam --broad -g hs --broad-cutoff 0.1 -n broad_ctr-ChIP_merge -f BAMPE
macs3 callpeak -t L1KD-ChIP_merge.bam -c filtered_s2.bam --broad -g hs --broad-cutoff 0.1 -n broad_L1KD-ChIP_merge -f BAMPE

echo "Finished chipseq MACS3 run with exit code $? at: `date`"
