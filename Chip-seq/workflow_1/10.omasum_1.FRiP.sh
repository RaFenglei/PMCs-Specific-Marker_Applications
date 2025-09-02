total_reads=$(samtools view -c omasum_1.output/omasum_1.rmdup.bam)
reads_in_peaks=$(bedtools intersect -u -a omasum_1.output/omasum_1.rmdup.bam -b omasum_1.callpeak/omasum_1_peaks.narrowPeak | samtools view -c)
frip=$(echo "scale=4; $reads_in_peaks/$total_reads" | bc)
echo "FRiP score for omasum_1: $frip" > omasum_1.qualiControl/omasum_1.FRiP.txt
echo "Total reads: $total_reads" >> omasum_1.qualiControl/omasum_1.FRiP.txt
echo "Reads in peaks: $reads_in_peaks" >> omasum_1.qualiControl/omasum_1.FRiP.txt
