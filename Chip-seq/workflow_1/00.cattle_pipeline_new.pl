use warnings;
use strict;

##Step1:初始设置和变量定义
my @samples = ("omasum_1");  # 定义所有样本名称，可以是一个也可以是多个
my $postfix = "fastq.gz"; #也可以是fq.gz
my $raw_data_dir = "/data02/zhangfenglei/project/01.omasum/03.ATAC/01.cattle/04.2025_07_03/00.raw_data"; #存放fq.gz的目录
my $refIdx = "/data02/zhangfenglei/project/01.omasum/03.ATAC/01.cattle/02.2025_06_05/00.rawdata/cattle_index"; # 修改为你的索引前缀,Path to reference index. U need to index reference first (bowtie2-build -f --threads 4 giraffe.fa giraffe_index/bwa index giraffe.fa)
my $tssref = "/data02/zhangfenglei/project/01.omasum/03.ATAC/01.cattle/02.2025_06_05/00.rawdata/cattle_tss_points.bed";

##Step2:软件路径定义
my $fastp = "/public/home/wangwen_lab/zhangfenglei/soft/fastp_0_23_4/fastp";
my $fastpCMD = "--detect_adapter_for_pe --trim_poly_g --length_required 30 --thread 16 --dedup --dup_calc_accuracy 5";
my $bowtie = "/public/software/bowtie2/bowtie2-2.3.4.1-linux-x86_64/bowtie2";
my $picard = "/public/home/wangwen_lab/zhoubotong/soft/jdk-11.0.15+10/bin/java -jar -Xms43152m -Xmx43152m /public/home/wangwen_lab/zhoubotong/soft/picard/picard_2.25.7.jar MarkDuplicates";
my $bamCoverage = "/public/home/wangwen_lab/zhoubotong/soft/conda3/bin/bamCoverage";
my $bedtools = "/public/home/wangwen_lab/zhoubotong/soft/bedtools/bedtools2-master/bin/bedtools";
my $macs3 = "/public/home/wangwen_lab/zhoubotong/soft/conda3/bin/macs2";
my $gatk = "/public/home/wangwen_lab/zhoubotong/soft/conda3/bin/gatk CollectInsertSizeMetrics";

## 为每个样本创建独立的分析流程
foreach my $sample (@samples) {
    ## Step3: 输出目录创建
    my $outdir = "$sample.output";
    my $qualiControl = "$sample.qualiControl";
    my $callpeak = "$sample.callpeak";
    `mkdir -p $outdir` if (!-d $outdir);
    `mkdir -p $qualiControl` if (!-d $qualiControl);
    `mkdir -p $callpeak` if (!-d $callpeak);
    
    ## Step4: 脚本文件打开
    open O, "> 01.$sample.fastp.sh" or die "Cannot open 01.$sample.fastp.sh: $!";
    open F, "> 02.$sample.bowtie2.sh" or die "Cannot open 02.$sample.bowtie2.sh: $!";
    open K, "> 03.$sample.picard.sh" or die "Cannot open 03.$sample.picard.sh: $!";
    open L, "> 04.$sample.bamCoverage.sh" or die "Cannot open 04.$sample.bamCoverage.sh: $!";
    open P, "> 05.$sample.bedtools.sh" or die "Cannot open 05.$sample.bedtools.sh: $!";
    open D, "> 06.$sample.callPeak.sh" or die "Cannot open 06.$sample.callPeak.sh: $!";
    open A, "> 07.$sample.fragmentsLength.sh" or die "Cannot open 07.$sample.fragmentsLength.sh: $!";
    open I, "> 08.$sample.computeMatrix.sh" or die "Cannot open 08.$sample.computeMatrix.sh: $!";
    open U, "> 09.$sample.plotHeatmap.sh" or die "Cannot open 09.$sample.plotHeatmap.sh: $!";
    open R, "> 10.$sample.FRiP.sh" or die "Cannot open 10.$sample.FRiP.sh: $!";
    
    ## Create master submission script
    open M, "> 00.$sample.submit_all.sh" or die "Cannot open 00.$sample.submit_all.sh: $!";
    print M "#!/bin/bash\n";
    print M "#PBS -N ATAC_$sample\n";
    print M "#PBS -q high1\n";
    print M "#PBS -l nodes=1:ppn=32\n";
    print M "#PBS -j oe\n";
    print M "#PBS -o $sample.pbs.log\n\n";
    
    print M "source /public/home/wangwen_lab/zhoubotong/.bashrc\n";
    print M "echo \"Starting analysis for $sample\"\n";
    print M "date\n\n";
    
    # Add all analysis steps to the master script
    print M "# Step 1: Quality control with fastp\n";
    print M "echo \"Running fastp\"\n";
    print M "bash 01.$sample.fastp.sh\n";
    print M "date\n\n";
    
    print M "# Step 2: Alignment with bowtie2\n";
    print M "echo \"Running bowtie2\"\n";
    print M "bash 02.$sample.bowtie2.sh\n";
    print M "date\n\n";
    
    print M "# Step 3: Remove duplicates with picard\n";
    print M "echo \"Running picard\"\n";
    print M "bash 03.$sample.picard.sh\n";
    print M "date\n\n";
    
    print M "# Step 4: Generate coverage with bamCoverage\n";
    print M "echo \"Running bamCoverage\"\n";
    print M "bash 04.$sample.bamCoverage.sh\n";
    print M "date\n\n";
    
    print M "# Step 5: Convert to BED format\n";
    print M "echo \"Running bedtools\"\n";
    print M "bash 05.$sample.bedtools.sh\n";
    print M "date\n\n";
    
    print M "# Step 6: Peak calling with macs3\n";
    print M "echo \"Running macs3\"\n";
    print M "bash 06.$sample.callPeak.sh\n";
    print M "date\n\n";
    
    print M "# Step 7: Fragment length analysis\n";
    print M "echo \"Running fragment length analysis\"\n";
    print M "bash 07.$sample.fragmentsLength.sh\n";
    print M "date\n\n";
    
    # Conditional inclusion of TSS-related steps (08 and 09)
    print M "# Step 8: Compute matrix around TSS (optional)\n";
    print M "echo \"Running computeMatrix (TSS)\"\n";
    print M "bash 08.$sample.computeMatrix.sh\n";
    print M "date\n\n";
    
    print M "# Step 9: Plot heatmap of TSS signal (optional)\n";
    print M "echo \"Running plotHeatmap (TSS)\"\n";
    print M "bash 09.$sample.plotHeatmap.sh\n";
    print M "date\n\n";
    
    print M "# Step 10: Calculate FRiP score\n";
    print M "echo \"Running FRiP calculation\"\n";
    print M "bash 10.$sample.FRiP.sh\n";
    print M "date\n\n";
    
    print M "echo \"Analysis completed for $sample\"\n";
    print M "date\n";
    close M;
    
    ## Step5: 主循环处理每个样本的fastq文件
    my $fq_1 = "$raw_data_dir/${sample}_R1.$postfix";
    my $fq_2 = "$raw_data_dir/${sample}_R2.$postfix";
    # 检查文件是否存在
    unless (-e $fq_1 && -e $fq_2) {
        warn "Warning: Missing fastq files for sample $sample\n";
        next;
    }
    # 创建样本特定的输出目录
    my $tmpoutdir = "$outdir";
    my $quaCondir = "$qualiControl";
    my $callpeakdir = "$callpeak";
    `mkdir -p $tmpoutdir` if (!-d $tmpoutdir);
    `mkdir -p $quaCondir` if (!-d $quaCondir);
    `mkdir -p $callpeakdir` if (!-d $callpeakdir);
    
    ## fastp质量控制
    print O "$fastp -i $fq_1 -I $fq_2 -o $tmpoutdir/${sample}_1.fq.gz -O $tmpoutdir/${sample}_2.fq.gz -h $tmpoutdir/${sample}.html $fastpCMD\n";

    
    ## bowtie2比对
    print F "$bowtie --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -X2000 -p 40 -x $refIdx -1 $tmpoutdir/${sample}_1.fq.gz -2 $tmpoutdir/${sample}_2.fq.gz 2> $tmpoutdir/${sample}.bowtie.log | samtools view -\@10 -F 1804 -h -f 2 -q 25 | samtools sort -\@10 -m 2G -O BAM -o $tmpoutdir/${sample}.bam\n";
    
    ## picard去重
    print K "$picard INPUT=$tmpoutdir/${sample}.bam OUTPUT=$tmpoutdir/${sample}.rmdup.bam METRICS_FILE=$tmpoutdir/${sample}.rmdup.bam.dupqc.txt VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=true; samtools index -\@10 $tmpoutdir/${sample}.rmdup.bam\n";
    
    ## bamCoverage生成覆盖度文件
    print L "$bamCoverage -b $tmpoutdir/${sample}.rmdup.bam -o $tmpoutdir/${sample}.rmdup.bw -of bigwig --binSize 10 --normalizeUsing RPKM --extendReads \n";
    
    ## bedtools转换格式
    print P "$bedtools bamtobed -i $tmpoutdir/${sample}.rmdup.bam > $tmpoutdir/${sample}.rmdup.bed \n";
    
    ## macs3 peak calling
    print D "$macs3 callpeak -t $tmpoutdir/${sample}.rmdup.bam -n $sample --outdir $callpeakdir -B -f BAMPE -q 0.01 -g 2715853792\n"; # -g 2715853792 为cattle基因组大小
    
    ## 片段长度分析
    print A "$gatk -H $quaCondir/${sample}.fragmentslength.pdf -M 0.5 -I $tmpoutdir/${sample}.rmdup.bam -O $quaCondir/${sample}.InsertSize.txt \n";
    
    ## 计算TSS附近信号矩阵 (optional)
    print I "computeMatrix reference-point --referencePoint TSS -S $tmpoutdir/${sample}.rmdup.bw -R $tssref -b 2000 -a 2000 -out $quaCondir/${sample}.scale_regions.tab.gz --skipZeros --missingDataAsZero\n";
    
    ## 绘制热图 (optional)
    print U "plotHeatmap -m $quaCondir/${sample}.scale_regions.tab.gz -out $quaCondir/${sample}.plotHeatmap.png --legendLocation none\n";
    
    ## 计算FRiP score
    print R "total_reads=\$(samtools view -c $tmpoutdir/${sample}.rmdup.bam)\n";
    print R "reads_in_peaks=\$(bedtools intersect -u -a $tmpoutdir/${sample}.rmdup.bam -b $callpeakdir/${sample}_peaks.narrowPeak | samtools view -c)\n";
    print R "frip=\$(echo \"scale=4; \$reads_in_peaks/\$total_reads\" | bc)\n";
    print R "echo \"FRiP score for $sample: \$frip\" > $quaCondir/${sample}.FRiP.txt\n";
    print R "echo \"Total reads: \$total_reads\" >> $quaCondir/${sample}.FRiP.txt\n";
    print R "echo \"Reads in peaks: \$reads_in_peaks\" >> $quaCondir/${sample}.FRiP.txt\n";
    
    ## 关闭文件句柄
    close O;
    close F;
    close K;
    close L;
    close P;
    close D;
    close A;
    close I;
    close U;
    close R;
    
    print "Generated analysis scripts for sample: $sample\n";
    print "Master submission script created: 00.$sample.submit_all.sh\n";
    print "You can submit the job with: qsub 00.$sample.submit_all.sh\n";
}