#!/bin/bash

source ~/lizihe/.bashrc
source ~/lizihe/soft/anaconda3/bin/activate agat2
export PATH="/public/home/wangwen_lab/lizihe/soft/anaconda3/envs/tobias/bin:$PATH"

conda info --envs
which python

# 获取命令行传入的第一个参数，即样本名称
sample=$1

# 切换到工作目录
cd /data02/zhangfenglei/project/01.omasum/03.ATAC/01.cattle/01.2025_03_15/05.tobias

# 使用 TOBIAS 的 ATACorrect 工具对 ATAC-seq 数据进行校正
# 参数说明：
# --bam: 输入的 BAM 文件路径
# --genome: 参考基因组文件路径
# --peaks: 输入的峰值文件路径（narrowPeak 格式）
# --outdir: 输出目录路径
# --cores: 使用的 CPU 核心数
TOBIAS ATACorrect  --bam /data02/zhangfenglei/project/01.omasum/03.ATAC/01.cattle/01.2025_03_15/05.tobias/01.bam/$sample.bam --genome /data01/zhangfenglei/project/master/01.encode/01.ATAC_seq/02.cattle/00.genome/cattle.ncbi.fa  --peaks /data02/zhangfenglei/project/01.omasum/03.ATAC/01.cattle/01.2025_03_15/05.tobias/02.narrowpeak/03.all.cluster.peak_score.bed --outdir /data02/zhangfenglei/project/01.omasum/03.ATAC/01.cattle/01.2025_03_15/05.tobias/04.tobias2/$sample --cores 20


# 使用 TOBIAS 的 FootprintScores 工具计算足迹评分
# 参数说明：
# --signal: 输入的校正后的信号文件路径（bigWig 格式）
# --regions: 输入的峰值文件路径（narrowPeak 格式）
# --cores: 使用的 CPU 核心数
# -o: 输出的足迹评分文件路径（bigWig 格式）
TOBIAS FootprintScores --signal /data02/zhangfenglei/project/01.omasum/03.ATAC/01.cattle/01.2025_03_15/05.tobias/04.tobias2/$sample/${sample}_corrected.bw --regions /data02/zhangfenglei/project/01.omasum/03.ATAC/01.cattle/01.2025_03_15/05.tobias/02.narrowpeak/03.all.cluster.peak_score.bed --cores 20 -o /data02/zhangfenglei/project/01.omasum/03.ATAC/01.cattle/01.2025_03_15/05.tobias/04.tobias2/$sample/footprint.bw


# 使用 TOBIAS 的 BINDetect 工具检测转录因子结合位点
# 参数说明：
# --motifs: 输入的转录因子 motif 文件路径（JASPAR 数据库格式）
# --signals: 输入的足迹评分文件路径（bigWig 格式）
# --genome: 参考基因组文件路径
# --cores: 使用的 CPU 核心数
# --peaks: 输入的峰值文件路径（BED 格式）
# --outdir: 输出目录路径
TOBIAS BINDetect --motifs /public/home/wangwen_lab/lizihe/data/motif/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt --signals /data02/zhangfenglei/project/01.omasum/03.ATAC/01.cattle/01.2025_03_15/05.tobias/04.tobias2/${sample}/footprint.bw --genome /data01/zhangfenglei/project/master/01.encode/01.ATAC_seq/02.cattle/00.genome/cattle.ncbi.fa --cores 20 --peaks /data02/zhangfenglei/project/01.omasum/03.ATAC/01.cattle/01.2025_03_15/05.tobias/02.narrowpeak/03.all.cluster.peak_score.bed --outdir /data02/zhangfenglei/project/01.omasum/03.ATAC/01.cattle/01.2025_03_15/05.tobias/04.tobias2/$sample/footprint --cond_names ${sample}
