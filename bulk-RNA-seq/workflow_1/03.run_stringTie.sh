#/public/home/liyongxin/lihaorong/soft/stringtie-2.1.4.Linux_x86_64/stringtie RB2.sort.bam -G XTJD_hic.gtf_chr -o RB2.gtf -p 10 -e

cd /data02/zhangfenglei/project/08.telomere/01.ABPC/01.RNA-seq/03.stringTie
mkdir -p 00.gtf 00.gtfqs


for i in `cat /data02/zhangfenglei/project/08.telomere/01.ABPC/01.RNA-seq/00.raw_data/all.name.list`;
do
	echo "/public/home/liyongxin/lihaorong/soft/stringtie-2.1.4.Linux_x86_64/stringtie /data02/zhangfenglei/project/08.telomere/01.ABPC/01.RNA-seq/02.hisat/01.hisat/${i}.mhl.sorted.bam -G /public/home/wangwen_lab/wangbin/project/myf/zfl/01.mhl/01.index/mhl.gtf -o /data02/zhangfenglei/project/08.telomere/01.ABPC/01.RNA-seq/03.stringTie/00.gtf/${i}.mhl.gtf -p 10 -e;" > /data02/zhangfenglei/project/08.telomere/01.ABPC/01.RNA-seq/03.stringTie/00.gtfqs/${i}.mhl.sh
done








