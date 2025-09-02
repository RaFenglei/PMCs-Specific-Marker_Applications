# 激活 Conda 环境
source /public/home/wangwen_lab/lizihe/soft/anaconda3/bin/activate atac

# 基因组 fasta 文件
genome="/data01/zhangfenglei/project/master/01.encode/01.ATAC_seq/02.cattle/00.genome/cattle.ncbi.fa"

# HOMER motif 集合设置
lineage="vertebrates"
cores=20

# 工作路径设置
input_dir="/public/home/wangwen_lab/wangxiaoyang/workspace/ATAC/07.homer"
output_dir="/public/home/wangwen_lab/wangxiaoyang/workspace/ATAC/07.homer"

# 分析的样本名列表
samples=("only_ABPC_v2.bed" "only_PP_v2.bed" "only_SFRP2_v2.bed")

# 创建输出路径
mkdir -p "$output_dir"

# 循环执行 motif 分析
for bed_file in "${samples[@]}"; do
    name=$(basename "$bed_file" .bed)
    peak_path="$input_dir/$bed_file"
    out_path="$output_dir/${name}_motif"

    echo "Running HOMER motif analysis for: $bed_file"
    findMotifsGenome.pl "$peak_path" "$genome" "$out_path" -size given -p "$cores" -mset "$lineage"
done
