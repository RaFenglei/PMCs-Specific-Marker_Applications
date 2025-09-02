import pandas as pd
import numpy as np

# 读取完整的CSV文件
df = pd.read_csv('/data02/zhangfenglei/project/01.omasum/01.RNA-seq/02.cattle_ref/01.cattle/01.LZH_2025_05_27/02.stringtie/00.gtf/01.gene_E53_TPM.csv')

# 设置gene_id为索引
df.set_index('gene_id', inplace=True)

# 找到表达量最高的organ
df['max_organ'] = df.idxmax(axis=1)

# 筛选出至少有一个organ表达量>=1的基因，只使用数值列计算最大值
expression_columns = df.columns.drop('max_organ')
df = df[df[expression_columns].max(axis=1) >= 1]

# 计算tau值函数
def calculate_tau(row):
    max_val = row.max()
    if max_val == 0:
        return 0
    tau = ((1 - (row / max_val)).sum()) / (len(row) - 1)
    return tau

# 计算每个基因的tau值，排除max_organ列
df['tau'] = df[expression_columns].apply(calculate_tau, axis=1)

# 筛选出tau>0.8且在特定器官中表达量最高的基因
# 注意：您需要将'E53_omasum_2'替换为您数据中实际的器官名称
specific_organ_df = df[(df['tau'] > 0.8) & (df['max_organ'] == 'E53_omasum_2')]

# 保存结果到TSV文件（制表符分隔）
df.to_csv('E53_combined_expression_with_tau.tsv', sep='\t')
specific_organ_df.to_csv('E53_specific_organ_gene.tsv', sep='\t')

# 可选：同时保存CSV文件（逗号分隔）
df.to_csv('E53_combined_expression_with_tau.csv', sep=',')
specific_organ_df.to_csv('E53_specific_organ_gene.csv', sep=',')

# 提取特异性基因列表并处理重复ID
def extract_processed_gene_list(gene_id_list):
    processed_genes = []
    for gene_id in gene_id_list:
        # 按"|"分割并取第一个元素（处理形如"CCDC181|CCDC181"的ID）
        processed_id = gene_id.split('|')[0]
        processed_genes.append(processed_id)
    return processed_genes

# 提取特异性基因列表并保存为文本文件
gene_ids = specific_organ_df.index.tolist()
processed_gene_ids = extract_processed_gene_list(gene_ids)

# 保存到TXT文件
output_file = 'E53_specific_organ_gene_list.txt'
with open(output_file, 'w') as f:
    for gene_id in processed_gene_ids:
        f.write(f"{gene_id}\n")

print(f"成功分析基因表达并提取了 {len(processed_gene_ids)} 个特异性基因ID")
print(f"结果已保存到:")
print(f"  - E53_combined_expression_with_tau.tsv/csv (完整数据)")
print(f"  - E53_specific_organ_gene.tsv/csv (特异性基因数据)")
print(f"  - {output_file} (特异性基因列表)")
