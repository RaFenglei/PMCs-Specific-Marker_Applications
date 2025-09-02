#!/bin/bash

# 定义样本列表
samples=("esophagus-1" "rumen-1" "omasum-1" "reticulum-1")

# 定义主脚本路径
main_script="/data02/zhangfenglei/project/01.omasum/03.ATAC/01.cattle/01.2025_03_15/05.tobias/00.tobias.sh"

# 定义输出目录，用于存放生成的投递脚本
output_dir="/data02/zhangfenglei/project/01.omasum/03.ATAC/01.cattle/01.2025_03_15/05.tobias/"

# 创建输出目录（如果不存在）
mkdir -p "$output_dir"

# 为每个样本生成一个投递脚本
for sample in "${samples[@]}"; do
    # 定义生成的脚本路径
    submit_script="${output_dir}/01.submit_${sample}.sh"

    # 写入内容到生成的脚本
    cat <<EOF > "$submit_script"
#!/bin/bash
bash $main_script $sample
EOF

    # 赋予生成的脚本可执行权限
    chmod +x "$submit_script"

    echo "Generated submit script: $submit_script"
done

echo "All submit scripts have been generated in: $output_dir"