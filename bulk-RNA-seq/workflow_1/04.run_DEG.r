#install.packages("tidyverse")
#install.packages("ggrepel")
#install.packages("ggplot2"
# 导入必要的数据操作和可视化库
library(tidyverse)
library(ggrepel)
library(ggplot2)

# 从.txt文件读取数据，假设存在标题并且字段由制表符分隔
DEGs <- read.csv('mhl_gumo_lncRNA_transcripts_line1-new.txt', header=TRUE, sep="\t")

# 设置基于log2倍数变化和调整后P值的筛选阈值
threshold_log2FoldChange <- log2(2)  # 2倍变化的log2值
threshold_padj <- 0.05               # 调整后的P值阈值

# 清理数据，移除没有padj或log2FoldChange值的项，并初始化一个空的标签列
DEGs <- DEGs[!is.na(DEGs$padj) & !is.na(DEGs$log2FoldChange),]

# 定义显著性标准
DEGs$significance <- DEGs$padj <= threshold_padj & abs(DEGs$log2FoldChange) >= threshold_log2FoldChange

# 对DEGs根据调整后P值进行排序
DEGs <- DEGs %>%
  arrange(padj) %>%
  mutate(rank = row_number())

# 初始化标签列
DEGs$label <- NA

# 只标记最显著的30个基因
DEGs$label[DEGs$rank <= 30] <- DEGs$gene[DEGs$rank <= 30]

# 统计上调和下调的数值
DEGs$change <- ifelse(DEGs$log2FoldChange >= threshold_log2FoldChange, "Up", 
                      ifelse(DEGs$log2FoldChange <= -threshold_log2FoldChange, "Down", "Neutral"))
up_down_genes <- DEGs[DEGs$significance & DEGs$change %in% c("Up", "Down"),]

# 将所有上调和下调的基因输出到一个文本文件中
write.table(up_down_genes, 'up_down_genes.txt', sep="\t", row.names=FALSE, quote=FALSE)

# 统计上调和下调的数值
up_count <- sum(DEGs$change == "Up" & DEGs$significance)
down_count <- sum(DEGs$change == "Down" & DEGs$significance)

# 使用ggplot2包准备火山图
volcano <- ggplot(data = DEGs, aes(x = log2FoldChange, y = -log10(padj), color = log2FoldChange, size = -log10(padj))) +
  geom_point(alpha=0.5) +
  geom_vline(xintercept = c(-threshold_log2FoldChange, threshold_log2FoldChange), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(threshold_padj), linetype = "dashed", color = "black") +
  geom_text_repel(aes(label = label), 
                  box.padding = unit(0.35, "lines"), 
                  point.padding = unit(0.3, "lines"), size = 2.5, segment.color = 'grey50') +
 scale_color_gradientn(colors = c("#3288bd", "#9e0142"), 
                        values = seq(0, 1, length.out = 5)) +
  labs(title = paste0("Volcano plot\nUp-regulated: ", up_count, " Down-regulated: ", down_count),
       x = "log2 Fold Change", y = "-log10 Adjusted P-Value") +
  theme_minimal() +
  theme(legend.position = "right", 
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        plot.margin = margin(5.5, 40, 5.5, 5.5))

# 将图表保存为PNG
ggsave('volcano_plot.png', volcano, width = 12, height = 8)

# 将图表保存为PDF
ggsave('volcano_plot.pdf', volcano, width = 12, height = 8)

# 渲染图表
volcano