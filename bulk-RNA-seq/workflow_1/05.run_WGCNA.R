library(WGCNA)
setwd("/data02/zhangfenglei/project/01.omasum/01.RNA-seq/02.cattle_ref/01.cattle/01.LZH_2025_05_27/04.WGCNA/04.all")

data = read.csv("00.cattle_all_gene_TPM.csv",header = T,row.names = 1,check.names=F)
dim(data)
head(data[, 1:4])
class <- read.delim("00.class.txt3", header=TRUE, 
                   colClasses=c("tissue"="character", "class"="character"))
head(data[, 1:4])
RNAseq_voom <- data
RNAseq_voom <- log2(RNAseq_voom+1)
WGCNA_matrix = t(RNAseq_voom[order(apply(RNAseq_voom,1,mad), decreasing = T)[1:20000],])

# apply() 函数用于对数组或矩阵的指定维度应用一个函数。
# 参数 RNAseq_voom 是输入的数据框或矩阵。
# 参数 1 表示按行操作（如果为 2 则表示按列操作）。
# 参数 mad 是一个函数，计算每行数据的中位数绝对偏差（Median Absolute Deviation）。mad 是一种衡量数据离散程度的方法，对异常值不敏感。

datExpr0 <- WGCNA_matrix
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
datExpr <- datExpr0
sampleNames = rownames(datExpr)
head(sampleNames)
# 这一段是为了调整矩阵的顺序
head(class)
datTraits = class
traitRows = match(sampleNames, datTraits$tissue)
rownames(datTraits) = datTraits[traitRows, 1]
save(datExpr,datTraits,file = 'wgcna-input.RData')
# 绘制软阈值和连接值
pdf("softthreshold.pdf")
powers = c(seq(1, 10, by=1), seq(12, 30, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
  ##sizeGrWindow(9, 5)
  par(mfrow = c(1,2));
  cex1 = 0.9;
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.90,col="red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
#自动算阈值
sft$powerEstimate

# 最重要的聚类计算，计算相关性（tom）
net = blockwiseModules(
                 datExpr,
                 power = sft$powerEstimate,
                 maxBlockSize = 6000,
                 TOMType = "signed", minModuleSize = 30,
                 reassignThreshold = 0, mergeCutHeight = 0.25,
                 numericLabels = FALSE, pamRespectsDendro = FALSE,
                 saveTOMs = F, 
                 verbose = 3
 )


#把颜色向量转化为颜色
table(net$colors)
mergedColors = net$colors
table(mergedColors)


# 层次聚类树
pdf("dendrogram_with_module_colors.pdf", width=10, height=8)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
## assign all of the gene to their corresponding module 
## hclust for the genes.


#明确样本数和基因数
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
#首先针对样本做个系统聚类树
pdf("datExpr_tree.pdf")
datExpr_tree<-hclust(dist(datExpr), method = "average")
par(mar = c(0,5,2,0))
plot(datExpr_tree, main = "Sample clustering", sub="", xlab="", cex.lab = 2, 
     cex.axis = 1, cex.main = 1,cex.lab=1)
dev.off()

# 生成相关性矩阵
# 首先要把分类做成因子型
head(datTraits)
table(datTraits$class)
new_levels <- c("rumen", "esophagus", "reticulum", "omasum", "abomasum")
datTraits$class <- factor(datTraits$class, levels = new_levels)
head(datTraits$class)
# design实际是提供一个比较方案
  design=model.matrix(~0+ datTraits$class)
  colnames(design)=levels(datTraits$class)
  moduleColors = net$colors 
  head(design)
  # Recalculate MEs with color labels
  MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
  MEs = orderMEs(MEs0); ##不同颜色的模块的ME值矩 (样本vs模块)
  moduleTraitCor = cor(MEs, design , use = "p");
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
# Will display correlations and their p-values
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)


# 创建模块-特性关联热图
png("c-Module-trait-relationships.png", width = 2400, height = 1600, res = 120)
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = colnames(design),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = TRUE,  # 尝试设置为 TRUE
  cex.text = 0.7,
  cex.lab = 1,
  yColorWidth = 0.01,
  xColorOffset = 0.01,
  font.lab.y = 1,
  yLabelsPosition = "left",
  zlim = c(-1, 1),
  main = "Module-trait relationships"
)
dev.off()

#提取出所有模块对应的基因
as.data.frame(table(moduleColors))
colnames(datExpr)
gene<-colnames(datExpr)   #转换处理后的datExpr0数据
genemodule<-cbind(gene,moduleColors)
write.csv(genemodule,"genemodule.csv")

#=======================各组织的高表达模块============================
# 定义输出目录（可选）
output_dir <- "module_genes_by_tissue"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# 遍历每个组织
for(tissue in colnames(moduleTraitCor)){
  # 1. 找到该组织相关性最高的模块（正相关优先）
  max_cor <- max(moduleTraitCor[, tissue])
  max_index <- which.max(moduleTraitCor[, tissue])
  module_name <- rownames(moduleTraitCor)[max_index]
  
  # 2. 验证模块名称格式（移除ME前缀）
  clean_module <- gsub("^ME", "", module_name)
  
  # 3. 提取模块对应的基因
  module_genes <- colnames(datExpr)[moduleColors == clean_module]
  
  # 4. 生成文件名（组织名+模块名）
  safe_tissue <- gsub("[^[:alnum:]]", "_", tissue)
  filename <- paste0(safe_tissue, "_", clean_module, "_genes.txt")
  filepath <- file.path(output_dir, filename)
  
  # 5. 保存基因列表
  if(length(module_genes) > 0){
    writeLines(module_genes, filepath)
    message(paste("Saved", length(module_genes), "genes for", tissue, "in module", clean_module))
  } else {
    warning(paste("No genes found for", tissue, "module", clean_module))
  }
}

# 生成汇总表格
library(dplyr)
top_modules <- data.frame(
  Tissue = colnames(moduleTraitCor),
  Module = sapply(colnames(moduleTraitCor), function(x) {
    gsub("^ME", "", rownames(moduleTraitCor)[which.max(moduleTraitCor[, x])])
  }),
  Correlation = apply(moduleTraitCor, 2, max),
  GeneCount = sapply(colnames(moduleTraitCor), function(x) {
    sum(moduleColors == gsub("^ME", "", rownames(moduleTraitCor)[which.max(moduleTraitCor[, x])]))
  })
) %>% arrange(desc(Correlation))

write.csv(top_modules, "top_modules_summary.csv", row.names = FALSE)