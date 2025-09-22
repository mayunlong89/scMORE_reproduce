
##----Figure C



##----OpenTargets--five autoimmune diseases

library(pheatmap)
library(RColorBrewer)

# 手动输入数据（或从 CSV 文件读取）
# 读取并准备数据
df <- read.csv("../01_opentarget_genes_heatmap.csv", header = TRUE, stringsAsFactors = FALSE)

# 设置行名为 Genes 列
rownames(df) <- df$Genes
df$Genes <- NULL

# 分离数值矩阵和 Regulon_TFs 注释列
# 重新排列疾病列为：IBD, UC, SLE, RA, PBC
expr_matrix <- df[, c("lymp_count", "lymp_percent")]
tf_anno <- data.frame(Regulon_TFs = df$Regulon_TFs)
rownames(tf_anno) <- rownames(df)

# 转换 Regulon_TFs 为因子型以便设定颜色
tf_anno$Regulon_TFs <- factor(tf_anno$Regulon_TFs, levels = c("Yes", "No"))


library(pheatmap)

# 自定义颜色（5段）
my_colors <- c("#FFFFFF", "#D0E1F2", "#91C8E4", "#3690C0", "#034E7B")  # 白→浅蓝→深蓝
my_colors <- c("#FFFFFF", "#E0ECF8", "#C2DFF3", "#91C8E4", "#56A3D9", "#2C7BB6", "#034E7B")
my_colors <- c("#FFFFFF",  "#3690C0")

# 设置断点：长度必须 = length(colors) + 1
my_breaks <- c(-0.1, 0.5, max(expr_matrix, na.rm = TRUE) + 1)


# 绘图
pheatmap(expr_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = my_colors,
         breaks = my_breaks,
         annotation_row = tf_anno,
         annotation_colors = list(Regulon_TFs = c(Yes = "#f38181", No = "#fae3d9")),
         fontsize_row = 6,
         fontsize_col = 10,
         main = "Gene–Disease Heatmap with 5 Color Breaks")

