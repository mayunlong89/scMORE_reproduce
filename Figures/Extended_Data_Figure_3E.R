
##----Figure E



##----S-MultiXcan--five autoimmune diseases
# 读取并准备数据
df <- read.csv("../02_smultixcan_genes_heatmap.csv", header = TRUE, stringsAsFactors = FALSE)
rownames(df) <- df$Genes
df$Genes <- NULL

# 提取表达矩阵
expr_matrix <- df[, c("lymp_count", "lymp_percent")]
log_expr_matrix <- -log10(expr_matrix)

# 可选：构造 TF 注释（目前转置后不使用）
tf_anno <- data.frame(Regulon_TFs = df$Regulon_TFs)
rownames(tf_anno) <- rownames(df)
tf_anno$Regulon_TFs <- factor(tf_anno$Regulon_TFs, levels = c("Yes", "No"))

# 构造星号标注矩阵
stars_matrix <- ifelse(log_expr_matrix >= 8, "***",
                       ifelse(log_expr_matrix >= 5, "**",
                              ifelse(log_expr_matrix >= 3, "*", "")))

# 自定义渐变颜色
my_colors <- colorRampPalette(c("#FFFFFF", "#E7D5E5", "#E9AEC1", "#A794C1", "#7A6090"))(100)

# 绘图（不加 annotation_row，因为已经转置）
pheatmap(log_expr_matrix,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         color = my_colors,
         annotation_row = tf_anno,
         display_numbers = stars_matrix,
         number_color = "black",
         annotation_colors = list(Regulon_TFs = c(Yes = "#f38181", No = "#fae3d9")),
         fontsize_row = 8,
         fontsize_col = 8,
         main = "Gene-Disease Heatmap (-log10 P-values ≥ 5 shown as **)")

