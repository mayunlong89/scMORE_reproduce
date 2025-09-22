#0---S-MultiXcan-results
#0---S-MultiXcan-results
#0---S-MultiXcan-results
df <- read.csv("02_smultixcan_genes_heatmap.csv", header = TRUE, stringsAsFactors = FALSE)


# 设置行名为 Genes 列
rownames(df) <- df$Genes
df$Genes <- NULL


# 提取表达矩阵列（确保存在）
expr_matrix <- df[, c("ASD","SCZ","AN", "ADHD", "BIP", "MDD", "OCD","TS")]
tf_anno <- data.frame(Regulon_TFs = df$Regulon_TFs)
rownames(tf_anno) <- rownames(df)
tf_anno$Regulon_TFs <- factor(tf_anno$Regulon_TFs, levels = c("Yes", "No"))

# 避免列名异常
colnames(tf_anno) <- trimws(colnames(tf_anno))

print(colnames(tf_anno))  # 必须是： "Regulon_TFs"


# 断点适配 p-value 范围（假设 log10转化后用于可视化）
#my_colors <- c("#FFFFFF", "#E0ECF8", "#C2DFF3", "#91C8E4", "#56A3D9", "#2C7BB6", "#034E7B")
#my_colors <- colorRampPalette(c("#7A6090", "#A794C1", "#E7D5E5","#E9AEC1","#CE7692"))(100)
#my_colors <- colorRampPalette(c("#FFFFFF", "#E7D5E5","#A794C1","#E9AEC1","#CE7692"))(100)
my_colors <- colorRampPalette(c("#FFFFFF", "#E7D5E5","#E9AEC1","#A794C1","#7A6090"))(100)

# 去除列名空格并强制设置正确列名
colnames(tf_anno) <- "Regulon_TFs"

# 注释颜色键
annotation_colors <- list(Regulon_TFs = c(Yes = "#f38181", No = "#fae3d9"))

# 建议转换为 -log10(p) 增强可视化
log_expr_matrix <- -log10(expr_matrix)
# 构建与 log_expr_matrix 同形状的矩阵
stars_matrix <- ifelse(log_expr_matrix >= 5, "**",
                       ifelse(log_expr_matrix >= 3, "*", ""))


pheatmap(log_expr_matrix,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         color = my_colors,
         annotation_row = tf_anno,
         annotation_colors = list(Regulon_TFs = c(Yes = "#f38181", No = "#fae3d9")),
         fontsize_row = 6,
         fontsize_col = 10,
         display_numbers = stars_matrix,  # ✅ 显示星号
         number_color = "black",          # 星号颜色
         main = "Gene-Disease Heatmap (-log10 P-values ≥ 5 shown as *)")





# 读取并准备数据
df <- read.csv("02_smultixcan_genes_heatmap.csv", header = TRUE, stringsAsFactors = FALSE)
rownames(df) <- df$Genes
df$Genes <- NULL

# 提取表达矩阵
expr_matrix <- df[, c("ASD", "SCZ", "AN", "ADHD", "BIP", "MDD", "OCD", "TS")]
log_expr_matrix <- -log10(expr_matrix)

# 可选：构造 TF 注释（目前转置后不使用）
tf_anno <- data.frame(Regulon_TFs = df$Regulon_TFs)
rownames(tf_anno) <- rownames(df)
tf_anno$Regulon_TFs <- factor(tf_anno$Regulon_TFs, levels = c("Yes", "No"))

# 构造星号标注矩阵
stars_matrix <- ifelse(log_expr_matrix >= 5, "**",
                       ifelse(log_expr_matrix >= 3, "*", ""))

# 转置矩阵（现在行为疾病，列为基因）
log_expr_matrix_T <- t(log_expr_matrix)
stars_matrix_T <- t(stars_matrix)

# 自定义渐变颜色
my_colors <- colorRampPalette(c("#FFFFFF", "#E7D5E5", "#E9AEC1", "#A794C1", "#7A6090"))(100)

# 绘图（不加 annotation_row，因为已经转置）
pheatmap(log_expr_matrix_T,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         color = my_colors,
         display_numbers = stars_matrix_T,
         number_color = "black",
         fontsize_row = 10,
         fontsize_col = 6,
         main = "Gene-Disease Heatmap (-log10 P-values ≥ 5 shown as **)")

