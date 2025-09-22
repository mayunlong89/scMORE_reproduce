
25-07-18

grn_pbmc8900 <- readRDS("grn_outputs_PBMC8900.rds")


grn_pbmc8900$grn

write.csv(grn_pbmc8900$grn, file="grn_outputs_PBMC8900_regulons.csv",quote=F, row.names = F)




# 假设你已经把 GWAS overlap 的表格读取为矩阵
library(corrplot)

# 构建 GWAS overlap 数量矩阵
overlap_mat <- matrix(c(
  0, 3, 8, 6, 26,
  3, 0, 9, 6, 1,
  8, 9, 0, 13, 4,
  6, 6, 13, 0, 3,
  26,1, 4, 3, 0
), nrow = 5, byrow = TRUE)

colnames(overlap_mat) <- rownames(overlap_mat) <- c("IBD", "PBC", "RA", "SLE", "UC")

# 可视化
corrplot(overlap_mat, is.corr = FALSE, method = "circle",
         type = "lower", tl.col = "black", tl.srt = 45,
         col = colorRampPalette(c("white", "red"))(200),
         number.cex = 0.8, addCoef.col = "black")





# 构建significant regulons数量矩阵
overlap_mat <- matrix(c(
  0, 21, 43, 37, 61,
  21, 0, 21, 23, 21,
  43, 21, 0, 39, 50,
  37, 23, 39, 0, 41,
  61,21, 50, 41, 0
), nrow = 5, byrow = TRUE)

colnames(overlap_mat) <- rownames(overlap_mat) <- c("IBD", "PBC", "RA", "SLE", "UC")

# 可视化
corrplot(overlap_mat, is.corr = FALSE, method = "circle",
         type = "lower", tl.col = "black", tl.srt = 45,
         col = colorRampPalette(c("white", "red"))(200),
         number.cex = 0.8, addCoef.col = "black")








# 构建OpenTarget-prioritized GWAS associations数量矩阵
overlap_mat <- matrix(c(
  0, 70, 507, 252, 863,
  70, 0, 69, 66, 46,
  507, 69, 0, 285, 407,
  252, 66, 285, 0, 178,
  863,46, 407, 178, 0
), nrow = 5, byrow = TRUE)

colnames(overlap_mat) <- rownames(overlap_mat) <- c("IBD", "PBC", "RA", "SLE", "UC")

# 可视化
corrplot(overlap_mat, is.corr = FALSE, method = "circle",
         type = "lower", tl.col = "black", tl.srt = 45,
         col = colorRampPalette(c("white", "red"))(200),
         number.cex = 0.8, addCoef.col = "black")





