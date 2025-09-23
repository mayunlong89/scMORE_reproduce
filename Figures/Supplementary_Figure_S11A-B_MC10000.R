
#2025-06-27


##load packages
suppressPackageStartupMessages({
  library(Pando)
  library(Seurat)
  library(Signac)
  library(ArchR)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(tidyr)
  library(dplyr)
  library(dotgen)
  library(org.Hs.eg.db)
  library(COSG)
  library(fitdistrplus)
  library(AnnotationDbi)
  library(GenomicRanges)
  library(IRanges)
})

library(Seurat)

###required load R packages

#devtools::install_github("mayunlong89/scMORE")

library(scMORE)
library(GenomicRanges)
library(IRanges)
library(Seurat)
library(Signac)

##-----PD------------------------------
##-----PD------------------------------
##-----PD------------------------------
##-----PD------------------------------
##-----PD------------------------------
human_PD_brain_10x <- readRDS("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/06-PD/PD_anno.rds")

human_PD_brain_10x <- JoinLayers(human_PD_brain_10x)

length(unique(human_PD_brain_10x$cell_type))

single_cell_58949  <- human_PD_brain_10x
Idents(single_cell_58949) <- human_PD_brain_10x$cell_type


grn_outputsBrain58949_phastCons_withoutExons <- createRegulon(single_cell_58949, n_targets=5,
                                                              peak2gene_method = 'Signac',
                                                              infer_method = 'glm',
                                                              exclude_exon_regions = TRUE,
                                                              conserved_regions = phastConsElements20Mammals.UCSC.hg38)


#grn_outputsBrain58949_phastCons_withoutExons


gene_info_PD <- read.table("~/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/06-PD/PD_processed_magma_results.genes.out",header = TRUE)
#snp_info
snp_info_PD <- read.csv("~/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/06-PD/PD_maf0.01.txt",header=T,sep="\t",stringsAsFactors = FALSE)

gene_info<- gene_info_PD
snp_info <- snp_info_PD

message("Step 3: Retrieving gene-level association scores...")
geneRiskScores <- getGeneScore(gene_info)

##----Celltype-specificity for each genes
target_scores <- suppressWarnings(getSpecificity(single_cell_58949))


#1.phastCons_withoutExons
regulon2disease_results_alpha1_PD_perm10000 <- regulon2disease(
  grn_outputs = grn_outputsBrain58949_phastCons_withoutExons,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 10000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  infer_method = 'glm',
  top_n = 5
)

write.csv(regulon2disease_results_alpha1_PD_perm10000, file="results_PD_scMore_trait_results_phastCons_withoutExons_perm10000.csv",quote=F, row.names = F)





#1.phastCons_withoutExons
regulon2disease_results_alpha1_PD_perm1000 <- regulon2disease(
  grn_outputs = grn_outputsBrain58949_phastCons_withoutExons,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  infer_method = 'glm',
  top_n = 5
)


##----------Response Figure for 10000 permutations
setwd("~/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2025-05-20-scMORE_revision/10_permutation_10000times")


###_-----------------------PD------------------------------------
# 加载必要包
library(ggplot2)
library(readr)

# 读取数据（请根据实际文件格式修改路径）
df1 <- read.csv("results_PD_scMore_trait_results_perm10000_vs_perm1000.csv")  # 或者 read_excel("your_file.xlsx") 如果是 Excel 文件


#df <- df1 %>% filter(df1$Celltype=="AS")
#df <- df1 %>% filter(df1$Celltype=="OPC")
df <- df1

# 清除 P=0 情况
df <- df[df$TRS_P_value_1000 > 0 & df$TRS_P_value_10000 > 0, ]

# 计算 -log10(P)
df$log10_p_1000 <- -log10(df$TRS_P_value_1000)
df$log10_p_10000 <- -log10(df$TRS_P_value_10000)

# 定义显著性颜色
df$color_class <- "Nonsignificant"
df$color_class[df$Significance_1000 == "Significant" & df$Significance_10000 == "Significant"] <- "Both"
df$color_class[df$Significance_1000 == "Significant" & df$Significance_10000 != "Significant"] <- "Only_1000"
df$color_class[df$Significance_1000 != "Significant" & df$Significance_10000 == "Significant"] <- "Only_10000"

# 转换为因子以控制图例顺序
df$color_class <- factor(df$color_class, levels = c("Both", "Only_10000", "Only_1000", "Nonsignificant"))

# 自定义颜色
color_map <- c(
  "Both" = "red",
  "Only_10000" = "gold",
  "Only_1000" = "green",
  "Nonsignificant" = "grey70"
)

# 绘图
ggplot(df, aes(x = log10_p_1000, y = log10_p_10000, color = color_class)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = color_map, name = "Significance Status") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(
    x = "-log10(TRS_P_value_1000)",
    y = "-log10(TRS_P_value_10000)",
    title = "Comparison of -log10 P-values (1000 vs. 10000 permutations)"
  ) +
  theme_bw()

# Pearson 和 Spearman 相关性
cat("Pearson:", cor(df$log10_p_1000, df$log10_p_10000, method = "pearson"), "\n")
cat("Spearman:", cor(df$log10_p_1000, df$log10_p_10000, method = "spearman"), "\n")




###_----------------------IBD---------------------------------
# 加载必要包
library(ggplot2)
library(readr)

# 读取数据（请根据实际文件格式修改路径）
df <- read.csv("results_IBD_scMore_trait_results_phastCons_withoutExons_perm10000.csv")  # 或者 read_excel("your_file.xlsx") 如果是 Excel 文件


# 清除 P=0 情况
df <- df[df$TRS_P_value_1000 > 0 & df$TRS_P_value_10000 > 0, ]

# 计算 -log10(P)
df$log10_p_1000 <- -log10(df$TRS_P_value_1000)
df$log10_p_10000 <- -log10(df$TRS_P_value_10000)

# 定义显著性颜色
df$color_class <- "Nonsignificant"
df$color_class[df$Significance_1000 == "Significant" & df$Significance_10000 == "Significant"] <- "Both"
df$color_class[df$Significance_1000 == "Significant" & df$Significance_10000 != "Significant"] <- "Only_1000"
df$color_class[df$Significance_1000 != "Significant" & df$Significance_10000 == "Significant"] <- "Only_10000"

# 转换为因子以控制图例顺序
df$color_class <- factor(df$color_class, levels = c("Both", "Only_10000", "Only_1000", "Nonsignificant"))

# 自定义颜色
color_map <- c(
  "Both" = "red",
  "Only_10000" = "gold",
  "Only_1000" = "green",
  "Nonsignificant" = "grey70"
)

# 绘图
ggplot(df, aes(x = log10_p_1000, y = log10_p_10000, color = color_class)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = color_map, name = "Significance Status") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(
    x = "-log10(TRS_P_value_1000)",
    y = "-log10(TRS_P_value_10000)",
    title = "Comparison of -log10 P-values (1000 vs. 10000 permutations)"
  ) +
  theme_bw()

# Pearson 和 Spearman 相关性
cat("Pearson:", cor(df$log10_p_1000, df$log10_p_10000, method = "pearson"), "\n")
cat("Spearman:", cor(df$log10_p_1000, df$log10_p_10000, method = "spearman"), "\n")





###_----------------------PBC---------------------------------
# 加载必要包
library(ggplot2)
library(readr)

# 读取数据（请根据实际文件格式修改路径）
df <- read.csv("results_PBC_scMore_trait_results_phastCons_withExons_perm10000.csv") #或者 read_excel("your_file.xlsx") 如果是 Excel 文件


# 清除 P=0 情况
df <- df[df$TRS_P_value_1000 > 0 & df$TRS_P_value_10000 > 0, ]

# 计算 -log10(P)
df$log10_p_1000 <- -log10(df$TRS_P_value_1000)
df$log10_p_10000 <- -log10(df$TRS_P_value_10000)

# 定义显著性颜色
df$color_class <- "Nonsignificant"
df$color_class[df$Significance_1000 == "Significant" & df$Significance_10000 == "Significant"] <- "Both"
df$color_class[df$Significance_1000 == "Significant" & df$Significance_10000 != "Significant"] <- "Only_1000"
df$color_class[df$Significance_1000 != "Significant" & df$Significance_10000 == "Significant"] <- "Only_10000"

# 转换为因子以控制图例顺序
df$color_class <- factor(df$color_class, levels = c("Both", "Only_10000", "Only_1000", "Nonsignificant"))

# 自定义颜色
color_map <- c(
  "Both" = "red",
  "Only_10000" = "gold",
  "Only_1000" = "green",
  "Nonsignificant" = "grey70"
)

# 绘图
ggplot(df, aes(x = log10_p_1000, y = log10_p_10000, color = color_class)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = color_map, name = "Significance Status") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(
    x = "-log10(TRS_P_value_1000)",
    y = "-log10(TRS_P_value_10000)",
    title = "Comparison of -log10 P-values (1000 vs. 10000 permutations)"
  ) +
  theme_bw()

# Pearson 和 Spearman 相关性
cat("Pearson:", cor(df$log10_p_1000, df$log10_p_10000, method = "pearson"), "\n")
cat("Spearman:", cor(df$log10_p_1000, df$log10_p_10000, method = "spearman"), "\n")



###_----------------------RA---------------------------------
# 加载必要包
library(ggplot2)
library(readr)

# 读取数据（请根据实际文件格式修改路径）
df <- read.csv("results_RA_scMore_trait_results_phastCons_withoutExons_perm10000.csv") #或者 read_excel("your_file.xlsx") 如果是 Excel 文件


# 清除 P=0 情况
df <- df[df$TRS_P_value_1000 > 0 & df$TRS_P_value_10000 > 0, ]

# 计算 -log10(P)
df$log10_p_1000 <- -log10(df$TRS_P_value_1000)
df$log10_p_10000 <- -log10(df$TRS_P_value_10000)

# 定义显著性颜色
df$color_class <- "Nonsignificant"
df$color_class[df$Significance_1000 == "Significant" & df$Significance_10000 == "Significant"] <- "Both"
df$color_class[df$Significance_1000 == "Significant" & df$Significance_10000 != "Significant"] <- "Only_1000"
df$color_class[df$Significance_1000 != "Significant" & df$Significance_10000 == "Significant"] <- "Only_10000"

# 转换为因子以控制图例顺序
df$color_class <- factor(df$color_class, levels = c("Both", "Only_10000", "Only_1000", "Nonsignificant"))

# 自定义颜色
color_map <- c(
  "Both" = "red",
  "Only_10000" = "gold",
  "Only_1000" = "green",
  "Nonsignificant" = "grey70"
)

# 绘图
ggplot(df, aes(x = log10_p_1000, y = log10_p_10000, color = color_class)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = color_map, name = "Significance Status") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(
    x = "-log10(TRS_P_value_1000)",
    y = "-log10(TRS_P_value_10000)",
    title = "Comparison of -log10 P-values (1000 vs. 10000 permutations)"
  ) +
  theme_bw()

# Pearson 和 Spearman 相关性
cat("Pearson:", cor(df$log10_p_1000, df$log10_p_10000, method = "pearson"), "\n")
cat("Spearman:", cor(df$log10_p_1000, df$log10_p_10000, method = "spearman"), "\n")







###_----------------------SLE--------------------------------
# 加载必要包
library(ggplot2)
library(readr)

# 读取数据（请根据实际文件格式修改路径）
df <- read.csv("results_SLE_scMore_trait_results_phastCons_withoutExons_perm10000.csv") #或者 read_excel("your_file.xlsx") 如果是 Excel 文件


# 清除 P=0 情况
df <- df[df$TRS_P_value_1000 > 0 & df$TRS_P_value_10000 > 0, ]

# 计算 -log10(P)
df$log10_p_1000 <- -log10(df$TRS_P_value_1000)
df$log10_p_10000 <- -log10(df$TRS_P_value_10000)

# 定义显著性颜色
df$color_class <- "Nonsignificant"
df$color_class[df$Significance_1000 == "Significant" & df$Significance_10000 == "Significant"] <- "Both"
df$color_class[df$Significance_1000 == "Significant" & df$Significance_10000 != "Significant"] <- "Only_1000"
df$color_class[df$Significance_1000 != "Significant" & df$Significance_10000 == "Significant"] <- "Only_10000"

# 转换为因子以控制图例顺序
df$color_class <- factor(df$color_class, levels = c("Both", "Only_10000", "Only_1000", "Nonsignificant"))

# 自定义颜色
color_map <- c(
  "Both" = "red",
  "Only_10000" = "gold",
  "Only_1000" = "green",
  "Nonsignificant" = "grey70"
)

# 绘图
ggplot(df, aes(x = log10_p_1000, y = log10_p_10000, color = color_class)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = color_map, name = "Significance Status") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(
    x = "-log10(TRS_P_value_1000)",
    y = "-log10(TRS_P_value_10000)",
    title = "Comparison of -log10 P-values (1000 vs. 10000 permutations)"
  ) +
  theme_bw()

# Pearson 和 Spearman 相关性
cat("Pearson:", cor(df$log10_p_1000, df$log10_p_10000, method = "pearson"), "\n")
cat("Spearman:", cor(df$log10_p_1000, df$log10_p_10000, method = "spearman"), "\n")







###_----------------------UC--------------------------------
# 加载必要包
library(ggplot2)
library(readr)

# 读取数据（请根据实际文件格式修改路径）
df <- read.csv("results_UC_scMore_trait_results_phastCons_withoutExons_perm10000.csv") #或者 read_excel("your_file.xlsx") 如果是 Excel 文件


# 清除 P=0 情况
df <- df[df$TRS_P_value_1000 > 0 & df$TRS_P_value_10000 > 0, ]

# 计算 -log10(P)
df$log10_p_1000 <- -log10(df$TRS_P_value_1000)
df$log10_p_10000 <- -log10(df$TRS_P_value_10000)

# 定义显著性颜色
df$color_class <- "Nonsignificant"
df$color_class[df$Significance_1000 == "Significant" & df$Significance_10000 == "Significant"] <- "Both"
df$color_class[df$Significance_1000 == "Significant" & df$Significance_10000 != "Significant"] <- "Only_1000"
df$color_class[df$Significance_1000 != "Significant" & df$Significance_10000 == "Significant"] <- "Only_10000"

# 转换为因子以控制图例顺序
df$color_class <- factor(df$color_class, levels = c("Both", "Only_10000", "Only_1000", "Nonsignificant"))

# 自定义颜色
color_map <- c(
  "Both" = "red",
  "Only_10000" = "gold",
  "Only_1000" = "green",
  "Nonsignificant" = "grey70"
)

# 绘图
ggplot(df, aes(x = log10_p_1000, y = log10_p_10000, color = color_class)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = color_map, name = "Significance Status") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(
    x = "-log10(TRS_P_value_1000)",
    y = "-log10(TRS_P_value_10000)",
    title = "Comparison of -log10 P-values (1000 vs. 10000 permutations)"
  ) +
  theme_bw()

# Pearson 和 Spearman 相关性
cat("Pearson:", cor(df$log10_p_1000, df$log10_p_10000, method = "pearson"), "\n")
cat("Spearman:", cor(df$log10_p_1000, df$log10_p_10000, method = "spearman"), "\n")



###_-----boxplot
# 原始数据（秒）
df <- data.frame(
  Disease = c("IBD", "PBC", "RA", "SLE", "UC", "PD"),
  sim1K_sec = c(3958.16, 4590.41, 4249.51, 5510.58, 4084.54, 8408.13),
  sim10K_sec = c(77238.96, 85894.28, 99918.45, 74072.93, 100497.90, 88278.68)
)

# 转为分钟
df$sim1K_min <- df$sim1K_sec / 60
df$sim10K_min <- df$sim10K_sec / 60

# 计算 fold change
df$fold_change <- df$sim10K_min / df$sim1K_min

# 输出 fold change 表
print(df[, c("Disease", "sim1K_min", "sim10K_min", "fold_change")])

# 平均 fold
mean_fc <- mean(df$fold_change)
cat("Average fold increase in runtime (sim10K vs sim1K):", round(mean_fc, 2), "x\n")

# 如果你想绘图用 long 格式
library(tidyr)
library(ggplot2)
library(ggpubr)

df_long <- df |>
  dplyr::select(Disease, sim1K_min, sim10K_min) |>
  pivot_longer(cols = c("sim1K_min", "sim10K_min"), names_to = "Simulation", values_to = "Time_min")

df_long$Simulation <- factor(df_long$Simulation, levels = c("sim1K_min", "sim10K_min"))

# 画图
ggplot(df_long, aes(x = Simulation, y = Time_min)) +
  geom_boxplot(outlier.shape = NA, fill = "lightgray", color = "black") +
  geom_point(aes(group = Disease, color = Disease), size = 3) +
  geom_line(aes(group = Disease, color = Disease), size = 1) +
  stat_compare_means(method = "wilcox.test", paired = TRUE, label.y = max(df_long$Time_min) * 1.05) +
  labs(
    title = "Runtime Comparison (Minutes): sim1K vs sim10K",
    y = "Time (minutes)",
    x = "Simulation",
    subtitle = paste0("Average fold increase: ", round(mean_fc, 2), "x")
  ) +
  theme_classic()







##----------Response Figure for 100000 permutations
##----------Response Figure for 100000 permutations
##----------Response Figure for 100000 permutations
setwd("~/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2025-05-20-scMORE_revision/10_permutation_10000times")


###_-----------------------PD------------------------------------
# 加载必要包
library(ggplot2)
library(readr)

# 读取数据（请根据实际文件格式修改路径）
df1 <- read.csv("results_PD_scMore_trait_results_perm100000_vs_perm1000.csv")  # 或者 read_excel("your_file.xlsx") 如果是 Excel 文件


#df <- df1 %>% filter(df1$Celltype=="AS")
#df <- df1 %>% filter(df1$Celltype=="OPC")
df <- df1

# 清除 P=0 情况
df <- df[df$TRS_P_value_1000 > 0 & df$TRS_P_value_100000 > 0, ]

# 计算 -log10(P)
df$log10_p_1000 <- -log10(df$TRS_P_value_1000)
df$log10_p_100000 <- -log10(df$TRS_P_value_100000)

# 定义显著性颜色
df$color_class <- "Nonsignificant"
df$color_class[df$Significance_1000 == "Significant" & df$Significance_100000 == "Significant"] <- "Both"
df$color_class[df$Significance_1000 == "Significant" & df$Significance_100000 != "Significant"] <- "Only_1000"
df$color_class[df$Significance_1000 != "Significant" & df$Significance_100000 == "Significant"] <- "Only_100000"

# 转换为因子以控制图例顺序
df$color_class <- factor(df$color_class, levels = c("Both", "Only_10000", "Only_1000", "Nonsignificant"))

# 自定义颜色
color_map <- c(
  "Both" = "red",
  "Only_10000" = "gold",
  "Only_1000" = "green",
  "Nonsignificant" = "grey70"
)

# 绘图
ggplot(df, aes(x = log10_p_1000, y = log10_p_100000, color = color_class)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = color_map, name = "Significance Status") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(
    x = "-log10(TRS_P_value_1000)",
    y = "-log10(TRS_P_value_100000)",
    title = "Comparison of -log10 P-values (1000 vs. 100000 permutations)"
  ) +
  theme_bw()

# Pearson 和 Spearman 相关性
cat("Pearson:", cor(df$log10_p_1000, df$log10_p_100000, method = "pearson"), "\n")
cat("Spearman:", cor(df$log10_p_1000, df$log10_p_100000, method = "spearman"), "\n")




###_----------------------IBD---------------------------------
# 加载必要包
library(ggplot2)
library(readr)

# 读取数据（请根据实际文件格式修改路径）
df <- read.csv("results_IBD_scMore_trait_results_phastCons_withoutExons_perm100000.csv")  # 或者 read_excel("your_file.xlsx") 如果是 Excel 文件


# 清除 P=0 情况
df <- df[df$TRS_P_value_1000 > 0 & df$TRS_P_value_100000 > 0, ]

# 计算 -log10(P)
df$log10_p_1000 <- -log10(df$TRS_P_value_1000)
df$log10_p_100000 <- -log10(df$TRS_P_value_100000)

# 定义显著性颜色
df$color_class <- "Nonsignificant"
df$color_class[df$Significance_1000 == "Significant" & df$Significance_100000 == "Significant"] <- "Both"
df$color_class[df$Significance_1000 == "Significant" & df$Significance_100000 != "Significant"] <- "Only_1000"
df$color_class[df$Significance_1000 != "Significant" & df$Significance_100000 == "Significant"] <- "Only_100000"

# 转换为因子以控制图例顺序
df$color_class <- factor(df$color_class, levels = c("Both", "Only_100000", "Only_1000", "Nonsignificant"))

# 自定义颜色
color_map <- c(
  "Both" = "red",
  "Only_100000" = "gold",
  "Only_1000" = "green",
  "Nonsignificant" = "grey70"
)

# 绘图
ggplot(df, aes(x = log10_p_1000, y = log10_p_100000, color = color_class)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = color_map, name = "Significance Status") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(
    x = "-log10(TRS_P_value_1000)",
    y = "-log10(TRS_P_value_100000)",
    title = "Comparison of -log10 P-values (1000 vs. 100000 permutations)"
  ) +
  theme_bw()

# Pearson 和 Spearman 相关性
cat("Pearson:", cor(df$log10_p_1000, df$log10_p_100000, method = "pearson"), "\n")
cat("Spearman:", cor(df$log10_p_1000, df$log10_p_100000, method = "spearman"), "\n")





###_----------------------PBC---------------------------------
# 加载必要包
library(ggplot2)
library(readr)

# 读取数据（请根据实际文件格式修改路径）
df <- read.csv("results_PBC_scMore_trait_results_phastCons_withExons_perm100000.csv") #或者 read_excel("your_file.xlsx") 如果是 Excel 文件


# 清除 P=0 情况
df <- df[df$TRS_P_value_1000 > 0 & df$TRS_P_value_100000 > 0, ]

# 计算 -log10(P)
df$log10_p_1000 <- -log10(df$TRS_P_value_1000)
df$log10_p_100000 <- -log10(df$TRS_P_value_100000)

# 定义显著性颜色
df$color_class <- "Nonsignificant"
df$color_class[df$Significance_1000 == "Significant" & df$Significance_100000 == "Significant"] <- "Both"
df$color_class[df$Significance_1000 == "Significant" & df$Significance_100000 != "Significant"] <- "Only_1000"
df$color_class[df$Significance_1000 != "Significant" & df$Significance_100000 == "Significant"] <- "Only_100000"

# 转换为因子以控制图例顺序
df$color_class <- factor(df$color_class, levels = c("Both", "Only_100000", "Only_1000", "Nonsignificant"))

# 自定义颜色
color_map <- c(
  "Both" = "red",
  "Only_100000" = "gold",
  "Only_1000" = "green",
  "Nonsignificant" = "grey70"
)

# 绘图
ggplot(df, aes(x = log10_p_1000, y = log10_p_100000, color = color_class)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = color_map, name = "Significance Status") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(
    x = "-log10(TRS_P_value_1000)",
    y = "-log10(TRS_P_value_100000)",
    title = "Comparison of -log10 P-values (1000 vs. 100000 permutations)"
  ) +
  theme_bw()

# Pearson 和 Spearman 相关性
cat("Pearson:", cor(df$log10_p_1000, df$log10_p_100000, method = "pearson"), "\n")
cat("Spearman:", cor(df$log10_p_1000, df$log10_p_100000, method = "spearman"), "\n")



###_----------------------RA---------------------------------(Predicted)
# 加载必要包
library(ggplot2)
library(readr)

# 读取数据（请根据实际文件格式修改路径）
df <- read.csv("results_RA_scMore_trait_results_phastCons_withoutExons_perm100000.csv") #或者 read_excel("your_file.xlsx") 如果是 Excel 文件


# 清除 P=0 情况
df <- df[df$TRS_P_value_1000 > 0 & df$TRS_P_value_100000 > 0, ]

# 计算 -log10(P)
df$log10_p_1000 <- -log10(df$TRS_P_value_1000)
df$log10_p_100000 <- -log10(df$TRS_P_value_100000)

# 定义显著性颜色
df$color_class <- "Nonsignificant"
df$color_class[df$Significance_1000 == "Significant" & df$Significance_100000 == "Significant"] <- "Both"
df$color_class[df$Significance_1000 == "Significant" & df$Significance_100000 != "Significant"] <- "Only_1000"
df$color_class[df$Significance_1000 != "Significant" & df$Significance_100000 == "Significant"] <- "Only_100000"

# 转换为因子以控制图例顺序
df$color_class <- factor(df$color_class, levels = c("Both", "Only_100000", "Only_1000", "Nonsignificant"))

# 自定义颜色
color_map <- c(
  "Both" = "red",
  "Only_100000" = "gold",
  "Only_1000" = "green",
  "Nonsignificant" = "grey70"
)

# 绘图
ggplot(df, aes(x = log10_p_1000, y = log10_p_100000, color = color_class)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = color_map, name = "Significance Status") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(
    x = "-log10(TRS_P_value_1000)",
    y = "-log10(TRS_P_value_100000)",
    title = "Comparison of -log10 P-values (1000 vs. 100000 permutations)"
  ) +
  theme_bw()

# Pearson 和 Spearman 相关性
cat("Pearson:", cor(df$log10_p_1000, df$log10_p_100000, method = "pearson"), "\n")
cat("Spearman:", cor(df$log10_p_1000, df$log10_p_100000, method = "spearman"), "\n")







###_----------------------SLE--------------------------------
# 加载必要包
library(ggplot2)
library(readr)

# 读取数据（请根据实际文件格式修改路径）
df <- read.csv("results_SLE_scMore_trait_results_phastCons_withoutExons_perm100000.csv") #或者 read_excel("your_file.xlsx") 如果是 Excel 文件


# 清除 P=0 情况
df <- df[df$TRS_P_value_1000 > 0 & df$TRS_P_value_100000 > 0, ]

# 计算 -log10(P)
df$log10_p_1000 <- -log10(df$TRS_P_value_1000)
df$log10_p_100000 <- -log10(df$TRS_P_value_100000)

# 定义显著性颜色
df$color_class <- "Nonsignificant"
df$color_class[df$Significance_1000 == "Significant" & df$Significance_100000 == "Significant"] <- "Both"
df$color_class[df$Significance_1000 == "Significant" & df$Significance_100000 != "Significant"] <- "Only_1000"
df$color_class[df$Significance_1000 != "Significant" & df$Significance_100000 == "Significant"] <- "Only_100000"

# 转换为因子以控制图例顺序
df$color_class <- factor(df$color_class, levels = c("Both", "Only_100000", "Only_1000", "Nonsignificant"))

# 自定义颜色
color_map <- c(
  "Both" = "red",
  "Only_100000" = "gold",
  "Only_1000" = "green",
  "Nonsignificant" = "grey70"
)

# 绘图
ggplot(df, aes(x = log10_p_1000, y = log10_p_100000, color = color_class)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = color_map, name = "Significance Status") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(
    x = "-log10(TRS_P_value_1000)",
    y = "-log10(TRS_P_value_100000)",
    title = "Comparison of -log10 P-values (1000 vs. 100000 permutations)"
  ) +
  theme_bw()

# Pearson 和 Spearman 相关性
cat("Pearson:", cor(df$log10_p_1000, df$log10_p_100000, method = "pearson"), "\n")
cat("Spearman:", cor(df$log10_p_1000, df$log10_p_100000, method = "spearman"), "\n")







###_----------------------UC--------------------------------
# 加载必要包
library(ggplot2)
library(readr)

# 读取数据（请根据实际文件格式修改路径）
df <- read.csv("results_UC_scMore_trait_results_phastCons_withoutExons_perm100000.csv") #或者 read_excel("your_file.xlsx") 如果是 Excel 文件


# 清除 P=0 情况
df <- df[df$TRS_P_value_1000 > 0 & df$TRS_P_value_100000 > 0, ]

# 计算 -log10(P)
df$log10_p_1000 <- -log10(df$TRS_P_value_1000)
df$log10_p_100000 <- -log10(df$TRS_P_value_100000)

# 定义显著性颜色
df$color_class <- "Nonsignificant"
df$color_class[df$Significance_1000 == "Significant" & df$Significance_100000 == "Significant"] <- "Both"
df$color_class[df$Significance_1000 == "Significant" & df$Significance_100000 != "Significant"] <- "Only_1000"
df$color_class[df$Significance_1000 != "Significant" & df$Significance_100000 == "Significant"] <- "Only_100000"

# 转换为因子以控制图例顺序
df$color_class <- factor(df$color_class, levels = c("Both", "Only_100000", "Only_1000", "Nonsignificant"))

# 自定义颜色
color_map <- c(
  "Both" = "red",
  "Only_100000" = "gold",
  "Only_1000" = "green",
  "Nonsignificant" = "grey70"
)

# 绘图
ggplot(df, aes(x = log10_p_1000, y = log10_p_100000, color = color_class)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = color_map, name = "Significance Status") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(
    x = "-log10(TRS_P_value_1000)",
    y = "-log10(TRS_P_value_100000)",
    title = "Comparison of -log10 P-values (1000 vs. 100000 permutations)"
  ) +
  theme_bw()

# Pearson 和 Spearman 相关性
cat("Pearson:", cor(df$log10_p_1000, df$log10_p_100000, method = "pearson"), "\n")
cat("Spearman:", cor(df$log10_p_1000, df$log10_p_100000, method = "spearman"), "\n")






###_-----boxplot
# 原始数据（秒）
df <- data.frame(
  Disease = c("IBD", "PBC", "RA", "SLE", "UC", "PD"),
  sim1K_sec = c(3958.16, 4590.41, 4249.51, 5510.58, 4084.54, 8408.13),
  sim10K_sec = c(77238.96, 85894.28, 99918.45, 74072.93, 100497.90, 88278.68)
)

# 转为分钟
df$sim1K_min <- df$sim1K_sec / 60
df$sim10K_min <- df$sim10K_sec / 60

# 计算 fold change
df$fold_change <- df$sim10K_min / df$sim1K_min

# 输出 fold change 表
print(df[, c("Disease", "sim1K_min", "sim10K_min", "fold_change")])

# 平均 fold
mean_fc <- mean(df$fold_change)
cat("Average fold increase in runtime (sim10K vs sim1K):", round(mean_fc, 2), "x\n")

# 如果你想绘图用 long 格式
library(tidyr)
library(ggplot2)
library(ggpubr)

df_long <- df |>
  select(Disease, sim1K_min, sim10K_min) |>
  pivot_longer(cols = c("sim1K_min", "sim10K_min"), names_to = "Simulation", values_to = "Time_min")

df_long$Simulation <- factor(df_long$Simulation, levels = c("sim1K_min", "sim10K_min"))

# 画图
ggplot(df_long, aes(x = Simulation, y = Time_min)) +
  geom_boxplot(outlier.shape = NA, fill = "lightgray", color = "black") +
  geom_point(aes(group = Disease, color = Disease), size = 3) +
  geom_line(aes(group = Disease, color = Disease), size = 1) +
  stat_compare_means(method = "wilcox.test", paired = TRUE, label.y = max(df_long$Time_min) * 1.05) +
  labs(
    title = "Runtime Comparison (Minutes): sim1K vs sim10K",
    y = "Time (minutes)",
    x = "Simulation",
    subtitle = paste0("Average fold increase: ", round(mean_fc, 2), "x")
  ) +
  theme_classic()







###_-----boxplot
# 原始数据（秒）
df <- data.frame(
  Disease = c("IBD", "PBC", "RA", "SLE", "UC", "PD"),
  sim1K_sec = c(3958.16, 4590.41, 4249.51, 5510.58, 4084.54, 8408.13),
  sim100K_sec = c(802054.9, 860585.2, 1062377.5, 816250.4, 929903, 1261219.5)
)

# 转为分钟
df$sim1K_min <- df$sim1K_sec / 60
df$sim100K_min <- df$sim100K_sec / 60

# 计算 fold change
df$fold_change <- df$sim100K_min / df$sim1K_min

# 输出 fold change 表
print(df[, c("Disease", "sim1K_min", "sim100K_min", "fold_change")])

# 平均 fold
mean_fc <- mean(df$fold_change)
cat("Average fold increase in runtime (sim100K vs sim1K):", round(mean_fc, 2), "x\n")

# 如果你想绘图用 long 格式
library(tidyr)
library(ggplot2)
library(ggpubr)
library(dplyr)

df_long <- df |>
  dplyr::select(Disease, sim1K_min, sim100K_min) |>
  pivot_longer(cols = c("sim1K_min", "sim100K_min"), names_to = "Simulation", values_to = "Time_min")

df_long$Simulation <- factor(df_long$Simulation, levels = c("sim1K_min", "sim100K_min"))

# 画图
ggplot(df_long, aes(x = Simulation, y = Time_min)) +
  geom_boxplot(outlier.shape = NA, fill = "lightgray", color = "black") +
  geom_point(aes(group = Disease, color = Disease), size = 3) +
  geom_line(aes(group = Disease, color = Disease), size = 1) +
  stat_compare_means(method = "wilcox.test", paired = TRUE, label.y = max(df_long$Time_min) * 1.05) +
  labs(
    title = "Runtime Comparison (Minutes): sim1K vs sim100K",
    y = "Time (minutes)",
    x = "Simulation",
    subtitle = paste0("Average fold increase: ", round(mean_fc, 2), "x")
  ) +
  theme_classic()











###--------------------------------AS cells-----------
PD_AS <- read.table("01_PD_AS_regulons.txt",header = TRUE,sep = "\t")
#PD_AS$x <- -log10(PD_AS$SpecificityScore_p)
#PD_AS$y <- -log10(PD_AS$ImportanceWeightScore_p)
PD_AS$x <- -log10(PD_AS$CTS_P_value + runif(nrow(PD_AS), min = 1e-4, max = 1e-3))
PD_AS$y <- -log10(PD_AS$GRS_P_value + runif(nrow(PD_AS), min = 1e-4, max = 1e-3))

threshold <- -log10(0.05)

# Highlight subset (e.g., x > 2 and y > 2)
PD_AS$highlight <- ifelse(PD_AS$x >= threshold & PD_AS$y >= threshold, "Significant", "Nonsignificant")

# Plot using ggplot2
library(ggplot2)
# 修改highlight颜色
ggplot(PD_AS, aes(x = x, y = y)) +
  # Scatter plot points
  geom_point(aes(color = highlight), size = 2, alpha = 0.7) +
  # Ellipse for entire data
  stat_ellipse(aes(x = x, y = y), level = 0.95, color = "lightblue", size = 0.5) +
  # Ellipse for highlighted points
  stat_ellipse(data = subset(PD_AS, highlight == "Significant"), aes(x = x, y = y),
               level = 0.95, color = "red", size = 0.5) +
  # Dashed vertical and horizontal lines
  geom_vline(xintercept = threshold, linetype = "dashed", color = "blue", size = 1) +
  geom_hline(yintercept = threshold, linetype = "dashed", color = "blue", size = 1) +
  # Title and labels
  labs(
    title = "PD_AS",
    x = "CTS",
    y = "GRS",
    color = "Group"
  ) +
  # Custom colors for highlight
  scale_color_manual(values = c("Significant" = "red", "Non-significant" = "gray")) +
  # Theme adjustments
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "top"
  )




###-------------------------------OPC cells-----------
PD_OPC <- read.table("01_PD_OPC_regulons.txt",header = TRUE,sep = "\t")
#PD_OPC$x <- -log10(PD_OPC$SpecificityScore_p)
#PD_OPC$y <- -log10(PD_OPC$ImportanceWeightScore_p)
set.seed(123)  # 设置随机种子以确保结果可复现
PD_OPC$x <- -log10(PD_OPC$CTS_P_value + runif(nrow(PD_OPC), min = 1e-4, max = 1e-3))
PD_OPC$y <- -log10(PD_OPC$GRS_P_value + runif(nrow(PD_OPC), min = 1e-4, max = 1e-3))

threshold <- -log10(0.05)

# Highlight subset (e.g., x > 2 and y > 2)
PD_OPC$highlight <- ifelse(PD_OPC$x >= threshold & PD_OPC$y >= threshold, "Significant", "Nonsignificant")

# Plot using ggplot2
library(ggplot2)
# 修改highlight颜色
ggplot(PD_OPC, aes(x = x, y = y)) +
  # Scatter plot points
  geom_point(aes(color = highlight), size = 2, alpha = 0.7) +
  # Ellipse for entire data
  stat_ellipse(aes(x = x, y = y), level = 0.95, color = "lightblue", size = 0.5) +
  # Ellipse for highlighted points
  stat_ellipse(data = subset(PD_OPC, highlight == "Significant"), aes(x = x, y = y),
               level = 0.95, color = "red", size = 0.5) +
  # Dashed vertical and horizontal lines
  geom_vline(xintercept = threshold, linetype = "dashed", color = "blue", size = 1) +
  geom_hline(yintercept = threshold, linetype = "dashed", color = "blue", size = 1) +
  # Title and labels
  labs(
    title = "PD_OPC",
    x = "CTS",
    y = "GRS",
    color = "Group"
  ) +
  # Custom colors for highlight
  scale_color_manual(values = c("Significant" = "red", "Non-significant" = "gray")) +
  # Theme adjustments
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "top"
  )




###-------------------------------ODC cells-----------
PD_ODC <- read.table("01_PD_ODC_regulons.txt",header = TRUE,sep = "\t")
set.seed(123)  # 设置随机种子以确保结果可复现
PD_ODC$x <- -log10(PD_ODC$CTS_P_value + runif(nrow(PD_ODC), min = 1e-4, max = 1e-3))
PD_ODC$y <- -log10(PD_ODC$GRS_P_value + runif(nrow(PD_ODC), min = 1e-4, max = 1e-3))

threshold <- -log10(0.05)

# Highlight subset (e.g., x > 2 and y > 2)
PD_ODC$highlight <- ifelse(PD_ODC$x >= threshold & PD_ODC$y >= threshold, "Significant", "Nonsignificant")

# Plot using ggplot2
library(ggplot2)
# 修改highlight颜色
ggplot(PD_ODC, aes(x = x, y = y)) +
  # Scatter plot points
  geom_point(aes(color = highlight), size = 2, alpha = 0.7) +
  # Ellipse for entire data
  stat_ellipse(aes(x = x, y = y), level = 0.95, color = "lightblue", size = 0.5) +
  # Ellipse for highlighted points
  stat_ellipse(data = subset(PD_ODC, highlight == "Significant"), aes(x = x, y = y),
               level = 0.95, color = "red", size = 0.5) +
  # Dashed vertical and horizontal lines
  geom_vline(xintercept = threshold, linetype = "dashed", color = "blue", size = 1) +
  geom_hline(yintercept = threshold, linetype = "dashed", color = "blue", size = 1) +
  # Title and labels
  labs(
    title = "PD_ODC",
    x = "CTS",
    y = "GRS",
    color = "Group"
  ) +
  # Custom colors for highlight
  scale_color_manual(values = c("Significant" = "red", "Non-significant" = "gray")) +
  # Theme adjustments
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "top"
  )




###-------------------------------MG cells-----------
PD_MG <- read.table("01_PD_MG_regulons.txt",header = TRUE,sep = "\t")
set.seed(123)  # 设置随机种子以确保结果可复现
PD_MG$x <- -log10(PD_MG$CTS_P_value + runif(nrow(PD_MG), min = 1e-4, max = 1e-3))
PD_MG$y <- -log10(PD_MG$GRS_P_value + runif(nrow(PD_MG), min = 1e-4, max = 1e-3))

threshold <- -log10(0.05)

# Highlight subset (e.g., x > 2 and y > 2)
PD_MG$highlight <- ifelse(PD_MG$x >= threshold & PD_MG$y >= threshold, "Significant", "Nonsignificant")

# Plot using ggplot2
library(ggplot2)
# 修改highlight颜色
ggplot(PD_MG, aes(x = x, y = y)) +
  # Scatter plot points
  geom_point(aes(color = highlight), size = 2, alpha = 0.7) +
  # Ellipse for entire data
  stat_ellipse(aes(x = x, y = y), level = 0.95, color = "lightblue", size = 0.5) +
  # Ellipse for highlighted points
  stat_ellipse(data = subset(PD_MG, highlight == "Significant"), aes(x = x, y = y),
               level = 0.95, color = "red", size = 0.5) +
  # Dashed vertical and horizontal lines
  geom_vline(xintercept = threshold, linetype = "dashed", color = "blue", size = 1) +
  geom_hline(yintercept = threshold, linetype = "dashed", color = "blue", size = 1) +
  # Title and labels
  labs(
    title = "PD_MG",
    x = "CTS",
    y = "GRS",
    color = "Group"
  ) +
  # Custom colors for highlight
  scale_color_manual(values = c("Significant" = "red", "Non-significant" = "gray")) +
  # Theme adjustments
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "top"
  )




###-------------------------------T cells-----------
PD_T <- read.table("01_PD_T_regulons.txt",header = TRUE,sep = "\t")
set.seed(123)  # 设置随机种子以确保结果可复现
PD_T$x <- -log10(PD_T$CTS_P_value + runif(nrow(PD_T), min = 1e-4, max = 1e-3))
PD_T$y <- -log10(PD_T$GRS_P_value + runif(nrow(PD_T), min = 1e-4, max = 1e-3))

threshold <- -log10(0.05)

# Highlight subset (e.g., x > 2 and y > 2)
PD_T$highlight <- ifelse(PD_T$x >= threshold & PD_T$y >= threshold, "Significant", "Nonsignificant")

# Plot using ggplot2
library(ggplot2)
# 修改highlight颜色
ggplot(PD_T, aes(x = x, y = y)) +
  # Scatter plot points
  geom_point(aes(color = highlight), size = 2, alpha = 0.7) +
  # Ellipse for entire data
  stat_ellipse(aes(x = x, y = y), level = 0.95, color = "lightblue", size = 0.5) +
  # Ellipse for highlighted points
  stat_ellipse(data = subset(PD_T, highlight == "Significant"), aes(x = x, y = y),
               level = 0.95, color = "red", size = 0.5) +
  # Dashed vertical and horizontal lines
  geom_vline(xintercept = threshold, linetype = "dashed", color = "blue", size = 1) +
  geom_hline(yintercept = threshold, linetype = "dashed", color = "blue", size = 1) +
  # Title and labels
  labs(
    title = "PD_T",
    x = "CTS",
    y = "GRS",
    color = "Group"
  ) +
  # Custom colors for highlight
  scale_color_manual(values = c("Significant" = "red", "Non-significant" = "gray")) +
  # Theme adjustments
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "top"
  )






###-------------------------------EC cells-----------
PD_EC <- read.table("01_PD_EC_regulons.txt",header = TRUE,sep = "\t")
set.seed(123)  # 设置随机种子以确保结果可复现
PD_EC$x <- -log10(PD_EC$CTS_P_value + runif(nrow(PD_EC), min = 1e-4, max = 1e-3))
PD_EC$y <- -log10(PD_EC$GRS_P_value + runif(nrow(PD_EC), min = 1e-4, max = 1e-3))

threshold <- -log10(0.05)

# Highlight subset (e.g., x > 2 and y > 2)
PD_EC$highlight <- ifelse(PD_EC$x >= threshold & PD_EC$y >= threshold, "Significant", "Nonsignificant")

# Plot using ggplot2
library(ggplot2)
# 修改highlight颜色
ggplot(PD_EC, aes(x = x, y = y)) +
  # Scatter plot points
  geom_point(aes(color = highlight), size = 2, alpha = 0.7) +
  # Ellipse for entire data
  stat_ellipse(aes(x = x, y = y), level = 0.95, color = "lightblue", size = 0.5) +
  # Ellipse for highlighted points
  stat_ellipse(data = subset(PD_EC, highlight == "Significant"), aes(x = x, y = y),
               level = 0.95, color = "red", size = 0.5) +
  # Dashed vertical and horizontal lines
  geom_vline(xintercept = threshold, linetype = "dashed", color = "blue", size = 1) +
  geom_hline(yintercept = threshold, linetype = "dashed", color = "blue", size = 1) +
  # Title and labels
  labs(
    title = "PD_EC",
    x = "CTS",
    y = "GRS",
    color = "Group"
  ) +
  # Custom colors for highlight
  scale_color_manual(values = c("Significant" = "red", "Non-significant" = "gray")) +
  # Theme adjustments
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "top"
  )





###-------------------------------N cells-----------
PD_N <- read.table("01_PD_N_regulons.txt",header = TRUE,sep = "\t")
set.seed(123)  # 设置随机种子以确保结果可复现
PD_N$x <- -log10(PD_N$CTS_P_value + runif(nrow(PD_N), min = 1e-4, max = 1e-3))
PD_N$y <- -log10(PD_N$GRS_P_value + runif(nrow(PD_N), min = 1e-4, max = 1e-3))

threshold <- -log10(0.05)

# Highlight subset (e.g., x > 2 and y > 2)
PD_N$highlight <- ifelse(PD_N$x >= threshold & PD_N$y >= threshold, "Significant", "Nonsignificant")

# Plot using ggplot2
library(ggplot2)
# 修改highlight颜色
ggplot(PD_N, aes(x = x, y = y)) +
  # Scatter plot points
  geom_point(aes(color = highlight), size = 2, alpha = 0.7) +
  # Ellipse for entire data
  stat_ellipse(aes(x = x, y = y), level = 0.95, color = "lightblue", size = 0.5) +
  # Ellipse for highlighted points
  stat_ellipse(data = subset(PD_N, highlight == "Significant"), aes(x = x, y = y),
               level = 0.95, color = "red", size = 0.5) +
  # Dashed vertical and horizontal lines
  geom_vline(xintercept = threshold, linetype = "dashed", color = "blue", size = 1) +
  geom_hline(yintercept = threshold, linetype = "dashed", color = "blue", size = 1) +
  # Title and labels
  labs(
    title = "PD_N",
    x = "CTS",
    y = "GRS",
    color = "Group"
  ) +
  # Custom colors for highlight
  scale_color_manual(values = c("Significant" = "red", "Non-significant" = "gray")) +
  # Theme adjustments
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "top"
  )


