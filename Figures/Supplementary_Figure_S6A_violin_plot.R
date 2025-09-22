#Visualization of plot------------------------------------------------------------------
# Load necessary libraries
library(dplyr)
library(ggplot2)


cell2regulon_Number <- read.table("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/05_glm_glmnet_xgb_tissues_cell_number_regulons_count.txt",header = T)


head(cell2regulon_Number)

ggplot(cell2regulon_Number, aes(x = method, y = Regulon_counts, fill = method)) +
  geom_violin(trim = FALSE, position = position_dodge(0.8)) + # 小提琴图
  geom_point(aes(color = Tissues), position = position_dodge(0.8), size = 1.5, alpha = 0.8) + # 数据点，用Tissues区分颜色
  scale_fill_manual(values = c( "#BFEFFF", "#EEDFCC","#FFDAB9")) + # 小提琴填充色
  scale_color_manual(values = c("orange", "blue", "#668B8B", "purple", "brown")) + # 数据点颜色映射
    labs(x = "", y = "Regulon counts", fill = "Method", color = "Tissues") + # 添加图例标签
  theme_minimal() +
  theme(legend.position = "right") # 图例放置在右侧

