

##correlation of lymp count and lymp percent
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
consine_average10traits_heatmap <- read.table("01_10blood_cell_traits_heatmap.txt",header = TRUE,sep = "\t")
cor(consine_average10traits_heatmap$lymp_count,consine_average10traits_heatmap$lymp_percent)
cor.test(consine_average10traits_heatmap$lymp_count,consine_average10traits_heatmap$lymp_percent)

data <- read.table("lymp_count_percent.txt",header = TRUE,sep = "\t")
head(data)
cor.test(data$lymp_count,data$lymp_percent)

# 加载 ggplot2
library(ggplot2)
# 绘制散点相关性图
ggplot(data, aes(x = lymp_count, y = lymp_percent,color=Celltype)) +
  geom_point(size = 2, alpha = 0.7) +  # 散点
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # 添加拟合线
  labs(title = "Correlation by Cell Type",
       x = "Lymphocyte count",
       y = "Lymphocyte percent") +
  theme_classic() +  # 经典主题
  facet_wrap(~ Celltype)  # 按细胞类型分面

# Perform correlation test for each cell type
cor_results <- data %>%
  dplyr::group_by(Celltype) %>%
  dplyr::summarise(
    correlation = cor(lymp_count, lymp_percent, use = "complete.obs"),  # Pearson correlation
    p_value = cor.test(lymp_count, lymp_percent)$p.value,  # Extract p-value
    conf_low = cor.test(lymp_count, lymp_percent)$conf.int[1],  # Lower bound of confidence interval
    conf_high = cor.test(lymp_count, lymp_percent)$conf.int[2]   # Upper bound of confidence interval
  )

print(cor_results)
