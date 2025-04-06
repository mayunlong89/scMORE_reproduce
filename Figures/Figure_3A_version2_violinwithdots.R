#---Signac method
#theta--visualization

library(ggplot2)
library(dplyr)

setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")


glm_glmnet_xgb_data <- read.table("05_glm_glmnet_xgb_performance.txt",header = T)

head(glm_glmnet_xgb_data)



library(ggplot2)
library(dplyr)

# 自定义颜色
custom_colors <- c("average" = "#F5D9B0", "cosine" = "#9A4A32")

# 创建 interaction 列用于分组显示
glm_glmnet_xgb_data$group <- interaction(glm_glmnet_xgb_data$method, glm_glmnet_xgb_data$specificity_method)

# 画图
ggplot(glm_glmnet_xgb_data, aes(x = group, y = Escore, fill = specificity_method)) +
  geom_violin(trim = FALSE, color = "black") +
  geom_boxplot(width = 0.15, fill = "white", color = "black", outlier.shape = NA, alpha = 0.9) +
  geom_jitter(size = 1.5, alpha = 0.6, color = "gray30", width = 0.1) +
  scale_fill_manual(values = custom_colors) +
  labs(
    x = "",
    y = "E statistics",
    fill = "Method"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
    text = element_text(size = 12)
  )


# Perform Wilcoxon test and calculate OR for each group
results <- glm_glmnet_xgb_data %>%
  group_by(method) %>%
  summarise(
    wilcox_p_value = wilcox.test(Escore[specificity_method == "cosine"], Escore[specificity_method == "average"])$p.value,
    OR = mean(Escore[specificity_method == "cosine"]) / mean(Escore[specificity_method == "average"]),
    cosine_variance = sd(Escore[specificity_method == "cosine"]),
    average_variance = sd(Escore[specificity_method == "average"]),
    OR_variance = sd(Escore[specificity_method == "cosine"])/sd(Escore[specificity_method == "average"])
  )

# Print the results
print(results)


