#---Signac method
#theta--visualization

library(ggplot2)
library(dplyr)

setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")


glm_glmnet_xgb_data <- read.table("05_glm_glmnet_xgb_performance.txt",header = T)

head(glm_glmnet_xgb_data)


# 自定义颜色：与上传图一致（浅米色 + 深砖红）
custom_colors <- c("average" = "#F5D9B0", "cosine" = "#9A4A32")

ggplot(glm_glmnet_xgb_data, aes(x = method, y = Escore, fill = specificity_method)) +
  geom_violin(trim = FALSE, position = position_dodge(0.8), color = "black") +
  geom_jitter(aes(color = specificity_method),  # 添加散点
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
              size = 1.5, alpha = 0.5) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = c("gray30","gray30")) +  # 散点颜色也匹配
  labs(
    x = "",
    y = "E statistics",
    fill = "Method",
    color = "Method"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
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
