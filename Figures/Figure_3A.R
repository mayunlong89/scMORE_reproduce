#---Signac method
#theta--visualization

setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")


glm_glmnet_xgb_data <- read.table("05_glm_glmnet_xgb_performance.txt",header = T)

head(glm_glmnet_xgb_data)

# Generate the violin plot
ggplot(glm_glmnet_xgb_data, aes(x = method, y = Escore, fill = specificity_method)) +
  geom_violin(trim = FALSE, position = position_dodge(0.8)) + # Adjust position for side-by-side violins
  scale_fill_manual(values = c("blue", "orange")) + # Match colors
  labs(x = "", y = "E score", fill = "method") +
  theme_minimal() +
  theme(legend.position = "right") # Place legend on the right


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
