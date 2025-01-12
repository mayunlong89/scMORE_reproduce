# 2025-01-06
#visualization
###-------use 10 blood cell traits to assess the difference of different parameters
#' buffer = 0bp, 50bp, 100bp, 200bp, 500bp, and 1000bp
#' window size = 0kb, 5kb, 10kb, 20kb, 50kb, 100kb.
#' theta = 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1
#' peak2gene_method = Signac and GREAT
#' regression_method = glm, glmnet, cv.glmnet, brms

# Load necessary libraries
library(dplyr)
library(ggplot2)


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



# Generate the violin plot
ggplot(glm_glmnet_xgb_data, aes(x = method, y = Escore, fill = specificity_method)) +
  geom_violin(trim = FALSE, position = position_dodge(0.8)) + # Adjust position for side-by-side violins
  scale_fill_manual(values = c("blue", "orange")) + # Match colors
  labs(x = "", y = "E score", fill = "method") +
  theme_classic() +
  theme(legend.position = "right") # Place legend on the right




results2 <- glm_glmnet_xgb_data %>%
  group_by(specificity_method) %>%
  summarise(
    AC_p_value = wilcox.test(Escore[method == "A_glm"], Escore[method == "C_glmnet"])$p.value,
    AC_OR = mean(Escore[method == "A_glm"]) / mean(Escore[method == "C_glmnet"]),
    AB_p_value = wilcox.test(Escore[method == "A_glm"], Escore[method == "B_xgb"])$p.value,
    AB_OR = mean(Escore[method == "A_glm"]) / mean(Escore[method == "B_xgb"]),
    BC_p_value = wilcox.test(Escore[method == "B_xgb"], Escore[method == "C_glmnet"])$p.value,
    BC_OR = mean(Escore[method == "B_xgb"]) / mean(Escore[method == "C_glmnet"])
  )
# Print the results
print(results2)





