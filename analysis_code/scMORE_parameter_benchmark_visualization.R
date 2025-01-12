
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


theta_data <- read.table("03_scMORE_paramter_bechmark_theta_cosine_average.txt",header = T)

head(theta_data)

# Generate the violin plot
ggplot(theta_data, aes(x = theta, y = E_score, fill = method)) +
  geom_violin(trim = FALSE, position = position_dodge(0.8)) + # Adjust position for side-by-side violins
  scale_fill_manual(values = c("blue", "orange")) + # Match colors
  labs(x = "theta parameter", y = "E score", fill = "method") +
  theme_minimal() +
  theme(legend.position = "right") # Place legend on the right


# Perform Wilcoxon test and calculate OR for each group
results <- theta_data %>%
  group_by(theta) %>%
  summarise(
    wilcox_p_value = wilcox.test(E_score[method == "cosine"], E_score[method == "average"])$p.value,
    OR = mean(E_score[method == "cosine"]) / mean(E_score[method == "average"]),
    cosine_variance = sd(E_score[method == "cosine"]),
    average_variance = sd(E_score[method == "average"]),
    OR_variance = sd(E_score[method == "cosine"])/sd(E_score[method == "average"])
  )

# Print the results
print(results)



#buffer--visualization

buffer_data <- read.table("03_scMORE_paramter_bechmark_buffer_cosine_average.txt",header = T)

head(buffer_data)

# Generate the violin plot
ggplot(buffer_data, aes(x = buffer, y = E_score, fill = method)) +
  geom_violin(trim = FALSE, position = position_dodge(0.8)) + # Adjust position for side-by-side violins
  scale_fill_manual(values = c("blue", "orange")) + # Match colors
  labs(x = "buffer parameter", y = "E score", fill = "method") +
  theme_minimal() +
  theme(legend.position = "right") # Place legend on the right

# Perform Wilcoxon test and calculate OR for each group
results <- buffer_data %>%
  group_by(buffer) %>%
  summarise(
    wilcox_p_value = wilcox.test(E_score[method == "cosine"], E_score[method == "average"])$p.value,
    OR = mean(E_score[method == "cosine"]) / mean(E_score[method == "average"]),
    cosine_variance = sd(E_score[method == "cosine"]),
    average_variance = sd(E_score[method == "average"]),
  )

# Print the results
print(results)



#---GREAT method
#theta--visualization

setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")

theta_data_GREAT <- read.table("03_scMORE_paramter_bechmark_theta_cosine_average_GREAT.txt",header = T)

head(theta_data_GREAT)

# Generate the violin plot
ggplot(theta_data_GREAT, aes(x = theta, y = E_score, fill = method)) +
  geom_violin(trim = FALSE, position = position_dodge(0.8)) + # Adjust position for side-by-side violins
  scale_fill_manual(values = c("blue", "orange")) + # Match colors
  labs(x = "theta parameter", y = "E score", fill = "method") +
  theme_minimal() +
  theme(legend.position = "right") # Place legend on the right


# Perform Wilcoxon test and calculate OR for each group
results <- theta_data_GREAT %>%
  group_by(theta) %>%
  summarise(
    wilcox_p_value = wilcox.test(E_score[method == "cosine"], E_score[method == "average"])$p.value,
    OR = mean(E_score[method == "cosine"]) / mean(E_score[method == "average"])
  )

# Print the results
print(results)



#buffer--visualization

buffer_data_GREAT <- read.table("03_scMORE_paramter_bechmark_buffer_cosine_average_GREAT.txt",header = T)

head(buffer_data_GREAT)

# Generate the violin plot
ggplot(buffer_data_GREAT, aes(x = buffer, y = E_score, fill = method)) +
  geom_violin(trim = FALSE, position = position_dodge(0.8)) + # Adjust position for side-by-side violins
  scale_fill_manual(values = c("blue", "orange")) + # Match colors
  labs(x = "buffer parameter", y = "E statistics", fill = "method") +
  theme_minimal() +
  theme(legend.position = "right") # Place legend on the right

# Perform Wilcoxon test and calculate OR for each group
results <- buffer_data_GREAT %>%
  group_by(buffer) %>%
  summarise(
    wilcox_p_value = wilcox.test(E_score[method == "cosine"], E_score[method == "average"])$p.value,
    OR = mean(E_score[method == "cosine"]) / mean(E_score[method == "average"])
  )

# Print the results
print(results)




#---GREAT vs Signac method benchmark
#theta--visualization

setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")

theta_data_GREAT_Signac <- read.table("04_scMORE_Signac_CREAT_bechmark_theta_cosine_average_violin.txt",header = T)

head(theta_data_GREAT_Signac)

#cosine
theta_data_GREAT_Signac1 <- theta_data_GREAT_Signac[which(theta_data_GREAT_Signac$method=="cosine"),]
# Generate the violin plot
ggplot(theta_data_GREAT_Signac1, aes(x = traits, y = E_score, fill = methd_group)) +
  geom_violin(trim = FALSE, position = position_dodge(0.8)) + # Adjust position for side-by-side violins
  scale_fill_manual(values = c("blue", "orange","grey","red")) + # Match colors
  labs(x = "theta parameter", y = "E score", fill = "method") +
  theme_minimal() +
  theme(legend.position = "right") # Place legend on the right


# Perform Wilcoxon test and calculate OR for each group
results <- theta_data_GREAT %>%
  group_by(theta) %>%
  summarise(
    wilcox_p_value = wilcox.test(E_score[method == "cosine"], E_score[method == "average"])$p.value,
    OR = mean(E_score[method == "cosine"]) / mean(E_score[method == "average"])
  )

# Print the results
print(results)


#average
theta_data_GREAT_Signac1 <- theta_data_GREAT_Signac[which(theta_data_GREAT_Signac$method=="average"),]
# Generate the violin plot
ggplot(theta_data_GREAT_Signac1, aes(x = traits, y = E_score, fill = methd_group)) +
  geom_violin(trim = FALSE, position = position_dodge(0.8)) + # Adjust position for side-by-side violins
  scale_fill_manual(values = c("blue", "orange","grey","red")) + # Match colors
  labs(x = "theta parameter", y = "E score", fill = "method") +
  theme_minimal() +
  theme(legend.position = "right") # Place legend on the right


# Perform Wilcoxon test and calculate OR for each group
results <- theta_data_GREAT %>%
  group_by(theta) %>%
  summarise(
    wilcox_p_value = wilcox.test(E_score[method == "cosine"], E_score[method == "average"])$p.value,
    OR = mean(E_score[method == "cosine"]) / mean(E_score[method == "average"])
  )

# Print the results
print(results)


###---------correlation---plot
data <- read.table("04_scMORE_Signac_CREAT_bechmark_theta_cosine_average_violin_correlation.txt",header = TRUE,sep = "\t")
head(data)
cor.test(data$E_score_CREAT,data$E_score_signac)

# 加载 ggplot2
library(ggplot2)
# 绘制散点相关性图
ggplot(data, aes(x = E_score_CREAT, y = E_score_signac,color=method)) +
  geom_point(size = 2, alpha = 0.7) +  # 散点
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # 添加拟合线
  labs(title = "Correlation by Cell Type",
       x = "GREAT-based E statistics",
       y = "Signac-based E statistics") +
  theme_classic() +  # 经典主题
  facet_wrap(~ method)  # 按细胞类型分面

# Perform correlation test for each cell type
cor_results <- data %>%
  dplyr::group_by(method) %>%
  dplyr::summarise(
    correlation = cor(E_score_CREAT, E_score_signac, use = "complete.obs"),  # Pearson correlation
    p_value = cor.test(E_score_CREAT, E_score_signac)$p.value,  # Extract p-value
    conf_low = cor.test(E_score_CREAT, E_score_signac)$conf.int[1],  # Lower bound of confidence interval
    conf_high = cor.test(E_score_CREAT, E_score_signac)$conf.int[2]   # Upper bound of confidence interval
  )

print(cor_results)


# 加载 ggplot2
library(ggplot2)

# 绘制散点相关性图
ggplot(data, aes(x = E_score_CREAT, y = E_score_signac)) +
  geom_point(size = 2, alpha = 0.7, color = "blue") +  # 修改点的颜色为蓝色
  geom_smooth(method = "lm", se = TRUE, linetype = "dashed", color = "red") +  # 添加拟合线并显示95%置信区间
  labs(
    title = "Correlation by Cell Type",
    x = "GREAT-based E statistics",
    y = "Signac-based E statistics"
  ) +
  theme_classic()

# 计算相关性和95%置信区间
cor_results <- data %>%
  dplyr::summarise(
    correlation = cor(E_score_CREAT, E_score_signac, use = "complete.obs"),  # Pearson correlation
    p_value = cor.test(E_score_CREAT, E_score_signac)$p.value,  # 提取p值
    conf_low = cor.test(E_score_CREAT, E_score_signac)$conf.int[1],  # 置信区间下界
    conf_high = cor.test(E_score_CREAT, E_score_signac)$conf.int[2]  # 置信区间上界
  )

# 打印相关性结果
print(cor_results)



#buffer--visualization

buffer_data_GREAT_Signac <- read.table("04_scMORE_Signac_CREAT_bechmark_buffer_cosine_average_violin.txt",header = T)

head(buffer_data_GREAT_Signac)

# Generate the violin plot
ggplot(buffer_data_GREAT_Signac, aes(x = traits, y = E_score1, fill = methd_group)) +
  geom_violin(trim = FALSE, position = position_dodge(0.8)) + # Adjust position for side-by-side violins
  scale_fill_manual(values = c("blue", "orange","red", "gray")) + # Match colors
  labs(x = "buffer parameter", y = "E statistics", fill = "method") +
  theme_minimal() +
  theme(legend.position = "right") # Place legend on the right


###---------correlation---plot
data <- read.table("04_scMORE_Signac_CREAT_bechmark_buffer_cosine_average_violin_correlation.txt",header = TRUE,sep = "\t")
head(data)
cor.test(data$E_score_CREAT,data$E_score_signac)

# 加载 ggplot2
library(ggplot2)
# 绘制散点相关性图
ggplot(data, aes(x = E_score_CREAT, y = E_score_signac,color=method)) +
  geom_point(size = 2, alpha = 0.7) +  # 散点
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # 添加拟合线
  labs(title = "Correlation by Cell Type",
       x = "GREAT-based E statistics",
       y = "Signac-based E statistics") +
  theme_classic() +  # 经典主题
  facet_wrap(~ method)  # 按细胞类型分面

# Perform correlation test for each cell type
cor_results <- data %>%
  dplyr::group_by(method) %>%
  dplyr::summarise(
    correlation = cor(E_score_CREAT, E_score_signac, use = "complete.obs"),  # Pearson correlation
    p_value = cor.test(E_score_CREAT, E_score_signac)$p.value,  # Extract p-value
    conf_low = cor.test(E_score_CREAT, E_score_signac)$conf.int[1],  # Lower bound of confidence interval
    conf_high = cor.test(E_score_CREAT, E_score_signac)$conf.int[2]   # Upper bound of confidence interval
  )

print(cor_results)


# 加载 ggplot2
library(ggplot2)

# 绘制散点相关性图
ggplot(data, aes(x = E_score_CREAT, y = E_score_signac)) +
  geom_point(size = 2, alpha = 0.7, color = "blue") +  # 修改点的颜色为蓝色
  geom_smooth(method = "lm", se = TRUE, linetype = "dashed", color = "red") +  # 添加拟合线并显示95%置信区间
  labs(
    title = "Correlation by Cell Type",
    x = "GREAT-based E statistics",
    y = "Signac-based E statistics"
  ) +
  theme_classic()

# 计算相关性和95%置信区间
cor_results <- data %>%
  dplyr::summarise(
    correlation = cor(E_score_CREAT, E_score_signac, use = "complete.obs"),  # Pearson correlation
    p_value = cor.test(E_score_CREAT, E_score_signac)$p.value,  # 提取p值
    conf_low = cor.test(E_score_CREAT, E_score_signac)$conf.int[1],  # 置信区间下界
    conf_high = cor.test(E_score_CREAT, E_score_signac)$conf.int[2]  # 置信区间上界
  )

# 打印相关性结果
print(cor_results)






