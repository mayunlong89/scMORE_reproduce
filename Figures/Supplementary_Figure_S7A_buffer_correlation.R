
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



