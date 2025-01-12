
##------glm_glmnet_xgb data transformed for correlation-----------------------

library(reshape2)

consine_average10traits_heatmap_glm <- read.table("01_10blood_cell_traits_sig_heatmap.txt",header = TRUE,sep = "\t")

consine_average10traits_heatmap_xgb <- read.table("01_10blood_cell_traits_sig_heatmap_xgb.txt",header = TRUE,sep = "\t")

consine_average10traits_heatmap_glmnet <- read.table("01_10blood_cell_traits_sig_heatmap_glmnet.txt",header = TRUE,sep = "\t")


md_glm <- melt(consine_average10traits_heatmap_glm, id=c("regulons"))
md_xgb <- melt(consine_average10traits_heatmap_xgb, id=c("regulons"))
md_glmnet <- melt(consine_average10traits_heatmap_glmnet, id=c("regulons"))

data <- md_glm

colnames(data) <- c("regulons","traits","glm_TRS")
head(data)
data$xgb_TRS <- md_xgb$value
data$glmnet_TRS <- md_glmnet$value


cor.test(data$glm_TRS,data$xgb_TRS)
cor.test(data$glm_TRS,data$glmnet_TRS)
cor.test(data$glmnet_TRS,data$xgb_TRS)


#Figure 3C-1
ggplot(data, aes(x = glm_TRS, y = xgb_TRS,color="red")) +
  geom_point(size = 2, alpha = 0.7) +  # 散点
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed",color="blue") +  # 添加拟合线
  labs(title = "Correlation by Cell Type",
       x = "glm_TRS",
       y = "xgb_TRS") +
  theme_classic()

cor_results2 <- data %>%
  dplyr::summarise(
    correlation = cor(glm_TRS, xgb_TRS, use = "complete.obs"),  # Pearson correlation
    p_value = cor.test(glm_TRS, xgb_TRS)$p.value,  # Extract p-value
    conf_low = cor.test(glm_TRS, xgb_TRS)$conf.int[1],  # Lower bound of confidence interval
    conf_high = cor.test(glm_TRS, xgb_TRS)$conf.int[2]   # Upper bound of confidence interval
  )
print(cor_results2)



#Figure 3C-2
ggplot(data, aes(x = glm_TRS, y = glmnet_TRS,color="red")) +
  geom_point(size = 2, alpha = 0.7) +  # 散点
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed",color="blue") +  # 添加拟合线
  labs(title = "Correlation by Cell Type",
       x = "glm_TRS",
       y = "glmnet_TRS") +
  theme_classic()

cor_results3 <- data %>%
  dplyr::summarise(
    correlation = cor(glm_TRS, glmnet_TRS, use = "complete.obs"),  # Pearson correlation
    p_value = cor.test(glm_TRS, glmnet_TRS)$p.value,  # Extract p-value
    conf_low = cor.test(glm_TRS, glmnet_TRS)$conf.int[1],  # Lower bound of confidence interval
    conf_high = cor.test(glm_TRS, glmnet_TRS)$conf.int[2]   # Upper bound of confidence interval
  )
print(cor_results3)

