
##--------------Figure 5e-----Regulon size
# 加载必要的库
library(ggplot2)
library(dplyr)

#统计regulon sizes-------
data2 <- grn_outputs$grn

# Group by TF and count the number of Targets, adding 1
result <- data2 %>%
  group_by(TF) %>%
  summarise(Target_Count = n() + 1)

write.csv(result, file="regulon_sizes_single_cellBrain_PD.csv",quote=F, row.names = F)

regulon_data <-  result
colnames(regulon_data) <- c("Regulon","Regulon_size")

# 绘制 density plot
ggplot(regulon_data, aes(x = Regulon_size)) +
  geom_density(fill = "blue", alpha = 0.4) +  # 绘制密度图
  geom_vline(xintercept = 20, color = "red", linetype = "dashed") +  # 标注 >= 20 的分界线
  geom_text(
    data = regulon_data %>% filter(Regulon_size >= 20),  # 筛选 >= 20 的数据点
    aes(x = Regulon_size, y = 0, label = Regulon),  # 添加文本标注
    color = "red",
    hjust = -0.2,
    vjust = 1.5
  ) +
  labs(
    title = "Density Plot of Regulon Size",
    x = "Regulon Size",
    y = "Density"
  ) +
  theme_minimal()



