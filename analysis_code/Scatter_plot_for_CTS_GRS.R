

###--------------------------------IBD_CD8+T--------
PD_AS <-  read.csv("~/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2025-05-20-scMORE_revision/05_parameter_assessment/results_IBD_scMore_trait_results_phastCons_withoutExons.csv")

# 读取数据
PD_AS <- results_IBD_phastCons_withoutExons_alpha0$scMore_trait_results

#PD_AS <-  PD_AS %>% filter(Celltype=="CD8+ T cells")
# 使用原始分数作为横纵坐标
PD_AS$x <- PD_AS$CTS
PD_AS$y <- PD_AS$GRS

# p-value 显著性阈值
threshold_p <- 0.05

# 标记是否显著（仍基于 p-value）
PD_AS$highlight <- ifelse(PD_AS$CTS_P_value < threshold_p & PD_AS$GRS_P_value < threshold_p,
                          "Significant", "Nonsignificant")

# 画图
library(ggplot2)

ggplot(PD_AS, aes(x = x, y = y)) +
  # 所有点
  geom_point(aes(color = highlight), size = 2, alpha = 0.7) +

  # 整体分布的椭圆
  stat_ellipse(level = 0.95, color = "lightblue", size = 0.5) +

  # 显著点的椭圆
  stat_ellipse(data = subset(PD_AS, highlight == "Significant"), level = 0.95,
               color = "red", size = 0.5) +

  # 虚线参考线（此时基于原始得分，可以设置或省略）
  # geom_vline(xintercept = some_value)
  # geom_hline(yintercept = some_value)

  # 图例和标签
  labs(
    title = "IBD_CD8+T Regulons",
    x = "CTS score",
    y = "GRS score",
    color = "Group"
  ) +

  # 自定义颜色
  scale_color_manual(values = c("Significant" = "red", "Nonsignificant" = "gray")) +

  # 美化主题
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "top"
  )

