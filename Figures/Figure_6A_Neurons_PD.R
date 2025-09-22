

###-------------------------------N cells-----------
PD_N <- read.table("01_PD_N_regulons.txt",header = TRUE,sep = "\t")
set.seed(123)  # 设置随机种子以确保结果可复现
PD_N$x <- -log10(PD_N$SpecificityScore_p + runif(nrow(PD_N), min = 1e-4, max = 1e-3))
PD_N$y <- -log10(PD_N$ImportanceWeightScore_p + runif(nrow(PD_N), min = 1e-4, max = 1e-3))

threshold <- -log10(0.05)

# Highlight subset (e.g., x > 2 and y > 2)
PD_N$highlight <- ifelse(PD_N$x >= threshold & PD_N$y >= threshold, "Significant", "Nonsignificant")

# Plot using ggplot2
library(ggplot2)
# 修改highlight颜色
ggplot(PD_N, aes(x = x, y = y)) +
  # Scatter plot points
  geom_point(aes(color = highlight), size = 2, alpha = 0.7) +
  # Ellipse for entire data
  stat_ellipse(aes(x = x, y = y), level = 0.95, color = "lightblue", size = 0.5) +
  # Ellipse for highlighted points
  stat_ellipse(data = subset(PD_N, highlight == "Significant"), aes(x = x, y = y),
               level = 0.95, color = "red", size = 0.5) +
  # Dashed vertical and horizontal lines
  geom_vline(xintercept = threshold, linetype = "dashed", color = "blue", size = 1) +
  geom_hline(yintercept = threshold, linetype = "dashed", color = "blue", size = 1) +
  # Title and labels
  labs(
    title = "PD_N",
    x = "CTS",
    y = "GRS",
    color = "Group"
  ) +
  # Custom colors for highlight
  scale_color_manual(values = c("Significant" = "red", "Non-significant" = "gray")) +
  # Theme adjustments
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "top"
  )


