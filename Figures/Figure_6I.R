

##----------Figure 6I
setwd("~/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/06-PD/")


###--------------------------------AS cells-----------
PD_AS <- read.table("01_PD_AS_regulons.txt",header = TRUE,sep = "\t")
#PD_AS$x <- -log10(PD_AS$SpecificityScore_p)
#PD_AS$y <- -log10(PD_AS$ImportanceWeightScore_p)
PD_AS$x <- -log10(PD_AS$SpecificityScore_p + runif(nrow(PD_AS), min = 1e-4, max = 1e-3))
PD_AS$y <- -log10(PD_AS$ImportanceWeightScore_p + runif(nrow(PD_AS), min = 1e-4, max = 1e-3))

threshold <- -log10(0.05)

# Highlight subset (e.g., x > 2 and y > 2)
PD_AS$highlight <- ifelse(PD_AS$x >= threshold & PD_AS$y >= threshold, "Significant", "Nonsignificant")

# Plot using ggplot2
library(ggplot2)
# 修改highlight颜色
ggplot(PD_AS, aes(x = x, y = y)) +
  # Scatter plot points
  geom_point(aes(color = highlight), size = 2, alpha = 0.7) +
  # Ellipse for entire data
  stat_ellipse(aes(x = x, y = y), level = 0.95, color = "lightblue", size = 0.5) +
  # Ellipse for highlighted points
  stat_ellipse(data = subset(PD_AS, highlight == "Significant"), aes(x = x, y = y),
               level = 0.95, color = "red", size = 0.5) +
  # Dashed vertical and horizontal lines
  geom_vline(xintercept = threshold, linetype = "dashed", color = "blue", size = 1) +
  geom_hline(yintercept = threshold, linetype = "dashed", color = "blue", size = 1) +
  # Title and labels
  labs(
    title = "PD_AS",
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



###-------------------------------OPC cells-----------
PD_OPC <- read.table("01_PD_OPC_regulons.txt",header = TRUE,sep = "\t")
#PD_OPC$x <- -log10(PD_OPC$SpecificityScore_p)
#PD_OPC$y <- -log10(PD_OPC$ImportanceWeightScore_p)
set.seed(123)  # 设置随机种子以确保结果可复现
PD_OPC$x <- -log10(PD_OPC$SpecificityScore_p + runif(nrow(PD_OPC), min = 1e-4, max = 1e-3))
PD_OPC$y <- -log10(PD_OPC$ImportanceWeightScore_p + runif(nrow(PD_OPC), min = 1e-4, max = 1e-3))

threshold <- -log10(0.05)

# Highlight subset (e.g., x > 2 and y > 2)
PD_OPC$highlight <- ifelse(PD_OPC$x >= threshold & PD_OPC$y >= threshold, "Significant", "Nonsignificant")

# Plot using ggplot2
library(ggplot2)
# 修改highlight颜色
ggplot(PD_OPC, aes(x = x, y = y)) +
  # Scatter plot points
  geom_point(aes(color = highlight), size = 2, alpha = 0.7) +
  # Ellipse for entire data
  stat_ellipse(aes(x = x, y = y), level = 0.95, color = "lightblue", size = 0.5) +
  # Ellipse for highlighted points
  stat_ellipse(data = subset(PD_OPC, highlight == "Significant"), aes(x = x, y = y),
               level = 0.95, color = "red", size = 0.5) +
  # Dashed vertical and horizontal lines
  geom_vline(xintercept = threshold, linetype = "dashed", color = "blue", size = 1) +
  geom_hline(yintercept = threshold, linetype = "dashed", color = "blue", size = 1) +
  # Title and labels
  labs(
    title = "PD_OPC",
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




