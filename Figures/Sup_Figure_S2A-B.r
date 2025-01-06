
###lymphocyte percent----dot plots
###--------------------------------CD8T cells-----------
data_CD8T <- read.table("lymp_percent_regulon_CD8T.txt",header = TRUE,sep = "\t")
data_CD8T$x <- -log10(data_CD8T$SpecificityScore_p)
data_CD8T$y <- -log10(data_CD8T$ImportanceWeightScore_p)
threshold <- -log10(0.05)

# Highlight subset (e.g., x > 2 and y > 2)
data_CD8T$highlight <- ifelse(data_CD8T$x >= threshold & data_CD8T$y >= threshold, "Significant", "Non-significant")

# Plot using ggplot2
library(ggplot2)
# 修改highlight颜色
ggplot(data_CD8T, aes(x = x, y = y)) +
  # Scatter plot points
  geom_point(aes(color = highlight), size = 2, alpha = 0.7) +
  # Ellipse for entire data
  stat_ellipse(aes(x = x, y = y), level = 0.95, color = "lightblue", size = 0.5) +
  # Ellipse for highlighted points
  stat_ellipse(data = subset(data_CD8T, highlight == "Significant"), aes(x = x, y = y),
               level = 0.95, color = "red", size = 0.5) +
  # Dashed vertical and horizontal lines
  geom_vline(xintercept = threshold, linetype = "dashed", color = "blue", size = 1) +
  geom_hline(yintercept = threshold, linetype = "dashed", color = "blue", size = 1) +
  # Title and labels
  labs(
    title = "CD8+T cells",
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


###--------------------------------Monocytes-----------
data_CD8T <- read.table("lymp_percent_regulon_monocyte.txt",header = TRUE,sep = "\t")
data_CD8T$x <- -log10(data_CD8T$SpecificityScore_p)
data_CD8T$y <- -log10(data_CD8T$ImportanceWeightScore_p)
threshold <- -log10(0.05)

# Highlight subset (e.g., x > 2 and y > 2)
data_CD8T$highlight <- ifelse(data_CD8T$x >= threshold & data_CD8T$y >= threshold, "Significant", "Non-significant")

# Plot using ggplot2
library(ggplot2)
# 修改highlight颜色
ggplot(data_CD8T, aes(x = x, y = y)) +
  # Scatter plot points
  geom_point(aes(color = highlight), size = 2, alpha = 0.7) +
  # Ellipse for entire data
  stat_ellipse(aes(x = x, y = y), level = 0.95, color = "lightblue", size = 0.5) +
  # Ellipse for highlighted points
  stat_ellipse(data = subset(data_CD8T, highlight == "Significant"), aes(x = x, y = y),
               level = 0.95, color = "orange", size = 1) +
  # Dashed vertical and horizontal lines
  geom_vline(xintercept = threshold, linetype = "dashed", color = "blue", size = 1) +
  geom_hline(yintercept = threshold, linetype = "dashed", color = "blue", size = 1) +
  # Title and labels
  labs(
    title = "Monocytes",
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





