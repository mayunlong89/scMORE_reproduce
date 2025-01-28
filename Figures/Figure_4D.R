##----------Figure 4D
setwd("~/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/08_five_immunediseases/")


###--------------------------------CD8T cells-----------
IBD_CD8T <- read.table("IBD_regulons_CD8T.txt",header = TRUE,sep = "\t")
IBD_CD8T$x <- -log10(IBD_CD8T$SpecificityScore_p)
IBD_CD8T$y <- -log10(IBD_CD8T$ImportanceWeightScore_p)
threshold <- -log10(0.05)

# Highlight subset (e.g., x > 2 and y > 2)
IBD_CD8T$highlight <- ifelse(IBD_CD8T$x >= threshold & IBD_CD8T$y >= threshold, "Significant", "Nonsignificant")

# Plot using ggplot2
library(ggplot2)
# 修改highlight颜色
ggplot(IBD_CD8T, aes(x = x, y = y)) +
  # Scatter plot points
  geom_point(aes(color = highlight), size = 2, alpha = 0.7) +
  # Ellipse for entire data
  stat_ellipse(aes(x = x, y = y), level = 0.95, color = "lightblue", size = 0.5) +
  # Ellipse for highlighted points
  stat_ellipse(data = subset(IBD_CD8T, highlight == "Significant"), aes(x = x, y = y),
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



###--------------------------------CD4T cells-----------
IBD_CD8T <- read.table("RA_regulons_CD4T.txt",header = TRUE,sep = "\t")
IBD_CD8T$x <- -log10(IBD_CD8T$SpecificityScore_p)
IBD_CD8T$y <- -log10(IBD_CD8T$ImportanceWeightScore_p)
threshold <- -log10(0.05)

# Highlight subset (e.g., x > 2 and y > 2)
IBD_CD8T$highlight <- ifelse(IBD_CD8T$x >= threshold & IBD_CD8T$y >= threshold, "Significant", "Nonsignificant")

# Plot using ggplot2
library(ggplot2)
# 修改highlight颜色
ggplot(IBD_CD8T, aes(x = x, y = y)) +
  # Scatter plot points
  geom_point(aes(color = highlight), size = 2, alpha = 0.7) +
  # Ellipse for entire data
  stat_ellipse(aes(x = x, y = y), level = 0.95, color = "lightblue", size = 0.5) +
  # Ellipse for highlighted points
  stat_ellipse(data = subset(IBD_CD8T, highlight == "Significant"), aes(x = x, y = y),
               level = 0.95, color = "red", size = 0.5) +
  # Dashed vertical and horizontal lines
  geom_vline(xintercept = threshold, linetype = "dashed", color = "blue", size = 1) +
  geom_hline(yintercept = threshold, linetype = "dashed", color = "blue", size = 1) +
  # Title and labels
  labs(
    title = "CD4+T cells",
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



###--------------------------------SLE monocyte cells-----------
IBD_CD8T <- read.table("SLE_regulons_monocytes.txt",header = TRUE,sep = "\t")
IBD_CD8T$x <- -log10(IBD_CD8T$SpecificityScore_p)
IBD_CD8T$y <- -log10(IBD_CD8T$ImportanceWeightScore_p)
threshold <- -log10(0.05)

# Highlight subset (e.g., x > 2 and y > 2)
IBD_CD8T$highlight <- ifelse(IBD_CD8T$x >= threshold & IBD_CD8T$y >= threshold, "Significant", "Nonsignificant")

# Plot using ggplot2
library(ggplot2)
# 修改highlight颜色
ggplot(IBD_CD8T, aes(x = x, y = y)) +
  # Scatter plot points
  geom_point(aes(color = highlight), size = 2, alpha = 0.7) +
  # Ellipse for entire data
  stat_ellipse(aes(x = x, y = y), level = 0.95, color = "lightblue", size = 0.5) +
  # Ellipse for highlighted points
  stat_ellipse(data = subset(IBD_CD8T, highlight == "Significant"), aes(x = x, y = y),
               level = 0.95, color = "red", size = 0.5) +
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




###--------------------------------UC mDC-----------
IBD_CD8T <- read.table("UC_regulons_mDC.txt",header = TRUE,sep = "\t")
IBD_CD8T$x <- -log10(IBD_CD8T$SpecificityScore_p)
IBD_CD8T$y <- -log10(IBD_CD8T$ImportanceWeightScore_p)
threshold <- -log10(0.05)

# Highlight subset (e.g., x > 2 and y > 2)
IBD_CD8T$highlight <- ifelse(IBD_CD8T$x >= threshold & IBD_CD8T$y >= threshold, "Significant", "Nonsignificant")

# Plot using ggplot2
library(ggplot2)
# 修改highlight颜色
ggplot(IBD_CD8T, aes(x = x, y = y)) +
  # Scatter plot points
  geom_point(aes(color = highlight), size = 2, alpha = 0.7) +
  # Ellipse for entire data
  stat_ellipse(aes(x = x, y = y), level = 0.95, color = "lightblue", size = 0.5) +
  # Ellipse for highlighted points
  stat_ellipse(data = subset(IBD_CD8T, highlight == "Significant"), aes(x = x, y = y),
               level = 0.95, color = "red", size = 0.5) +
  # Dashed vertical and horizontal lines
  geom_vline(xintercept = threshold, linetype = "dashed", color = "blue", size = 1) +
  geom_hline(yintercept = threshold, linetype = "dashed", color = "blue", size = 1) +
  # Title and labels
  labs(
    title = "mDC",
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

