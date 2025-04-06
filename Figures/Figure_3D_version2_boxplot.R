# Required Libraries
library(ggplot2)
library(dplyr)
library(tidyr)


setwd("~/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")

# Example Data Structure
# Load your data (modify the file path and data structure based on your data)
data <- read.csv("01_figure3D_data.csv", header = TRUE)

# Ensure that 'Method' and 'Cellcount' are factors for proper grouping
data$Method <- factor(data$Method, levels = c( "Average", "MAGMA_Celltyping","scMORE"))
data$Cellcount <- factor(data$Cellcount, levels = c("1K", "3K", "5K", "7K", "9K"))




##---------meam 10 blood cell trait

# Calculate Mean and Standard Error for Each Method and Cellcount
summary_data <- data %>%
  group_by(Method, Cellcount) %>%
  summarise(
    Mean_Score = mean(Score),
    SD_Score = sd(Score),
    SE_Score = SD_Score / sqrt(n()),  # Standard Error
    .groups = 'drop'
  )


library(ggplot2)
library(dplyr)

# Boxplot + jitter 展示分布（每个 Cellcount 中不同 Method）
ggplot(data, aes(x = Cellcount, y = Score, fill = Method)) +
  geom_boxplot(width = 0.5, alpha = 0.8, color = "black", outlier.shape = NA,
               position = position_dodge(width = 0.8)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8),
              size = 1.5, alpha = 0.5, color = "gray30") +
  scale_fill_manual(values = c("#E9E4A6", "#E9B78A", "#F07B52")) +
  labs(
    title = "Scores Across 10 Traits by Method and Cell Count",
    x = "Cell Count",
    y = "E-statistics",
    fill = "Method"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 12)
  )+
  coord_cartesian(ylim = c(0.8, NA))





