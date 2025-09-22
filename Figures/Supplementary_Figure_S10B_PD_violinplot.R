

setwd("~/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/06-PD/")



# Required Libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Example Data Structure
# Load your data (modify the file path and data structure based on your data)
data <- read.csv("01_benchmarking_PD_aging.csv", header = TRUE)

# Ensure that 'Method' and 'Cellcount' are factors for proper grouping
data$Method <- factor(data$Method, levels = c( "Average", "MAGMA_Celltyping","scMORE"))


##----violin plot

# Calculate Mean and Standard Error for Each Method and Cellcount
summary_data <- data %>%
  group_by(Method) %>%
  summarise(
    Mean_Score = mean(Score),
    SD_Score = sd(Score),
    SE_Score = SD_Score / sqrt(n()),  # Standard Error
    .groups = 'drop'
  )


library(ggplot2)
library(dplyr)

## non-classic theme
ggplot(data, aes(x = Method, y = Score, fill = Method)) +
  geom_violin(trim = FALSE, width=1,alpha = 0.7, color = "black") +
  geom_boxplot(width = 0.1, alpha = 0.7, color = "black", fill="white",outlier.shape = NA) +
  #stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "black", fatten = 0.5) +  # 添加中位数线
  geom_jitter(width = 0.1, size = 1.5, alpha = 0.6, color = "gray30") +
  scale_fill_manual(values = c("#E9E4A6", "#E9B78A", "#F07B52")) +
  labs(
    title = "Scores Across 5 Autoimmune Diseases by Method",
    x = "Method",
    y = "E-statistics",
    fill = "Method"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 12)
  ) +
  coord_cartesian(ylim = c(0.8, NA))

