

#####---------------Figure 5E---(old Sup_Figure_S26 codes were deprecated)

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
#data$Cellcount <- factor(data$Cellcount, levels = c("1K", "3K", "5K", "7K", "9K"))

# Create the Bar Plot
ggplot(data, aes(x = Method, y = Score, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  facet_wrap(~ Trait, scales = "free_y") +  # Separate plots for each trait
  scale_fill_manual(values = c("#E9E4A6","#E9B78A","#F07B52")) +  # Custom colors
  labs(
    title = "Method Comparison Across Methods and Traits",
    x = "",
    y = "E-statistics",
    fill = "Method"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 10, face = "bold")
  )



##---------mean 8 aging-related traits

# Calculate Mean and Standard Error for Each Method and Cellcount
summary_data <- data %>%
  group_by(Method) %>%
  summarise(
    Mean_Score = mean(Score),
    SD_Score = sd(Score),
    SE_Score = SD_Score / sqrt(n()),  # Standard Error
    .groups = 'drop'
  )



# Create Bar Plot with Error Bars
ggplot(summary_data, aes(x = Method, y = Mean_Score, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = Mean_Score - SE_Score, ymax = Mean_Score + SE_Score),
                position = position_dodge(0.7), width = 0.1) +
  scale_fill_manual(values = c("#E9E4A6","#E9B78A","#F07B52")) +  # Custom colors
  labs(
    title = "Average Scores Across 8 aging-related diseases by Method and Cell Count",
    x = "",
    y = "E-statistics",
    fill = "Method"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 12)
  )




