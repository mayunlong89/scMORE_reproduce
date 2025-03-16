
# Required Libraries
library(ggplot2)
library(dplyr)
library(tidyr)


setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")


# Example Data Structure
# Load your data (modify the file path and data structure based on your data)
data <- read.csv("01_Supplemental_Figure_S10.csv", header = TRUE)

# Ensure that 'Method' and 'Cellcount' are factors for proper grouping
data$Method <- factor(data$Method, levels = c( "fastBAT", "MAGMA","SMR","FUSION","S-PrediXcan","S-MultiXcan"))
data$Cellcount <- factor(data$Cellcount, levels = c("1K", "3K", "5K", "7K", "9K"))

# Create the Bar Plot
ggplot(data, aes(x = Cellcount, y = Score, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  facet_wrap(~ Trait, scales = "free_y") +  # Separate plots for each trait
  scale_fill_manual(values = c("lightblue","aquamarine4","chartreuse3","#E9E4A6","#E9B78A","#F07B52")) +  # Custom colors
  labs(
    title = "Method Comparison Across Cell Counts and Traits",
    x = "Cell Count",
    y = "Score",
    fill = "Method"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 10, face = "bold")
  )



##---------meam 10 blood cell trait bar

# Calculate Mean and Standard Error for Each Method and Cellcount
summary_data <- data %>%
  group_by(Method, Cellcount) %>%
  summarise(
    Mean_Score = mean(Score),
    SD_Score = sd(Score),
    SE_Score = SD_Score / sqrt(n()),  # Standard Error
    .groups = 'drop'
  )



# Create Bar Plot with Error Bars
ggplot(summary_data, aes(x = Cellcount, y = Mean_Score, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = Mean_Score - SE_Score, ymax = Mean_Score + SE_Score),
                position = position_dodge(0.7), width = 0.2) +
  scale_fill_manual(values = c("lightblue","aquamarine4","chartreuse3","#E9E4A6","#E9B78A","#F07B52"))+
  labs(
    title = "Average Scores Across 10 Traits by Method and Cell Count",
    x = "Cell Count",
    y = "Mean Score ± SE",
    fill = "Method"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 12)
  )
