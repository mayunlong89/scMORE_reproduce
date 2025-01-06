

###--------------Figure 2G boxplot with lines

#install.packages("ggpubr")
library(ggpubr)  #------mi
### 1) mi
#setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/ctDRTF_assessment_results/Monocyte_assessment/magma_different_windows_monocytes_results_plots/")
#mode0_vs_1 <- read.table("mode_1_vs_0_mi.txt",header = TRUE)
#ggpaired(mode0_vs_1,cond1 = "mode0",cond2 = "mode1",fill=c("grey","lightblue"),palette = "jco")
#t.test(mode0_vs_1$mode1,mode0_vs_1$mode0,paired = TRUE)


setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
consine_average10traits <- read.table("01_10blood_cell_traits_cosine_average_benchamrk.txt",header = TRUE,sep = "\t")

###----Boxplot for mode0 vs mode1
ggpaired(consine_average10traits,cond1 = "average",cond2 = "cosine",fill=c("grey","lightblue"),palette = "jco")

t.test(consine_average10traits$average,consine_average10traits$cosine,paired = TRUE)

##----boxplot with line---
# Enhanced ggpaired plot with improved aesthetics and flexibility
ggpaired(
  data = consine_average10traits,  # Specify the dataset
  cond1 = "average",              # Column name for the first condition
  cond2 = "cosine",               # Column name for the second condition
  fill = c("grey", "lightblue"),  # Colors for the paired conditions
  palette = "jco",                # Color palette for the plot
  line.color = "black",           # Color for connecting lines
  line.size = 0.5,                # Thickness of the connecting lines
  point.size = 3,                 # Size of the points
  title = "",  # Title
  xlab = "",            # X-axis label
  ylab = "E-statistics"                  # Y-axis label
) +
  theme_minimal() +  # Apply minimal theme for cleaner appearance
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  # Style for title
    axis.title = element_text(size = 14),  # Style for axis titles
    axis.text = element_text(size = 12)    # Style for axis text
  )
