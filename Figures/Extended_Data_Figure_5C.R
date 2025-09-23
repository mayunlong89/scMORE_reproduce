library(ggplot2)
library(ggpubr)
library(dplyr)
library(forcats)
library(patchwork)
library(rstatix)

setwd('/share2/pub/zhouyj/zhouyj/ctDRTF_test/10x_pbmc/simu_cellcounts/scMORE')

# Set phenotypes and cell counts
phenotypes <- c('Baso_count', 'Eosin_count', 'HL_count', 'Lymphocyte_count', 
               'Lymphocyte_percent', 'MCHC_count', 'MCV_count', 'Monocyte_count', 
               'Neutr_count', 'WBC_count')
cell_counts <- c(1000,5000,10000,15000,20000,25000)

results_df <- data.frame()

# Load scMORE results
for (pheno in phenotypes) {
  for (sz in cell_counts) {
      file_path <- glue::glue('{pheno}/scMore_{sz}.rds')
      if (file.exists(file_path)) {
        data <- readRDS(file_path)
        data <- data$scMore_trait_results
        print(head(data$Significance))
        # Check file 
          # Calculate significant eRegulon counts
          sig_count <- sum(data$Significance == "Significant", na.rm = TRUE)
          print(sig_count)
          
          # Bind results
          results_df <- rbind(results_df, data.frame(
            Phenotype = pheno,
            CellCount = sz,
            SignificantCount = sig_count
          ))
      } else {
        warning(paste("File not found:", file_path))
      }
    }
}

# Convert factor
results_df$CellCount <- factor(results_df$CellCount, levels = unique(results_df$CellCount))
results_df$Phenotype <- factor(results_df$Phenotype, levels = phenotypes)


# 1. Polt
p1 <- ggplot(results_df, aes(x = CellCount, y = SignificantCount)) +
  geom_boxplot(fill = "lightblue", alpha = 0.7,outlier.shape=NA) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.6) +
  stat_compare_means(
    comparisons = list(
      c("1000", "5000"),
      c("1000", "10000"),
      c("1000", "15000"),
      c("1000", "20000"),
      c("1000", "25000")
    ),
    method = "wilcox.test",          # Wilcox test
    label = "p.format",              # Add P value
    aes(label = ifelse(..p.. < 0.001, 
                      paste0(..p.signif.., "\n(p = ", scales::scientific(..p.., digits = 2), ")"),
                      paste0(..p.signif.., "\n(p = ", round(..p.., 3), ")"))),
    bracket.size = 0.3,              
    tip.length = 0.02,               
    vjust = -0.5                     
  ) +
  labs(
    title = "Significant Regulons by Cell Counts",
    x = "Number of Cell Counts",
    y = "Count of Significant Regulons"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))


ggsave('/share2/pub/zhouyj/zhouyj/ctDRTF_test/10x_pbmc/figures/var_cell_cellcounts_all_pheno.pdf',p1,width = 8)

