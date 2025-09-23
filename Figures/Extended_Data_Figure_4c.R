library(ggplot2)
library(ggpubr)
library(dplyr)
library(forcats)
library(patchwork)
library(rstatix)
library(scMORE)

setwd('/share2/pub/zhouyj/zhouyj/ctDRTF_test/10x_pbmc/test_Omega/scMORE')

# Set phenotypes and Omega
phenotypes <- c('Baso_count', 'Eosin_count', 'HL_count', 'Lymphocyte_count', 
               'Lymphocyte_percent', 'MCHC_count', 'MCV_count', 'Monocyte_count', 
               'Neutr_count', 'WBC_count')
Omega <- c(0.2,0.4,0.6,0.8,1,2,5,10)

results_df <- data.frame()

# Load scMORE results
for (pheno in phenotypes) {
      print(pheno)
  for (sz in Omega) {
      file_path <- glue::glue('{pheno}/scMore_Omega{sz}.rds')
      print(sz)
      #print(file_path)
      if (file.exists(file_path)) {
        data <- readRDS(file_path)
        data <- data$scMore_trait_results
        if( 'RegulonScore' %in% colnames(data))
        data$TRS <- data$RegulonScore
        # Calculate E-score
          for(i in unique(data$Celltype)){
          escore <- getEnergyScore(data, targetCelltype = i)
          print(escore)
          
          # Bind results
          results_df <- rbind(results_df, data.frame(
            Phenotype = pheno,
            Celltype = i,
            Omega = sz,
            E_statistic = escore
          ))
          }
      } else {
        warning(paste("File not found:", file_path))
      }
    }
}

# Convert factor
results_df$Omega <- factor(results_df$Omega, levels = unique(results_df$Omega))
results_df$Phenotype <- factor(results_df$Phenotype, levels = phenotypes)


# Plot and test
p1 <- ggplot(results_df, aes(x = Omega, y = E_statistic)) +
  geom_boxplot(fill = "lightblue", alpha = 0.7,outlier.shape=NA) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.6) +
  stat_compare_means(
    comparisons = list(
      c("0.2", "1"),
      c("0.4", "1"),
      c("0.6", "1"),
      c("0.8", "1"),
      c("1", "2"),
      c("1", "5"),
      c("1", "10")
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
    title = "E-statistic Score by Omega",
    x = "Omega",
    y = "E-statistics"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))


ggsave('/share2/pub/zhouyj/zhouyj/ctDRTF_test/10x_pbmc/figures/var_cell_Omega_all_pheno.pdf',p1,width = 12)
write.csv(p1$data,'/share2/pub/zhouyj/zhouyj/ctDRTF_test/10x_pbmc/Extended_Data_Figure_4_c.csv',quote=F,row.names=F)


