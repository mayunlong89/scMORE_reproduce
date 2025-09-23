library(ggplot2)
library(ggpubr)
library(dplyr)
library(forcats)
library(patchwork)
library(rstatix)

setwd('/share2/pub/zhouyj/zhouyj/ctDRTF_test/10x_pbmc/simu_celltype/scMORE')

# Set phenotypes and cell_types
phenotypes <- c('Baso_count', 'Eosin_count', 'HL_count', 'Lymphocyte_count', 
               'Lympnocyte_percent', 'MCHC_count', 'MCV_count', 'Monocyte_count', 
               'Neutr_count', 'WBC_count')
cell_types <- c(2,4,6,8)
replicates <- 1:5

results_df <- data.frame()

# Load scMORE results
for (pheno in phenotypes) {
  for (ct in cell_types) {
    for (rep in replicates) {
      file_path <- glue::glue('{pheno}/scMore_ct{ct}_{rep}.rds')
      #print(file_path)
      if (file.exists(file_path)) {
        data <- readRDS(file_path)
        data <- data$scMore_trait_results
        print(head(data$Significance))
          # Calculate significant eRegulon counts
          sig_count <- sum(data$Significance == "Significant", na.rm = TRUE)
          print(sig_count)
          
          # Bind results
          results_df <- rbind(results_df, data.frame(
            Phenotype = pheno,
            CellTypeNum = paste('ct',ct,sep=''),
            Replicate = rep,
            SignificantCount = sig_count,
            Avg_SignificantCount = sig_count/ct
          ))
      } else {
        warning(paste("File not found:", file_path))
      }
    }
  }
}

# Convert fatcor
results_df$CellTypeNum <- factor(results_df$CellTypeNum, levels = unique(results_df$CellTypeNum))
results_df$Phenotype <- factor(results_df$Phenotype, levels = phenotypes)


# 1. Plot all spots
p1 <- ggplot(results_df, aes(x = CellTypeNum, y = SignificantCount)) +
  geom_boxplot(fill = "lightblue", alpha = 0.7,outlier.shape=NA) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.6) +
  stat_compare_means(
    comparisons = list(
      c("ct2", "ct4"),
      c("ct2", "ct6"),
      c("ct2", "ct8")
    ),
    method = "wilcox.test",          # Wilcox test
    label = "p.format",              # Add pvalue
    aes(label = ifelse(..p.. < 0.001, 
                      paste0(..p.signif.., "\n(p = ", scales::scientific(..p.., digits = 2), ")"),
                      paste0(..p.signif.., "\n(p = ", round(..p.., 3), ")"))),
    bracket.size = 0.3,              
    tip.length = 0.02,               
    vjust = -0.5                     
  ) +
  labs(
    title = "Significant Regulons by Cell Type Number",
    x = "Number of Cell Types",
    y = "Count of Significant Regulons"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))


ggsave('/share2/pub/zhouyj/zhouyj/ctDRTF_test/10x_pbmc/figures/var_cell_types_all_pheno.pdf',p1,width = 8)
write.csv(p1$data,'/share2/pub/zhouyj/zhouyj/ctDRTF_test/10x_pbmc/ Extended_Data_Figure_5_a.csv',quote=F,row.names=F) 

p2 <- ggplot(results_df, aes(x = CellTypeNum, y = Avg_SignificantCount)) +
  geom_boxplot(fill = "lightblue", alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.6) +
  stat_compare_means(
    comparisons = list(
      c("ct2", "ct4"),
      c("ct2", "ct6"),
      c("ct2", "ct8")
    ),
    method = "wilcox.test",          # 或 "t.test"
    label = "p.format",              # 显示 p 值（默认格式）
    aes(label = ifelse(..p.. < 0.001, 
                      paste0(..p.signif.., "\n(p = ", scales::scientific(..p.., digits = 2), ")"),
                      paste0(..p.signif.., "\n(p = ", round(..p.., 3), ")"))),
    bracket.size = 0.3,              # 调整括号粗细
    tip.length = 0.02,               # 调整指示线长度
    vjust = -0.5                     # 调整标签垂直位置
  ) +
  labs(
    title = "Average Significant Regulons by Cell Type Number",
    x = "Number of Cell Types",
    y = "Count of Significant Regulons"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave('/share2/pub/zhouyj/zhouyj/ctDRTF_test/10x_pbmc/figures/var_cell_types_avg_all_pheno.pdf',p2,width = 8)
write.csv(p2$data,'/share2/pub/zhouyj/zhouyj/ctDRTF_test/10x_pbmc/ Extended_Data_Figure_5_b.csv',quote=F,row.names=F)




