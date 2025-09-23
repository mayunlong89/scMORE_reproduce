library(ggplot2)
library(ggpubr)
library(dplyr)
library(forcats)
library(patchwork)
library(rstatix)

setwd('/share2/pub/zhouyj/zhouyj/ctDRTF_test/10x_pbmc/simu_GWAS/scMORE')

# Set phenotypes, logFC and significant regulons
phenotypes <- c( 'LogFC')
logfc <- c(-0.6,-0.5,-0.4,-0.3,-0.2,-0.1)
sigregulon <- c('ETS1','FOXO1','FOXP1','IKZF1','LEF1','ZEB1','ZNF721')

results_df <- data.frame()
data_df <- data.frame()

# Load scMORE results
for (pheno in phenotypes) {
for(i in c(1:100)){
  for (sz in logfc) {
      file_path <- glue::glue('{pheno}/scMore_logfc{sz}_{i}.rds')
      #print(file_path)
      if (file.exists(file_path)) {
        data <- readRDS(file_path)
        data <- data$scMore_trait_results 
        data <- subset(data, RegulonName %in% sigregulon)
        data <- subset(data, Celltype == 'CD8+ T cells')
        if('RegulonScore_p' %in% colnames(data)){
         data$TRS_P_value <- data$RegulonScore_p
        }
        print(head(data))
          ## Calculate significant regulon counts
          sig_count <- sum(data$TRS_P_value < 0.05)
                    print(sig_count)
          
          # Bind results and set logFC -0.6 as baseline
          results_df <- rbind(results_df, data.frame(
            Phenotype = 'Lymphocyte_count',
            LogFC = sz+0.6,
            Repeat = i,
            RegulonCount = sig_count
          ))
      } else {
        warning(paste("File not found:", file_path))
      }
    }
}
}

# Convert factor
results_df$LogFC <- factor(results_df$LogFC, levels =unique(results_df$LogFC))
results_df$Phenotype <- factor(results_df$Phenotype, levels = 'Lymphocyte_count')

# Set RegulonCount/length(sigregulon) as power
results_df <- results_df %>% mutate(Power = RegulonCount/7)

results_df <- results_df %>%
  mutate(
    LogFC_group = factor(LogFC),  
    x_index = as.numeric(LogFC_group)
  )

p_combined <- ggplot(results_df, aes(x = LogFC, y = Power, group = LogFC)) +
  geom_boxplot(fill = "lightblue", alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.6) +
  ylim(0,1)+
  stat_summary(
    aes(y = Power), 
    fun = mean, 
    geom = "point", 
    shape = 19, 
    size = 2, 
    color = "darkorange"
  ) +
  stat_summary(
    aes(y = Power, group = 1), 
    fun = mean, 
    geom = "line", 
    color = "grey70", 
    linewidth = 0.8
  )+
  theme_minimal()

ggsave('/share2/pub/zhouyj/zhouyj/ctDRTF_test/10x_pbmc/figures/var_logfc_lyc_combined.pdf',p_combined)
write.csv(p_combined$data,'/share2/pub/zhouyj/zhouyj/ctDRTF_test/10x_pbmc/Extended_Data_Figure_4_e.csv',quote=F,row.names=F)



