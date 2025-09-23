library(ggplot2)
library(ggpubr)
library(dplyr)
library(forcats)
library(patchwork)
library(rstatix)

setwd('/share2/pub/zhouyj/zhouyj/ctDRTF_test/10x_pbmc/simu_GWAS/scMORE')

# Set phonotypes and simulation data
phenotypes <- c( 'Lymphocyte_count_simu')
simu <- c('data_LC_New_100W_h0.5_p_0.01_a_0.28')

#Load real data results
raw_scmore <- readRDS('/share2/pub/zhouyj/zhouyj/ctDRTF_test/10x_pbmc/simu_GWAS/scMORE/LogFC/scMore_logfc0.rds')
raw_scmore <- raw_scmore$scMore_trait_results

results_df <- data.frame()

# Load sumulation results
for (pheno in phenotypes) {
  for (sz in simu) {
      file_path <- glue::glue('{pheno}/scMore_{sz}.rds')
      #print(file_path)
      if (file.exists(file_path)) {
        data <- readRDS(file_path)
        data <- data$scMore_trait_results
        if('RegulonScore_p' %in% colnames(data)){
        data$TRS_P_value <- data$RegulonScore_p
        }
        #print(head(data$Significance))
          ## Calculate significant eRegulon counts
          #sig_count <- sum(data$Significance == "Significant", na.rm = TRUE)
          
          # Bind results
          inner_results_df <- inner_join(
                                    raw_scmore %>% select(RegulonName, Celltype, TRS_P_value),
                                    data %>% select(RegulonName, Celltype, TRS_P_value),
                                    by = c("RegulonName", "Celltype"),  
                                    suffix = c("_Exp", "_Obs"))  
          inner_results_df$simu <- sz
          inner_results_df$Phenotype <- pheno
          results_df <- rbind(results_df, inner_results_df )
      } else {
        warning(paste("File not found:", file_path))
      }
    }
}

#Prepare qq plot data
qq_data <- results_df %>%
  arrange(TRS_P_value_Exp) %>%
  mutate(
    rank_exp = row_number(),
    quantile_exp = rank_exp / (n() + 1),
    TRS_P_value_Obs_sorted = sort(TRS_P_value_Obs)
  )

# Split simulation data
qq_data2 <- results_df %>%
  group_by(simu) %>%
  arrange(TRS_P_value_Exp) %>%
  mutate(
    rank_exp = row_number(),
    quantile_exp = rank_exp / (n() + 1),
    TRS_P_value_Obs_sorted = sort(TRS_P_value_Obs)
  ) %>%
  ungroup()

# QQ plot
p_logqq2 <- ggplot(qq_data2, aes(x = -log10(TRS_P_value_Exp), y = -log10(TRS_P_value_Obs_sorted))) +
  geom_point(aes(color = Phenotype), alpha = 0.7, size = 2) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  facet_wrap(~ simu, scales = "free") +
  labs(
    title = "GWAS simulation under the null",
    x = "Expected -log10(P-value)",
    y = "Observed -log10(P-value)",
    color = "Phenotype"
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))


ggsave('/share2/pub/zhouyj/zhouyj/ctDRTF_test/10x_pbmc/figures/var_null_lyc_split_simu.pdf',p_logqq2,width = 8)
write.csv(p_logqq2$data,'/share2/pub/zhouyj/zhouyj/ctDRTF_test/10x_pbmc/ Extended_Data_Figure_4_d.csv',quote=F,row.names=F)


