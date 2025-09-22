
#2025-06-27
setwd("~/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2025-05-20-scMORE_revision/07-羿君-修稿/simu_GWAS/PD_sex_stratified/")


#--5691
ieu_a_812 <- readRDS("../ieu-a-812/scMore_ieu-a-812.rds")

write.csv(ieu_a_812$scMore_trait_results, file="results_ieu_a_812_PD_sample5691.csv",quote=F, row.names = F)




#--1672
ieu_a_818 <- readRDS("../ieu-a-818/scMore_ieu-a-818.rds")

write.csv(ieu_a_818$scMore_trait_results, file="results_ieu_a_818_PD_sample1672.csv",quote=F, row.names = F)



#--482730
ieu_b_7_maf_not_filter <- readRDS("../ieu-b-7_没有maf_filter/scMore_ieu-b-7.rds")

write.csv(ieu_b_7_maf_not_filter$scMore_trait_results, file="results_ieu_b_7_PD_sample482730_maf_not_filter.csv",quote=F, row.names = F)


#--482730---final version
ieu_b_7_maf_filter <- read.csv("../results_PD_scMore_trait_results_phastCons_withoutExons_perm1000_alpha=1.csv")

length(ieu_b_7_maf_filter$Significance[which(ieu_b_7_maf_filter$Significance=="Significant")])



###------analysis for sample sizes:
###------analysis for sample sizes:
###------analysis for sample sizes:

#setwd("~/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2025-05-20-scMORE_revision/01_PD_GWAS/01_PD_sex/")
setwd("~/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2025-05-20-scMORE_revision/07-羿君-修稿/simu_GWAS/PD_sex_stratified/")


female100K <-  readRDS("scMore_PD_female_100k.rds")
male100K <-  readRDS("scMore_PD_male_100k.rds")


write.csv(male100K$scMore_trait_results,file="scMORE_PD_male_100k.csv",quote = F,row.names = F)
write.csv(female100K$scMore_trait_results,file="scMORE_PD_female_100k.csv",quote = F,row.names = F)



setwd("~/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2025-05-20-scMORE_revision/07-羿君-修稿/simu_GWAS/PD_sex_stratified/")


female20K <-  readRDS("scMore_PD_female_no_ukb_20k.rds")
male20K <-  readRDS("scMore_PD_male_no_ukb_20K.rds")


write.csv(male20K$scMore_trait_results,file="scMORE_PD_male_20k.csv",quote = F,row.names = F)
write.csv(female20K$scMore_trait_results,file="scMORE_PD_female_20k.csv",quote = F,row.names = F)



# 假设你的数据框是 df
library(dplyr)

ieu_b_7_maf_filter %>%
  filter(Significance == "Significant") %>%
  group_by(Celltype) %>%
  summarise(Significant_Count = n()) %>%
  arrange(desc(Significant_Count))


ieu_a_818$scMore_trait_results %>%
  filter(Significance == "Significant") %>%
  group_by(Celltype) %>%
  summarise(Significant_Count = n()) %>%
  arrange(desc(Significant_Count))


ieu_a_812$scMore_trait_results %>%
  filter(Significance == "Significant") %>%
  group_by(Celltype) %>%
  summarise(Significant_Count = n()) %>%
  arrange(desc(Significant_Count))


library(dplyr)
library(tidyr)  # 解决 replace_na 的错误

# Female summary (N = 24,054 samples)
female_summary <- female20K$scMore_trait_results %>%
  filter(Significance == "Significant") %>%
  group_by(Celltype) %>%
  summarise(Female = n())

# Male summary (N = 19,773 males)
male_summary <- male20K$scMore_trait_results %>%
  filter(Significance == "Significant") %>%
  group_by(Celltype) %>%
  summarise(Male = n())

# 合并并计算平均值
combined_summary_20k <- full_join(female_summary, male_summary, by = "Celltype") %>%
  mutate(
    Female = replace_na(Female, 0),
    Male = replace_na(Male, 0),
    Mean = (Female + Male) / 2
  ) %>%
  arrange(desc(Mean))

print(combined_summary_20k)
(19773+24054)/2

library(dplyr)
library(tidyr)  # 解决 replace_na 的错误

# Female summary (N = 104,082 samples)
female_summary <- female100K$scMore_trait_results %>%
  filter(Significance == "Significant") %>%
  group_by(Celltype) %>%
  summarise(Female = n())

# Male summary (N = 102,680 samples)
male_summary <- male100K$scMore_trait_results %>%
  filter(Significance == "Significant") %>%
  group_by(Celltype) %>%
  summarise(Male = n())

# 合并并计算平均值
combined_summary <- full_join(female_summary, male_summary, by = "Celltype") %>%
  mutate(
    Female = replace_na(Female, 0),
    Male = replace_na(Male, 0),
    Mean = (Female + Male) / 2
  ) %>%
  arrange(desc(Mean))

print(combined_summary)
(104082+ 102680)/2



library(ggplot2)
library(dplyr)

# 三组数据手动输入（也可从文件读取或转换）
#PD_sample_1672
df1 <- data.frame(Celltype = c("N", "OPC", "AS", "EC", "T", "MG", "ODC"),
                  Significant_Count = c(13, 11, 11, 8, 7, 6, 3),
                  Dataset = "PD_sample_2K")

#PD_sample_5691
df2 <- data.frame(Celltype = c("AS", "N", "OPC", "T", "EC", "MG", "ODC"),
                  Significant_Count = c(14, 14, 14, 9, 8, 8, 6),
                  Dataset = "PD_sample_6K")

#PD_sample_21914
df3 <- data.frame(Celltype = c("AS", "N", "OPC", "T", "EC", "MG", "ODC"),
                  Significant_Count = c(56, 33, 49, 35, 39, 31, 27),
                  Dataset = "PD_sample_20K")

#PD_sample_103381
df4 <- data.frame(Celltype = c("AS", "N", "OPC", "T", "EC", "MG", "ODC"),
                  Significant_Count = c(53, 32, 50, 35, 37, 32, 26),
                  Dataset = "PD_sample_100K")

#PD_sample_482730
df5 <- data.frame(Celltype = c("AS", "OPC", "EC", "T", "MG", "N", "ODC"),
                  Significant_Count = c(61, 55, 41, 37, 35, 34, 29),
                  Dataset = "PD_sample_500K")



# 合并数据
df_all <- bind_rows(df1, df2, df3,df4,df5)

#Save data
#write.csv(df_all,file="response_to_PD_GWAS_Sample_sizes.csv",quote = F, row.names = F)

#df_all$Dataset <- factor(df_all$Dataset, level=c("PD_sample_1672","PD_sample_5691","PD_sample_21914","PD_sample_103381","PD_sample_482730"))
df_all$Dataset <- factor(df_all$Dataset, level=c("PD_sample_2K","PD_sample_6K","PD_sample_20K","PD_sample_100K","PD_sample_500K"))

df_all$Celltype <- factor(df_all$Celltype, level=c("AS","OPC","EC","T","MG","N","ODC"))

# 画图：每个Celltype为x轴，显著count为y轴，点为原始值，箱线图为整体分布
ggplot(df_all, aes(x = Celltype, y = Significant_Count)) +
  geom_boxplot(outlier.shape = NA, fill = "lightblue", color = "black") +
  geom_jitter(aes(color = Dataset), width = 0.2, size = 2, alpha = 0.8) +
  scale_color_manual(values = c("PD_sample_2K" = "#1f77b4", "PD_sample_6K" = "#97EC71","PD_sample_20K"="#DE9DD6", "PD_sample_100K" = "orange","PD_sample_500K" = "#d62728")) +  # 可手动指定颜色
  theme_bw() +
  labs(
    title = "Significant Regulon Count per Celltype Across Datasets",
    y = "Significant Regulon Count",
    x = "Cell Type",
    color = "Dataset"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


write.csv(df_all,file="scMORE_PD_GWAS_sample_sizes_boxplot_Extended_Data_Figure_5d.csv",quote = F,row.names = F)



