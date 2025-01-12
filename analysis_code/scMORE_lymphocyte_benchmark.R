###------------------------------------------------------------------------------------
#2025-01-02
###------------------------------------------------------------------------------------

##load packages
suppressPackageStartupMessages({
  library(Pando)
  library(Seurat)
  library(Signac)
  library(ArchR)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(tidyr)
  library(dplyr)
  library(dotgen)
  library(org.Hs.eg.db)
  library(COSG)
  library(fitdistrplus)
  library(AnnotationDbi)
  library(GenomicRanges)
  library(IRanges)
})

library(Seurat)

###required load R packages

#devtools::install_github("mayunlong89/scMORE")

library(scMORE)
library(GenomicRanges)
library(IRanges)
library(Seurat)
library(Signac)


##---2000Cells
#saveRDS(pbmc_10x_real_data_downsampled_2000,file = "/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/10X_PBMC/10X_PBMC_downsample_2000cells.rds")

#---load pbmc_10x_real_data_downsampled_2000
pbmc_10x_real_data_downsampled_2000 <- readRDS("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/10X_PBMC/10X_PBMC_downsample_2000cells.rds")
DimPlot(pbmc_10x_real_data_downsampled_2000, reduction = "umap.new",cols = c("#67ADB7","#E77A77"))
single_cell <-pbmc_10x_real_data_downsampled_2000


#----generate regulons from single-cell multiomics data
n_targets=5
grn_outputs <- createRegulon(single_cell, n_targets)


##-------------------------------------------------------------------------------------------
##-------------------------------------------------------------------------------------------
##-------------------------------------------------------------------------------------------






##--scMORE for lymphocyte count-----------------------------------------------------------------------------------------
##MAGMA-results of lymphocyte count
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
magma_lymp <- read.table("lymp_count_processed_magma_results.genes.out",header = TRUE)
#magma_results <- getGeneScore(magma_mono)
#head(magma_results)
gene_info_lymp <- magma_lymp

snp_info_lymp <- read.csv("lymphocyte_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)


MORE_results_lymp1  <- scMore(single_cell = pbmc_10x_real_data_downsampled_2000,
                              gene_info = gene_info_lymp,
                              snp_info = snp_info_lymp,
                              n_targets= 5,
                              perm_n = 1000,
                              theta=0.5,
                              alpha=1,
                              top_n=5,
                              buffer=500,
                              method = 'cosine',
                              nSeed=1234)

MORE_results_lymp2  <- scMore(single_cell = pbmc_10x_real_data_downsampled_2000,
                        gene_info = gene_info_lymp,
                        snp_info = snp_info_lymp,
                        n_targets= 5,
                        perm_n = 1000,
                        theta=0.5,
                        alpha=1,
                        top_n=5,
                        buffer=500,
                        method = 'average',
                        nSeed=1234)

write.csv(MORE_results_lymp2$scMore_trait_results,file="01_lymp_count_scMORE_withJSI_2000_v5_scale_twoCelltypes_separate_nonScale_TF_average_pow1_notAddModule_monocyte_CD8T_V3_1000_average.csv")
write.csv(MORE_results_lymp1$scMore_trait_results,file="01_lymp_count_scMORE_withJSI_2000_v5_scale_twoCelltypes_separate_nonScale_TF_average_pow1_notAddModule_monocyte_CD8T_V3_1000_consine.csv")

lym_count_average <- getEnergyScore(MORE_results_lymp2$scMore_trait_results,targetCelltype = 2)
log2(lym_count_average+1)+1
lym_count_cosine <- getEnergyScore(MORE_results_lymp1$scMore_trait_results,targetCelltype = 2)
log2(lym_count_cosine+1)+1



###--------------------------------CD8T cells-----------
data_CD8T <- read.table("lymp_count_regulon_CD8T.txt",header = TRUE,sep = "\t")
data_CD8T$x <- -log10(data_CD8T$SpecificityScore_p)
data_CD8T$y <- -log10(data_CD8T$ImportanceWeightScore_p)
threshold <- -log10(0.05)

# Highlight subset (e.g., x > 2 and y > 2)
data_CD8T$highlight <- ifelse(data_CD8T$x >= threshold & data_CD8T$y >= threshold, "Significant", "Non-significant")

# Plot using ggplot2
library(ggplot2)
# 修改highlight颜色
ggplot(data_CD8T, aes(x = x, y = y)) +
  # Scatter plot points
  geom_point(aes(color = highlight), size = 2, alpha = 0.7) +
  # Ellipse for entire data
  stat_ellipse(aes(x = x, y = y), level = 0.95, color = "lightblue", size = 0.5) +
  # Ellipse for highlighted points
  stat_ellipse(data = subset(data_CD8T, highlight == "Significant"), aes(x = x, y = y),
               level = 0.95, color = "red", size = 0.5) +
  # Dashed vertical and horizontal lines
  geom_vline(xintercept = threshold, linetype = "dashed", color = "blue", size = 1) +
  geom_hline(yintercept = threshold, linetype = "dashed", color = "blue", size = 1) +
  # Title and labels
  labs(
    title = "CD8+T cells",
    x = "CTS",
    y = "GRS",
    color = "Group"
  ) +
  # Custom colors for highlight
  scale_color_manual(values = c("Significant" = "red", "Non-significant" = "gray")) +
  # Theme adjustments
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "top"
  )


###--------------------------------Monocytes-----------
data_CD8T <- read.table("lymp_count_regulon_monocyte.txt",header = TRUE,sep = "\t")
data_CD8T$x <- -log10(data_CD8T$SpecificityScore_p)
data_CD8T$y <- -log10(data_CD8T$ImportanceWeightScore_p)
threshold <- -log10(0.05)

# Highlight subset (e.g., x > 2 and y > 2)
data_CD8T$highlight <- ifelse(data_CD8T$x >= threshold & data_CD8T$y >= threshold, "Significant", "Non-significant")

# Plot using ggplot2
library(ggplot2)
# 修改highlight颜色
ggplot(data_CD8T, aes(x = x, y = y)) +
  # Scatter plot points
  geom_point(aes(color = highlight), size = 2, alpha = 0.7) +
  # Ellipse for entire data
  stat_ellipse(aes(x = x, y = y), level = 0.95, color = "lightblue", size = 0.5) +
  # Ellipse for highlighted points
  stat_ellipse(data = subset(data_CD8T, highlight == "Significant"), aes(x = x, y = y),
               level = 0.95, color = "orange", size = 1) +
  # Dashed vertical and horizontal lines
  geom_vline(xintercept = threshold, linetype = "dashed", color = "blue", size = 1) +
  geom_hline(yintercept = threshold, linetype = "dashed", color = "blue", size = 1) +
  # Title and labels
  labs(
    title = "Monocytes",
    x = "CTS",
    y = "GRS",
    color = "Group"
  ) +
  # Custom colors for highlight
  scale_color_manual(values = c("Significant" = "red", "Non-significant" = "gray")) +
  # Theme adjustments
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "top"
  )




##-------------------------------------------------------------------------------------------
##-------------------------------------------------------------------------------------------
##-------------------------------------------------------------------------------------------
##-------------------------------------------------------------------------------------------








##-------------------------------------------------------------------------------------------
###---Lymphocyte percent
gene_info_lymp_percent <- read.table("lymp_percent_processed_magma_results.genes.out",header = TRUE)
#snp_info
snp_info_lymp_percent <- read.csv("lymphocyte_percent_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)


MORE_results_lymp_percent1  <- scMore(single_cell = pbmc_10x_real_data_downsampled_2000,
                              gene_info = gene_info_lymp_percent,
                              snp_info = snp_info_lymp_percent,
                              n_targets= 5,
                              perm_n = 1000,
                              theta=0.5,
                              alpha=1,
                              top_n=5,
                              buffer=500,
                              method = 'cosine',
                              nSeed=1234)

MORE_results_lymp_percent2  <- scMore(single_cell = pbmc_10x_real_data_downsampled_2000,
                                      gene_info = gene_info_lymp_percent,
                                      snp_info = snp_info_lymp_percent,
                                      n_targets= 5,
                                      perm_n = 1000,
                                      theta=0.5,
                                      alpha=1,
                                      top_n=5,
                                      buffer=500,
                                      method = 'average',
                                      nSeed=1234)


write.csv(MORE_results_lymp_percent2$scMore_trait_results,file="01_lymp_percent_scMORE_withJSI_2000_v5_scale_twoCelltypes_separate_nonScale_TF_average_pow1_notAddModule_monocyte_CD8T_V3_1000_average.csv")
write.csv(MORE_results_lymp_percent1$scMore_trait_results,file="01_lymp_percent_scMORE_withJSI_2000_v5_scale_twoCelltypes_separate_nonScale_TF_average_pow1_notAddModule_monocyte_CD8T_V3_1000_consine.csv")

lym_percent_average <- getEnergyScore(MORE_results_lymp_percent2$scMore_trait_results,targetCelltype = 2)
log2(lym_percent_average+1)+1
lym_percent_cosine <- getEnergyScore(MORE_results_lymp_percent1$scMore_trait_results,targetCelltype = 2)
log2(lym_percent_cosine+1)+1






##correlation of lymp count and lymp percent
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
consine_average10traits_heatmap <- read.table("01_10blood_cell_traits_heatmap.txt",header = TRUE,sep = "\t")
cor(consine_average10traits_heatmap$lymp_count,consine_average10traits_heatmap$lymp_percent)
cor.test(consine_average10traits_heatmap$lymp_count,consine_average10traits_heatmap$lymp_percent)

data <- read.table("lymp_count_percent.txt",header = TRUE,sep = "\t")
head(data)
cor.test(data$lymp_count,data$lymp_percent)

# 加载 ggplot2
library(ggplot2)
# 绘制散点相关性图
ggplot(data, aes(x = lymp_count, y = lymp_percent,color=Celltype)) +
  geom_point(size = 2, alpha = 0.7) +  # 散点
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # 添加拟合线
  labs(title = "Correlation by Cell Type",
       x = "Lymphocyte count",
       y = "Lymphocyte percent") +
  theme_classic() +  # 经典主题
  facet_wrap(~ Celltype)  # 按细胞类型分面

# Perform correlation test for each cell type
cor_results <- data %>%
  dplyr::group_by(Celltype) %>%
  dplyr::summarise(
    correlation = cor(lymp_count, lymp_percent, use = "complete.obs"),  # Pearson correlation
    p_value = cor.test(lymp_count, lymp_percent)$p.value,  # Extract p-value
    conf_low = cor.test(lymp_count, lymp_percent)$conf.int[1],  # Lower bound of confidence interval
    conf_high = cor.test(lymp_count, lymp_percent)$conf.int[2]   # Upper bound of confidence interval
  )

print(cor_results)


custom_colors <- c("Lymphocyte count" = "darkblue", "Lymphocyte percent" = "orange")
# Plot with updated colors
ggplot(data, aes(x = lymp_count, y = lymp_percent,color=Celltype)) +
  geom_point(size = 1.5, alpha = 0.8,color="#33c4ff") +  # Adjust point size and transparency
  geom_smooth(method = "lm", color = "red", linetype = "dashed", size = 0.8) +  # Add regression line
  labs(
    title = "Correlation of two lymphocyte traits",
    x = "Lymphocyte count",
    y = "Lymphocyte percent",
    color = "Cell types"
  ) +
  scale_color_manual(values = custom_colors) +  # Apply custom colors
  theme_classic() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )



##correlation of lymphocyte count and monocytes count-----Figure 2F

data <- read.table("lymp_count_mono.txt",header = TRUE,sep = "\t")
head(data)
cor.test(data$lymp_count,data$mono_count)

# 加载 ggplot2
library(ggplot2)
# 绘制散点相关性图
ggplot(data, aes(x = lymp_count, y = mono_count,color=Celltype)) +
  geom_point(size = 2, alpha = 0.7) +  # 散点
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # 添加拟合线
  labs(title = "Correlation by Cell Type",
       x = "Lymphocyte count",
       y = "Lymphocyte percent") +
  theme_classic() +  # 经典主题
  facet_wrap(~ Celltype)  # 按细胞类型分面

# Perform correlation test for each cell type
cor_results <- data %>%
  dplyr::group_by(Celltype) %>%
  dplyr::summarise(
    correlation = cor(lymp_count, lymp_percent, use = "complete.obs"),  # Pearson correlation
    p_value = cor.test(lymp_count, lymp_percent)$p.value,  # Extract p-value
    conf_low = cor.test(lymp_count, lymp_percent)$conf.int[1],  # Lower bound of confidence interval
    conf_high = cor.test(lymp_count, lymp_percent)$conf.int[2]   # Upper bound of confidence interval
  )

print(cor_results)

custom_colors <- c("Lymphocyte count" = "darkblue", "Lymphocyte percent" = "orange")
# Plot with updated colors
ggplot(data, aes(x = lymp_count, y = mono_count,color=Celltype)) +
  geom_point(size = 1.5, alpha = 0.8,color="#33c4ff") +  # Adjust point size and transparency
  geom_smooth(method = "lm", color = "red", linetype = "dashed", size = 0.8) +  # Add regression line
  labs(
    title = "Correlation of two lymphocyte traits",
    x = "Lymphocyte count",
    y = "Monocyte count",
    color = "Cell types"
  ) +
  scale_color_manual(values = custom_colors) +  # Apply custom colors
  theme_classic() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )




##----------Heatmap plots------------------------------------Figure 2E
#@ Heatmap
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")

library(pheatmap)


lym_heatmap <- read.table("lymp_count_percent_heatmap.txt",header = TRUE,sep = "\t")


data_h1 <- lym_heatmap[,-1]
rownames(data_h1) <- lym_heatmap[,1]


data_h1_score <- data_h1[,c(1,3,5,7)]

data_h1_p <- data_h1[,c(2,4,6,8)]

data_h1_p2 <- filter(data_h1_p,Lymp_count_CD8T_P < 0.05 | Lymp_percent_CD8T_P <0.05 | Lymp_count_monocyte_P < 0.05 | Lymp_percent_monocyte_P < 0.05 )
data_h1_score2 <- data_h1_score[which(rownames(data_h1_score) %in% rownames(data_h1_p2)),]

p_vals<-t(data_h1_p2)
heatmap_data <- t(data_h1_score2)

p_vals<- data_h1_p2
heatmap_data <- data_h1_score2

pearson_heatmap<-pheatmap(heatmap_data,cellwidth =15, cellheight =15,fontsize=8,
                          border_color="black",scale="none",
                          fontsize_row=8,
                          color=colorRampPalette(rev(c("#E77A77", "#ECBA84","white","lightblue",'#336699')))(102),
                          cluster_rows = F,
                          cluster_cols = F,
                          display_numbers = matrix(ifelse(p_vals <0.01 & heatmap_data > 1.96, "**",
                                                          ifelse(p_vals<0.05 & heatmap_data>1.96,"*","")), nrow(p_vals)))


##-------------------------------------------------------------------------------------------
##-------------------------------------------------------------------------------------------
##-------------------------------------------------------------------------------------------
##-------------------------------------------------------------------------------------------




##-------------------------------------------------------------------------------------------
##MAGMA-results
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")

##---Test results from MAGMA--------------------
#gene_info
magma_mono <- read.table("mono_count_processed_magma_results.genes.out",header = TRUE)
#magma_results <- getGeneScore(magma_mono)
#head(magma_results)
gene_info_mono <- magma_mono
geneRiskScores <- getGeneScore(gene_info_mono)

#snp_info
snp_info_mono <- read.csv("monocyte_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

#geneRiskScores <- getGeneScore(gene_info_lymp)
#snp_info<- snp_info_lymp

n_targets = 5
perm_n = 1000
theta = 0.5
alpha = 1
buffer = 500
top_n = 5
p1 = 0.05
p2 = 0.05
p3 = 0.05
nSeed = 1234
MORE_results_mono1  <- scMore(single_cell = pbmc_10x_real_data_downsampled_2000,
                              gene_info = gene_info_mono,
                              snp_info = snp_info_mono,
                              n_targets= 5,
                              perm_n = 1000,
                              theta=0.5,
                              alpha=1,
                              top_n=5,
                              buffer=500,
                              method ='cosine',
                              nSeed=1234)

MORE_results_mono2  <- scMore(single_cell = pbmc_10x_real_data_downsampled_2000,
                             gene_info = gene_info_mono,
                             snp_info = snp_info_mono,
                             n_targets= 5,
                             perm_n = 1000,
                             theta=0.5,
                             alpha=1,
                             top_n=5,
                             buffer=500,
                             method ='average',
                             nSeed=1234)

write.csv(MORE_results_mono$scMore_trait_results,file="01_mono_count_scMORE_withJSI_2000_v5_scale_twoCelltypes_separate_nonScale_TF_average_pow1_notAddModule_monocyte_CD8T_100000.csv")
write.csv(MORE_results_mono2$scMore_trait_results,file="01_mono_count_scMORE_withJSI_2000_v5_scale_twoCelltypes_separate_nonScale_TF_average_pow1_notAddModule_monocyte_CD8T_1000_average.csv")
write.csv(MORE_results_mono1$scMore_trait_results,file="01_mono_count_scMORE_withJSI_2000_v5_scale_twoCelltypes_separate_nonScale_TF_average_pow1_notAddModule_monocyte_CD8T_1000_cosine.csv")

mono_count_average <- getEnergyScore(MORE_results_mono2$scMore_trait_results,targetCelltype = 2)
log2(mono_count_average+1)+1
mono_count_cosine <- getEnergyScore(MORE_results_mono$scMore_trait_results,targetCelltype = 2)
log2(mono_count_cosine+1)+1




###-----------------------------
###---baso count
gene_info_baso <- read.table("baso_count_processed_magma_results.genes.out",header = TRUE)
#snp_info
snp_info_baso <- read.csv("baso_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)


MORE_results_baso1  <- scMore(single_cell = pbmc_10x_real_data_downsampled_2000,
                                      gene_info = gene_info_baso,
                                      snp_info = snp_info_baso,
                                      n_targets= 5,
                                      perm_n = 1000,
                                      theta=0.5,
                                      alpha=1,
                                      top_n=5,
                                      buffer=500,
                                      method = 'cosine',
                                      nSeed=1234)

MORE_results_baso2  <- scMore(single_cell = pbmc_10x_real_data_downsampled_2000,
                                      gene_info = gene_info_baso,
                                      snp_info = snp_info_baso,
                                      n_targets= 5,
                                      perm_n = 1000,
                                      theta=0.5,
                                      alpha=1,
                                      top_n=5,
                                      buffer=500,
                                      method = 'average',
                                      nSeed=1234)


write.csv(MORE_results_baso2$scMore_trait_results,file="01_baso_count_1000_average.csv")
write.csv(MORE_results_baso1$scMore_trait_results,file="01_baso_count_1000_cosine.csv")

baso_average <- getEnergyScore(MORE_results_baso2$scMore_trait_results,targetCelltype = 2)
log2(baso_average+1)+1
baso_cosine <- getEnergyScore(MORE_results_baso1$scMore_trait_results,targetCelltype = 2)
log2(baso_cosine+1)+1


log2(baso_cosine+2)/log2(baso_average+2)
(baso_cosine+1)/(baso_average+1)


###---eosin count
gene_info_eosin <- read.table("eosin_count_processed_magma_results.genes.out",header = TRUE)
#snp_info
snp_info_eosin <- read.csv("eosin_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)


# Step 3: Get gene-level association scores
message("Step 3: Retrieving gene-level association scores...")
geneRiskScores <- getGeneScore(gene_info_eosin)
snp_info <-snp_info_eosin


target_scores <- suppressWarnings(getSpecificity(single_cell,method = "cosine"))
# Step 4: Identify cell type-specific regulons relevant to disease
# Integrate GRN outputs, target specificity scores, and gene association data
message("Step 4: Identifying cell type-specific regulons relevant to disease...")
regulon2disease_results_eosin1 <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  top_n = 5
)

target_scores <- suppressWarnings(getSpecificity(single_cell,method = "average"))
# Step 4: Identify cell type-specific regulons relevant to disease
# Integrate GRN outputs, target specificity scores, and gene association data
message("Step 4: Identifying cell type-specific regulons relevant to disease...")
regulon2disease_results_eosin2 <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  top_n = 5
)

write.csv(regulon2disease_results_eosin1,file="01_eosin_count_1000_cosine.csv")
write.csv(regulon2disease_results_eosin2,file="01_eosin_count_1000_average.csv")

eosin_average <- getEnergyScore(regulon2disease_results_eosin2,targetCelltype = 2)
log2(eosin_average+1)+1
eosin_cosine <- getEnergyScore(regulon2disease_results_eosin1,targetCelltype = 2)
log2(eosin_cosine+1)+1





###---HL count
gene_info_HL <- read.table("HL_count_processed_magma_results.genes.out",header = TRUE)
#snp_info
snp_info_HL <- read.csv("HL_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)



# Step 3: Get gene-level association scores
message("Step 3: Retrieving gene-level association scores...")
geneRiskScores <- getGeneScore(gene_info_HL)
snp_info <-snp_info_HL


target_scores <- suppressWarnings(getSpecificity(single_cell,method = "cosine"))
# Step 4: Identify cell type-specific regulons relevant to disease
# Integrate GRN outputs, target specificity scores, and gene association data
message("Step 4: Identifying cell type-specific regulons relevant to disease...")
regulon2disease_results_HL1 <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  top_n = 5
)

target_scores <- suppressWarnings(getSpecificity(single_cell,method = "average"))
# Step 4: Identify cell type-specific regulons relevant to disease
# Integrate GRN outputs, target specificity scores, and gene association data
message("Step 4: Identifying cell type-specific regulons relevant to disease...")
regulon2disease_results_HL2 <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  top_n = 5
)

write.csv(regulon2disease_results_HL1,file="01_HL_count_1000_cosine.csv")
write.csv(regulon2disease_results_HL2,file="01_HL_count_1000_average.csv")

HL_average <- getEnergyScore(regulon2disease_results_HL2,targetCelltype = 2)
log2(HL_average+1)+1
HL_cosine <- getEnergyScore(regulon2disease_results_HL1,targetCelltype = 2)
log2(HL_cosine+1)+1






###---MCHC count
gene_info_MCHC <- read.table("MCHC_count_processed_magma_results.genes.out",header = TRUE)
#snp_info
snp_info_MCHC <- read.csv("MCHC_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)



# Step 3: Get gene-level association scores
message("Step 3: Retrieving gene-level association scores...")
geneRiskScores <- getGeneScore(gene_info_MCHC)
snp_info <-snp_info_MCHC


target_scores <- suppressWarnings(getSpecificity(single_cell,method = "cosine"))
# Step 4: Identify cell type-specific regulons relevant to disease
# Integrate GRN outputs, target specificity scores, and gene association data
message("Step 4: Identifying cell type-specific regulons relevant to disease...")
regulon2disease_results_MCHC1 <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  top_n = 5
)

target_scores <- suppressWarnings(getSpecificity(single_cell,method = "average"))
# Step 4: Identify cell type-specific regulons relevant to disease
# Integrate GRN outputs, target specificity scores, and gene association data
message("Step 4: Identifying cell type-specific regulons relevant to disease...")
regulon2disease_results_MCHC2 <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  top_n = 5
)

write.csv(regulon2disease_results_MCHC1,file="01_MCHC_count_1000_cosine.csv")
write.csv(regulon2disease_results_MCHC2,file="01_MCHC_count_1000_average.csv")

MCHC_average <- getEnergyScore(regulon2disease_results_MCHC2,targetCelltype = 2)
log2(MCHC_average+1)+1
MCHC_cosine <- getEnergyScore(regulon2disease_results_MCHC1,targetCelltype = 2)
log2(MCHC_cosine+1)+1





###---MCV count
gene_info_MCV <- read.table("MCV_count_processed_magma_results.genes.out",header = TRUE)
#snp_info
snp_info_MCV <- read.csv("MCV_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)


# Step 3: Get gene-level association scores
message("Step 3: Retrieving gene-level association scores...")
geneRiskScores <- getGeneScore(gene_info_MCV)
snp_info <-snp_info_MCV


target_scores <- suppressWarnings(getSpecificity(single_cell,method = "cosine"))
# Step 4: Identify cell type-specific regulons relevant to disease
# Integrate GRN outputs, target specificity scores, and gene association data
message("Step 4: Identifying cell type-specific regulons relevant to disease...")
regulon2disease_results_MCV1 <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  top_n = 5
)

target_scores <- suppressWarnings(getSpecificity(single_cell,method = "average"))
# Step 4: Identify cell type-specific regulons relevant to disease
# Integrate GRN outputs, target specificity scores, and gene association data
message("Step 4: Identifying cell type-specific regulons relevant to disease...")
regulon2disease_results_MCV2 <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  top_n = 5
)

write.csv(regulon2disease_results_MCV1,file="01_MCV_count_1000_cosine.csv")
write.csv(regulon2disease_results_MCV2,file="01_MCV_count_1000_average.csv")

MCV_average <- getEnergyScore(regulon2disease_results_MCV2,targetCelltype = 2)
log2(MCV_average+1)+1
MCV_cosine <- getEnergyScore(regulon2disease_results_MCV1,targetCelltype = 2)
log2(MCV_cosine+1)+1






###---neutr count
gene_info_neutr <- read.table("neutr_count_processed_magma_results.genes.out",header = TRUE)
#snp_info
snp_info_neutr <- read.csv("neutr_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)


# Step 3: Get gene-level association scores
message("Step 3: Retrieving gene-level association scores...")
geneRiskScores <- getGeneScore(gene_info_neutr)
snp_info <-snp_info_neutr


target_scores <- suppressWarnings(getSpecificity(single_cell,method = "cosine"))
# Step 4: Identify cell type-specific regulons relevant to disease
# Integrate GRN outputs, target specificity scores, and gene association data
message("Step 4: Identifying cell type-specific regulons relevant to disease...")
regulon2disease_results_neutr1 <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  top_n = 5
)

target_scores <- suppressWarnings(getSpecificity(single_cell,method = "average"))
# Step 4: Identify cell type-specific regulons relevant to disease
# Integrate GRN outputs, target specificity scores, and gene association data
message("Step 4: Identifying cell type-specific regulons relevant to disease...")
regulon2disease_results_neutr2 <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  top_n = 5
)

write.csv(regulon2disease_results_neutr1,file="01_neutr_count_1000_cosine.csv")
write.csv(regulon2disease_results_neutr2,file="01_neutr_count_1000_average.csv")

neutr_average <- getEnergyScore(regulon2disease_results_neutr2,targetCelltype = 2)
log2(neutr_average+1)+1
neutr_cosine <- getEnergyScore(regulon2disease_results_neutr1,targetCelltype = 2)
log2(neutr_cosine+1)+1




###---wbc count
gene_info_wbc <- read.table("wbc_count_processed_magma_results.genes.out",header = TRUE)
#snp_info
snp_info_wbc <- read.csv("wbc_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)


# Step 3: Get gene-level association scores
message("Step 3: Retrieving gene-level association scores...")
geneRiskScores <- getGeneScore(gene_info_wbc)
snp_info <-snp_info_wbc


target_scores <- suppressWarnings(getSpecificity(single_cell,method = "cosine"))
# Step 4: Identify cell type-specific regulons relevant to disease
# Integrate GRN outputs, target specificity scores, and gene association data
message("Step 4: Identifying cell type-specific regulons relevant to disease...")
regulon2disease_results_wbc1 <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  top_n = 5
)

target_scores <- suppressWarnings(getSpecificity(single_cell,method = "average"))
# Step 4: Identify cell type-specific regulons relevant to disease
# Integrate GRN outputs, target specificity scores, and gene association data
message("Step 4: Identifying cell type-specific regulons relevant to disease...")
regulon2disease_results_wbc2 <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  top_n = 5
)

write.csv(regulon2disease_results_wbc1,file="01_wbc_count_1000_cosine.csv")
write.csv(regulon2disease_results_wbc2,file="01_wbc_count_1000_average.csv")

wbc_average <- getEnergyScore(regulon2disease_results_wbc2,targetCelltype = 2)
log2(wbc_average+1)+1
wbc_cosine <- getEnergyScore(regulon2disease_results_wbc1,targetCelltype = 2)
log2(wbc_cosine+1)+1






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


###----------Figure 2相关-----heatmap--correlation of 10 blood cell traits

setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")

consine_average10traits_heatmap <- read.table("01_10blood_cell_traits_heatmap.txt",header = TRUE,sep = "\t")
cor(consine_average10traits_heatmap$lymp_count,consine_average10traits_heatmap$lymp_percent)
cor.test(consine_average10traits_heatmap$lymp_count,consine_average10traits_heatmap$lymp_percent)


library(pheatmap)


data_h1 <- consine_average10traits_heatmap[,c(-1,-2,-3,-4)]
rownames(data_h1) <- consine_average10traits_heatmap[,1]

da <- cor(data_h1)

pearson_heatmap<-pheatmap(da,cellwidth =15, cellheight =15,fontsize=8,
                          border_color="black",scale="none",
                          fontsize_row=8,
                          color=colorRampPalette(rev(c("#E77A77", "#ECBA84", "#ECBA84","white","lightblue",'#336699')))(102),
                          cluster_rows = T,
                          cluster_cols = T)



cor(consine_average10traits_heatmap$lymp_count,consine_average10traits_heatmap$lymp_percent)
cor.test(consine_average10traits_heatmap$lymp_count,consine_average10traits_heatmap$lymp_percent)


##-2
consine_average10traits_heatmap <- read.table("01_10blood_cell_traits_sig_heatmap.txt",header = TRUE,sep = "\t")

library(pheatmap)

data_h1 <- consine_average10traits_heatmap[,c(-1)]
rownames(data_h1) <- consine_average10traits_heatmap[,1]
pearson_heatmap<-pheatmap(data_h1,cellwidth =15, cellheight =15,fontsize=8,
                          border_color="black",scale="row",
                          fontsize_row=8,
                          color=colorRampPalette(rev(c("red","#E77A77", "#ECBA84","white","lightblue",'#336699')))(102),
                          cluster_rows = T,
                          cluster_cols = T)

consine_average10traits_heatmap <- read.table("01_10blood_cell_traits_sig_heatmap.txt",header = TRUE,sep = "\t")
data_h1 <- consine_average10traits_heatmap[,c(2,4,6,8,10,12,14,16,18,20)]
rownames(data_h1) <- consine_average10traits_heatmap[,1]
pearson_heatmap<-pheatmap(data_h1,cellwidth =15, cellheight =15,fontsize=8,
                          border_color="black",scale="row",
                          fontsize_row=8,
                          color=colorRampPalette(rev(c("red","#E77A77", "#ECBA84","white","lightblue",'#336699')))(102),
                          cluster_rows = T,
                          cluster_cols = T)

##-3
consine_average10traits_heatmap <- read.table("01_10blood_cell_traits_sig_heatmap2_transform.txt",header = TRUE,sep = "\t")

library(pheatmap)

data_h1 <- consine_average10traits_heatmap[,c(-1,-2,-3)]
rownames(data_h1) <- consine_average10traits_heatmap[,1]
data_h1<-t(data_h1)
pearson_heatmap<-pheatmap(data_h1,cellwidth =15, cellheight =15,fontsize=8,
                          border_color="black",scale="none",
                          fontsize_row=8,
                          color=colorRampPalette(rev(c("red","#E77A77", "#ECBA84","white","lightblue",'#336699')))(102),
                          cluster_rows = T,
                          cluster_cols = T)







###----------significant TFs---------regulon feature plot---------Figure 2D and Sup Figure S2C

# Load necessary libraries
library(Seurat)
library(ggplot2)
library(patchwork) # For combining plots
#---load pbmc_10x_real_data_downsampled_2000
pbmc_10x_real_data_downsampled_2000 <- readRDS("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/10X_PBMC/10X_PBMC_downsample_2000cells.rds")

DimPlot(pbmc_10x_real_data_downsampled_2000, reduction = "umap.new",cols = c("#67ADB7","#E77A77"))


single_cell <-pbmc_10x_real_data_downsampled_2000


# Load necessary libraries
library(Seurat)
library(ggplot2)
library(patchwork) # For combining plots

#Default use 'Signac' method for analysis
grn_outputs <- createRegulon(single_cell,
                             n_targets = 5,
                             peak2gene_method="Signac",
                             infer_method = 'glm')


regulons <- grn_outputs$grn

#1
Module <- regulons[regulons$TF == "ZEB1", ]
ZEB1_regulon <- c(unique(Module$TF), unique(Module$Target))
# Add module score for ZEB1 regulon
single_cell <- AddModuleScore(
  object = single_cell,
  features = list(ZEB1_regulon),  # Use only the genes present in the object
  name = "ZEB1_regulon"
)


FeaturePlot(object = single_cell, reduction = "umap.new", features = "ZEB1_regulon1")
library(RColorBrewer)
FeaturePlot(object = single_cell, reduction = "umap.new", features = "ZEB1_regulon1") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu")))

# 使用 FeaturePlot 并自定义颜色梯度
FeaturePlot(
  object = single_cell,
  reduction = "umap.new",
  features = "ZEB1_regulon1"
) +
  scale_colour_gradientn(
    colours = c("gray", "white", "red")  # 自定义颜色梯度
  ) +
  theme_classic()  # 可选：更改图形主题

# 使用 FeaturePlot 并自定义颜色梯度
FeaturePlot(
  object = single_cell,
  reduction = "umap.new",
  features = "ZEB1"
) +
  scale_colour_gradientn(
    colours = c("gray", "white", "red")  # 自定义颜色梯度
  ) +
  theme_classic()  # 可选：更改图形主题


#2
Module <- regulons[regulons$TF == "FOXO1", ]
FOXO1_regulon <- c(unique(Module$TF), unique(Module$Target))
# Add module score for FOXO1 regulon
single_cell <- AddModuleScore(
  object = single_cell,
  features = list(FOXO1_regulon),  # Use only the genes present in the object
  name = "FOXO1_regulon"
)

library(RColorBrewer)
FeaturePlot(object = single_cell, reduction = "umap.new", features = "FOXO1_regulon1") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu")))

# 使用 FeaturePlot 并自定义颜色梯度
FeaturePlot(
  object = single_cell,
  reduction = "umap.new",
  features = "FOXO1_regulon1"
) +
  scale_colour_gradientn(
    colours = c("gray", "white", "red")  # 自定义颜色梯度
  ) +
  theme_classic()  # 可选：更改图形主题

# 使用 FeaturePlot 并自定义颜色梯度
FeaturePlot(
  object = single_cell,
  reduction = "umap.new",
  features = "FOXO1"
) +
  scale_colour_gradientn(
    colours = c("gray", "white", "red")  # 自定义颜色梯度
  ) +
  theme_classic()  # 可选：更改图形主题


#3
Module <- regulons[regulons$TF == "IKZF1", ]
IKZF1_regulon <- c(unique(Module$TF), unique(Module$Target))
# Add module score for IKZF1 regulon
single_cell <- AddModuleScore(
  object = single_cell,
  features = list(IKZF1_regulon),  # Use only the genes present in the object
  name = "IKZF1_regulon"
)

library(RColorBrewer)
FeaturePlot(object = single_cell, reduction = "umap.new", features = "IKZF1_regulon1") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu")))

FeaturePlot(object = single_cell, reduction = "umap.new", features = "IKZF1") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu")))

# 使用 FeaturePlot 并自定义颜色梯度
FeaturePlot(
  object = single_cell,
  reduction = "umap.new",
  features = "IKZF1"
) +
  scale_colour_gradientn(
    colours = c("gray", "white", "red")  # 自定义颜色梯度
  ) +
  theme_classic()  # 可选：更改图形主题

#4
Module <- regulons[regulons$TF == "ZNF721", ]
ZNF721_regulon <- c(unique(Module$TF), unique(Module$Target))
# Add module score for ZNF721 regulon
single_cell <- AddModuleScore(
  object = single_cell,
  features = list(ZNF721_regulon),  # Use only the genes present in the object
  name = "ZNF721_regulon"
)

library(RColorBrewer)
FeaturePlot(object = single_cell, reduction = "umap.new", features = "ZNF721_regulon1") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu")))

FeaturePlot(object = single_cell, reduction = "umap.new", features = "ZNF721") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu")))

# 使用 FeaturePlot 并自定义颜色梯度
FeaturePlot(
  object = single_cell,
  reduction = "umap.new",
  features = "ZNF721"
) +
  scale_colour_gradientn(
    colours = c("gray", "white", "red")  # 自定义颜色梯度
  ) +
  theme_classic()  # 可选：更改图形主题

# 使用 FeaturePlot 并自定义颜色梯度
FeaturePlot(
  object = single_cell,
  reduction = "umap.new",
  features = "ZNF721_regulon1"
) +
  scale_colour_gradientn(
    colours = c("gray", "white", "red")  # 自定义颜色梯度
  ) +
  theme_classic()  # 可选：更改图形主题


#5
Module <- regulons[regulons$TF == "ETS1", ]
ETS1_regulon <- c(unique(Module$TF), unique(Module$Target))
# Add module score for FOXO1 regulon
single_cell <- AddModuleScore(
  object = single_cell,
  features = list(ETS1_regulon),  # Use only the genes present in the object
  name = "ETS1_regulon"
)

library(RColorBrewer)
FeaturePlot(object = single_cell, reduction = "umap.new", features = "ETS1_regulon1") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu")))

# 使用 FeaturePlot 并自定义颜色梯度
FeaturePlot(
  object = single_cell,
  reduction = "umap.new",
  features = "ETS1_regulon1"
) +
  scale_colour_gradientn(
    colours = c("gray", "white", "red")  # 自定义颜色梯度
  ) +
  theme_classic()  # 可选：更改图形主题

# 使用 FeaturePlot 并自定义颜色梯度
FeaturePlot(
  object = single_cell,
  reduction = "umap.new",
  features = "ETS1"
) +
  scale_colour_gradientn(
    colours = c("gray", "white", "red")  # 自定义颜色梯度
  ) +
  theme_classic()  # 可选：更改图形主题


#6
Module <- regulons[regulons$TF == "FOXP1", ]
FOXP1_regulon <- c(unique(Module$TF), unique(Module$Target))
# Add module score for FOXO1 regulon
single_cell <- AddModuleScore(
  object = single_cell,
  features = list(FOXP1_regulon),  # Use only the genes present in the object
  name = "FOXP1_regulon"
)

library(RColorBrewer)
FeaturePlot(object = single_cell, reduction = "umap.new", features = "FOXP1_regulon1") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu")))

# 使用 FeaturePlot 并自定义颜色梯度
FeaturePlot(
  object = single_cell,
  reduction = "umap.new",
  features = "FOXP1_regulon1"
) +
  scale_colour_gradientn(
    colours = c("gray", "white", "red")  # 自定义颜色梯度
  ) +
  theme_classic()  # 可选：更改图形主题

# 使用 FeaturePlot 并自定义颜色梯度
FeaturePlot(
  object = single_cell,
  reduction = "umap.new",
  features = "FOXP1"
) +
  scale_colour_gradientn(
    colours = c("gray", "white", "red")  # 自定义颜色梯度
  ) +
  theme_classic()  # 可选：更改图形主题





#7
Module <- regulons[regulons$TF == "LEF1", ]
LEF1_regulon <- c(unique(Module$TF), unique(Module$Target))
# Add module score for FOXO1 regulon
single_cell <- AddModuleScore(
  object = single_cell,
  features = list(LEF1_regulon),  # Use only the genes present in the object
  name = "LEF1_regulon"
)

library(RColorBrewer)
FeaturePlot(object = single_cell, reduction = "umap.new", features = "LEF1_regulon1") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu")))

# 使用 FeaturePlot 并自定义颜色梯度
FeaturePlot(
  object = single_cell,
  reduction = "umap.new",
  features = "LEF1_regulon1"
) +
  scale_colour_gradientn(
    colours = c("gray", "white", "red")  # 自定义颜色梯度
  ) +
  theme_classic()  # 可选：更改图形主题

# 使用 FeaturePlot 并自定义颜色梯度
FeaturePlot(
  object = single_cell,
  reduction = "umap.new",
  features = "LEF1"
) +
  scale_colour_gradientn(
    colours = c("gray", "white", "red")  # 自定义颜色梯度
  ) +
  theme_classic()  # 可选：更改图形主题





###------------------------------------------------------------------------------------
#2024-11-24
###------------------------------------------------------------------------------------

library(COSG)

#load data on 10x pbmc example data
#Load scMultiomic data
pbmc_10x <- readRDS("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/10X_PBMC/10X_PBMC.rds")
single_cell_10x <- pbmc_10x
Idents(single_cell_10x) <- single_cell_10x$cell_type

#single_cell$cell_type2[which(single_cell$cell_type2 %in% c("CD14+ Monocytes","FCGR3A+ Monocytes"))] <- "Monocytes"

cell_type2<- as.character(single_cell_10x$cell_type)
cell_type2[which(cell_type2 %in% c("CD14+ Monocytes","FCGR3A+ Monocytes"))] <- "Monocytes"
single_cell_10x$cell_type2 <- cell_type2

Idents(single_cell_10x) <- single_cell_10x$cell_type2


#subset cells of two cell types
pbmc_10x_real_data <- subset(single_cell,idents=c("Monocytes","NK cells"))
#downsample cell list
cell.list <- WhichCells(pbmc_10x_real_data,idents = c("Monocytes","NK cells"),downsample = 1000)
pbmc_10x_real_data_downsampled_2000 <- pbmc_10x_real_data[,cell.list]
table(pbmc_10x_real_data_downsampled_2000$cell_type2)
pbmc_10x_real_data_downsampled_2000 #this real dataset used for assessing ctDRTF performance



#subset cells of two cell types
pbmc_10x_real_data <- subset(single_cell_10x,idents=c("Monocytes","CD8+ T cells"))
#downsample cell list
cell.list <- WhichCells(pbmc_10x_real_data,idents = c("Monocytes","CD8+ T cells"),downsample = 1000)
pbmc_10x_real_data_downsampled_2000 <- pbmc_10x_real_data[,cell.list]
table(pbmc_10x_real_data_downsampled_2000$cell_type2)
pbmc_10x_real_data_downsampled_2000 #this real dataset used for assessing ctDRTF performance


#DimPlot(pbmc_10x_real_data_downsampled_2000)

#re-clustering the UMAP
pbmc_10x_real_data_downsampled_2000 <- NormalizeData(pbmc_10x_real_data_downsampled_2000)
pbmc_10x_real_data_downsampled_2000 <- FindVariableFeatures(pbmc_10x_real_data_downsampled_2000)
pbmc_10x_real_data_downsampled_2000 <- ScaleData(pbmc_10x_real_data_downsampled_2000)
pbmc_10x_real_data_downsampled_2000 <- RunPCA(pbmc_10x_real_data_downsampled_2000,reduction.name = "pca")
pbmc_10x_real_data_downsampled_2000 <- FindNeighbors(pbmc_10x_real_data_downsampled_2000,
                                                     dims = 1:30,reduction = "pca")
pbmc_10x_real_data_downsampled_2000 <- FindClusters(pbmc_10x_real_data_downsampled_2000,
                                                    resolution = 0.5)
#Run UMAP
pbmc_10x_real_data_downsampled_2000 <- RunUMAP(pbmc_10x_real_data_downsampled_2000, dims = 1:30,
                                               reduction = "pca",reduction.name = "umap.new")

#Visualization
Idents(pbmc_10x_real_data_downsampled_2000) <- pbmc_10x_real_data_downsampled_2000$cell_type2

DimPlot(pbmc_10x_real_data_downsampled_2000, reduction = "umap.new",cols = c("#67ADB7","#E77A77"))

single_cell<-pbmc_10x_real_data_downsampled_2000

#saveRDS(pbmc_10x_real_data_downsampled_2000,file = "/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/10X_PBMC/10X_PBMC_downsample_2000cells.rds")

#---load pbmc_10x_real_data_downsampled_2000

pbmc_10x_real_data_downsampled_2000 <- readRDS("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/10X_PBMC/10X_PBMC_downsample_2000cells.rds")
single_cell <-pbmc_10x_real_data_downsampled_2000
# run 'createRegulon()' function
#grn_outputs <- createRegulon(single_cell = pbmc_10x_real_data_downsampled_2000)
grn_outputs_2000_2 <- createRegulon(single_cell = pbmc_10x_real_data_downsampled_2000,
                                      peak2gene_method ='Signac',
                                      infer_method = 'glm',
                                      n_targets = 5)
length(unique(grn_outputs_2000_2$grn$TF))

grn_outputs_2000_2$grn$Target[which(grn_outputs_2000_2$grn$TF=="ARID3A")]

#Final verion for 2000 cells: grn_outputs_2000_2
grn_outputs<-grn_outputs_2000_2




#load single-cell data---1600 cells
pbmc_10x_real_data_downsampled_1600 <- readRDS("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/ctDRTF_assessment_results/pbmc_10x_real_data_downsampled_1600.rds")

Idents(pbmc_10x_real_data_downsampled_1600) <- pbmc_10x_real_data_downsampled_1600$cell_type2
DimPlot(pbmc_10x_real_data_downsampled_1600, reduction = "umap.new",cols = c("#67ADB7","#E77A77"))

single_cell<-pbmc_10x_real_data_downsampled_1600

target_scores <- getSpecificity(single_cell = pbmc_10x_real_data_downsampled_1600)


grn_outputs_1600 <- createRegulon(single_cell = pbmc_10x_real_data_downsampled_1600,
                                  peak2gene_method ='GREAT',
                                  infer_method = 'glm')

grn_outputs_1600_2_g <- createRegulon(single_cell = pbmc_10x_real_data_downsampled_1600,
                                  peak2gene_method ='GREAT',
                                  infer_method = 'glm',
                                  n_targets = 5)

grn_outputs_4443 <- createRegulon(single_cell = pbmc_10x_real_data,
                                  peak2gene_method ='GREAT',
                                  infer_method = 'glm')


grn_outputs_1600_2_2 <- createRegulon(single_cell = pbmc_10x_real_data_downsampled_1600,
                                      peak2gene_method ='Signac',
                                      infer_method = 'glm',
                                      n_targets = 10)

grn_outputs_1600_2_1 <- createRegulon(single_cell = pbmc_10x_real_data_downsampled_1600,
                                      peak2gene_method ='Signac',
                                      infer_method = 'glm',
                                      n_targets = 5)
#final version:

grn_outputs_1600_2_1$grn$Target[which(grn_outputs_1600_2_1$grn$TF=="ARID3A")]

grn_outputs<-grn_outputs_1600_2_1


#
grn_outputs1600 <- readRDS("/users/mayunlong/Desktop/grn_outputs_1600_monocyte_NK_cells.rds")





regulons <- grn_outputs_1600$grn
tf_list <- grn_outputs_1600$tf_names

regulons <- regulons[,c(1,2)]



#magma_GWAS_monocyte
#Monocyte count GWAS



##MAGMA-results
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")

##---Test results from MAGMA--------------------
#gene_info
magma_mono <- read.table("mono_count_processed_magma_results.genes.out",header = TRUE)
#magma_results <- getGeneScore(magma_mono)
#head(magma_results)
gene_info <- magma_mono
geneRiskScores <- getGeneScore(gene_info)

#snp_info
snp_info <- read.csv("monocyte_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

target_scores <- getSpecificity(pbmc_10x_real_data_downsampled_2000)


final_results <- regulon2disease(grn_outputs_2000_2,
                                      target_scores,
                                      snp_info,
                                      geneRiskScores,
                                      perm_n = 100,
                                      pow=1)
write.csv(final_results,file="mono_count_scMORE_withJSI_1600_v5_scale_twoCelltypes_separate_Scale_zscores.csv")
write.csv(final_results,file="mono_count_scMORE_withJSI_1600_v5_scale_twoCelltypes_separate_nonScale_zscores.csv")
write.csv(final_results,file="mono_count_scMORE_withJSI_1600_v5_scale_twoCelltypes_separate_Scale_zscores_use_realzscore_and_specificity.csv")
write.csv(final_results,file="mono_count_scMORE_withJSI_1600_v5_scale_twoCelltypes_separate_Scale_TF_average.csv")
write.csv(final_results,file="mono_count_scMORE_withJSI_1600_v5_scale_twoCelltypes_separate_nonScale_TF_average.csv")
write.csv(final_results,file="mono_count_scMORE_withJSI_1600_v5_scale_twoCelltypes_separate_nonScale_TF_average_pow0.0000001.csv")
write.csv(final_results,file="mono_count_scMORE_withJSI_2000_v5_scale_twoCelltypes_separate_nonScale_TF_average_pow0.0000001.csv")
write.csv(final_results,file="mono_count_scMORE_withJSI_2000_v5_scale_twoCelltypes_separate_nonScale_TF_average.csv")
write.csv(final_results,file="mono_count_scMORE_withJSI_2000_v5_scale_twoCelltypes_separate_nonScale_TF_sum_addmodel.csv")
write.csv(final_results,file="mono_count_scMORE_withJSI_2000_v5_scale_twoCelltypes_separate_nonScale_TF_average_addmodel.csv")
write.csv(final_results,file="mono_count_scMORE_withJSI_2000_v5_scale_twoCelltypes_separate_nonScale_TF_average_addmodel_pow0.0000001.csv")
write.csv(final_results,file="mono_count_scMORE_withJSI_2000_v5_scale_twoCelltypes_separate_nonScale_TF_average_pow1_notAddModule.csv")

#----
write.csv(final_results,file="mono_count_scMORE_withJSI_2000_v5_scale_twoCelltypes_separate_nonScale_TF_average_pow1_notAddModule_monocyte_CD8T.csv")



write.csv(final_results,file="mono_count_scMORE_withJSI_1600_3.csv")

write.csv(final_results,file="mono_count_scMORE_withoutJSI_1600_3.csv")


magma_lymp <- read.table("lymp_count_processed_magma_results.genes.out",header = TRUE)
#magma_results <- getGeneScore(magma_mono)
#head(magma_results)
geneRiskScores <- getGeneScore(magma_lymp)

#snp_info
snp_info_lymp <- read.csv("lymphocyte_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)
snp_info<- snp_info_lymp


final_results_lymp <- regulon2disease(grn_outputs_2000_2,
                                 target_scores,
                                 snp_info,
                                 geneRiskScores,
                                 perm_n = 100,
                                 pow=1)

write.csv(final_results_lymp,file="lymp_count_scMORE_withJSI_2000_v5_scale_twoCelltypes_separate_nonScale_TF_average_pow1_notAddModule.csv")

write.csv(final_results_lymp,file="lymp_count_scMORE_withJSI_1600_v1.csv")
write.csv(final_results_lymp,file="lymp_count_scMORE_withJSI_1600_v2_scale_twoCelltypes_separate.csv")
write.csv(final_results_lymp,file="lymp_count_scMORE_withJSI_1600_v2_scale_twoCelltypes_separate_nonScale_zscores.csv")
write.csv(final_results_lymp,file="lymp_count_scMORE_withJSI_1600_v3_scale_twoCelltypes_separate_nonScale_zscores.csv")
write.csv(final_results_lymp,file="lymp_count_scMORE_withJSI_1600_v4_scale_twoCelltypes_separate_nonScale_zscores.csv")
write.csv(final_results_lymp,file="lymp_count_scMORE_withJSI_1600_v5_scale_twoCelltypes_separate_nonScale_zscores.csv")
write.csv(final_results_lymp,file="lymp_count_scMORE_withJSI_1600_v5_scale_twoCelltypes_separate_nonScale_TF_average.csv")
#----
write.csv(final_results_lymp,file="lymp_count_scMORE_withJSI_2000_v5_scale_twoCelltypes_separate_nonScale_TF_average_pow1_notAddModule_monocyte_CD8T.csv")



library(ctDRTF)

results <- ctdrtf(single_cell = pbmc_10x_real_data_downsampled_1600,
                  MAGMA_GWAS_data = geneRiskScores)



gene_info_lymp <- read.table("baso_count_processed_magma_results.genes.out",header = TRUE)
#snp_info
snp_info_lymp <- read.csv("baso_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)
geneRiskScores <- getGeneScore(gene_info_lymp)
snp_info<- snp_info_lymp


gene_info_eosin <- read.table("eosin_count_processed_magma_results.genes.out",header = TRUE)
#snp_info
snp_info_eosin <- read.csv("eosin_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

geneRiskScores <- getGeneScore(gene_info_eosin)
snp_info<- snp_info_eosin


gene_info_lymp_percent <- read.table("lymp_percent_processed_magma_results.genes.out",header = TRUE)
#snp_info
snp_info_lymp_percent <- read.csv("lymphocyte_percent_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

geneRiskScores <- getGeneScore(gene_info_lymp_percent)
snp_info<- snp_info_lymp_percent


final_results <- regulon2disease(grn_outputs_2000_2,
                                      target_scores,
                                      snp_info,
                                      geneRiskScores,
                                      perm_n = 100,
                                      pow=1,
                                      mo=1)
write.csv(final_results,file="lymp_percent_scMORE_withJSI_1600_v5_scale_twoCelltypes_separate_nonScale_zscores.csv")
write.csv(final_results,file="lymp_percent_scMORE_withJSI_1600_v5_scale_twoCelltypes_separate_nonScale_TF_average_pow1_v2.csv")
write.csv(final_results,file="eosin_count_scMORE_withJSI_1600_v5_scale_twoCelltypes_separate_nonScale_TF_average_pow1_v2.csv")
write.csv(final_results,file="baso_count_scMORE_withJSI_1600_v5_scale_twoCelltypes_separate_nonScale_TF_average_pow1_v2.csv")


write.csv(final_results,file="eosin_count_scMORE_withJSI_2000_v5_scale_twoCelltypes_separate_nonScale_TF_average_pow1_v2--------1.csv")
write.csv(final_results,file="lymp_percent_scMORE_withJSI_2000_v5_scale_twoCelltypes_separate_nonScale_TF_average_pow1_v2--------1.csv")

write.csv(final_results,file="baso_percent_scMORE_withJSI_2000_v5_scale_twoCelltypes_separate_nonScale_TF_average_pow1_v2--------1.csv")

write.csv(final_results,file="eosin_count_scMORE_withJSI_2000_v5_scale_twoCelltypes_separate_nonScale_TF_average_pow1_notAddModule_monocyte_CD8T.csv")



n_targets= 5
top_genes = 500
perm_n = 1000
theta=0.5
pow=1
mo=1
buffer=50
nSeed=1234

final_results <- scMore(single_cell = pbmc_10x_real_data_downsampled_2000,snp_info,gene_info)

final_results_withoutJSI <- scMore(single_cell = pbmc_10x_real_data_downsampled_1600,snp_info,gene_info)
write.csv(final_results_withoutJSI,file="mono_count_scMORE_withoutJSI_1600.csv")

final_results_withJSI <- scMore(single_cell = pbmc_10x_real_data_downsampled_1600,snp_info,gene_info)
write.csv(final_results_withJSI,file="mono_count_scMORE_withJSI_1600.csv")




##--------------提取----GRN all information---##-------------##-------------##-------------##-------------
regulons <- modules

# Step 7: Extract regulatory network data
data_regulons <- data.frame(
  TF = regulons@meta$tf,                # Transcription factors
  Target = regulons@meta$target,        # Target genes
  Regions = regulons@meta$regions,      # Regulatory regions
  Pval = regulons@meta$pval             # Peak-gene association significance
)

# Step 8: Filter regulons with at least `n_genes` target genes
tf_gene_counts <- data.frame(table(data_regulons$TF))
valid_tfs <- tf_gene_counts$Var1[tf_gene_counts$Freq >= 10]


# Filter GRN to include only valid TFs
filtered_regulons <- data_regulons[data_regulons$TF %in% valid_tfs, ]

# link peaks to target


# Step 9: Prepare outputs
grn_outputs <- list(
  grn = filtered_regulons,              # Filtered GRN data
  tf_names = as.vector(valid_tfs)       # List of TFs with valid regulons

)

##-------------##-------------##-------------##-------------##-------------##-------------##-------------





#' @param gene_info Gene-based genetic association results based on MAGMA or FUMA
#'                  scMORE support the format of results from MAGMA and FUMA

##-----getGeneScore() 函数---逐步分析
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")

##---Test results from MAGMA--------------------
#gene_info
magma_mono <- read.table("mono_count_processed_magma_results.genes.out",header = TRUE)

magma_lymp_c <- read.table("lymp_count_processed_magma_results.genes.out",header = TRUE)
magma_lymp_p <- read.table("lymp_percent_processed_magma_results.genes.out",header = TRUE)

magma_baso <- read.table("baso_count_processed_magma_results.genes.out",header = TRUE)

magma_results <- getGeneScore(magma_mono)
head(magma_results)


#----Test results from FUMA
data_monocyte <- read.table("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/GWAS-data/MAGMA_result_fuma/ieu-b-31.vcf/magma.genes.out",header = TRUE)


fuma_results <- getGeneScore(data_monocyte)
head(fuma_results)


magmatop500 <- magma_results$SYMBOL[1:500]
fumatop500 <- fuma_results$SYMBOL[1:500]

overlap_size <- length(intersect(magmatop500,fumatop500))

write.csv(fuma_results,file="fuma_monocyte_results.csv")
write.csv(magma_results,file="magma_monocyte_results.csv")


##-----previous version
#BiocManager::install("org.Hs.eg.db")

library(org.Hs.eg.db)

# 示例 Entrez Gene ID 列表
entrez_ids <- c(79501, 100996442, 729759)

# 获取基因符号
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = as.character(entrez_ids),
                       column = "SYMBOL",
                       keytype = "ENTREZID",
                       multiVals = "first")

gene_symbols <- as.data.frame(gene_symbols)


#@' 1) Calculating the P value for each regulon-disease link

setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")




if ('SYMBOL' %in% colnames(data_monocyte)){

  data_monocyte <- data_monocyte %>% mutate(logP = -log10(P)) %>% arrange(desc(logP))
  magma_GWAS_monocyte <- data_monocyte[,c("SYMBOL","logP","ZSTAT")]


} else {

  magma <- read.table("mono_count_processed_magma_results.genes.out",header = TRUE)

  entrez_ids <- magma$GENE

  # 获取基因符号
  gene_symbols <- mapIds(org.Hs.eg.db,
                         keys = as.character(entrez_ids),
                         column = "SYMBOL",
                         keytype = "ENTREZID",
                         multiVals = "first")

  gene_symbols <- as.data.frame(gene_symbols)

  magma$SYMBOL <- gene_symbols$gene_symbols

  #remove NA
  magma_clean <- magma %>%
    filter(!is.na(SYMBOL))

  magma_clean <- magma_clean %>% mutate(logP = -log10(P)) %>% arrange(desc(logP))
  magma_clean_final <- magma_clean[,c("SYMBOL","logP","ZSTAT")]

}



snp_data <- read.table("snp.bed", header = FALSE, stringsAsFactors = FALSE)
colnames(snp_data) <- c("chr", "start", "end", "snp_id","pval","beta")





##Pando-generated peaks
peak_ids <- modules@meta$regions

expanded_peaks <- unlist(strsplit(peak_ids, split = ";"))

peak_data <- do.call(rbind, lapply(expanded_peaks, function(peak) {
  parts <- unlist(strsplit(peak, split = "-"))
  data.frame(chr = parts[1], start = as.integer(parts[2]), end = as.integer(parts[3]), stringsAsFactors = FALSE)
}))

peak_data$peak_id <- expanded_peaks
head(peak_data)


# create GenomicRanges project
snp_ranges <- GRanges(seqnames = snp_data$chr,
                      ranges = IRanges(start = snp_data$start, end = snp_data$end),
                      snp_id = snp_data$snp_id)

peak_ranges <- GRanges(seqnames = peak_data$chr,
                       ranges = IRanges(start = peak_data$start, end = peak_data$end),
                       peak_id = peak_data$peak_id)


#Optionally, we use 'buffer' parameter to extend the peak range
#extend peak range (upstream and downstream + 1000bp)
buffer=2000
extended_peak_ranges <- resize(peak_ranges, width = width(peak_ranges) + buffer, fix = "center")

# make sure peak range not surpass chromosome border (i.e., start < 1)
start(extended_peak_ranges) <- pmax(start(extended_peak_ranges), 1)





# Overlap SNPs to peaks
overlaps <- findOverlaps(snp_ranges, extended_peak_ranges)

# extract overlap info
snp_peak_map <- data.frame(
  snp_id = mcols(snp_ranges)$snp_id[queryHits(overlaps)],
  peak_id = mcols(extended_peak_ranges)$peak_id[subjectHits(overlaps)],
  stringsAsFactors = FALSE
)


# 安装 igraph 包
if (!requireNamespace("igraph", quietly = TRUE)) {
  install.packages("igraph")
}

library(igraph)

# 创建调控网络（TF-gene）
edges <- data.frame(
  TF = c("TF1", "TF1", "TF2", "TF3", "TF3", "TF3", "TF3"),
  Gene = c("GeneA", "GeneB", "GeneD", "GeneF", "GeneG", "GeneH", "GeneI")
)

# 创建网络图
g <- graph_from_data_frame(edges, directed = TRUE)

# 计算 Degree 中心性
degree_centrality <- degree(g, mode = "out")

# 计算 Betweenness 中心性
betweenness_centrality <- betweenness(g)

# 查看结果
cat("Degree Centrality:\n")
print(degree_centrality)

cat("\nBetweenness Centrality:\n")
print(betweenness_centrality)








####------------regulon2disease() 函数逐步分析

#Name all regulons
regulon_names <- data.frame(paste(tf_list,"_regulon",sep = ""))
#tf_list: the vector of all TFs
#Specificity score of each regulon in a particular cell type

#Collecting specificity scores
Final_regulon_score <- data.frame(matrix(ncol = 1, nrow = length(tf_list)))
names(Final_regulon_score)<-"ID_regulon"
Final_regulon_score$ID_regulon <- tf_list

#Collecting specificity*JSI for each regulon-disease link
Final_regulon_MORE_score<- data.frame(matrix(ncol = 1, nrow = length(tf_list)))
names(Final_regulon_MORE_score)<-"ID_regulon"
Final_regulon_MORE_score$ID_regulon <- tf_list

#Collecting MC P value for each regulon-disease link
Final_regulon_MORE_perm_p<- data.frame(matrix(ncol = 1, nrow = length(tf_list)))
names(Final_regulon_MORE_perm_p)<-"ID_regulon"
Final_regulon_MORE_perm_p$ID_regulon <- tf_list

#Obtain all names of the cell types
all_celltype_names <- unique(target_scores[,3])


#max-min scaling the specificity score

for (celltype in all_celltype_names ){

  target_scores$scores[which(target_scores$celltypes == celltype)] <- max_min_scale(target_scores$scores[which(target_scores$celltypes == celltype)])

}

"
   Using the COSR method to calculate the specificity scores (S) of TF-regulons;

   Using disease-specific genes and regulons weighted by S to identify cell type-specific TFs relevant to disease

   "
#Open progress bar
pb <- txtProgressBar(style=3)
start_time <- Sys.time() ##record start time


total_run <- length(all_celltype_names)*length(tf_list)
count <- 0
#COSR and JSI interaction analysis
for (i in 1:length(all_celltype_names)){
  regulon_MORE_score <- c()
  regulon_MORE_perm_p <- c()
  regulon_score_all <-c()

  for (j in 1:length(tf_list)){

    #extracting the gene list of each regulon
    Module <- regulons[which(regulons$TF == tf_list[j]),]
    Module_regulon <- c(unique(Module$TF),Module$Target)


    #all specificty score and z score of all regulon genes
    target_scores_sub <- target_scores[which(target_scores[,3] == all_celltype_names[i]),]

    #annotating magma z-score
    temp <-genetics_data[which(genetics_data$SYMBOL %in% target_scores_sub$genes),]
    temp <- temp[!duplicated(temp[,c("SYMBOL")]),]

    #overlap target_scores regulon genes with magma genes
    target_scores_sub<- target_scores_sub[which(target_scores_sub$genes %in%temp$SYMBOL),]

    #match() function
    #data_set for all regulon genes specificity and z scores
    target_scores_sub$magma_zscore <- temp$ZSTAT[match(target_scores_sub$genes, temp$SYMBOL )]


    #extracting the specificity score of each regulon
    each_module_score <- target_scores_sub[!is.na(match(target_scores_sub[,1], Module_regulon)),]
    #each_module_score$anno <- rep("Gene",n_num,length(each_module_score[,1]))
    each_module_score$anno <- rep("Gene", length(each_module_score[,1]))
    each_module_score$anno[which(each_module_score[,1] == Module_regulon[1])] <- "TF"

    #annotation MAGMA z-score
    #each_module_score<- each_module_score[which(each_module_score$genes %in%genetics_data$SYMBOL),]
    #each_module_score$magma_zscore <-genetics_data$ZSTAT[which(genetics_data$SYMBOL %in% each_module_score$genes)]


    #Calculating the module specificity score for TF in each regulon
    #tf_s_z <- each_module_score[,c("adj_score","magma_zscore")][which(each_module_score[,1] == Module_regulon[1]),]
    tf_Score_Zscore <- each_module_score[,c("scores","magma_zscore")][which(each_module_score[,1] == Module_regulon[1]),]
    tf_combined_score <- as.numeric((tf_Score_Zscore[1])^pow*(tf_Score_Zscore[2])^mo)

    if(is.na(tf_combined_score)){

      gene_Score_Zscore <- each_module_score[,c("scores","magma_zscore")][which(each_module_score[,1] != Module_regulon[1]),]
      gene_combined_score <- as.numeric((gene_Score_Zscore[,1])^pow*(gene_Score_Zscore[,2])^mo)
      average_score <- sum(gene_combined_score)/(length(gene_combined_score)+1)

      tf_combined_score <- 0

      #theta = 0.5  #theta range from 0.1 ~ 1, default set to 0.5
      regulon_score <- as.numeric(tf_combined_score) + as.numeric(theta*average_score) #regulon-specific score for each cell type


    } else{

      #Calculating the module specificity score for genes in each regulon
      #gene_s_z <- each_module_score[,c("adj_score","magma_zscore")][which(each_module_score[,1] != Module_regulon[1]),]
      gene_Score_Zscore <- each_module_score[,c("scores","magma_zscore")][which(each_module_score[,1] != Module_regulon[1]),]
      gene_combined_score <- as.numeric((gene_Score_Zscore[,1])^pow*(gene_Score_Zscore[,2])^mo)
      average_score <- mean(gene_combined_score)

      #theta = 0.5  #theta range from 0.1 ~ 1, default set to 0.5
      regulon_score <- as.numeric(tf_combined_score) + as.numeric(theta*average_score) #regulon-specific score for each cell type

    }

    #Sum
    regulon_score_all <- c(regulon_score_all,regulon_score)

    #Calculating the Jaccard Similarity Index (JSI)
    top_ranked_genes <-genetics_data$SYMBOL[1:top_genes]
    inter_genes_n <- length(intersect(top_ranked_genes,Module_regulon))
    union_genes_n <- length(union(top_ranked_genes,Module_regulon))
    JSI_score <- inter_genes_n/union_genes_n # Jaccard similarity index

    #Interaction: specificity*JSI for each regulon-disease link
    MORE_score <- regulon_score*JSI_score

    #print(paste0("Regulon ",tf_list[j]," ctDRTF score is: ", MORE_score, sep=""))


    #MC simulation for random specificity*JSI scores
    #Function: permutation()
    #perm_n = 1000
    tf_list_1 <- tf_list[which(tf_list!=tf_list[j])] #removing the targeted TF as controls
    len_of_regulon = length(Module_regulon)
    all_genes <-genetics_data$SYMBOL

    perm_results <- replicate(perm_n,permutation(target_scores_sub,
                                               tf_list_1,
                                               len_of_regulon,
                                               all_genes,
                                               top_genes,
                                               theta,
                                               pow,
                                               mo))

    #Calculating the MC p-values
    perm_p <- (1+length(perm_results[perm_results> MORE_score]))/(1+length(perm_results))


    #Running
    print(paste("Running the regulon of ",tf_list[j], " for the cell type of ",all_celltype_names[i],sep = ""))

    #Running percent:
    count=count+1
    completed_percent <- count/total_run
    print(sprintf('Completed percent: %1.2f%%',100*completed_percent))
    #Real-time progress bar
    #print(paste("Runing percent: ",percent((i+j)/(length(all_celltype_names)*length(tf_list))),sep = ""))
    setTxtProgressBar(pb,(count)/total_run)


    #Saving results
    regulon_MORE_score <- c(regulon_MORE_score,MORE_score)
    regulon_MORE_perm_p <- c(regulon_MORE_perm_p,perm_p)


  }

  #Collecting specificity scores
  regulon_score_all<- as.data.frame(regulon_score_all)
  names(regulon_score_all) <- all_celltype_names[i]
  Final_regulon_score <- cbind(Final_regulon_score,regulon_score_all)

  #Collecting MC P values
  regulon_MORE_perm_p<- as.data.frame(regulon_MORE_perm_p)
  names(regulon_MORE_perm_p) <- all_celltype_names[i]
  Final_regulon_MORE_perm_p <- cbind(Final_regulon_MORE_perm_p,regulon_MORE_perm_p)

  #Normalization

  regulon_MORE_score_norm <- (regulon_MORE_score - mean(regulon_MORE_score))/sd(regulon_MORE_score)

  #Alternative normalized method
  #max-min normalization
  #regulon_ct_score_norm <- max_min_scale(regulon_ct_score)

  #Collecting specificity*JSI for each regulon-disease link
  regulon_MORE_score_norm <- as.data.frame(regulon_MORE_score_norm)
  names(regulon_MORE_score_norm) <- all_celltype_names[i]
  Final_regulon_MORE_score <- cbind(Final_regulon_MORE_score,regulon_MORE_score_norm)

}

##Record end time
end_time <- Sys.time()

#Close progress bar
close(pb)

#Calculating the running time
run_time <- end_time - start_time
print(paste("Running time: ",run_time),sep = "")

#Outputs
#out_results <- list(ctDRTF_score = Final_regulon_ct_score,
#                     MC_p = Final_regulon_ct_mc_p,
#                      regulon_specificity_s=Final_regulon_s)
#
MORE_results <- list(MORE_score = Final_regulon_MORE_score,
                    MORE_perm_Pval = Final_regulon_MORE_perm_p)













###-----createRegulon
createRegulon <- function(single_cell, n_genes = 10) {

  # Step 1: Select variable features from single-cell data
  single_cell <- Seurat::FindVariableFeatures(single_cell, assay = "RNA")

  # Step 2: Initiate GRN object and select candidate regions
  single_cell <- Pando::initiate_grn(single_cell)

  # Step 3: Scan candidate regions for TF binding motifs
  single_cell <- Pando::find_motifs(
    single_cell,
    pfm = motifs,                          # Position Frequency Matrix (PFM) for motifs
    genome = BSgenome.Hsapiens.UCSC.hg38   # Human genome reference
  )

  # Step 4: Infer the gene regulatory network
  grn_object <- Pando::infer_grn(single_cell)

  # Step 5: Identify regulatory modules
  grn_object <- Pando::find_modules(grn_object)

  # Step 6: Extract GRN modules (regulons)
  regulons <- Pando::NetworkModules(grn_object)

  # Step 7: Extract regulatory network data
  data_regulons <- data.frame(
    TF = regulons@meta$tf,                # Transcription factors
    Target = regulons@meta$target,        # Target genes
    Regions = regulons@meta$regions       # Regulatory regions
  )

  # Step 8: Filter regulons with at least `n_genes` target genes
  tf_gene_counts <- data.frame(table(data_regulons$TF))
  valid_tfs <- tf_gene_counts$Var1[tf_gene_counts$Freq >= n_genes]


  # Filter GRN to include only valid TFs
  filtered_regulons <- data_regulons[data_regulons$TF %in% valid_tfs, ]

  # link peaks to target


  # Step 9: Prepare outputs
  grn_outputs <- list(
    grn = filtered_regulons,              # Filtered GRN data
    tf_names = as.vector(valid_tfs)       # List of TFs with valid regulons

  )

  return(grn_outputs)
}











