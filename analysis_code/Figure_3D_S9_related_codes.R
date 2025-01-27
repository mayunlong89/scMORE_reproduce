#2025-01-26

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
  library(scMORE)
})




#load data on 10x pbmc example data
#Load scMultiomic data
pbmc_10x <- readRDS("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/10X_PBMC/10X_PBMC.rds")
single_cell_8900<- pbmc_10x
Idents(single_cell_8900) <- single_cell_8900$cell_type

#single_cell$cell_type2[which(single_cell$cell_type2 %in% c("CD14+ Monocytes","FCGR3A+ Monocytes"))] <- "Monocytes"

cell_type2<- as.character(single_cell_8900$cell_type)
cell_type2[which(cell_type2 %in% c("CD14+ Monocytes","FCGR3A+ Monocytes"))] <- "Monocytes"
single_cell_8900$cell_type2 <- cell_type2

Idents(single_cell_8900) <- single_cell_8900$cell_type2 


##-----------------------------7000 cells
#subset cells of two cell types
#pbmc_10x_real_data <- subset(single_cell_8900,idents=c("Monocytes","CD8+ T cells"))
#downsample cell list
cell.list <- WhichCells(single_cell_8900,downsample = 2000)
pbmc_10x_real_data_downsampled_7000 <- single_cell_8900[,cell.list]
table(pbmc_10x_real_data_downsampled_7000$cell_type2)
pbmc_10x_real_data_downsampled_7000 #this real dataset used for assessing ctDRTF performance

DimPlot(pbmc_10x_real_data_downsampled_7000)



##-----------------------------5000 cells
#subset cells of two cell types
#pbmc_10x_real_data <- subset(single_cell_8900,idents=c("Monocytes","CD8+ T cells"))
#downsample cell list
cell.list <- WhichCells(single_cell_8900,downsample = 1000)
pbmc_10x_real_data_downsampled_5000 <- single_cell_8900[,cell.list]
table(pbmc_10x_real_data_downsampled_5000$cell_type2)
pbmc_10x_real_data_downsampled_5000 #this real dataset used for assessing ctDRTF performance

DimPlot(pbmc_10x_real_data_downsampled_5000)




##-----------------------------3000 cells
#subset cells of two cell types
#pbmc_10x_real_data <- subset(single_cell_8900,idents=c("Monocytes","CD8+ T cells"))
#downsample cell list
cell.list <- WhichCells(single_cell_8900,downsample = 500)
pbmc_10x_real_data_downsampled_3000 <- single_cell_8900[,cell.list]
table(pbmc_10x_real_data_downsampled_3000$cell_type2)
pbmc_10x_real_data_downsampled_3000 #this real dataset used for assessing ctDRTF performance

DimPlot(pbmc_10x_real_data_downsampled_3000)




##-----------------------------1000 cells
#subset cells of two cell types
cell.list <- WhichCells(single_cell_8900,downsample = 200)
pbmc_10x_real_data_downsampled_1000 <- single_cell_8900[,cell.list]
table(pbmc_10x_real_data_downsampled_1000$cell_type2)
pbmc_10x_real_data_downsampled_1000 #this real dataset used for assessing ctDRTF performance

DimPlot(pbmc_10x_real_data_downsampled_1000)






##--scMORE for lymphocyte count-----------------------------------------------------------------------------------------
##MAGMA-results of lymphocyte count
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
gene_info_lymp_count  <- read.table("lymp_count_processed_magma_results.genes.out",header = TRUE)
snp_info_lymp_count <- read.csv("lymphocyte_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

geneRiskScores <- getGeneScore(gene_info_lymp_count)
snp_info <-snp_info_lymp_count


###---Lymphocyte percent
gene_info_lymp_percent <- read.table("lymp_percent_processed_magma_results.genes.out",header = TRUE)
#snp_info
snp_info_lymp_percent <- read.csv("lymphocyte_percent_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

geneRiskScores <- getGeneScore(gene_info_lymp_percent)
snp_info <-snp_info_lymp_percent


###-----Monocyte count
gene_info_mono <- read.table("mono_count_processed_magma_results.genes.out",header = TRUE)
#snp_info
snp_info_mono <- read.csv("monocyte_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

geneRiskScores <- getGeneScore(gene_info_mono)
snp_info <-snp_info_mono

###---baso count
gene_info_baso <- read.table("baso_count_processed_magma_results.genes.out",header = TRUE)
#snp_info
snp_info_baso <- read.csv("baso_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

geneRiskScores <- getGeneScore(gene_info_baso)
snp_info <-snp_info_baso


###---eosin count
gene_info_eosin <- read.table("eosin_count_processed_magma_results.genes.out",header = TRUE)
#snp_info
snp_info_eosin <- read.csv("eosin_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)
geneRiskScores <- getGeneScore(gene_info_eosin)
snp_info <-snp_info_eosin 



###---HL count
gene_info_HL <- read.table("HL_count_processed_magma_results.genes.out",header = TRUE)
#snp_info
snp_info_HL <- read.csv("HL_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

geneRiskScores <- getGeneScore(gene_info_HL)
snp_info <-snp_info_HL



###---MCHC count
gene_info_MCHC <- read.table("MCHC_count_processed_magma_results.genes.out",header = TRUE)
#snp_info
snp_info_MCHC <- read.csv("MCHC_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

geneRiskScores <- getGeneScore(gene_info_MCHC)
snp_info <-snp_info_MCHC


###---MCV count
gene_info_MCV <- read.table("MCV_count_processed_magma_results.genes.out",header = TRUE)
#snp_info
snp_info_MCV <- read.csv("MCV_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

geneRiskScores <- getGeneScore(gene_info_MCV)
snp_info <-snp_info_MCV



###---neutr count
gene_info_neutr <- read.table("neutr_count_processed_magma_results.genes.out",header = TRUE)
#snp_info
snp_info_neutr <- read.csv("neutr_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

geneRiskScores <- getGeneScore(gene_info_neutr)
snp_info <-snp_info_neutr


###---wbc count
gene_info_wbc <- read.table("wbc_count_processed_magma_results.genes.out",header = TRUE)
#snp_info
snp_info_wbc <- read.csv("wbc_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

geneRiskScores <- getGeneScore(gene_info_wbc)
snp_info <-snp_info_wbc





##----1000cells
#----generate regulons from single-cell multiomics data
n_targets=5
#grn_outputs1 <- createRegulon(pbmc_10x_real_data_downsampled_1000, n_targets)

target_scores1 <- suppressWarnings(getSpecificity(pbmc_10x_real_data_downsampled_1000,method = "cosine"))

target_scores1 <- suppressWarnings(getSpecificity(pbmc_10x_real_data_downsampled_1000,method = "average"))

target_scores1 <- suppressWarnings(getSpecificity(pbmc_10x_real_data_downsampled_1000,method = "mean_specificity"))


regulon2disease_results1 <- regulon2disease(
  grn_outputs = grn_outputs1,
  target_scores = target_scores1,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 10,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  #infer_method = 'glmnet',
  #infer_method = 'xgb',
  infer_method = 'glm',
  top_n = 5
)

#E-statistics
lymp_count_cosine1 <- getEnergyScore(regulon2disease_results1,targetCelltype = 2)
log2(lymp_count_cosine1+1)+1





##----3000 cells
#######------------
#----generate regulons from single-cell multiomics data
n_targets=5
#grn_outputs3 <- createRegulon(pbmc_10x_real_data_downsampled_3000, n_targets)

target_scores3 <- suppressWarnings(getSpecificity(pbmc_10x_real_data_downsampled_3000,method = "cosine"))

target_scores3 <- suppressWarnings(getSpecificity(pbmc_10x_real_data_downsampled_3000,method = "average"))

target_scores3 <- suppressWarnings(getSpecificity(pbmc_10x_real_data_downsampled_3000,method = "mean_specificity"))


regulon2disease_results3 <- regulon2disease(
  grn_outputs = grn_outputs3,
  target_scores = target_scores3,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 10,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  #infer_method = 'glmnet',
  #infer_method = 'xgb',
  infer_method = 'glm',
  top_n = 5
)

#E-statistics
lymp_count_cosine3 <- getEnergyScore(regulon2disease_results3,targetCelltype = 2)
log2(lymp_count_cosine3+1)+1






##----5000cells
#----generate regulons from single-cell multiomics data
n_targets=5
#grn_outputs5 <- createRegulon(pbmc_10x_real_data_downsampled_5000, n_targets)

target_scores5 <- suppressWarnings(getSpecificity(pbmc_10x_real_data_downsampled_5000,method = "cosine"))

target_scores5 <- suppressWarnings(getSpecificity(pbmc_10x_real_data_downsampled_5000,method = "average"))

target_scores5 <- suppressWarnings(getSpecificity(pbmc_10x_real_data_downsampled_5000,method = "mean_specificity"))


regulon2disease_results5 <- regulon2disease(
  grn_outputs = grn_outputs5,
  target_scores = target_scores5,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 10,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  #infer_method = 'glmnet',
  #infer_method = 'xgb',
  infer_method = 'glm',
  top_n = 5
)

#E-statistics
lymp_count_cosine5 <- getEnergyScore(regulon2disease_results5,targetCelltype = 2)
log2(lymp_count_cosine5+1)+1







##----7000cells
#----generate regulons from single-cell multiomics data
n_targets=5
#grn_outputs7 <- createRegulon(pbmc_10x_real_data_downsampled_7000, n_targets)

target_scores7 <- suppressWarnings(getSpecificity(pbmc_10x_real_data_downsampled_7000,method = "cosine"))

target_scores7 <- suppressWarnings(getSpecificity(pbmc_10x_real_data_downsampled_7000,method = "average"))

target_scores7 <- suppressWarnings(getSpecificity(pbmc_10x_real_data_downsampled_7000,method = "mean_specificity"))


regulon2disease_results7 <- regulon2disease(
  grn_outputs = grn_outputs7,
  target_scores = target_scores7,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 10,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  #infer_method = 'glmnet',
  #infer_method = 'xgb',
  infer_method = 'glm',
  top_n = 5
)

#E-statistics
lymp_count_cosine7 <- getEnergyScore(regulon2disease_results7,targetCelltype = 2)
log2(lymp_count_cosine7+1)+1






##----9000cells
#----generate regulons from single-cell multiomics data
n_targets=5
#grn_outputs9 <- createRegulon(single_cell_8900, n_targets)

target_scores9 <- suppressWarnings(getSpecificity(single_cell_8900,method = "cosine"))

target_scores9 <- suppressWarnings(getSpecificity(single_cell_8900,method = "average"))

target_scores9 <- suppressWarnings(getSpecificity(single_cell_8900,method = "mean_specificity"))


regulon2disease_results9 <- regulon2disease(
  grn_outputs = grn_outputs9,
  target_scores = target_scores9,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 10,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  #infer_method = 'glmnet',
  #infer_method = 'xgb',
  infer_method = 'glm',
  top_n = 5
)

#E-statistics
lymp_count_cosine9 <- getEnergyScore(regulon2disease_results9,targetCelltype = 2)
log2(lymp_count_cosine9+1)+1






# Required Libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Example Data Structure
# Load your data (modify the file path and data structure based on your data)
data <- read.csv("01_figure3D_data.csv", header = TRUE)

# Ensure that 'Method' and 'Cellcount' are factors for proper grouping
data$Method <- factor(data$Method, levels = c( "Average", "MAGMA_Celltyping","scMORE"))
data$Cellcount <- factor(data$Cellcount, levels = c("1K", "3K", "5K", "7K", "9K"))

# Create the Bar Plot
ggplot(data, aes(x = Cellcount, y = Score, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  facet_wrap(~ Trait, scales = "free_y") +  # Separate plots for each trait
  scale_fill_manual(values = c("#E9E4A6","#E9B78A","#F07B52")) +  # Custom colors
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
  scale_fill_manual(values = c("#E9E4A6","#E9B78A","#F07B52")) +  # Custom colors
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




