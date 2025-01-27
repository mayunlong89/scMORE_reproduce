#2025-01-27


##----five autoimmune diseases


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

pal1 <- c("#85C89C","#67ADB7","#36600E","#6A8473","#C0BFDF",
          "#E77A77","#A13B46","#7E6148FF","#ECBA84","#CA8C74")


pal <- c("#A13B46","#67ADB7","#36600E","#6A8473","#C0BFDF","#7E6148FF",
         "#85C89C")

DimPlot(single_cell_8900,reduction = "umap.biMod",label = T, cols=pal, pt.size = 0.1, repel = T)




##-----IBD
gene_info_IBD <- read.table("IBD_processed_magma_results.genes.out",header = TRUE)
#snp_info
snp_info_IBD <- read.csv("IBD_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

gene_info<- gene_info_IBD 
snp_info <- snp_info_IBD


grn_outputs_8900 <- createRegulon(single_cell = single_cell_8900,
                                      peak2gene_method ='Signac',
                                      infer_method = 'glm',
                                      n_targets = 5)


target_scores9 <- suppressWarnings(getSpecificity(single_cell_8900,method = "cosine"))

target_scores9 <- suppressWarnings(getSpecificity(single_cell_8900,method = "average"))

target_scores9 <- suppressWarnings(getSpecificity(single_cell_8900,method = "mean_specificity"))


regulon2disease_results9 <- regulon2disease(
  grn_outputs = grn_outputs_8900 ,
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





##-----PBC
gene_info_PBC <- read.table("PBC_count_GWAS_final.hg19_SNP_Gene_Analysis_results.genes.out",header = TRUE)
#snp_info
snp_info_PBC <- read.table("PBC_maf0.01.txt",header=T,stringsAsFactors = FALSE)


gene_info<- gene_info_PBC
snp_info <- snp_info_PBC

target_scores9 <- suppressWarnings(getSpecificity(single_cell_8900,method = "cosine"))

target_scores9 <- suppressWarnings(getSpecificity(single_cell_8900,method = "average"))

target_scores9 <- suppressWarnings(getSpecificity(single_cell_8900,method = "mean_specificity"))


regulon2disease_results9 <- regulon2disease(
  grn_outputs = grn_outputs_8900 ,
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





##-----RA
gene_info_RA <- read.table("RA_count_GWAS_final.hg19_SNP_Gene_Analysis_results.genes.out",header = TRUE)
#snp_info
snp_info_RA <- read.table("RA_maf0.01.txt",header=T,stringsAsFactors = FALSE)

gene_info<- gene_info_RA
snp_info <- snp_info_RA


target_scores9 <- suppressWarnings(getSpecificity(single_cell_8900,method = "cosine"))

target_scores9 <- suppressWarnings(getSpecificity(single_cell_8900,method = "average"))

target_scores9 <- suppressWarnings(getSpecificity(single_cell_8900,method = "mean_specificity"))


regulon2disease_results9 <- regulon2disease(
  grn_outputs = grn_outputs_8900 ,
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



##-----SLE
gene_info_SLE <- read.table("SLE_count_GWAS_final.hg19_SNP_Gene_Analysis_results.genes.out",header = TRUE)
#snp_info
snp_info_SLE <- read.table("SLE_maf0.01.txt",header=T,stringsAsFactors = FALSE)

gene_info<- gene_info_SLE
snp_info <- snp_info_SLE



target_scores9 <- suppressWarnings(getSpecificity(single_cell_8900,method = "cosine"))

target_scores9 <- suppressWarnings(getSpecificity(single_cell_8900,method = "average"))

target_scores9 <- suppressWarnings(getSpecificity(single_cell_8900,method = "mean_specificity"))


regulon2disease_results9 <- regulon2disease(
  grn_outputs = grn_outputs_8900 ,
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




##-----UC
gene_info_UC <- read.table("UC_processed_magma_results.genes.out",header = TRUE)
#snp_info
snp_info_UC <- read.table("UC_maf0.01.txt",header=T,stringsAsFactors = FALSE)

gene_info<- gene_info_UC
snp_info <- snp_info_UC



target_scores9 <- suppressWarnings(getSpecificity(single_cell_8900,method = "cosine"))

target_scores9 <- suppressWarnings(getSpecificity(single_cell_8900,method = "average"))

target_scores9 <- suppressWarnings(getSpecificity(single_cell_8900,method = "mean_specificity"))


regulon2disease_results9 <- regulon2disease(
  grn_outputs = grn_outputs_8900 ,
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
data <- read.csv("01_figure4B_data.csv", header = TRUE)

# Ensure that 'Method' and 'Cellcount' are factors for proper grouping
data$Method <- factor(data$Method, levels = c( "Average", "MAGMA_Celltyping","scMORE"))
#data$Cellcount <- factor(data$Cellcount, levels = c("1K", "3K", "5K", "7K", "9K"))

# Create the Bar Plot
ggplot(data, aes(x = Method, y = Score, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  facet_wrap(~ Trait, scales = "free_y") +  # Separate plots for each trait
  scale_fill_manual(values = c("#E9E4A6","#E9B78A","#F07B52")) +  # Custom colors
  labs(
    title = "Method Comparison Across Methods and Traits",
    x = "",
    y = "E-statistics",
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
  group_by(Method) %>%
  summarise(
    Mean_Score = mean(Score),
    SD_Score = sd(Score),
    SE_Score = SD_Score / sqrt(n()),  # Standard Error
    .groups = 'drop'
  )



# Create Bar Plot with Error Bars
ggplot(summary_data, aes(x = Method, y = Mean_Score, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = Mean_Score - SE_Score, ymax = Mean_Score + SE_Score),
                position = position_dodge(0.7), width = 0.2) +
  scale_fill_manual(values = c("#E9E4A6","#E9B78A","#F07B52")) +  # Custom colors
  labs(
    title = "Average Scores Across 5 autoimmune diseases by Method and Cell Count",
    x = "",
    y = "E-statistics",
    fill = "Method"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 12)
  )


#统计regulon sizes
data <- grn_outputs_8900$grn

# Group by TF and count the number of Targets, adding 1
result <- data %>%
  group_by(TF) %>%
  summarise(Target_Count = n() + 1)

write.csv(result, file="regulon_size_single_cell8900.csv",quote=F, row.names = F)
