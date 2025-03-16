
#---2025-03-12


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
library(scMORE)

##--------------------

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

n_targets=5
grn_output_9K <- createRegulon(single_cell_8900, n_targets)
target_scores1_9K <- suppressWarnings(getSpecificity(single_cell_8900,method = "cosine"))






##-----------------------------7000 cells
#subset cells of two cell types
#pbmc_10x_real_data <- subset(single_cell_8900,idents=c("Monocytes","CD8+ T cells"))
#downsample cell list
cell.list <- WhichCells(single_cell_8900,downsample = 2000)
pbmc_10x_real_data_downsampled_7000 <- single_cell_8900[,cell.list]
table(pbmc_10x_real_data_downsampled_7000$cell_type2)
pbmc_10x_real_data_downsampled_7000 #this real dataset used for assessing ctDRTF performance

DimPlot(pbmc_10x_real_data_downsampled_7000)


n_targets=5
grn_output_7K <- createRegulon(pbmc_10x_real_data_downsampled_7000, n_targets)
target_scores1_7K <- suppressWarnings(getSpecificity(pbmc_10x_real_data_downsampled_7000,method = "cosine"))



##-----------------------------5000 cells
#subset cells of two cell types
#pbmc_10x_real_data <- subset(single_cell_8900,idents=c("Monocytes","CD8+ T cells"))
#downsample cell list
cell.list <- WhichCells(single_cell_8900,downsample = 1000)
pbmc_10x_real_data_downsampled_5000 <- single_cell_8900[,cell.list]
table(pbmc_10x_real_data_downsampled_5000$cell_type2)
pbmc_10x_real_data_downsampled_5000 #this real dataset used for assessing ctDRTF performance

DimPlot(pbmc_10x_real_data_downsampled_5000)

n_targets=5
grn_output_5K <- createRegulon(pbmc_10x_real_data_downsampled_5000, n_targets)
target_scores1_5K <- suppressWarnings(getSpecificity(pbmc_10x_real_data_downsampled_5000,method = "cosine"))



##-----------------------------3000 cells
#subset cells of two cell types
#pbmc_10x_real_data <- subset(single_cell_8900,idents=c("Monocytes","CD8+ T cells"))
#downsample cell list
cell.list <- WhichCells(single_cell_8900,downsample = 500)
pbmc_10x_real_data_downsampled_3000 <- single_cell_8900[,cell.list]
table(pbmc_10x_real_data_downsampled_3000$cell_type2)
pbmc_10x_real_data_downsampled_3000 #this real dataset used for assessing ctDRTF performance

DimPlot(pbmc_10x_real_data_downsampled_3000)

n_targets=5
grn_output_3K <- createRegulon(pbmc_10x_real_data_downsampled_3000, n_targets)
target_scores1_3K <- suppressWarnings(getSpecificity(pbmc_10x_real_data_downsampled_3000,method = "cosine"))


##-----------------------------1000 cells
#subset cells of two cell types
cell.list <- WhichCells(single_cell_8900,downsample = 200)
pbmc_10x_real_data_downsampled_1000 <- single_cell_8900[,cell.list]
table(pbmc_10x_real_data_downsampled_1000$cell_type2)
pbmc_10x_real_data_downsampled_1000 #this real dataset used for assessing ctDRTF performance

DimPlot(pbmc_10x_real_data_downsampled_1000)


n_targets=5
grn_output_1K <- createRegulon(pbmc_10x_real_data_downsampled_1000, n_targets)
target_scores1_1K <- suppressWarnings(getSpecificity(pbmc_10x_real_data_downsampled_1000,method = "cosine"))





########-------------------------------------------------------------------------------------
########-------------------------------------------------------------------------------------
########-------------------------------------------------------------------------------------
########-------------------------------------------------------------------------------------
########-------------------------------------------------------------------------------------
########-------------------------------------------------------------------------------------
########-------------------------------------------------------------------------------------
########-------------------------------------------------------------------------------------
########-------------------------------------------------------------------------------------
########-------------------------------------------------------------------------------------



####------fastBAT assessment-------9K

###baso_fastBAT
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/Method_envalue/")

baso_fastBAT <- read.table("./fastBAT/baso_count_maf0.01.txt.gene.fastbat",header = T)

gene_info <- baso_fastBAT[,c(1,8,9)]
#head(gene_info)
colnames(gene_info)<- c("SYMBOL","ZSTAT","P")

geneRiskScores_9K <- getGeneScore(gene_info)

#snp_info
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
snp_info <- read.csv("baso_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

regulon2disease_results1 <- regulon2disease(
  grn_outputs = grn_output_9K,
  target_scores = target_scores1_9K,
  geneRiskScores = geneRiskScores_9K,
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

#E-statistics---9K

baso_count_cosine1 <- getEnergyScore(regulon2disease_results1,targetCelltype = 2)
log2(baso_count_cosine1+1)+1



####------lymph_count_fastBAT

setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/Method_envalue/")

lymph_count_fastBAT <- read.table("./fastBAT/lymphocyte_count_maf0.01.txt.SS.gene.fastbat",header = T)

gene_info <- lymph_count_fastBAT[,c(1,8,9)]
#head(gene_info)
colnames(gene_info)<- c("SYMBOL","ZSTAT","P")

geneRiskScores_9K <- getGeneScore(gene_info)

#snp_info
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
snp_info <- read.csv("lymphocyte_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

regulon2disease_results1 <- regulon2disease(
  grn_outputs = grn_output_9K,
  target_scores = target_scores1_9K,
  geneRiskScores = geneRiskScores_9K,
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

#E-statistics---9K

lymph_count_cosine1 <- getEnergyScore(regulon2disease_results1,targetCelltype = 2)
log2(lymph_count_cosine1+1)+1




####------lymph_percent_fastBAT

setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/Method_envalue/")

lymph_percent_fastBAT <- read.table("./fastBAT/lymphocyte_percent_maf0.01.txt.SS.gene.fastbat",header = T)

gene_info <- lymph_percent_fastBAT[,c(1,8,9)]
#head(gene_info)
colnames(gene_info)<- c("SYMBOL","ZSTAT","P")

geneRiskScores_9K <- getGeneScore(gene_info)


#snp_info
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
snp_info <- read.csv("lymphocyte_percent_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

regulon2disease_results1 <- regulon2disease(
  grn_outputs = grn_output_9K,
  target_scores = target_scores1_9K,
  geneRiskScores = geneRiskScores_9K,
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

#E-statistics---9K

lymph_percent_cosine1 <- getEnergyScore(regulon2disease_results1,targetCelltype = 2)
log2(lymph_percent_cosine1+1)+1






####------esoin_fastBAT

setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/Method_envalue/")

esoin_fastBAT <- read.table("./fastBAT/eosin_count_maf0.01.txt.gene.fastbat",header = T)

gene_info <- esoin_fastBAT[,c(1,8,9)]
#head(gene_info)
colnames(gene_info)<- c("SYMBOL","ZSTAT","P")

geneRiskScores_9K <- getGeneScore(gene_info)


#snp_info
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
snp_info <- read.csv("eosin_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

regulon2disease_results1 <- regulon2disease(
  grn_outputs = grn_output_9K,
  target_scores = target_scores1_9K,
  geneRiskScores = geneRiskScores_9K,
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

#E-statistics---9K

esoin_cosine1 <- getEnergyScore(regulon2disease_results1,targetCelltype = 2)
log2(esoin_cosine1+1)+1






####------HL_fastBAT

setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/Method_envalue/")

HL_fastBAT <- read.table("./fastBAT/HL_count_maf0.01.txt.SS.gene.fastbat",header = T)

gene_info <- HL_fastBAT[,c(1,8,9)]
#head(gene_info)
colnames(gene_info)<- c("SYMBOL","ZSTAT","P")

geneRiskScores_9K <- getGeneScore(gene_info)


#snp_info
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
snp_info <- read.csv("HL_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

regulon2disease_results1 <- regulon2disease(
  grn_outputs = grn_output_9K,
  target_scores = target_scores1_9K,
  geneRiskScores = geneRiskScores_9K,
  snp_info = snp_info,
  perm_n = 10,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  infer_method = 'glm',
  top_n = 5
)

#E-statistics---9K

HL_cosine1 <- getEnergyScore(regulon2disease_results1,targetCelltype = 2)
log2(HL_cosine1+1)+1






####------MCHC_fastBAT

setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/Method_envalue/")

MCHC_fastBAT <- read.table("./fastBAT/MCHC_count_maf0.01.txt.SS.gene.fastbat",header = T)

gene_info <- MCHC_fastBAT[,c(1,8,9)]
#head(gene_info)
colnames(gene_info)<- c("SYMBOL","ZSTAT","P")

geneRiskScores_9K <- getGeneScore(gene_info)


#snp_info
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
snp_info <- read.csv("MCHC_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

regulon2disease_results1 <- regulon2disease(
  grn_outputs = grn_output_9K,
  target_scores = target_scores1_9K,
  geneRiskScores = geneRiskScores_9K,
  snp_info = snp_info,
  perm_n = 10,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  infer_method = 'glm',
  top_n = 5
)

#E-statistics---9K

MCHC_cosine1 <- getEnergyScore(regulon2disease_results1,targetCelltype = 2)
log2(MCHC_cosine1+1)+1






####------WBC_fastBAT

setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/Method_envalue/")

xx_fastBAT <- read.table("./fastBAT/wbc_count_maf0.01.txt.gene.fastbat",header = T)

gene_info <- xx_fastBAT[,c(1,8,9)]
#head(gene_info)
colnames(gene_info)<- c("SYMBOL","ZSTAT","P")

geneRiskScores_9K <- getGeneScore(gene_info)


#snp_info
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
snp_info <- read.csv("WBC_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

regulon2disease_results1 <- regulon2disease(
  grn_outputs = grn_output_9K,
  target_scores = target_scores1_9K,
  geneRiskScores = geneRiskScores_9K,
  snp_info = snp_info,
  perm_n = 10,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  infer_method = 'glm',
  top_n = 5
)

#E-statistics---9K

WBC_cosine1 <- getEnergyScore(regulon2disease_results1,targetCelltype = 2)
log2(WBC_cosine1+1)+1





####------Neutr_fastBAT

setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/Method_envalue/")

xx_fastBAT <- read.table("./fastBAT/neutr_count_maf0.01.txt.gene.fastbat",header = T)

gene_info <- xx_fastBAT[,c(1,8,9)]
#head(gene_info)
colnames(gene_info)<- c("SYMBOL","ZSTAT","P")

geneRiskScores_9K <- getGeneScore(gene_info)


#snp_info
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
snp_info <- read.csv("neutr_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

regulon2disease_results1 <- regulon2disease(
  grn_outputs = grn_output_9K,
  target_scores = target_scores1_9K,
  geneRiskScores = geneRiskScores_9K,
  snp_info = snp_info,
  perm_n = 10,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  infer_method = 'glm',
  top_n = 5
)

#E-statistics---9K

neutro_cosine1 <- getEnergyScore(regulon2disease_results1,targetCelltype = 2)
log2(neutro_cosine1+1)+1





####------MCV_fastBAT

setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/Method_envalue/")

xx_fastBAT <- read.table("./fastBAT/MCV_count_maf0.01.txt.SS.gene.fastbat",header = T)

gene_info <- xx_fastBAT[,c(1,8,9)]
#head(gene_info)
colnames(gene_info)<- c("SYMBOL","ZSTAT","P")

geneRiskScores_9K <- getGeneScore(gene_info)


#snp_info
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
snp_info <- read.csv("MCV_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

regulon2disease_results1 <- regulon2disease(
  grn_outputs = grn_output_9K,
  target_scores = target_scores1_9K,
  geneRiskScores = geneRiskScores_9K,
  snp_info = snp_info,
  perm_n = 10,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  infer_method = 'glm',
  top_n = 5
)

#E-statistics---9K

MCV_cosine1 <- getEnergyScore(regulon2disease_results1,targetCelltype = 2)
log2(MCV_cosine1+1)+1





####------Mono_fastBAT

setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/Method_envalue/")

xx_fastBAT <- read.table("./fastBAT/monocyte_count_maf0.01.txt.gene.fastbat",header = T)

gene_info <- xx_fastBAT[,c(1,8,9)]
#head(gene_info)
colnames(gene_info)<- c("SYMBOL","ZSTAT","P")

geneRiskScores_9K <- getGeneScore(gene_info)


#snp_info
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
snp_info <- read.csv("monocyte_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

regulon2disease_results1 <- regulon2disease(
  grn_outputs = grn_output_9K,
  target_scores = target_scores1_9K,
  geneRiskScores = geneRiskScores_9K,
  snp_info = snp_info,
  perm_n = 10,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  infer_method = 'glm',
  top_n = 5
)

#E-statistics---9K

mono_cosine1 <- getEnergyScore(regulon2disease_results1,targetCelltype = 2)
log2(mono_cosine1+1)+1




########-------------------------------------------------------------------------------------
########-------------------------------------------------------------------------------------
########-------------------------------------------------------------------------------------
########-------------------------------------------------------------------------------------
########-------------------------------------------------------------------------------------
########-------------------------------------------------------------------------------------
########-------------------------------------------------------------------------------------
########-------------------------------------------------------------------------------------
########-------------------------------------------------------------------------------------
########-------------------------------------------------------------------------------------





####------fastBAT assessment-------7K

###-----baso_fastBAT
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/Method_envalue/")

baso_fastBAT <- read.table("./fastBAT/baso_count_maf0.01.txt.gene.fastbat",header = T)

gene_info <- baso_fastBAT[,c(1,8,9)]
#head(gene_info)
colnames(gene_info)<- c("SYMBOL","ZSTAT","P")

geneRiskScores_7K <- getGeneScore(gene_info)

#snp_info
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
snp_info <- read.csv("baso_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

regulon2disease_results1 <- regulon2disease(
  grn_outputs = grn_output_7K,
  target_scores = target_scores1_7K,
  geneRiskScores = geneRiskScores_7K,
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

#E-statistics---

baso_count_cosine1 <- getEnergyScore(regulon2disease_results1,targetCelltype = 2)
log2(baso_count_cosine1+1)+1



####------lymph_count_fastBAT

setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/Method_envalue/")

lymph_count_fastBAT <- read.table("./fastBAT/lymphocyte_count_maf0.01.txt.SS.gene.fastbat",header = T)

gene_info <- lymph_count_fastBAT[,c(1,8,9)]
#head(gene_info)
colnames(gene_info)<- c("SYMBOL","ZSTAT","P")

geneRiskScores_7K <- getGeneScore(gene_info)

#snp_info
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
snp_info <- read.csv("lymphocyte_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

regulon2disease_results1 <- regulon2disease(
  grn_outputs = grn_output_7K,
  target_scores = target_scores1_7K,
  geneRiskScores = geneRiskScores_7K,
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

#E-statistics---

lymph_count_cosine1 <- getEnergyScore(regulon2disease_results1,targetCelltype = 2)
log2(lymph_count_cosine1+1)+1




####------lymph_percent_fastBAT

setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/Method_envalue/")

lymph_percent_fastBAT <- read.table("./fastBAT/lymphocyte_percent_maf0.01.txt.SS.gene.fastbat",header = T)

gene_info <- lymph_percent_fastBAT[,c(1,8,9)]
#head(gene_info)
colnames(gene_info)<- c("SYMBOL","ZSTAT","P")

geneRiskScores_7K <- getGeneScore(gene_info)


#snp_info
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
snp_info <- read.csv("lymphocyte_percent_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

regulon2disease_results1 <- regulon2disease(
  grn_outputs = grn_output_7K,
  target_scores = target_scores1_7K,
  geneRiskScores = geneRiskScores_7K,
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

#E-statistics---

lymph_percent_cosine1 <- getEnergyScore(regulon2disease_results1,targetCelltype = 2)
log2(lymph_percent_cosine1+1)+1






####------esoin_fastBAT

setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/Method_envalue/")

esoin_fastBAT <- read.table("./fastBAT/eosin_count_maf0.01.txt.gene.fastbat",header = T)

gene_info <- esoin_fastBAT[,c(1,8,9)]
#head(gene_info)
colnames(gene_info)<- c("SYMBOL","ZSTAT","P")

geneRiskScores_7K <- getGeneScore(gene_info)


#snp_info
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
snp_info <- read.csv("eosin_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

regulon2disease_results1 <- regulon2disease(
  grn_outputs = grn_output_7K,
  target_scores = target_scores1_7K,
  geneRiskScores = geneRiskScores_7K,
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

#E-statistics--

esoin_cosine1 <- getEnergyScore(regulon2disease_results1,targetCelltype = 2)
log2(esoin_cosine1+1)+1






####------HL_fastBAT

setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/Method_envalue/")

HL_fastBAT <- read.table("./fastBAT/HL_count_maf0.01.txt.SS.gene.fastbat",header = T)

gene_info <- HL_fastBAT[,c(1,8,9)]
#head(gene_info)
colnames(gene_info)<- c("SYMBOL","ZSTAT","P")

geneRiskScores_7K <- getGeneScore(gene_info)


#snp_info
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
snp_info <- read.csv("HL_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

regulon2disease_results1 <- regulon2disease(
  grn_outputs = grn_output_7K,
  target_scores = target_scores1_7K,
  geneRiskScores = geneRiskScores_7K,
  snp_info = snp_info,
  perm_n = 10,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  infer_method = 'glm',
  top_n = 5
)

#E-statistics--

HL_cosine1 <- getEnergyScore(regulon2disease_results1,targetCelltype = 2)
log2(HL_cosine1+1)+1






####------MCHC_fastBAT

setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/Method_envalue/")

MCHC_fastBAT <- read.table("./fastBAT/MCHC_count_maf0.01.txt.SS.gene.fastbat",header = T)

gene_info <- MCHC_fastBAT[,c(1,8,9)]
#head(gene_info)
colnames(gene_info)<- c("SYMBOL","ZSTAT","P")

geneRiskScores_7K <- getGeneScore(gene_info)


#snp_info
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
snp_info <- read.csv("MCHC_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

regulon2disease_results1 <- regulon2disease(
  grn_outputs = grn_output_7K,
  target_scores = target_scores1_7K,
  geneRiskScores = geneRiskScores_7K,
  snp_info = snp_info,
  perm_n = 10,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  infer_method = 'glm',
  top_n = 5
)

#E-statistics---

MCHC_cosine1 <- getEnergyScore(regulon2disease_results1,targetCelltype = 2)
log2(MCHC_cosine1+1)+1






####------WBC_fastBAT

setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/Method_envalue/")

xx_fastBAT <- read.table("./fastBAT/wbc_count_maf0.01.txt.gene.fastbat",header = T)

gene_info <- xx_fastBAT[,c(1,8,9)]
#head(gene_info)
colnames(gene_info)<- c("SYMBOL","ZSTAT","P")

geneRiskScores_7K <- getGeneScore(gene_info)


#snp_info
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
snp_info <- read.csv("WBC_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

regulon2disease_results1 <- regulon2disease(
  grn_outputs = grn_output_7K,
  target_scores = target_scores1_7K,
  geneRiskScores = geneRiskScores_7K,
  snp_info = snp_info,
  perm_n = 10,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  infer_method = 'glm',
  top_n = 5
)

#E-statistics---

WBC_cosine1 <- getEnergyScore(regulon2disease_results1,targetCelltype = 2)
log2(WBC_cosine1+1)+1





####------Neutr_fastBAT

setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/Method_envalue/")

xx_fastBAT <- read.table("./fastBAT/neutr_count_maf0.01.txt.gene.fastbat",header = T)

gene_info <- xx_fastBAT[,c(1,8,9)]
#head(gene_info)
colnames(gene_info)<- c("SYMBOL","ZSTAT","P")

geneRiskScores_7K <- getGeneScore(gene_info)


#snp_info
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
snp_info <- read.csv("neutr_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

regulon2disease_results1 <- regulon2disease(
  grn_outputs = grn_output_7K,
  target_scores = target_scores1_7K,
  geneRiskScores = geneRiskScores_7K,
  snp_info = snp_info,
  perm_n = 10,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  infer_method = 'glm',
  top_n = 5
)

#E-statistics---7K

neutro_cosine1 <- getEnergyScore(regulon2disease_results1,targetCelltype = 2)
log2(neutro_cosine1+1)+1





####------MCV_fastBAT

setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/Method_envalue/")

xx_fastBAT <- read.table("./fastBAT/MCV_count_maf0.01.txt.SS.gene.fastbat",header = T)

gene_info <- xx_fastBAT[,c(1,8,9)]
#head(gene_info)
colnames(gene_info)<- c("SYMBOL","ZSTAT","P")

geneRiskScores_7K <- getGeneScore(gene_info)


#snp_info
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
snp_info <- read.csv("MCV_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

regulon2disease_results1 <- regulon2disease(
  grn_outputs = grn_output_7K,
  target_scores = target_scores1_7K,
  geneRiskScores = geneRiskScores_7K,
  snp_info = snp_info,
  perm_n = 10,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  infer_method = 'glm',
  top_n = 5
)

#E-statistics---7K

MCV_cosine1 <- getEnergyScore(regulon2disease_results1,targetCelltype = 2)
log2(MCV_cosine1+1)+1





####------Mono_fastBAT

setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/Method_envalue/")

xx_fastBAT <- read.table("./fastBAT/monocyte_count_maf0.01.txt.gene.fastbat",header = T)

gene_info <- xx_fastBAT[,c(1,8,9)]
#head(gene_info)
colnames(gene_info)<- c("SYMBOL","ZSTAT","P")

geneRiskScores_7K <- getGeneScore(gene_info)


#snp_info
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
snp_info <- read.csv("monocyte_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

regulon2disease_results1 <- regulon2disease(
  grn_outputs = grn_output_7K,
  target_scores = target_scores1_7K,
  geneRiskScores = geneRiskScores_7K,
  snp_info = snp_info,
  perm_n = 10,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  infer_method = 'glm',
  top_n = 5
)

#E-statistics---7K

mono_cosine1 <- getEnergyScore(regulon2disease_results1,targetCelltype = 2)
log2(mono_cosine1+1)+1






########-------------------------------------------------------------------------------------
########-------------------------------------------------------------------------------------
########-------------------------------------------------------------------------------------
########-------------------------------------------------------------------------------------
########-------------------------------------------------------------------------------------
########-------------------------------------------------------------------------------------
########-------------------------------------------------------------------------------------
########-------------------------------------------------------------------------------------
########-------------------------------------------------------------------------------------
########-------------------------------------------------------------------------------------





####------fastBAT assessment-------5K

###-----baso_fastBAT
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/Method_envalue/")

baso_fastBAT <- read.table("./fastBAT/baso_count_maf0.01.txt.gene.fastbat",header = T)

gene_info <- baso_fastBAT[,c(1,8,9)]
#head(gene_info)
colnames(gene_info)<- c("SYMBOL","ZSTAT","P")

geneRiskScores_5K <- getGeneScore(gene_info)

#snp_info
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
snp_info <- read.csv("baso_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

regulon2disease_results1 <- regulon2disease(
  grn_outputs = grn_output_5K,
  target_scores = target_scores1_5K,
  geneRiskScores = geneRiskScores_5K,
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

#E-statistics---

baso_count_cosine1 <- getEnergyScore(regulon2disease_results1,targetCelltype = 2)
log2(baso_count_cosine1+1)+1



####------lymph_count_fastBAT

setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/Method_envalue/")

lymph_count_fastBAT <- read.table("./fastBAT/lymphocyte_count_maf0.01.txt.SS.gene.fastbat",header = T)

gene_info <- lymph_count_fastBAT[,c(1,8,9)]
#head(gene_info)
colnames(gene_info)<- c("SYMBOL","ZSTAT","P")

geneRiskScores_5K <- getGeneScore(gene_info)

#snp_info
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
snp_info <- read.csv("lymphocyte_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

regulon2disease_results1 <- regulon2disease(
  grn_outputs = grn_output_5K,
  target_scores = target_scores1_5K,
  geneRiskScores = geneRiskScores_5K,
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

#E-statistics---

lymph_count_cosine1 <- getEnergyScore(regulon2disease_results1,targetCelltype = 2)
log2(lymph_count_cosine1+1)+1




####------lymph_percent_fastBAT

setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/Method_envalue/")

lymph_percent_fastBAT <- read.table("./fastBAT/lymphocyte_percent_maf0.01.txt.SS.gene.fastbat",header = T)

gene_info <- lymph_percent_fastBAT[,c(1,8,9)]
#head(gene_info)
colnames(gene_info)<- c("SYMBOL","ZSTAT","P")

geneRiskScores_5K <- getGeneScore(gene_info)


#snp_info
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
snp_info <- read.csv("lymphocyte_percent_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

regulon2disease_results1 <- regulon2disease(
  grn_outputs = grn_output_5K,
  target_scores = target_scores1_5K,
  geneRiskScores = geneRiskScores_5K,
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

#E-statistics---

lymph_percent_cosine1 <- getEnergyScore(regulon2disease_results1,targetCelltype = 2)
log2(lymph_percent_cosine1+1)+1






####------esoin_fastBAT

setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/Method_envalue/")

esoin_fastBAT <- read.table("./fastBAT/eosin_count_maf0.01.txt.gene.fastbat",header = T)

gene_info <- esoin_fastBAT[,c(1,8,9)]
#head(gene_info)
colnames(gene_info)<- c("SYMBOL","ZSTAT","P")

geneRiskScores_5K <- getGeneScore(gene_info)


#snp_info
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
snp_info <- read.csv("eosin_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

regulon2disease_results1 <- regulon2disease(
  grn_outputs = grn_output_5K,
  target_scores = target_scores1_5K,
  geneRiskScores = geneRiskScores_5K,
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

#E-statistics--

esoin_cosine1 <- getEnergyScore(regulon2disease_results1,targetCelltype = 2)
log2(esoin_cosine1+1)+1






####------HL_fastBAT

setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/Method_envalue/")

HL_fastBAT <- read.table("./fastBAT/HL_count_maf0.01.txt.SS.gene.fastbat",header = T)

gene_info <- HL_fastBAT[,c(1,8,9)]
#head(gene_info)
colnames(gene_info)<- c("SYMBOL","ZSTAT","P")

geneRiskScores_5K <- getGeneScore(gene_info)


#snp_info
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
snp_info <- read.csv("HL_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

regulon2disease_results1 <- regulon2disease(
  grn_outputs = grn_output_5K,
  target_scores = target_scores1_5K,
  geneRiskScores = geneRiskScores_5K,
  snp_info = snp_info,
  perm_n = 10,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  infer_method = 'glm',
  top_n = 5
)

#E-statistics--

HL_cosine1 <- getEnergyScore(regulon2disease_results1,targetCelltype = 2)
log2(HL_cosine1+1)+1






####------MCHC_fastBAT

setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/Method_envalue/")

MCHC_fastBAT <- read.table("./fastBAT/MCHC_count_maf0.01.txt.SS.gene.fastbat",header = T)

gene_info <- MCHC_fastBAT[,c(1,8,9)]
#head(gene_info)
colnames(gene_info)<- c("SYMBOL","ZSTAT","P")

geneRiskScores_5K <- getGeneScore(gene_info)


#snp_info
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
snp_info <- read.csv("MCHC_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

regulon2disease_results1 <- regulon2disease(
  grn_outputs = grn_output_5K,
  target_scores = target_scores1_5K,
  geneRiskScores = geneRiskScores_5K,
  snp_info = snp_info,
  perm_n = 10,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  infer_method = 'glm',
  top_n = 5
)

#E-statistics---

MCHC_cosine1 <- getEnergyScore(regulon2disease_results1,targetCelltype = 2)
log2(MCHC_cosine1+1)+1






####------WBC_fastBAT

setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/Method_envalue/")

xx_fastBAT <- read.table("./fastBAT/wbc_count_maf0.01.txt.gene.fastbat",header = T)

gene_info <- xx_fastBAT[,c(1,8,9)]
#head(gene_info)
colnames(gene_info)<- c("SYMBOL","ZSTAT","P")

geneRiskScores_5K <- getGeneScore(gene_info)


#snp_info
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
snp_info <- read.csv("WBC_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

regulon2disease_results1 <- regulon2disease(
  grn_outputs = grn_output_5K,
  target_scores = target_scores1_5K,
  geneRiskScores = geneRiskScores_5K,
  snp_info = snp_info,
  perm_n = 10,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  infer_method = 'glm',
  top_n = 5
)

#E-statistics---

WBC_cosine1 <- getEnergyScore(regulon2disease_results1,targetCelltype = 2)
log2(WBC_cosine1+1)+1





####------Neutr_fastBAT

setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/Method_envalue/")

xx_fastBAT <- read.table("./fastBAT/neutr_count_maf0.01.txt.gene.fastbat",header = T)

gene_info <- xx_fastBAT[,c(1,8,9)]
#head(gene_info)
colnames(gene_info)<- c("SYMBOL","ZSTAT","P")

geneRiskScores_5K <- getGeneScore(gene_info)


#snp_info
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
snp_info <- read.csv("neutr_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

regulon2disease_results1 <- regulon2disease(
  grn_outputs = grn_output_5K,
  target_scores = target_scores1_5K,
  geneRiskScores = geneRiskScores_5K,
  snp_info = snp_info,
  perm_n = 10,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  infer_method = 'glm',
  top_n = 5
)

#E-statistics---

neutro_cosine1 <- getEnergyScore(regulon2disease_results1,targetCelltype = 2)
log2(neutro_cosine1+1)+1





####------MCV_fastBAT

setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/Method_envalue/")

xx_fastBAT <- read.table("./fastBAT/MCV_count_maf0.01.txt.SS.gene.fastbat",header = T)

gene_info <- xx_fastBAT[,c(1,8,9)]
#head(gene_info)
colnames(gene_info)<- c("SYMBOL","ZSTAT","P")

geneRiskScores_5K <- getGeneScore(gene_info)


#snp_info
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
snp_info <- read.csv("MCV_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

regulon2disease_results1 <- regulon2disease(
  grn_outputs = grn_output_5K,
  target_scores = target_scores1_5K,
  geneRiskScores = geneRiskScores_5K,
  snp_info = snp_info,
  perm_n = 10,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  infer_method = 'glm',
  top_n = 5
)

#E-statistics--

MCV_cosine1 <- getEnergyScore(regulon2disease_results1,targetCelltype = 2)
log2(MCV_cosine1+1)+1





####------Mono_fastBAT

setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/Method_envalue/")

xx_fastBAT <- read.table("./fastBAT/monocyte_count_maf0.01.txt.gene.fastbat",header = T)

gene_info <- xx_fastBAT[,c(1,8,9)]
#head(gene_info)
colnames(gene_info)<- c("SYMBOL","ZSTAT","P")

geneRiskScores_5K <- getGeneScore(gene_info)


#snp_info
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
snp_info <- read.csv("monocyte_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

regulon2disease_results1 <- regulon2disease(
  grn_outputs = grn_output_5K,
  target_scores = target_scores1_5K,
  geneRiskScores = geneRiskScores_5K,
  snp_info = snp_info,
  perm_n = 10,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  infer_method = 'glm',
  top_n = 5
)

#E-statistics---

mono_cosine1 <- getEnergyScore(regulon2disease_results1,targetCelltype = 2)
log2(mono_cosine1+1)+1





########-------------------------------------------------------------------------------------
########-------------------------------------------------------------------------------------
########-------------------------------------------------------------------------------------
########-------------------------------------------------------------------------------------
########-------------------------------------------------------------------------------------
########-------------------------------------------------------------------------------------
########-------------------------------------------------------------------------------------
########-------------------------------------------------------------------------------------
########-------------------------------------------------------------------------------------
########-------------------------------------------------------------------------------------





####------fastBAT assessment-------3K

###-----baso_fastBAT
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/Method_envalue/")

baso_fastBAT <- read.table("./fastBAT/baso_count_maf0.01.txt.gene.fastbat",header = T)

gene_info <- baso_fastBAT[,c(1,8,9)]
#head(gene_info)
colnames(gene_info)<- c("SYMBOL","ZSTAT","P")

geneRiskScores_3K <- getGeneScore(gene_info)

#snp_info
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
snp_info <- read.csv("baso_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

regulon2disease_results1 <- regulon2disease(
  grn_outputs = grn_output_3K,
  target_scores = target_scores1_3K,
  geneRiskScores = geneRiskScores_3K,
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

#E-statistics---

baso_count_cosine1 <- getEnergyScore(regulon2disease_results1,targetCelltype = 2)
log2(baso_count_cosine1+1)+1



####------lymph_count_fastBAT

setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/Method_envalue/")

lymph_count_fastBAT <- read.table("./fastBAT/lymphocyte_count_maf0.01.txt.SS.gene.fastbat",header = T)

gene_info <- lymph_count_fastBAT[,c(1,8,9)]
#head(gene_info)
colnames(gene_info)<- c("SYMBOL","ZSTAT","P")

geneRiskScores_3K <- getGeneScore(gene_info)

#snp_info
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
snp_info <- read.csv("lymphocyte_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

regulon2disease_results1 <- regulon2disease(
  grn_outputs = grn_output_3K,
  target_scores = target_scores1_3K,
  geneRiskScores = geneRiskScores_3K,
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

#E-statistics---

lymph_count_cosine1 <- getEnergyScore(regulon2disease_results1,targetCelltype = 2)
log2(lymph_count_cosine1+1)+1




####------lymph_percent_fastBAT

setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/Method_envalue/")

lymph_percent_fastBAT <- read.table("./fastBAT/lymphocyte_percent_maf0.01.txt.SS.gene.fastbat",header = T)

gene_info <- lymph_percent_fastBAT[,c(1,8,9)]
#head(gene_info)
colnames(gene_info)<- c("SYMBOL","ZSTAT","P")

geneRiskScores_3K <- getGeneScore(gene_info)


#snp_info
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
snp_info <- read.csv("lymphocyte_percent_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

regulon2disease_results1 <- regulon2disease(
  grn_outputs = grn_output_3K,
  target_scores = target_scores1_3K,
  geneRiskScores = geneRiskScores_3K,
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

#E-statistics---

lymph_percent_cosine1 <- getEnergyScore(regulon2disease_results1,targetCelltype = 2)
log2(lymph_percent_cosine1+1)+1






####------esoin_fastBAT

setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/Method_envalue/")

esoin_fastBAT <- read.table("./fastBAT/eosin_count_maf0.01.txt.gene.fastbat",header = T)

gene_info <- esoin_fastBAT[,c(1,8,9)]
#head(gene_info)
colnames(gene_info)<- c("SYMBOL","ZSTAT","P")

geneRiskScores_3K <- getGeneScore(gene_info)


#snp_info
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
snp_info <- read.csv("eosin_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

regulon2disease_results1 <- regulon2disease(
  grn_outputs = grn_output_3K,
  target_scores = target_scores1_3K,
  geneRiskScores = geneRiskScores_3K,
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

#E-statistics--

esoin_cosine1 <- getEnergyScore(regulon2disease_results1,targetCelltype = 2)
log2(esoin_cosine1+1)+1






####------HL_fastBAT

setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/Method_envalue/")

HL_fastBAT <- read.table("./fastBAT/HL_count_maf0.01.txt.SS.gene.fastbat",header = T)

gene_info <- HL_fastBAT[,c(1,8,9)]
#head(gene_info)
colnames(gene_info)<- c("SYMBOL","ZSTAT","P")

geneRiskScores_3K <- getGeneScore(gene_info)


#snp_info
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
snp_info <- read.csv("HL_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

regulon2disease_results1 <- regulon2disease(
  grn_outputs = grn_output_3K,
  target_scores = target_scores1_3K,
  geneRiskScores = geneRiskScores_3K,
  snp_info = snp_info,
  perm_n = 10,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  infer_method = 'glm',
  top_n = 5
)

#E-statistics--

HL_cosine1 <- getEnergyScore(regulon2disease_results1,targetCelltype = 2)
log2(HL_cosine1+1)+1






####------MCHC_fastBAT

setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/Method_envalue/")

MCHC_fastBAT <- read.table("./fastBAT/MCHC_count_maf0.01.txt.SS.gene.fastbat",header = T)

gene_info <- MCHC_fastBAT[,c(1,8,9)]
#head(gene_info)
colnames(gene_info)<- c("SYMBOL","ZSTAT","P")

geneRiskScores_3K <- getGeneScore(gene_info)


#snp_info
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
snp_info <- read.csv("MCHC_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

regulon2disease_results1 <- regulon2disease(
  grn_outputs = grn_output_3K,
  target_scores = target_scores1_3K,
  geneRiskScores = geneRiskScores_3K,
  snp_info = snp_info,
  perm_n = 10,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  infer_method = 'glm',
  top_n = 5
)

#E-statistics---

MCHC_cosine1 <- getEnergyScore(regulon2disease_results1,targetCelltype = 2)
log2(MCHC_cosine1+1)+1






####------WBC_fastBAT

setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/Method_envalue/")

xx_fastBAT <- read.table("./fastBAT/wbc_count_maf0.01.txt.gene.fastbat",header = T)

gene_info <- xx_fastBAT[,c(1,8,9)]
#head(gene_info)
colnames(gene_info)<- c("SYMBOL","ZSTAT","P")

geneRiskScores_3K <- getGeneScore(gene_info)


#snp_info
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
snp_info <- read.csv("WBC_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

regulon2disease_results1 <- regulon2disease(
  grn_outputs = grn_output_3K,
  target_scores = target_scores1_3K,
  geneRiskScores = geneRiskScores_3K,
  snp_info = snp_info,
  perm_n = 10,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  infer_method = 'glm',
  top_n = 5
)

#E-statistics---

WBC_cosine1 <- getEnergyScore(regulon2disease_results1,targetCelltype = 2)
log2(WBC_cosine1+1)+1





####------Neutr_fastBAT

setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/Method_envalue/")

xx_fastBAT <- read.table("./fastBAT/neutr_count_maf0.01.txt.gene.fastbat",header = T)

gene_info <- xx_fastBAT[,c(1,8,9)]
#head(gene_info)
colnames(gene_info)<- c("SYMBOL","ZSTAT","P")

geneRiskScores_3K <- getGeneScore(gene_info)


#snp_info
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
snp_info <- read.csv("neutr_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

regulon2disease_results1 <- regulon2disease(
  grn_outputs = grn_output_3K,
  target_scores = target_scores1_3K,
  geneRiskScores = geneRiskScores_3K,
  snp_info = snp_info,
  perm_n = 10,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  infer_method = 'glm',
  top_n = 5
)

#E-statistics---

neutro_cosine1 <- getEnergyScore(regulon2disease_results1,targetCelltype = 2)
log2(neutro_cosine1+1)+1





####------MCV_fastBAT

setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/Method_envalue/")

xx_fastBAT <- read.table("./fastBAT/MCV_count_maf0.01.txt.SS.gene.fastbat",header = T)

gene_info <- xx_fastBAT[,c(1,8,9)]
#head(gene_info)
colnames(gene_info)<- c("SYMBOL","ZSTAT","P")

geneRiskScores_3K <- getGeneScore(gene_info)


#snp_info
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
snp_info <- read.csv("MCV_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

regulon2disease_results1 <- regulon2disease(
  grn_outputs = grn_output_3K,
  target_scores = target_scores1_3K,
  geneRiskScores = geneRiskScores_3K,
  snp_info = snp_info,
  perm_n = 10,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  infer_method = 'glm',
  top_n = 5
)

#E-statistics--

MCV_cosine1 <- getEnergyScore(regulon2disease_results1,targetCelltype = 2)
log2(MCV_cosine1+1)+1





####------Mono_fastBAT

setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/Method_envalue/")

xx_fastBAT <- read.table("./fastBAT/monocyte_count_maf0.01.txt.gene.fastbat",header = T)

gene_info <- xx_fastBAT[,c(1,8,9)]
#head(gene_info)
colnames(gene_info)<- c("SYMBOL","ZSTAT","P")

geneRiskScores_3K <- getGeneScore(gene_info)


#snp_info
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
snp_info <- read.csv("monocyte_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

regulon2disease_results1 <- regulon2disease(
  grn_outputs = grn_output_3K,
  target_scores = target_scores1_3K,
  geneRiskScores = geneRiskScores_3K,
  snp_info = snp_info,
  perm_n = 10,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  infer_method = 'glm',
  top_n = 5
)

#E-statistics---

mono_cosine1 <- getEnergyScore(regulon2disease_results1,targetCelltype = 2)
log2(mono_cosine1+1)+1


######---------------------------------------------------------------------------------------
######---------------------------------------------------------------------------------------
######---------------------------------------------------------------------------------------
######---------------------------------------------------------------------------------------
######---------------------------------------------------------------------------------------
######---------------------------------------------------------------------------------------
######---------------------------------------------------------------------------------------
######---------------------------------------------------------------------------------------


####------fastBAT assessment-------1K

###-----baso_fastBAT
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/Method_envalue/")

baso_fastBAT <- read.table("./fastBAT/baso_count_maf0.01.txt.gene.fastbat",header = T)

gene_info <- baso_fastBAT[,c(1,8,9)]
#head(gene_info)
colnames(gene_info)<- c("SYMBOL","ZSTAT","P")

geneRiskScores_1K <- getGeneScore(gene_info)

#snp_info
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
snp_info <- read.csv("baso_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

regulon2disease_results1 <- regulon2disease(
  grn_outputs = grn_output_1K,
  target_scores = target_scores1_1K,
  geneRiskScores = geneRiskScores_1K,
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

#E-statistics---

baso_count_cosine1 <- getEnergyScore(regulon2disease_results1,targetCelltype = 2)
log2(baso_count_cosine1+1)+1



####------lymph_count_fastBAT

setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/Method_envalue/")

lymph_count_fastBAT <- read.table("./fastBAT/lymphocyte_count_maf0.01.txt.SS.gene.fastbat",header = T)

gene_info <- lymph_count_fastBAT[,c(1,8,9)]
#head(gene_info)
colnames(gene_info)<- c("SYMBOL","ZSTAT","P")

geneRiskScores_1K <- getGeneScore(gene_info)

#snp_info
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
snp_info <- read.csv("lymphocyte_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

regulon2disease_results1 <- regulon2disease(
  grn_outputs = grn_output_1K,
  target_scores = target_scores1_1K,
  geneRiskScores = geneRiskScores_1K,
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

#E-statistics---

lymph_count_cosine1 <- getEnergyScore(regulon2disease_results1,targetCelltype = 2)
log2(lymph_count_cosine1+1)+1




####------lymph_percent_fastBAT

setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/Method_envalue/")

lymph_percent_fastBAT <- read.table("./fastBAT/lymphocyte_percent_maf0.01.txt.SS.gene.fastbat",header = T)

gene_info <- lymph_percent_fastBAT[,c(1,8,9)]
#head(gene_info)
colnames(gene_info)<- c("SYMBOL","ZSTAT","P")

geneRiskScores_1K <- getGeneScore(gene_info)


#snp_info
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
snp_info <- read.csv("lymphocyte_percent_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

regulon2disease_results1 <- regulon2disease(
  grn_outputs = grn_output_1K,
  target_scores = target_scores1_1K,
  geneRiskScores = geneRiskScores_1K,
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

#E-statistics---

lymph_percent_cosine1 <- getEnergyScore(regulon2disease_results1,targetCelltype = 2)
log2(lymph_percent_cosine1+1)+1






####------esoin_fastBAT

setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/Method_envalue/")

esoin_fastBAT <- read.table("./fastBAT/eosin_count_maf0.01.txt.gene.fastbat",header = T)

gene_info <- esoin_fastBAT[,c(1,8,9)]
#head(gene_info)
colnames(gene_info)<- c("SYMBOL","ZSTAT","P")

geneRiskScores_1K <- getGeneScore(gene_info)


#snp_info
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
snp_info <- read.csv("eosin_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

regulon2disease_results1 <- regulon2disease(
  grn_outputs = grn_output_1K,
  target_scores = target_scores1_1K,
  geneRiskScores = geneRiskScores_1K,
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

#E-statistics--

esoin_cosine1 <- getEnergyScore(regulon2disease_results1,targetCelltype = 2)
log2(esoin_cosine1+1)+1






####------HL_fastBAT

setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/Method_envalue/")

HL_fastBAT <- read.table("./fastBAT/HL_count_maf0.01.txt.SS.gene.fastbat",header = T)

gene_info <- HL_fastBAT[,c(1,8,9)]
#head(gene_info)
colnames(gene_info)<- c("SYMBOL","ZSTAT","P")

geneRiskScores_1K <- getGeneScore(gene_info)


#snp_info
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
snp_info <- read.csv("HL_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

regulon2disease_results1 <- regulon2disease(
  grn_outputs = grn_output_1K,
  target_scores = target_scores1_1K,
  geneRiskScores = geneRiskScores_1K,
  snp_info = snp_info,
  perm_n = 10,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  infer_method = 'glm',
  top_n = 5
)

#E-statistics--

HL_cosine1 <- getEnergyScore(regulon2disease_results1,targetCelltype = 2)
log2(HL_cosine1+1)+1






####------MCHC_fastBAT

setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/Method_envalue/")

MCHC_fastBAT <- read.table("./fastBAT/MCHC_count_maf0.01.txt.SS.gene.fastbat",header = T)

gene_info <- MCHC_fastBAT[,c(1,8,9)]
#head(gene_info)
colnames(gene_info)<- c("SYMBOL","ZSTAT","P")

geneRiskScores_1K <- getGeneScore(gene_info)


#snp_info
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
snp_info <- read.csv("MCHC_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

regulon2disease_results1 <- regulon2disease(
  grn_outputs = grn_output_1K,
  target_scores = target_scores1_1K,
  geneRiskScores = geneRiskScores_1K,
  snp_info = snp_info,
  perm_n = 10,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  infer_method = 'glm',
  top_n = 5
)

#E-statistics---

MCHC_cosine1 <- getEnergyScore(regulon2disease_results1,targetCelltype = 2)
log2(MCHC_cosine1+1)+1






####------WBC_fastBAT

setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/Method_envalue/")

xx_fastBAT <- read.table("./fastBAT/wbc_count_maf0.01.txt.gene.fastbat",header = T)

gene_info <- xx_fastBAT[,c(1,8,9)]
#head(gene_info)
colnames(gene_info)<- c("SYMBOL","ZSTAT","P")

geneRiskScores_1K <- getGeneScore(gene_info)


#snp_info
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
snp_info <- read.csv("WBC_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

regulon2disease_results1 <- regulon2disease(
  grn_outputs = grn_output_1K,
  target_scores = target_scores1_1K,
  geneRiskScores = geneRiskScores_1K,
  snp_info = snp_info,
  perm_n = 10,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  infer_method = 'glm',
  top_n = 5
)

#E-statistics---

WBC_cosine1 <- getEnergyScore(regulon2disease_results1,targetCelltype = 2)
log2(WBC_cosine1+1)+1





####------Neutr_fastBAT

setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/Method_envalue/")

xx_fastBAT <- read.table("./fastBAT/neutr_count_maf0.01.txt.gene.fastbat",header = T)

gene_info <- xx_fastBAT[,c(1,8,9)]
#head(gene_info)
colnames(gene_info)<- c("SYMBOL","ZSTAT","P")

geneRiskScores_1K <- getGeneScore(gene_info)


#snp_info
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
snp_info <- read.csv("neutr_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

regulon2disease_results1 <- regulon2disease(
  grn_outputs = grn_output_1K,
  target_scores = target_scores1_1K,
  geneRiskScores = geneRiskScores_1K,
  snp_info = snp_info,
  perm_n = 10,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  infer_method = 'glm',
  top_n = 5
)

#E-statistics---

neutro_cosine1 <- getEnergyScore(regulon2disease_results1,targetCelltype = 2)
log2(neutro_cosine1+1)+1





####------MCV_fastBAT

setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/Method_envalue/")

xx_fastBAT <- read.table("./fastBAT/MCV_count_maf0.01.txt.SS.gene.fastbat",header = T)

gene_info <- xx_fastBAT[,c(1,8,9)]
#head(gene_info)
colnames(gene_info)<- c("SYMBOL","ZSTAT","P")

geneRiskScores_1K <- getGeneScore(gene_info)


#snp_info
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
snp_info <- read.csv("MCV_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

regulon2disease_results1 <- regulon2disease(
  grn_outputs = grn_output_1K,
  target_scores = target_scores1_1K,
  geneRiskScores = geneRiskScores_1K,
  snp_info = snp_info,
  perm_n = 10,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  infer_method = 'glm',
  top_n = 5
)

#E-statistics--

MCV_cosine1 <- getEnergyScore(regulon2disease_results1,targetCelltype = 2)
log2(MCV_cosine1+1)+1





####------Mono_fastBAT

setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/Method_envalue/")

xx_fastBAT <- read.table("./fastBAT/monocyte_count_maf0.01.txt.gene.fastbat",header = T)

gene_info <- xx_fastBAT[,c(1,8,9)]
#head(gene_info)
colnames(gene_info)<- c("SYMBOL","ZSTAT","P")

geneRiskScores_1K <- getGeneScore(gene_info)


#snp_info
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
snp_info <- read.csv("monocyte_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

regulon2disease_results1 <- regulon2disease(
  grn_outputs = grn_output_1K,
  target_scores = target_scores1_1K,
  geneRiskScores = geneRiskScores_1K,
  snp_info = snp_info,
  perm_n = 10,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  infer_method = 'glm',
  top_n = 5
)

#E-statistics---

mono_cosine1 <- getEnergyScore(regulon2disease_results1,targetCelltype = 2)
log2(mono_cosine1+1)+1






























