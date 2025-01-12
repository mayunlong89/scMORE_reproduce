
#-2025-01-05
#analyses
###-------use 10 blood cell traits to assess the difference of different parameters
#' regression_method = glm, glmnet, cv.glmnet, xgb


#devtools::install_github("mayunlong89/scMORE")

library(scMORE)
library(GenomicRanges)
library(IRanges)
library(Seurat)
library(Signac)


##---2000Cells

#---load pbmc_10x_real_data_downsampled_2000

pbmc_10x_real_data_downsampled_2000 <- readRDS("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/10X_PBMC/10X_PBMC_downsample_2000cells.rds")

DimPlot(pbmc_10x_real_data_downsampled_2000, reduction = "umap.new",cols = c("#67ADB7","#E77A77"))

single_cell <-pbmc_10x_real_data_downsampled_2000

cosine_1 <- getSpecificity(single_cell, method = 'cosine')
head(cosine_1)

average_1 <- getSpecificity(single_cell, method = 'average')
head(average_1)

mean_specificity_1 <- getSpecificity(single_cell, method = 'mean_specificity')
head(mean_specificity_1)



set.seed(12356)

#Default use 'Signac' method for analysis
#condition 1
grn_outputs1 <- createRegulon(single_cell, n_targets = 5,
                             peak2gene_method="Signac",
                             infer_method = "glm")
target_scores <- suppressWarnings(getSpecificity(single_cell,method = "cosine"))

grn_outputs<-grn_outputs1

target_scores <- suppressWarnings(getSpecificity(single_cell,method = "average"))

target_scores <- suppressWarnings(getSpecificity(single_cell,method = "mean_specificity"))


#condition 3
grn_outputs2 <- createRegulon(single_cell, n_targets = 5,
                             peak2gene_method="Signac",
                             infer_method = "glmnet")
target_scores <- suppressWarnings(getSpecificity(single_cell,method = "cosine"))


grn_outputs<-grn_outputs2
target_scores <- suppressWarnings(getSpecificity(single_cell,method = "average"))

target_scores <- suppressWarnings(getSpecificity(single_cell,method = "mean_specificity"))


#condition 2
grn_outputs3 <- createRegulon(single_cell, n_targets = 5,
                              peak2gene_method="Signac",
                              infer_method = "xgb")

target_scores <- suppressWarnings(getSpecificity(single_cell,method = "cosine"))
grn_outputs<-grn_outputs3

target_scores <- suppressWarnings(getSpecificity(single_cell,method = "average"))
grn_outputs<-grn_outputs3

target_scores <- suppressWarnings(getSpecificity(single_cell,method = "mean_specificity"))
grn_outputs<-grn_outputs3



grn_outputs <- createRegulon(single_cell, n_targets = 5,
                             peak2gene_method="GREAT",
                             infer_method = "glm")


target_scores <- suppressWarnings(getSpecificity(single_cell,method = "cosine"))
target_scores <- suppressWarnings(getSpecificity(single_cell,method = "average"))
target_scores <- suppressWarnings(getSpecificity(single_cell,method = "mean_specificity"))


# Here, the theta and buffer parameters were used the default setting
# theta = 0.5
# buffer = 500bp



##--scMORE for lymphocyte count-----------------------------------------------------------------------------------------
##MAGMA-results of lymphocyte count
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
gene_info_lymp_count  <- read.table("lymp_count_processed_magma_results.genes.out",header = TRUE)
snp_info_lymp_count <- read.csv("lymphocyte_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)


geneRiskScores <- getGeneScore(gene_info_lymp_count)
snp_info <-snp_info_lymp_count

regulon2disease_results <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  #infer_method = 'glmnet',
  infer_method = 'xgb',
  #infer_method = 'glm',
  top_n = 5
)

#E-statistics
lymp_count_cosine <- getEnergyScore(regulon2disease_results,targetCelltype = 2)
log2(lymp_count_cosine+1)+1

#proportion
total_regulon_n <- length(regulon2disease_results$Significance)
sig_n <-length(regulon2disease_results$Significance[which(regulon2disease_results$Significance=="Significant")])
prop1 <- (sig_n+1)/(total_regulon_n+1)
sig_n
print(prop1)

write.csv(regulon2disease_results, file="./01_glm_glmnet_xgb_results/01_lymp_count_cosine_xgb.csv",quote = F)


###---Lymphocyte percent
gene_info_lymp_percent <- read.table("lymp_percent_processed_magma_results.genes.out",header = TRUE)
#snp_info
snp_info_lymp_percent <- read.csv("lymphocyte_percent_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)



geneRiskScores <- getGeneScore(gene_info_lymp_percent)
snp_info <-snp_info_lymp_percent


regulon2disease_results <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  #infer_method = 'glmnet',
  infer_method = 'xgb',
  #infer_method = 'glm',
  top_n = 5
)

#E-statistics
lymp_percent_cosine <- getEnergyScore(regulon2disease_results,targetCelltype = 2)
log2(lymp_percent_cosine+1)+1

#proportion
total_regulon_n <- length(regulon2disease_results$Significance)
sig_n <-length(regulon2disease_results$Significance[which(regulon2disease_results$Significance=="Significant")])
prop2 <- (sig_n+1)/(total_regulon_n+1)
print(prop2)


##---Test results from MAGMA--------------------
#gene_info
gene_info_mono <- read.table("mono_count_processed_magma_results.genes.out",header = TRUE)
#snp_info
snp_info_mono <- read.csv("monocyte_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)



geneRiskScores <- getGeneScore(gene_info_mono)
snp_info <-snp_info_mono

regulon2disease_results <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  #infer_method = 'glmnet',
  #infer_method = 'xgb',
  infer_method = 'glm',
  top_n = 5
)

#E-statistics
mono_count_cosine <- getEnergyScore(regulon2disease_results,targetCelltype = 2)
log2(mono_count_cosine+1)+1

#proportion
total_regulon_n <- length(regulon2disease_results$Significance)
sig_n <-length(regulon2disease_results$Significance[which(regulon2disease_results$Significance=="Significant")])
prop3 <- (sig_n+1)/(total_regulon_n+1)
print(prop3)



###-----------------------------
###---baso count
gene_info_baso <- read.table("baso_count_processed_magma_results.genes.out",header = TRUE)
#snp_info
snp_info_baso <- read.csv("baso_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)



geneRiskScores <- getGeneScore(gene_info_baso)
snp_info <-snp_info_baso


regulon2disease_results <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  #infer_method = 'glmnet',
  infer_method = 'xgb',
  #infer_method = 'glm',
  top_n = 5
)

#E-statistics
baso_count_cosine <- getEnergyScore(regulon2disease_results,targetCelltype = 2)
log2(baso_count_cosine+1)+1

#proportion
total_regulon_n <- length(regulon2disease_results$Significance)
sig_n <-length(regulon2disease_results$Significance[which(regulon2disease_results$Significance=="Significant")])
prop4 <- (sig_n+1)/(total_regulon_n+1)
print(prop4)



###---eosin count
gene_info_eosin <- read.table("eosin_count_processed_magma_results.genes.out",header = TRUE)
#snp_info
snp_info_eosin <- read.csv("eosin_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)



geneRiskScores <- getGeneScore(gene_info_eosin)
snp_info <-snp_info_eosin

regulon2disease_results <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  #infer_method = 'glmnet',
  infer_method = 'xgb',
  #infer_method = 'glm',
  top_n = 5
)

#E-statistics
eosin_count_cosine <- getEnergyScore(regulon2disease_results,targetCelltype = 2)
log2(eosin_count_cosine+1)+1


#proportion
total_regulon_n <- length(regulon2disease_results$Significance)
sig_n <-length(regulon2disease_results$Significance[which(regulon2disease_results$Significance=="Significant")])
prop5 <- (sig_n+1)/(total_regulon_n+1)
print(prop5)


###---HL count
gene_info_HL <- read.table("HL_count_processed_magma_results.genes.out",header = TRUE)
#snp_info
snp_info_HL <- read.csv("HL_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)



geneRiskScores <- getGeneScore(gene_info_HL)
snp_info <-snp_info_HL

regulon2disease_results <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  #infer_method = 'glmnet',
  infer_method = 'xgb',
  #infer_method = 'glm',
  top_n = 5
)

#E-statistics
HL_count_cosine <- getEnergyScore(regulon2disease_results,targetCelltype = 2)
log2(HL_count_cosine+1)+1

#proportion
total_regulon_n <- length(regulon2disease_results$Significance)
sig_n <-length(regulon2disease_results$Significance[which(regulon2disease_results$Significance=="Significant")])
prop6 <- (sig_n+1)/(total_regulon_n+1)
print(prop6)





###---MCHC count
gene_info_MCHC <- read.table("MCHC_count_processed_magma_results.genes.out",header = TRUE)
#snp_info
snp_info_MCHC <- read.csv("MCHC_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)



geneRiskScores <- getGeneScore(gene_info_MCHC)
snp_info <-snp_info_MCHC

regulon2disease_results <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  #infer_method = 'glmnet',
  infer_method = 'xgb',
  #infer_method = 'glm',
  top_n = 5
)

#E-statistics
MCHC_count_cosine <- getEnergyScore(regulon2disease_results,targetCelltype = 2)
log2(MCHC_count_cosine+1)+1

#proportion
total_regulon_n <- length(regulon2disease_results$Significance)
sig_n <-length(regulon2disease_results$Significance[which(regulon2disease_results$Significance=="Significant")])
prop7 <- (sig_n+1)/(total_regulon_n+1)
print(prop7)




###---MCV count
gene_info_MCV <- read.table("MCV_count_processed_magma_results.genes.out",header = TRUE)
#snp_info
snp_info_MCV <- read.csv("MCV_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)


geneRiskScores <- getGeneScore(gene_info_MCV)
snp_info <-snp_info_MCV


regulon2disease_results <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  #infer_method = 'glmnet',
  infer_method = 'xgb',
  #infer_method = 'glm',
  top_n = 5
)

#E-statistics
MCV_count_cosine <- getEnergyScore(regulon2disease_results,targetCelltype = 2)
log2(MCV_count_cosine+1)+1

#proportion
total_regulon_n <- length(regulon2disease_results$Significance)
sig_n <-length(regulon2disease_results$Significance[which(regulon2disease_results$Significance=="Significant")])
prop8 <- (sig_n+1)/(total_regulon_n+1)
print(prop8)



###---neutr count
gene_info_neutr <- read.table("neutr_count_processed_magma_results.genes.out",header = TRUE)
#snp_info
snp_info_neutr <- read.csv("neutr_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)



geneRiskScores <- getGeneScore(gene_info_neutr)
snp_info <-snp_info_neutr


regulon2disease_results <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  #infer_method = 'glmnet',
  infer_method = 'xgb',
  #infer_method = 'glm',
  top_n = 5
)

#E-statistics
neutr_count_cosine <- getEnergyScore(regulon2disease_results,targetCelltype = 2)
log2(neutr_count_cosine+1)+1

#proportion
total_regulon_n <- length(regulon2disease_results$Significance)
sig_n <-length(regulon2disease_results$Significance[which(regulon2disease_results$Significance=="Significant")])
prop9 <- (sig_n+1)/(total_regulon_n+1)
print(prop9)




###---wbc count
gene_info_wbc <- read.table("wbc_count_processed_magma_results.genes.out",header = TRUE)
#snp_info
snp_info_wbc <- read.csv("wbc_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)



geneRiskScores <- getGeneScore(gene_info_wbc)
snp_info <-snp_info_wbc

regulon2disease_results <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  #infer_method = 'glmnet',
  infer_method = 'xgb',
  #infer_method = 'glm',
  top_n = 5
)

#E-statistics
wbc_count_cosine <- getEnergyScore(regulon2disease_results,targetCelltype = 2)
log2(wbc_count_cosine+1)+1

#proportion
total_regulon_n <- length(regulon2disease_results$Significance)
sig_n <-length(regulon2disease_results$Significance[which(regulon2disease_results$Significance=="Significant")])
prop10 <- (sig_n+1)/(total_regulon_n+1)
print(prop10)

##-##--##-#####-####-####-####-###-####-####-####-####-####-####-####-####-##
##-##--##-#####-####-####-####-###-####-####-####-####-####-####-####-####-##





##-##--##-#####-####-####-####-###-####-####-####-####-####-####-####-####-##
##XGBoost
##---------cosine------------------------------
grn_outputs3 <- createRegulon(single_cell, n_targets = 5,
                              peak2gene_method="Signac",
                              infer_method = "xgb")

target_scores <- suppressWarnings(getSpecificity(single_cell,method = "cosine"))
grn_outputs<-grn_outputs3



#1
geneRiskScores <- getGeneScore(gene_info_lymp_count)
snp_info <-snp_info_lymp_count

regulon2disease_results <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  #infer_method = 'glmnet',
  infer_method = 'xgb',
  #infer_method = 'glm',
  top_n = 5
)

#E-statistics
lymp_count_cosine <- getEnergyScore(regulon2disease_results,targetCelltype = 2)
log2(lymp_count_cosine+1)+1

#proportion
total_regulon_n <- length(regulon2disease_results$Significance)
sig_n <-length(regulon2disease_results$Significance[which(regulon2disease_results$Significance=="Significant")])
prop1 <- (sig_n+1)/(total_regulon_n+1)
sig_n
print(prop1)

write.csv(regulon2disease_results, file="./01_glm_glmnet_xgb_results/01_lymp_count_cosine_xgb.csv",quote = F)


#2
geneRiskScores <- getGeneScore(gene_info_lymp_percent)
snp_info <-snp_info_lymp_percent


regulon2disease_results <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  #infer_method = 'glmnet',
  infer_method = 'xgb',
  #infer_method = 'glm',
  top_n = 5
)

#E-statistics
lymp_percent_cosine <- getEnergyScore(regulon2disease_results,targetCelltype = 2)
log2(lymp_percent_cosine+1)+1

#proportion
total_regulon_n <- length(regulon2disease_results$Significance)
sig_n <-length(regulon2disease_results$Significance[which(regulon2disease_results$Significance=="Significant")])
prop2 <- (sig_n+1)/(total_regulon_n+1)
print(prop2)


write.csv(regulon2disease_results, file="./01_glm_glmnet_xgb_results/01_lymp_percent_cosine_xgb.csv",quote = F)


#3
geneRiskScores <- getGeneScore(gene_info_mono)
snp_info <-snp_info_mono

regulon2disease_results <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  #infer_method = 'glmnet',
  infer_method = 'xgb',
  #infer_method = 'glm',
  top_n = 5
)

#E-statistics
mono_count_cosine <- getEnergyScore(regulon2disease_results,targetCelltype = 2)
log2(mono_count_cosine+1)+1

#proportion
total_regulon_n <- length(regulon2disease_results$Significance)
sig_n <-length(regulon2disease_results$Significance[which(regulon2disease_results$Significance=="Significant")])
prop3 <- (sig_n+1)/(total_regulon_n+1)
print(prop3)

write.csv(regulon2disease_results, file="./01_glm_glmnet_xgb_results/01_mono_count_cosine_xgb2.csv",quote = F)


###-----------------------------
#4
geneRiskScores <- getGeneScore(gene_info_baso)
snp_info <-snp_info_baso


regulon2disease_results <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  #infer_method = 'glmnet',
  infer_method = 'xgb',
  #infer_method = 'glm',
  top_n = 5
)

#E-statistics
baso_count_cosine <- getEnergyScore(regulon2disease_results,targetCelltype = 2)
log2(baso_count_cosine+1)+1

#proportion
total_regulon_n <- length(regulon2disease_results$Significance)
sig_n <-length(regulon2disease_results$Significance[which(regulon2disease_results$Significance=="Significant")])
prop4 <- (sig_n+1)/(total_regulon_n+1)
print(prop4)

write.csv(regulon2disease_results, file="./01_glm_glmnet_xgb_results/01_baso_count_cosine_xgb.csv",quote = F)




#5
geneRiskScores <- getGeneScore(gene_info_eosin)
snp_info <-snp_info_eosin

regulon2disease_results <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  #infer_method = 'glmnet',
  infer_method = 'xgb',
  #infer_method = 'glm',
  top_n = 5
)

#E-statistics
eosin_count_cosine <- getEnergyScore(regulon2disease_results,targetCelltype = 2)
log2(eosin_count_cosine+1)+1


#proportion
total_regulon_n <- length(regulon2disease_results$Significance)
sig_n <-length(regulon2disease_results$Significance[which(regulon2disease_results$Significance=="Significant")])
prop5 <- (sig_n+1)/(total_regulon_n+1)
print(prop5)

write.csv(regulon2disease_results, file="./01_glm_glmnet_xgb_results/01_eosin_count_cosine_xgb.csv",quote = F)


#6
geneRiskScores <- getGeneScore(gene_info_HL)
snp_info <-snp_info_HL

regulon2disease_results <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  #infer_method = 'glmnet',
  infer_method = 'xgb',
  #infer_method = 'glm',
  top_n = 5
)

#E-statistics
HL_count_cosine <- getEnergyScore(regulon2disease_results,targetCelltype = 2)
log2(HL_count_cosine+1)+1

#proportion
total_regulon_n <- length(regulon2disease_results$Significance)
sig_n <-length(regulon2disease_results$Significance[which(regulon2disease_results$Significance=="Significant")])
prop6 <- (sig_n+1)/(total_regulon_n+1)
print(prop6)

write.csv(regulon2disease_results, file="./01_glm_glmnet_xgb_results/01_HL_count_cosine_xgb.csv",quote = F)





#7
geneRiskScores <- getGeneScore(gene_info_MCHC)
snp_info <-snp_info_MCHC

regulon2disease_results <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  #infer_method = 'glmnet',
  infer_method = 'xgb',
  #infer_method = 'glm',
  top_n = 5
)

#E-statistics
MCHC_count_cosine <- getEnergyScore(regulon2disease_results,targetCelltype = 2)
log2(MCHC_count_cosine+1)+1

#proportion
total_regulon_n <- length(regulon2disease_results$Significance)
sig_n <-length(regulon2disease_results$Significance[which(regulon2disease_results$Significance=="Significant")])
prop7 <- (sig_n+1)/(total_regulon_n+1)
print(prop7)

write.csv(regulon2disease_results, file="./01_glm_glmnet_xgb_results/01_MCHC_count_cosine_xgb.csv",quote = F)




#8
geneRiskScores <- getGeneScore(gene_info_MCV)
snp_info <-snp_info_MCV


regulon2disease_results <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  #infer_method = 'glmnet',
  infer_method = 'xgb',
  #infer_method = 'glm',
  top_n = 5
)

#E-statistics
MCV_count_cosine <- getEnergyScore(regulon2disease_results,targetCelltype = 2)
log2(MCV_count_cosine+1)+1

#proportion
total_regulon_n <- length(regulon2disease_results$Significance)
sig_n <-length(regulon2disease_results$Significance[which(regulon2disease_results$Significance=="Significant")])
prop8 <- (sig_n+1)/(total_regulon_n+1)
print(prop8)

write.csv(regulon2disease_results, file="./01_glm_glmnet_xgb_results/01_MCV_count_cosine_xgb.csv",quote = F)




#9
geneRiskScores <- getGeneScore(gene_info_neutr)
snp_info <-snp_info_neutr


regulon2disease_results <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  #infer_method = 'glmnet',
  infer_method = 'xgb',
  #infer_method = 'glm',
  top_n = 5
)

#E-statistics
neutr_count_cosine <- getEnergyScore(regulon2disease_results,targetCelltype = 2)
log2(neutr_count_cosine+1)+1

#proportion
total_regulon_n <- length(regulon2disease_results$Significance)
sig_n <-length(regulon2disease_results$Significance[which(regulon2disease_results$Significance=="Significant")])
prop9 <- (sig_n+1)/(total_regulon_n+1)
print(prop9)

write.csv(regulon2disease_results, file="./01_glm_glmnet_xgb_results/01_neutr_count_cosine_xgb.csv",quote = F)



#10
geneRiskScores <- getGeneScore(gene_info_wbc)
snp_info <-snp_info_wbc

regulon2disease_results <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  #infer_method = 'glmnet',
  infer_method = 'xgb',
  #infer_method = 'glm',
  top_n = 5
)

#E-statistics
wbc_count_cosine <- getEnergyScore(regulon2disease_results,targetCelltype = 2)
log2(wbc_count_cosine+1)+1

#proportion
total_regulon_n <- length(regulon2disease_results$Significance)
sig_n <-length(regulon2disease_results$Significance[which(regulon2disease_results$Significance=="Significant")])
prop10 <- (sig_n+1)/(total_regulon_n+1)
print(prop10)

write.csv(regulon2disease_results, file="./01_glm_glmnet_xgb_results/01_wbc_count_cosine_xgb.csv",quote = F)

#######_######_######_######_######_######_######_######_######_######_######_######_




#######_######_######_######_######_######_######_######_######_######_######_######_
##XGBoost
##---------average------------------------------
#condition 2


target_scores <- suppressWarnings(getSpecificity(single_cell,method = "average"))
grn_outputs<-grn_outputs3


#1
geneRiskScores <- getGeneScore(gene_info_lymp_count)
snp_info <-snp_info_lymp_count

regulon2disease_results <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  #infer_method = 'glmnet',
  infer_method = 'xgb',
  #infer_method = 'glm',
  top_n = 5
)

#E-statistics
lymp_count_cosine <- getEnergyScore(regulon2disease_results,targetCelltype = 2)
log2(lymp_count_cosine+1)+1

#proportion
total_regulon_n <- length(regulon2disease_results$Significance)
sig_n <-length(regulon2disease_results$Significance[which(regulon2disease_results$Significance=="Significant")])
prop1 <- (sig_n+1)/(total_regulon_n+1)
sig_n
print(prop1)

write.csv(regulon2disease_results, file="./01_glm_glmnet_xgb_results/01_lymp_count_average_xgb.csv",quote = F)


#2
geneRiskScores <- getGeneScore(gene_info_lymp_percent)
snp_info <-snp_info_lymp_percent


regulon2disease_results <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  #infer_method = 'glmnet',
  infer_method = 'xgb',
  #infer_method = 'glm',
  top_n = 5
)

#E-statistics
lymp_percent_cosine <- getEnergyScore(regulon2disease_results,targetCelltype = 2)
log2(lymp_percent_cosine+1)+1

#proportion
total_regulon_n <- length(regulon2disease_results$Significance)
sig_n <-length(regulon2disease_results$Significance[which(regulon2disease_results$Significance=="Significant")])
prop2 <- (sig_n+1)/(total_regulon_n+1)
print(prop2)


write.csv(regulon2disease_results, file="./01_glm_glmnet_xgb_results/01_lymp_percent_average_xgb.csv",quote = F)


#3
geneRiskScores <- getGeneScore(gene_info_mono)
snp_info <-snp_info_mono

regulon2disease_results <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  #infer_method = 'glmnet',
  infer_method = 'xgb',
  #infer_method = 'glm',
  top_n = 5
)

#E-statistics
mono_count_cosine <- getEnergyScore(regulon2disease_results,targetCelltype = 2)
log2(mono_count_cosine+1)+1

#proportion
total_regulon_n <- length(regulon2disease_results$Significance)
sig_n <-length(regulon2disease_results$Significance[which(regulon2disease_results$Significance=="Significant")])
prop3 <- (sig_n+1)/(total_regulon_n+1)
print(prop3)

write.csv(regulon2disease_results, file="./01_glm_glmnet_xgb_results/01_mono_count_average_xgb2.csv",quote = F)



###-----------------------------
#4
geneRiskScores <- getGeneScore(gene_info_baso)
snp_info <-snp_info_baso


regulon2disease_results <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  #infer_method = 'glmnet',
  infer_method = 'xgb',
  #infer_method = 'glm',
  top_n = 5
)

#E-statistics
baso_count_cosine <- getEnergyScore(regulon2disease_results,targetCelltype = 2)
log2(baso_count_cosine+1)+1

#proportion
total_regulon_n <- length(regulon2disease_results$Significance)
sig_n <-length(regulon2disease_results$Significance[which(regulon2disease_results$Significance=="Significant")])
prop4 <- (sig_n+1)/(total_regulon_n+1)
print(prop4)

write.csv(regulon2disease_results, file="./01_glm_glmnet_xgb_results/01_baso_count_average_xgb.csv",quote = F)




#5
geneRiskScores <- getGeneScore(gene_info_eosin)
snp_info <-snp_info_eosin

regulon2disease_results <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  #infer_method = 'glmnet',
  infer_method = 'xgb',
  #infer_method = 'glm',
  top_n = 5
)

#E-statistics
eosin_count_cosine <- getEnergyScore(regulon2disease_results,targetCelltype = 2)
log2(eosin_count_cosine+1)+1


#proportion
total_regulon_n <- length(regulon2disease_results$Significance)
sig_n <-length(regulon2disease_results$Significance[which(regulon2disease_results$Significance=="Significant")])
prop5 <- (sig_n+1)/(total_regulon_n+1)
print(prop5)

write.csv(regulon2disease_results, file="./01_glm_glmnet_xgb_results/01_eosin_count_average_xgb.csv",quote = F)


#6
geneRiskScores <- getGeneScore(gene_info_HL)
snp_info <-snp_info_HL

regulon2disease_results <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  #infer_method = 'glmnet',
  infer_method = 'xgb',
  #infer_method = 'glm',
  top_n = 5
)

#E-statistics
HL_count_cosine <- getEnergyScore(regulon2disease_results,targetCelltype = 2)
log2(HL_count_cosine+1)+1

#proportion
total_regulon_n <- length(regulon2disease_results$Significance)
sig_n <-length(regulon2disease_results$Significance[which(regulon2disease_results$Significance=="Significant")])
prop6 <- (sig_n+1)/(total_regulon_n+1)
print(prop6)

write.csv(regulon2disease_results, file="./01_glm_glmnet_xgb_results/01_HL_count_average_xgb.csv",quote = F)





#7
geneRiskScores <- getGeneScore(gene_info_MCHC)
snp_info <-snp_info_MCHC

regulon2disease_results <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  #infer_method = 'glmnet',
  infer_method = 'xgb',
  #infer_method = 'glm',
  top_n = 5
)

#E-statistics
MCHC_count_cosine <- getEnergyScore(regulon2disease_results,targetCelltype = 2)
log2(MCHC_count_cosine+1)+1

#proportion
total_regulon_n <- length(regulon2disease_results$Significance)
sig_n <-length(regulon2disease_results$Significance[which(regulon2disease_results$Significance=="Significant")])
prop7 <- (sig_n+1)/(total_regulon_n+1)
print(prop7)

write.csv(regulon2disease_results, file="./01_glm_glmnet_xgb_results/01_MCHC_count_average_xgb.csv",quote = F)




#8
geneRiskScores <- getGeneScore(gene_info_MCV)
snp_info <-snp_info_MCV


regulon2disease_results <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  #infer_method = 'glmnet',
  infer_method = 'xgb',
  #infer_method = 'glm',
  top_n = 5
)

#E-statistics
MCV_count_cosine <- getEnergyScore(regulon2disease_results,targetCelltype = 2)
log2(MCV_count_cosine+1)+1

#proportion
total_regulon_n <- length(regulon2disease_results$Significance)
sig_n <-length(regulon2disease_results$Significance[which(regulon2disease_results$Significance=="Significant")])
prop8 <- (sig_n+1)/(total_regulon_n+1)
print(prop8)

write.csv(regulon2disease_results, file="./01_glm_glmnet_xgb_results/01_MCV_count_average_xgb.csv",quote = F)




#9
geneRiskScores <- getGeneScore(gene_info_neutr)
snp_info <-snp_info_neutr


regulon2disease_results <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  #infer_method = 'glmnet',
  infer_method = 'xgb',
  #infer_method = 'glm',
  top_n = 5
)

#E-statistics
neutr_count_cosine <- getEnergyScore(regulon2disease_results,targetCelltype = 2)
log2(neutr_count_cosine+1)+1

#proportion
total_regulon_n <- length(regulon2disease_results$Significance)
sig_n <-length(regulon2disease_results$Significance[which(regulon2disease_results$Significance=="Significant")])
prop9 <- (sig_n+1)/(total_regulon_n+1)
print(prop9)

write.csv(regulon2disease_results, file="./01_glm_glmnet_xgb_results/01_neutr_count_average_xgb.csv",quote = F)



#10
geneRiskScores <- getGeneScore(gene_info_wbc)
snp_info <-snp_info_wbc

regulon2disease_results <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  #infer_method = 'glmnet',
  infer_method = 'xgb',
  #infer_method = 'glm',
  top_n = 5
)

#E-statistics
wbc_count_cosine <- getEnergyScore(regulon2disease_results,targetCelltype = 2)
log2(wbc_count_cosine+1)+1

#proportion
total_regulon_n <- length(regulon2disease_results$Significance)
sig_n <-length(regulon2disease_results$Significance[which(regulon2disease_results$Significance=="Significant")])
prop10 <- (sig_n+1)/(total_regulon_n+1)
print(prop10)

write.csv(regulon2disease_results, file="./01_glm_glmnet_xgb_results/01_wbc_count_average_xgb.csv",quote = F)
##-##-####-##-####-##-####-##-####-##-####-##-####-##-####-##-####-##-####-##-####-##-##



##-##-####-##-####-##-####-##-####-##-####-##-####-##-####-##-####-##-####-##-####-##-##
#condition 3
#-----glmnet--cosine
grn_outputs<-grn_outputs2
target_scores <- suppressWarnings(getSpecificity(single_cell,method = "cosine"))



#1
geneRiskScores <- getGeneScore(gene_info_lymp_count)
snp_info <-snp_info_lymp_count

regulon2disease_results <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  infer_method = 'glmnet',
  #infer_method = 'xgb',
  #infer_method = 'glm',
  top_n = 5
)

#E-statistics
lymp_count_cosine <- getEnergyScore(regulon2disease_results,targetCelltype = 2)
log2(lymp_count_cosine+1)+1

#proportion
total_regulon_n <- length(regulon2disease_results$Significance)
sig_n <-length(regulon2disease_results$Significance[which(regulon2disease_results$Significance=="Significant")])
prop1 <- (sig_n+1)/(total_regulon_n+1)
sig_n
print(prop1)

write.csv(regulon2disease_results, file="./01_glm_glmnet_xgb_results/01_lymp_count_cosine_glmnet.csv",quote = F)


#2
geneRiskScores <- getGeneScore(gene_info_lymp_percent)
snp_info <-snp_info_lymp_percent


regulon2disease_results <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  infer_method = 'glmnet',
  #infer_method = 'xgb',
  #infer_method = 'glm',
  top_n = 5
)

#E-statistics
lymp_percent_cosine <- getEnergyScore(regulon2disease_results,targetCelltype = 2)
log2(lymp_percent_cosine+1)+1

#proportion
total_regulon_n <- length(regulon2disease_results$Significance)
sig_n <-length(regulon2disease_results$Significance[which(regulon2disease_results$Significance=="Significant")])
prop2 <- (sig_n+1)/(total_regulon_n+1)
print(prop2)


write.csv(regulon2disease_results, file="./01_glm_glmnet_xgb_results/01_lymp_percent_cosine_glmnet.csv",quote = F)


#3
geneRiskScores <- getGeneScore(gene_info_mono)
snp_info <-snp_info_mono

regulon2disease_results <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  infer_method = 'glmnet',
  #infer_method = 'xgb',
  #infer_method = 'glm',
  top_n = 5
)

#E-statistics
mono_count_cosine <- getEnergyScore(regulon2disease_results,targetCelltype = 2)
log2(mono_count_cosine+1)+1

#proportion
total_regulon_n <- length(regulon2disease_results$Significance)
sig_n <-length(regulon2disease_results$Significance[which(regulon2disease_results$Significance=="Significant")])
prop3 <- (sig_n+1)/(total_regulon_n+1)
print(prop3)

write.csv(regulon2disease_results, file="./01_glm_glmnet_xgb_results/01_mono_count_cosine_glmnet.csv",quote = F)


###-----------------------------
#4
geneRiskScores <- getGeneScore(gene_info_baso)
snp_info <-snp_info_baso


regulon2disease_results <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  infer_method = 'glmnet',
  #infer_method = 'xgb',
  #infer_method = 'glm',
  top_n = 5
)

#E-statistics
baso_count_cosine <- getEnergyScore(regulon2disease_results,targetCelltype = 2)
log2(baso_count_cosine+1)+1

#proportion
total_regulon_n <- length(regulon2disease_results$Significance)
sig_n <-length(regulon2disease_results$Significance[which(regulon2disease_results$Significance=="Significant")])
prop4 <- (sig_n+1)/(total_regulon_n+1)
print(prop4)

write.csv(regulon2disease_results, file="./01_glm_glmnet_xgb_results/01_baso_count_cosine_glmnet.csv",quote = F)




#5
geneRiskScores <- getGeneScore(gene_info_eosin)
snp_info <-snp_info_eosin

regulon2disease_results <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  infer_method = 'glmnet',
  #infer_method = 'xgb',
  #infer_method = 'glm',
  top_n = 5
)

#E-statistics
eosin_count_cosine <- getEnergyScore(regulon2disease_results,targetCelltype = 2)
log2(eosin_count_cosine+1)+1


#proportion
total_regulon_n <- length(regulon2disease_results$Significance)
sig_n <-length(regulon2disease_results$Significance[which(regulon2disease_results$Significance=="Significant")])
prop5 <- (sig_n+1)/(total_regulon_n+1)
print(prop5)

write.csv(regulon2disease_results, file="./01_glm_glmnet_xgb_results/01_eosin_count_cosine_glmnet.csv",quote = F)


#6
geneRiskScores <- getGeneScore(gene_info_HL)
snp_info <-snp_info_HL

regulon2disease_results <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  infer_method = 'glmnet',
  #infer_method = 'xgb',
  #infer_method = 'glm',
  top_n = 5
)

#E-statistics
HL_count_cosine <- getEnergyScore(regulon2disease_results,targetCelltype = 2)
log2(HL_count_cosine+1)+1

#proportion
total_regulon_n <- length(regulon2disease_results$Significance)
sig_n <-length(regulon2disease_results$Significance[which(regulon2disease_results$Significance=="Significant")])
prop6 <- (sig_n+1)/(total_regulon_n+1)
print(prop6)

write.csv(regulon2disease_results, file="./01_glm_glmnet_xgb_results/01_HL_count_cosine_glmnet.csv",quote = F)





#7
geneRiskScores <- getGeneScore(gene_info_MCHC)
snp_info <-snp_info_MCHC

regulon2disease_results <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  infer_method = 'glmnet',
  #infer_method = 'xgb',
  #infer_method = 'glm',
  top_n = 5
)

#E-statistics
MCHC_count_cosine <- getEnergyScore(regulon2disease_results,targetCelltype = 2)
log2(MCHC_count_cosine+1)+1

#proportion
total_regulon_n <- length(regulon2disease_results$Significance)
sig_n <-length(regulon2disease_results$Significance[which(regulon2disease_results$Significance=="Significant")])
prop7 <- (sig_n+1)/(total_regulon_n+1)
print(prop7)

write.csv(regulon2disease_results, file="./01_glm_glmnet_xgb_results/01_MCHC_count_cosine_glmnet.csv",quote = F)




#8
geneRiskScores <- getGeneScore(gene_info_MCV)
snp_info <-snp_info_MCV


regulon2disease_results <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  infer_method = 'glmnet',
  #infer_method = 'xgb',
  #infer_method = 'glm',
  top_n = 5
)

#E-statistics
MCV_count_cosine <- getEnergyScore(regulon2disease_results,targetCelltype = 2)
log2(MCV_count_cosine+1)+1

#proportion
total_regulon_n <- length(regulon2disease_results$Significance)
sig_n <-length(regulon2disease_results$Significance[which(regulon2disease_results$Significance=="Significant")])
prop8 <- (sig_n+1)/(total_regulon_n+1)
print(prop8)

write.csv(regulon2disease_results, file="./01_glm_glmnet_xgb_results/01_MCV_count_cosine_glmnet.csv",quote = F)




#9
geneRiskScores <- getGeneScore(gene_info_neutr)
snp_info <-snp_info_neutr


regulon2disease_results <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  infer_method = 'glmnet',
  #infer_method = 'xgb',
  #infer_method = 'glm',
  top_n = 5
)

#E-statistics
neutr_count_cosine <- getEnergyScore(regulon2disease_results,targetCelltype = 2)
log2(neutr_count_cosine+1)+1

#proportion
total_regulon_n <- length(regulon2disease_results$Significance)
sig_n <-length(regulon2disease_results$Significance[which(regulon2disease_results$Significance=="Significant")])
prop9 <- (sig_n+1)/(total_regulon_n+1)
print(prop9)

write.csv(regulon2disease_results, file="./01_glm_glmnet_xgb_results/01_neutr_count_cosine_glmnet.csv",quote = F)



#10
geneRiskScores <- getGeneScore(gene_info_wbc)
snp_info <-snp_info_wbc

regulon2disease_results <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  infer_method = 'glmnet',
  #infer_method = 'xgb',
  #infer_method = 'glm',
  top_n = 5
)

#E-statistics
wbc_count_cosine <- getEnergyScore(regulon2disease_results,targetCelltype = 2)
log2(wbc_count_cosine+1)+1

#proportion
total_regulon_n <- length(regulon2disease_results$Significance)
sig_n <-length(regulon2disease_results$Significance[which(regulon2disease_results$Significance=="Significant")])
prop10 <- (sig_n+1)/(total_regulon_n+1)
print(prop10)

write.csv(regulon2disease_results, file="./01_glm_glmnet_xgb_results/01_wbc_count_cosine_glmnet.csv",quote = F)


##-##-####-##-####-##-####-##-####-##-####-##-####-##-####-##-####-##-####-##-####-##-##
#condition 3
#-----glmnet--cosine

grn_outputs<-grn_outputs2
target_scores <- suppressWarnings(getSpecificity(single_cell,method = "average"))


#1
geneRiskScores <- getGeneScore(gene_info_lymp_count)
snp_info <-snp_info_lymp_count

regulon2disease_results <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  infer_method = 'glmnet',
  #infer_method = 'xgb',
  #infer_method = 'glm',
  top_n = 5
)

#E-statistics
lymp_count_cosine <- getEnergyScore(regulon2disease_results,targetCelltype = 2)
log2(lymp_count_cosine+1)+1

#proportion
total_regulon_n <- length(regulon2disease_results$Significance)
sig_n <-length(regulon2disease_results$Significance[which(regulon2disease_results$Significance=="Significant")])
prop1 <- (sig_n+1)/(total_regulon_n+1)
sig_n
print(prop1)

write.csv(regulon2disease_results, file="./01_glm_glmnet_xgb_results/01_lymp_count_average_glmnet.csv",quote = F)


#2
geneRiskScores <- getGeneScore(gene_info_lymp_percent)
snp_info <-snp_info_lymp_percent


regulon2disease_results <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  infer_method = 'glmnet',
  #infer_method = 'xgb',
  #infer_method = 'glm',
  top_n = 5
)

#E-statistics
lymp_percent_cosine <- getEnergyScore(regulon2disease_results,targetCelltype = 2)
log2(lymp_percent_cosine+1)+1

#proportion
total_regulon_n <- length(regulon2disease_results$Significance)
sig_n <-length(regulon2disease_results$Significance[which(regulon2disease_results$Significance=="Significant")])
prop2 <- (sig_n+1)/(total_regulon_n+1)
print(prop2)


write.csv(regulon2disease_results, file="./01_glm_glmnet_xgb_results/01_lymp_percent_average_glmnet.csv",quote = F)


#3
geneRiskScores <- getGeneScore(gene_info_mono)
snp_info <-snp_info_mono

regulon2disease_results <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  infer_method = 'glmnet',
  #infer_method = 'xgb',
  #infer_method = 'glm',
  top_n = 5
)

#E-statistics
mono_count_cosine <- getEnergyScore(regulon2disease_results,targetCelltype = 2)
log2(mono_count_cosine+1)+1

#proportion
total_regulon_n <- length(regulon2disease_results$Significance)
sig_n <-length(regulon2disease_results$Significance[which(regulon2disease_results$Significance=="Significant")])
prop3 <- (sig_n+1)/(total_regulon_n+1)
print(prop3)

write.csv(regulon2disease_results, file="./01_glm_glmnet_xgb_results/01_mono_count_average_glmnet.csv",quote = F)



###-----------------------------
#4
geneRiskScores <- getGeneScore(gene_info_baso)
snp_info <-snp_info_baso


regulon2disease_results <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  infer_method = 'glmnet',
  #infer_method = 'xgb',
  #infer_method = 'glm',
  top_n = 5
)

#E-statistics
baso_count_cosine <- getEnergyScore(regulon2disease_results,targetCelltype = 2)
log2(baso_count_cosine+1)+1

#proportion
total_regulon_n <- length(regulon2disease_results$Significance)
sig_n <-length(regulon2disease_results$Significance[which(regulon2disease_results$Significance=="Significant")])
prop4 <- (sig_n+1)/(total_regulon_n+1)
print(prop4)

write.csv(regulon2disease_results, file="./01_glm_glmnet_xgb_results/01_baso_count_average_glmnet.csv",quote = F)




#5
geneRiskScores <- getGeneScore(gene_info_eosin)
snp_info <-snp_info_eosin

regulon2disease_results <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  infer_method = 'glmnet',
  #infer_method = 'xgb',
  #infer_method = 'glm',
  top_n = 5
)

#E-statistics
eosin_count_cosine <- getEnergyScore(regulon2disease_results,targetCelltype = 2)
log2(eosin_count_cosine+1)+1


#proportion
total_regulon_n <- length(regulon2disease_results$Significance)
sig_n <-length(regulon2disease_results$Significance[which(regulon2disease_results$Significance=="Significant")])
prop5 <- (sig_n+1)/(total_regulon_n+1)
print(prop5)

write.csv(regulon2disease_results, file="./01_glm_glmnet_xgb_results/01_eosin_count_average_glmnet.csv",quote = F)


#6
geneRiskScores <- getGeneScore(gene_info_HL)
snp_info <-snp_info_HL

regulon2disease_results <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  infer_method = 'glmnet',
  #infer_method = 'xgb',
  #infer_method = 'glm',
  top_n = 5
)

#E-statistics
HL_count_cosine <- getEnergyScore(regulon2disease_results,targetCelltype = 2)
log2(HL_count_cosine+1)+1

#proportion
total_regulon_n <- length(regulon2disease_results$Significance)
sig_n <-length(regulon2disease_results$Significance[which(regulon2disease_results$Significance=="Significant")])
prop6 <- (sig_n+1)/(total_regulon_n+1)
print(prop6)

write.csv(regulon2disease_results, file="./01_glm_glmnet_xgb_results/01_HL_count_average_glmnet.csv",quote = F)





#7
geneRiskScores <- getGeneScore(gene_info_MCHC)
snp_info <-snp_info_MCHC

regulon2disease_results <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  infer_method = 'glmnet',
  #infer_method = 'xgb',
  #infer_method = 'glm',
  top_n = 5
)

#E-statistics
MCHC_count_cosine <- getEnergyScore(regulon2disease_results,targetCelltype = 2)
log2(MCHC_count_cosine+1)+1

#proportion
total_regulon_n <- length(regulon2disease_results$Significance)
sig_n <-length(regulon2disease_results$Significance[which(regulon2disease_results$Significance=="Significant")])
prop7 <- (sig_n+1)/(total_regulon_n+1)
print(prop7)

write.csv(regulon2disease_results, file="./01_glm_glmnet_xgb_results/01_MCHC_count_average_glmnet.csv",quote = F)




#8
geneRiskScores <- getGeneScore(gene_info_MCV)
snp_info <-snp_info_MCV


regulon2disease_results <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  infer_method = 'glmnet',
  #infer_method = 'xgb',
  #infer_method = 'glm',
  top_n = 5
)

#E-statistics
MCV_count_cosine <- getEnergyScore(regulon2disease_results,targetCelltype = 2)
log2(MCV_count_cosine+1)+1

#proportion
total_regulon_n <- length(regulon2disease_results$Significance)
sig_n <-length(regulon2disease_results$Significance[which(regulon2disease_results$Significance=="Significant")])
prop8 <- (sig_n+1)/(total_regulon_n+1)
print(prop8)

write.csv(regulon2disease_results, file="./01_glm_glmnet_xgb_results/01_MCV_count_average_glmnet.csv",quote = F)




#9
geneRiskScores <- getGeneScore(gene_info_neutr)
snp_info <-snp_info_neutr


regulon2disease_results <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  infer_method = 'glmnet',
  #infer_method = 'xgb',
  #infer_method = 'glm',
  top_n = 5
)

#E-statistics
neutr_count_cosine <- getEnergyScore(regulon2disease_results,targetCelltype = 2)
log2(neutr_count_cosine+1)+1

#proportion
total_regulon_n <- length(regulon2disease_results$Significance)
sig_n <-length(regulon2disease_results$Significance[which(regulon2disease_results$Significance=="Significant")])
prop9 <- (sig_n+1)/(total_regulon_n+1)
print(prop9)

write.csv(regulon2disease_results, file="./01_glm_glmnet_xgb_results/01_neutr_count_average_glmnet.csv",quote = F)



#10
geneRiskScores <- getGeneScore(gene_info_wbc)
snp_info <-snp_info_wbc

regulon2disease_results <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  infer_method = 'glmnet',
  #infer_method = 'xgb',
  #infer_method = 'glm',
  top_n = 5
)

#E-statistics
wbc_count_cosine <- getEnergyScore(regulon2disease_results,targetCelltype = 2)
log2(wbc_count_cosine+1)+1

#proportion
total_regulon_n <- length(regulon2disease_results$Significance)
sig_n <-length(regulon2disease_results$Significance[which(regulon2disease_results$Significance=="Significant")])
prop10 <- (sig_n+1)/(total_regulon_n+1)
print(prop10)

write.csv(regulon2disease_results, file="./01_glm_glmnet_xgb_results/01_wbc_count_average_glmnet.csv",quote = F)




