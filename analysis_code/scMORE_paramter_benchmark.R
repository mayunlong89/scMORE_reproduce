#-2025-01-05
#analyses
###-------use 10 blood cell traits to assess the difference of different parameters
#' buffer = 0bp, 50bp, 100bp, 200bp, 500bp, and 1000bp
#' window size = 0kb, 5kb, 10kb, 20kb, 50kb, 100kb.
#' theta = 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1
#' peak2gene_method = Signac and GREAT


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



set.seed(12356)

#Default use 'Signac' method for analysis
grn_outputs <- createRegulon(single_cell, n_targets = 5, peak2gene_method="Signac")

grn_outputs <- createRegulon(single_cell, n_targets = 5, peak2gene_method="GREAT")


target_scores <- suppressWarnings(getSpecificity(single_cell,method = "cosine"))
target_scores <- suppressWarnings(getSpecificity(single_cell,method = "average"))


##--scMORE for lymphocyte count-----------------------------------------------------------------------------------------
##MAGMA-results of lymphocyte count
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
magma_lymp <- read.table("lymp_count_processed_magma_results.genes.out",header = TRUE)
#magma_results <- getGeneScore(magma_mono)
#head(magma_results)
gene_info_lymp <- magma_lymp

snp_info_lymp <- read.csv("lymphocyte_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

geneRiskScores <- getGeneScore(gene_info_lymp)
snp_info <-snp_info_lymp

regulon2disease_results <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 10,
  theta = 1,
  alpha = 1,
  buffer = 500,
  top_n = 5
)

lym_count_cosine <- getEnergyScore(regulon2disease_results,targetCelltype = 2)
log2(lym_count_cosine+1)+1



###---Lymphocyte percent
gene_info_lymp_percent <- read.table("lymp_percent_processed_magma_results.genes.out",header = TRUE)
#snp_info
snp_info_lymp_percent <- read.csv("lymphocyte_percent_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

geneRiskScores <- getGeneScore(gene_info_lymp_percent)
snp_info <-snp_info_lymp_percent



##MAGMA-results
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")

##---Test results from MAGMA--------------------
#gene_info
gene_info_mono <- read.table("mono_count_processed_magma_results.genes.out",header = TRUE)
#snp_info
snp_info_mono <- read.csv("monocyte_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)


geneRiskScores <- getGeneScore(gene_info_mono)
snp_info <-snp_info_mono



###-----------------------------
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


# Step 3: Get gene-level association scores
message("Step 3: Retrieving gene-level association scores...")
geneRiskScores <- getGeneScore(gene_info_MCHC)
snp_info <-snp_info_MCHC






###---MCV count
gene_info_MCV <- read.table("MCV_count_processed_magma_results.genes.out",header = TRUE)
#snp_info
snp_info_MCV <- read.csv("MCV_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)


# Step 3: Get gene-level association scores
message("Step 3: Retrieving gene-level association scores...")
geneRiskScores <- getGeneScore(gene_info_MCV)
snp_info <-snp_info_MCV




###---neutr count
gene_info_neutr <- read.table("neutr_count_processed_magma_results.genes.out",header = TRUE)
#snp_info
snp_info_neutr <- read.csv("neutr_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)


# Step 3: Get gene-level association scores
message("Step 3: Retrieving gene-level association scores...")
geneRiskScores <- getGeneScore(gene_info_neutr)
snp_info <-snp_info_neutr



###---wbc count
gene_info_wbc <- read.table("wbc_count_processed_magma_results.genes.out",header = TRUE)
#snp_info
snp_info_wbc <- read.csv("wbc_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)


# Step 3: Get gene-level association scores
message("Step 3: Retrieving gene-level association scores...")
geneRiskScores <- getGeneScore(gene_info_wbc)
snp_info <-snp_info_wbc








##--theta----0.1 ~ 1.0---running

regulon2disease_results_0.1 <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 10,
  theta = 0.1,
  alpha = 1,
  buffer = 500,
  top_n = 5
)


regulon2disease_results_0.2 <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 10,
  theta = 0.2,
  alpha = 1,
  buffer = 500,
  top_n = 5
)

regulon2disease_results_0.3 <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 10,
  theta = 0.3,
  alpha = 1,
  buffer = 500,
  top_n = 5
)

regulon2disease_results_0.4 <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 10,
  theta = 0.4,
  alpha = 1,
  buffer = 500,
  top_n = 5
)

regulon2disease_results_0.5 <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 10,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  top_n = 5
)

regulon2disease_results_0.6 <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 10,
  theta = 0.6,
  alpha = 1,
  buffer = 500,
  top_n = 5
)

regulon2disease_results_0.7 <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 10,
  theta = 0.7,
  alpha = 1,
  buffer = 500,
  top_n = 5
)

regulon2disease_results_0.8 <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 10,
  theta = 0.8,
  alpha = 1,
  buffer = 500,
  top_n = 5
)

regulon2disease_results_0.9 <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 10,
  theta = 0.9,
  alpha = 1,
  buffer = 500,
  top_n = 5
)

regulon2disease_results_1 <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 10,
  theta = 1,
  alpha = 1,
  buffer = 500,
  top_n = 5
)






##--buffer------running

regulon2disease_results_0 <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 10,
  theta = 0.5,
  alpha = 1,
  buffer = 0,
  top_n = 5
)


regulon2disease_results_50 <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 10,
  theta = 0.5,
  alpha = 1,
  buffer = 50,
  top_n = 5
)

regulon2disease_results_100 <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 10,
  theta = 0.5,
  alpha = 1,
  buffer = 100,
  top_n = 5
)

regulon2disease_results_200 <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 10,
  theta = 0.5,
  alpha = 1,
  buffer = 200,
  top_n = 5
)

regulon2disease_results_500 <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 10,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  top_n = 5
)

regulon2disease_results_800 <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 10,
  theta = 0.5,
  alpha = 1,
  buffer = 800,
  top_n = 5
)

regulon2disease_results_1000 <- regulon2disease(
  grn_outputs = grn_outputs,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 10,
  theta = 0.5,
  alpha = 1,
  buffer = 1000,
  top_n = 5
)


##--theta
getEnergyScore_results_0.1 <- getEnergyScore(regulon2disease_results_0.1,targetCelltype = 2)
getEnergyScore_results_0.2 <- getEnergyScore(regulon2disease_results_0.2,targetCelltype = 2)
getEnergyScore_results_0.3 <- getEnergyScore(regulon2disease_results_0.3,targetCelltype = 2)
getEnergyScore_results_0.4 <- getEnergyScore(regulon2disease_results_0.4,targetCelltype = 2)
getEnergyScore_results_0.5 <- getEnergyScore(regulon2disease_results_0.5,targetCelltype = 2)
getEnergyScore_results_0.6 <- getEnergyScore(regulon2disease_results_0.6,targetCelltype = 2)
getEnergyScore_results_0.7 <- getEnergyScore(regulon2disease_results_0.7,targetCelltype = 2)
getEnergyScore_results_0.8 <- getEnergyScore(regulon2disease_results_0.8,targetCelltype = 2)
getEnergyScore_results_0.9 <- getEnergyScore(regulon2disease_results_0.9,targetCelltype = 2)
getEnergyScore_results_1 <- getEnergyScore(regulon2disease_results_1,targetCelltype = 2)
p01<- log2(getEnergyScore_results_0.1+1)+1
p02<-log2(getEnergyScore_results_0.2+1)+1
p03<-log2(getEnergyScore_results_0.3+1)+1
p04<-log2(getEnergyScore_results_0.4+1)+1
p05<-log2(getEnergyScore_results_0.5+1)+1
p06<-log2(getEnergyScore_results_0.6+1)+1
p07<-log2(getEnergyScore_results_0.7+1)+1
p08<-log2(getEnergyScore_results_0.8+1)+1
p09<-log2(getEnergyScore_results_0.9+1)+1
p010<-log2(getEnergyScore_results_1+1)+1

print(c(p01,p02,p03,p04,p05,p06,p07,p08,p09,p010))


c01<-cor(regulon2disease_results_0.1$RegulonScore,regulon2disease_results_1$RegulonScore)
c02<-cor(regulon2disease_results_0.2$RegulonScore,regulon2disease_results_1$RegulonScore)
c03<-cor(regulon2disease_results_0.3$RegulonScore,regulon2disease_results_1$RegulonScore)
c04<-cor(regulon2disease_results_0.4$RegulonScore,regulon2disease_results_1$RegulonScore)
c05<-cor(regulon2disease_results_0.5$RegulonScore,regulon2disease_results_1$RegulonScore)
c06<-cor(regulon2disease_results_0.6$RegulonScore,regulon2disease_results_1$RegulonScore)
c07<-cor(regulon2disease_results_0.7$RegulonScore,regulon2disease_results_1$RegulonScore)
c08<-cor(regulon2disease_results_0.8$RegulonScore,regulon2disease_results_1$RegulonScore)
c09<-cor(regulon2disease_results_0.9$RegulonScore,regulon2disease_results_1$RegulonScore)
c010<-cor(regulon2disease_results_1$RegulonScore,regulon2disease_results_1$RegulonScore)
print(c(c01,c02,c03,c04,c05,c06,c07,c08,c09,c010))


##buffer
getEnergyScore_results_0 <- getEnergyScore(regulon2disease_results_0,targetCelltype = 2)
p1<-log2(getEnergyScore_results_0+1)+1
getEnergyScore_results_50 <- getEnergyScore(regulon2disease_results_50,targetCelltype = 2)
p2<-log2(getEnergyScore_results_50+1)+1
getEnergyScore_results_100 <- getEnergyScore(regulon2disease_results_100,targetCelltype = 2)
p3<-log2(getEnergyScore_results_100+1)+1
getEnergyScore_results_200 <- getEnergyScore(regulon2disease_results_200,targetCelltype = 2)
p4<-log2(getEnergyScore_results_200+1)+1
getEnergyScore_results_500 <- getEnergyScore(regulon2disease_results_500,targetCelltype = 2)
p5<-log2(getEnergyScore_results_500+1)+1
getEnergyScore_results_800 <- getEnergyScore(regulon2disease_results_800,targetCelltype = 2)
p6<-log2(getEnergyScore_results_800+1)+1
getEnergyScore_results_1000 <- getEnergyScore(regulon2disease_results_1000,targetCelltype = 2)
p7<-log2(getEnergyScore_results_1000+1)+1

print(c(p1,p2,p3,p4,p5,p6,p7))



c1<-cor(regulon2disease_results_0$RegulonScore,regulon2disease_results_0$RegulonScore)
c2<-cor(regulon2disease_results_50$RegulonScore,regulon2disease_results_0$RegulonScore)
c3<-cor(regulon2disease_results_100$RegulonScore,regulon2disease_results_0$RegulonScore)
c4<-cor(regulon2disease_results_200$RegulonScore,regulon2disease_results_0$RegulonScore)
c5<-cor(regulon2disease_results_500$RegulonScore,regulon2disease_results_0$RegulonScore)
c6<-cor(regulon2disease_results_800$RegulonScore,regulon2disease_results_0$RegulonScore)
c7<-cor(regulon2disease_results_1000$RegulonScore,regulon2disease_results_0$RegulonScore)
print(c(c1,c2,c3,c4,c5,c6,c7))

