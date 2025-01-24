
#2025-01-24
#2025-01-12

#analyses
###-------use 10 blood cell traits to benchmark scMORE to VAM, AUCell, RSS, AddModuleScore, MAGMA_CellTyping

# using groundtruth data

#devtools::install_github("mayunlong89/scMORE")

library(scMORE)
library(GenomicRanges)
library(IRanges)
library(Seurat)
library(Signac)


set.seed(12356)

##---2000Cells

#---load pbmc_10x_real_data_downsampled_2000

pbmc_10x_real_data_downsampled_2000 <- readRDS("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/10X_PBMC/10X_PBMC_downsample_2000cells.rds")

DimPlot(pbmc_10x_real_data_downsampled_2000, reduction = "umap.new",cols = c("#67ADB7","#E77A77"))

single_cell <-pbmc_10x_real_data_downsampled_2000


#Default use 'Signac' method for GRN inference analysis
#condition 1
grn_outputs1 <- createRegulon(single_cell, n_targets = 5,
                              peak2gene_method="Signac",
                              infer_method = "glm")

grn_outputs<-grn_outputs1


#----benchmarking gene set specificity methods

#scMORE default cosine similarity method

target_scores <- suppressWarnings(getSpecificity(single_cell,method = "cosine"))
head(target_scores)

target_scores1 <- suppressWarnings(getSpecificity(single_cell,method = "average"))
head(target_scores1)

#scMORE using mean_specificity method implemented in MAGMA_CellTyping
target_scores2 <- suppressWarnings(getSpecificity(single_cell,method = "mean_specificity"))
head(target_scores2)


setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
gene_info_lymp  <- read.table("lymp_count_processed_magma_results.genes.out",header = TRUE)
snp_info_lymp <- read.csv("lymphocyte_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

gene_info<-gene_info_lymp
snp_info<-snp_info_lymp
geneRiskScores <- getGeneScore(gene_info)


##函数需要修改：
#getRegulonScore.R
#regulon2disease.R
#getSpecificity.R
#参考github: for AUCell, AddModuleScore, VAM, VISION,
#https://github.com/mayunlong89/windows_Rscript/blob/main/Rscripts/ctDRTF_vision_AUCell_addModuleScore.R
#RSS
#https://rdrr.io/github/aertslab/SCENIC/man/calcRSS.html

#构建不同genetic architecture GWAS
#https://github.com/dengchunyu/scPagwas_reproduce/blob/main/Analysis/random_1000test_celltypes.r

###Random cell types and permuted genetic architectures
#https://github.com/dengchunyu/scPagwas_reproduce/blob/main/Analysis/random_1000test_celltypes.r


###----regulon2disease parameters
perm_n = 1000
theta = 0.5
alpha = 1
top_n = 5
buffer = 500
infer_method = 'glm'
p1 = 0.05
p2 = 0.05
p3 = 0.05


## for a regulon extracted from all TF-regulons
Module_regulon

Idents(single_cell) <- single_cell$cell_type2


####--------AUCell------1---Mean
library(AUCell)

cells_rankings <- AUCell_buildRankings(single_cell@assays$RNA@data)

#Calculating AUCell for each regulon
#cells_rankings <- AUCell_buildRankings(organoid_brain_multi@assays$RNA@data)
cells_AUC <- AUCell_calcAUC(Module_regulon,cells_rankings,aucMaxRank = nrow(cells_rankings)*0.05)
#cells_assignment <- AUCell_exploreThresholds(cells_AUC,plotHist = TRUE, nCores = 1, assign =TRUE)

#Extracting AUCell score
AUCell_auc <- as.numeric(getAUC(cells_AUC)[1,])

#Assign AUCell score to each cell
auc_cell <- as.data.frame(Idents(single_cell))

auc_cell$AUCell_auc <- AUCell_auc

colnames(auc_cell) <- c("celltypes","scores")
results_auc <- aggregate(scores~celltypes,auc_cell,mean)

results_auc$regluons <- Module_regulon[1]

head(results_auc)



####--------AUCell------2---RSS----calcRSS()
library(AUCell)
source("~/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/scMORE_test_codes/calcRSS.R")

cells_rankings <- AUCell_buildRankings(single_cell@assays$RNA@data)

#Calculating AUCell for each regulon
#cells_rankings <- AUCell_buildRankings(organoid_brain_multi@assays$RNA@data)
cells_AUC <- AUCell_calcAUC(Module_regulon,cells_rankings,aucMaxRank = nrow(cells_rankings)*0.05)
#cells_assignment <- AUCell_exploreThresholds(cells_AUC,plotHist = TRUE, nCores = 1, assign =TRUE)

#Extracting AUCell score
AUCell_rss <- as.numeric(calcRSS(cells_AUC,cellAnnotation = Idents(single_cell))[1,])
names(AUCell_rss) <- levels(Idents(single_cell))
#Assign AUCell score to each cell
rss_cell <- as.data.frame(AUCell_rss)

rss_cell$regluons <- Module_regulon[1]
rss_cell$celltypes <- rownames(rss_cell)
results_rss <- rss_cell

colnames(results_rss) <- c("scores","regulons","celltypes")

head(results_rss)





##----2
library(VISION)
# 检查和安装依赖
#install.packages(c("irlba", "Matrix", "methods"))

Module_regulon

# 获取 RNA counts 并归一化
counts <- single_cell@assays$RNA@counts
n.umi <- colSums(counts)

# 使用矩阵操作进行归一化

library(Matrix)
# 确保输入是稀疏矩阵
scaled_counts <- as(Matrix(scaled_counts, sparse = TRUE), "dgCMatrix")
scaled_counts <- t(t(counts) / n.umi) * median(n.umi)

# Function to convert to GMT format and store in memory
output_gmt <- function(geneset) {
  # Convert each gene set to GMT format as a string
  gmt_data <- sapply(names(geneset), function(i) {
    paste(c(i, "temp", geneset[[i]]), collapse = "\t")
  })

  # Return as a character vector
  return(gmt_data)
}

#generating .gmt file
# Use the function to get GMT format in memory
geneset <- list(Module_regulon)
names(geneset) <- Module_regulon[1]
#print(names(geneset))
Module_regulon_gmt <- output_gmt(geneset)
gmt_file <- tempfile(fileext = ".gmt")
writeLines(Module_regulon_gmt, con = gmt_file)

# 使用 Vision 运行分析
vis <- Vision(scaled_counts,       # 归一化的表达矩阵
  signatures = gmt_file       # 临时生成的签名文件
)

# 执行 Vision 分析
vis <- analyze(vis)
#Viewing the results
#viewResults(vis)
getSignatureAutocorrelation(vis)




##---3
########-----------------VAM
library(VAM)  #install.packages("VAM")
########-----------------VAM
library(VAM)  #install.packages("VAM")
#BiocManager::install("qusage")
library(qusage)
#devtools::install_github("cellgeni/sceasy")
#install.packages('reticulate')
library(sceasy)
library(reticulate)

## Convert seurat format to sce format
#sceasy::convertFormat(single_cell,from = "seurat",to="sce",outFile = "organoids_brain_multi_sce.rds")
library(Seurat)
library(SeuratDisk)
library(SingleCellExperiment)

# 将 Seurat 对象转换为 SingleCellExperiment 对象
sce_multi <- as.SingleCellExperiment(single_cell)


# Function to convert to GMT format and store in memory
output_gmt <- function(geneset) {
  # Convert each gene set to GMT format as a string
  gmt_data <- sapply(names(geneset), function(i) {
    paste(c(i, "temp", geneset[[i]]), collapse = "\t")
  })

  # Return as a character vector
  return(gmt_data)
}

#generating .gmt file
# Use the function to get GMT format in memory
geneset <- list(Module_regulon)
names(geneset) <- Module_regulon[1]
#print(names(geneset))
Module_regulon_gmt <- output_gmt(geneset)
gmt_file <- tempfile(fileext = ".gmt")
writeLines(Module_regulon_gmt, con = gmt_file)


gs <- qusage::read.gmt(gmt_file)

sce_multi <- importGeneSetsFromList(inSCE = sce_multi,geneSetList = gs,
                                    by = "rownames")

sce_multi <- runVAM(inSCE = sce_multi,
                    geneSetCollectionName = "GeneSetCollection",
                    useAssay = "logcounts")

#VAM values
sce_multi@int_colData$reducedDims$VAM_GeneSetCollection_Distance
results_VAM <- sce_multi@int_colData$reducedDims$VAM_GeneSetCollection_CDF
VAM_cell<-c()
VAM_cell <- as.data.frame(Idents(single_cell))

VAM_cell$vam_score <- results_VAM[,1]

colnames(VAM_cell) <- c("celltypes","scores")
results_vam2 <- aggregate(scores~celltypes,VAM_cell,mean)

results_vam2$regluons <- Module_regulon[1]

head(results_vam2)
head(results)



###--------addmodulescore

  #Calculating addModuleScore for each regulon
  DefaultAssay(single_cell) <- "RNA"
  single_cell <- AddModuleScore(single_cell,
                                features = list(Module_regulon),
                                name="Module_regulon")

  add_cell <- as.data.frame(Idents(single_cell))

  add_cell$add_score <- single_cell$Module_regulon1
  colnames(add_cell) <- c("celltypes","regulon_score")
  results_add <- aggregate(regulon_score~celltypes,add_cell,mean)

  results_add$regluons <- Module_regulon[1]
  colnames(results_add) <- c("celltypes","scores","regulons")
  head(results_add)



#---------UCell
  if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
  BiocManager::install("UCell")


  library(UCell)

  seu_matrix <- single_cell@assays$RNA@data
  gene.sets <- list(Module_regulon)

  UCellscores <- ScoreSignatures_UCell(seu_matrix, features=gene.sets)
  head(UCellscores)

  results_ucell <- as.data.frame(Idents(single_cell))

  results_ucell$ucell_score <-  UCellscores
  colnames(results_ucell) <- c("celltypes","regulon_score")
  results_ucell2 <- aggregate(regulon_score~celltypes,results_ucell,mean)

  results_ucell2$regluons <- Module_regulon[1]
  colnames(results_ucell2) <- c("celltypes","scores","regulons")
  head(results_ucell2)















