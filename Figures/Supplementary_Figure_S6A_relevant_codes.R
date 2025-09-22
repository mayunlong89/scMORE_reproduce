###------------------------------------------------------------------------------------
#2025-01-09
###------------------------------------------------------------------------------------

###required load R packages
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


#devtools::install_github("mayunlong89/scMORE")

library(scMORE)
library(GenomicRanges)
library(IRanges)
library(Seurat)
library(Signac)


#load data on 10x pbmc example data
#Load scMultiomic data
pbmc_10x <- readRDS("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/10X_PBMC/10X_PBMC.rds")

length(unique(pbmc_10x$cell_type))

single_cell_8900<-pbmc_10x


grn_outputs8900 <- createRegulon(single_cell_8900, n_targets=5,
                             peak2gene_method = 'Signac',
                             infer_method = 'glm')

grn_outputs8900_glmnet <- createRegulon(single_cell_8900, n_targets = 5,
                              peak2gene_method="Signac",
                              infer_method = "glmnet")

grn_outputs8900_xgb <- createRegulon(single_cell_8900, n_targets=5,
                              peak2gene_method = 'Signac',
                              infer_method = 'xgb')



#load data on 10x human brain data
#Load scMultiomic data
human_brain_10x <- readRDS("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/10X_human_brain/10X_Human_Brain.rds")

length(unique(human_brain_10x$cell_type))

single_cell_2464   <- human_brain_10x

grn_outputs2464_glm <- createRegulon(single_cell_2464, n_targets=5,
                                 peak2gene_method = 'Signac',
                                 infer_method = 'glm')

grn_outputs8900_glmnet <- createRegulon(single_cell_2464, n_targets = 5,
                                        peak2gene_method="Signac",
                                        infer_method = "glmnet")

grn_outputs2464_xgb <- createRegulon(single_cell_2464, n_targets=5,
                                     peak2gene_method = 'Signac',
                                     infer_method = 'xgb')



#load data on 10x human B Cell lymphoma
#Load scMultiomic data
human_Bcell_10x <- readRDS("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/10X_B_Cell_Lymphoma/10X_B_Cell_Lymphoma.rds")

length(unique(human_Bcell_10x$cell_type))

single_cell_10478 <- human_Bcell_10x


grn_outputs10478_glm <- createRegulon(single_cell_10478, n_targets=5,
                                     peak2gene_method = 'Signac',
                                     infer_method = 'glm')

grn_outputs10478_glmnet <- createRegulon(single_cell_10478, n_targets = 5,
                                        peak2gene_method="Signac",
                                        infer_method = "glmnet")

grn_outputs10478_xgb <- createRegulon(single_cell_10478, n_targets=5,
                                     peak2gene_method = 'Signac',
                                     infer_method = 'xgb')




#load data on 10x human PD Brain data
#Load scMultiomic data
human_PD_brain_10x <- readRDS("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/06-PD/PD_anno.rds")

human_PD_brain_10x <- JoinLayers(human_PD_brain_10x)

length(unique(human_PD_brain_10x$cell_type))

single_cell_58949  <- human_PD_brain_10x

grn_outputs58949_glm <- createRegulon(single_cell_58949, n_targets=5,
                                      peak2gene_method = 'Signac',
                                      infer_method = 'glm')

grn_outputs58949_glmnet <- createRegulon(single_cell_58949, n_targets = 5,
                                         peak2gene_method="Signac",
                                         infer_method = "glmnet")

grn_outputs58949_xgb <- createRegulon(single_cell_58949, n_targets=5,
                                      peak2gene_method = 'Signac',
                                      infer_method = 'xgb')




#load data on 10x human organoid eye cells
#Load scMultiomic data
human_organoid_eye_10x <- readRDS("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/01_five_eye_diseases_and_organoids/Final_anno_organoids.rds")
Idents(human_organoid_eye_10x) <- human_organoid_eye_10x$Final_anno

length(unique(human_organoid_eye_10x$Final_anno))

single_cell_2827 <- human_organoid_eye_10x


grn_outputs2827_glm <- createRegulon(single_cell_2827, n_targets=5,
                                      peak2gene_method = 'Signac',
                                      infer_method = 'glm')

grn_outputs2827_glmnet <- createRegulon(single_cell_2827, n_targets = 5,
                                         peak2gene_method="Signac",
                                         infer_method = "glmnet")

grn_outputs2827_xgb <- createRegulon(single_cell_2827, n_targets=5,
                                      peak2gene_method = 'Signac',
                                      infer_method = 'xgb')
