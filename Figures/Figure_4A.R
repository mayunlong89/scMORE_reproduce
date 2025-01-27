#2025-01-26

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

