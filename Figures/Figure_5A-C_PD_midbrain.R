
setwd("~/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/06-PD/")

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


###---pando_results for brain organoids single-cell data

PD_seu <- readRDS("./PD_anno.rds")

PD_seu <- JoinLayers(PD_seu)


head(PD_seu@meta.data)

#layer ID
table(unique(PD_seu$orig.ident))
length(unique(PD_seu$orig.ident))



Idents(PD_seu) <- PD_seu$cell_type
pal <- c("#A13B46","#67ADB7","#36600E","#6A8473","#C0BFDF","#7E6148FF",
         "#85C89C")

#DimPlot---
DimPlot(PD_seu, reduction = "umap.cca",label = T, cols=pal, pt.size = 0.1, repel = T)


#pal <- c("#A13B46","#67ADB7","#36600E","#6A8473","#C0BFDF","#ECBA84", "#85C89C")

#pal <- c("#A13B46","#67ADB7","#36600E","#6A8473","#C0BFDF","darkorchid2","#85C89C")
#"#E77A77", "#ECBA84"


##read metadata

meta_anno <- read.csv("01_metadata_PD.csv", header=T)
head(meta_anno)

PD_seu$group <- meta_anno$group[match(PD_seu$orig.ident, meta_anno$ID)]
PD_seu$age <- meta_anno$age[match(PD_seu$orig.ident, meta_anno$ID)]
PD_seu$sex <- meta_anno$sex[match(PD_seu$orig.ident, meta_anno$ID)]

Idents(PD_seu) <- PD_seu$group
pal <- c("#6A8473","#A13B46","#C0BFDF")

Idents(PD_seu) <- PD_seu$sex
pal <- c("#A13B46","#67ADB7")

#DimPlot---
DimPlot(PD_seu, reduction = "umap.cca",label = F, cols=pal, pt.size = 0.1, repel = T)


Idents(PD_seu) <- PD_seu$age

#DimPlot---
DimPlot(PD_seu, reduction = "umap.cca", pt.size = 0.1, repel = T)


