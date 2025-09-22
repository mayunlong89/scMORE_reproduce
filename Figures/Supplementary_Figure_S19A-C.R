
#----2025-01-31


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

PD_seu1 <- readRDS("./PD_anno.rds")

PD_seu2 <- JoinLayers(PD_seu1)

PD_seu <-PD_seu2

head(PD_seu@meta.data)

table(PD_seu$cell_type)

#layer ID
table(unique(PD_seu$orig.ident))
length(unique(PD_seu$orig.ident))

PD_seu@meta.data$cell_type <- factor(PD_seu@meta.data$cell_type,levels=c("ODC","AS","MG","N","EC","OPC","T"))


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

setdiff(PD_seu$orig.ident, meta_anno$ID)  # 哪些 sample 在 Seurat 中但不在 metadata 中


PD_seu$group <- meta_anno$group[match(PD_seu$orig.ident, meta_anno$ID)]
PD_seu$age <- meta_anno$age[match(PD_seu$orig.ident, meta_anno$ID)]
PD_seu$sex <- meta_anno$sex[match(PD_seu$orig.ident, meta_anno$ID)]

Idents(PD_seu) <- PD_seu$group
pal <- c("#6A8473","#A13B46","#C0BFDF")

Idents(PD_seu) <- PD_seu$sex
pal <- c("#A13B46","#67ADB7")

#DimPlot---Figure 5C
DimPlot(PD_seu, reduction = "umap.cca",label = F, cols=pal, pt.size = 0.1, repel = T)


Idents(PD_seu) <- PD_seu$age

#DimPlot---
DimPlot(PD_seu, reduction = "umap.cca", pt.size = 0.1, repel = T)


##------Figure 5-G Proportion of significant trait-regulon pairs

# 加载必要的包
library(ggplot2)

# 示例数据
data <- data.frame(
  CellType = factor(c("ASs", "OPCs", "ECs", "Ts", "MGs", "Ns", "ODCs"),
                    levels = c("ASs", "OPCs", "ECs", "Ts", "MGs", "Ns", "ODCs")),
  Proportion = c(0.21978, 0.1904762, 0.14285714, 0.12820513, 0.11721612, 0.1025641, 0.0989011)
)

# 简洁柱状图
p <- ggplot(data, aes(x = CellType, y = Proportion)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black", width = 0.7) +
  theme_classic(base_size = 14) +
  labs(
    title = "Proportion of PD-related regulons",
    x = NULL,
    y = "Proportion"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 13)
  )

# 显示图形
print(p)






###_--Sup Figure S B

##----dotPlot for markers
markers <- c('SLC17A7','SLC17A6','SLC32A1','TBR1','STMN2','MOBP','MOG','SOX2','OLIG1','OLIG2','SLC1A2','GFAP','AQP4','FLT1','VWF','PDGFRB','TREM2','GPR34','CD3E','CD8A')
PD_seu@meta.data$cell_type <- factor(PD_seu@meta.data$cell_type,levels=c("T","MG","EC","AS","OPC","ODC","N"))

Idents(PD_seu) <- PD_seu@meta.data$cell_type
DotPlot(PD_seu,features = markers)+coord_flip()
