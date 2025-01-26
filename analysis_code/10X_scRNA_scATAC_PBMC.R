############### Seurat整合scRNA，scATAC数据 ######################
############### Seurat整合多模态数据 #############################
############### scRNA的数据是由Seurat处理的
############### scATAC的数据是由Signac处理的
library(Seurat)
library(Signac)
library(ggplot2)
library(hdf5r)
library(biovizBase)
library(patchwork)
library(EnsDb.Hsapiens.v86)
library(AnnotationHub)
library(dplyr)
set.seed(123)
setwd("F:/Project/scHBO/10X_PBMC/")

############## 读入h5文件 ##################
############## 文件中包含scRNA和scATAC的数据
counts <- Read10X_h5(filename ="F:/Project/scHBO/Datasets/10x_PBMC/scATAC/pbmc_unsorted_10k_filtered_feature_bc_matrix.h5")

############## 读取scRNA数据
single_data<-CreateSeuratObject(counts = counts$`Gene Expression`,project="PBMC",min.cells = 10,min.features = 0)
single_data[["percent.mt"]] <- PercentageFeatureSet(single_data, pattern = "^MT-")
single_data[["percent.rp"]] <- PercentageFeatureSet(single_data, pattern = "^RP[SL][[:digit:]]")

############## 加入metadata信息，该数据只有一个女性样本
single_data$sex<-"Female"
single_data$age<-"25"

############## Signac创建scATAC数据 ################
counts$Peaks<-counts$Peaks[grep("^chr",row.names(counts$Peaks)),] #### 这里是提取和参考基因注释信息一样的染色体编号
chrom_assay <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = 'F:/Project/scHBO/Datasets/10x_PBMC/scATAC/pbmc_unsorted_10k_atac_fragments.tsv.gz',
  min.cells = 10,
  min.features = 0
)

############## 构建基因注释信息
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"
Annotation(chrom_assay) <- annotations

############# 将ATAC加入到单细胞中
single_data[["peaks"]]<-chrom_assay

############ 在ATAC矩阵下，进行质控 ################
DefaultAssay(single_data)<-"peaks"
############ 核小体条带模式
single_data <- NucleosomeSignal(object = single_data)
############ 转录起始位点富集分数(TSS)
single_data <- TSSEnrichment(object = single_data, fast = FALSE)
############ 计算blacklist_fraction
single_data$blacklist_fraction <- FractionCountsInRegion(
  object = single_data, 
  assay = 'peaks',
  regions = blacklist_hg38
)

########### 结合scRNA和scATAC的数据进行质控 ################
VlnPlot(object = single_data,
        features = c("nCount_peaks","nFeature_peaks","nCount_RNA","nFeature_RNA","nucleosome_signal",
                     "TSS.enrichment","blacklist_fraction","percent.mt","percent.rp"),
        pt.size = 0.1,
        ncol = 5)

########### 剩下8900个细胞
single_data1<-subset(single_data,subset = nCount_peaks < 80000 &
                      nCount_peaks > 2000 &
                      nCount_RNA < 20000 &
                      nCount_RNA > 500 &
                      nFeature_RNA > 300 &
                      nucleosome_signal < 4 &
                      TSS.enrichment > 2 &
                      blacklist_fraction < 0.05 &
                      percent.mt < 15 &
                      percent.rp < 15)

############ 分别在RNA和ATAC做不同组学的标准流程 ##################
############ scRNA
DefaultAssay(single_data1)<-"RNA"
single_data1 <- NormalizeData(single_data1, normalization.method = "LogNormalize", scale.factor = 10000)
single_data1 <- FindVariableFeatures(single_data1, selection.method = "vst", nfeatures = 3000)
single_data1 <- ScaleData(single_data1, verbose = FALSE)
single_data1 <- RunPCA(single_data1, npcs = 20, verbose = FALSE)
single_data1 <- RunUMAP(single_data1, reduction = "pca", dims = 1:20,reduction.name = 'umap.rna', reduction.key = 'UMAP.RNA_')

############ scATAC
DefaultAssay(single_data1)<-"peaks"
single_data1 <- RunTFIDF(single_data1)
single_data1 <- FindTopFeatures(single_data1, min.cutoff = 'q0')
single_data1 <- RunSVD(single_data1)
single_data1 <- RunUMAP(single_data1, reduction = 'lsi', dims = 2:30, reduction.name = "umap.atac", reduction.key = "UMAP.ATAC_")

############ 整合两种组学聚类分群 ##################
single_data1 <- FindMultiModalNeighbors(object = single_data1,
                                        reduction.list = list("pca", "lsi"), 
                                        dims.list = list(1:20, 2:30),
                                        modality.weight.name = "RNA.weight",
                                        verbose = TRUE)
single_data1 <- FindClusters(object = single_data1, resolution = 0.2, graph.name="wsnn")

############ umap
single_data1 <- RunUMAP(object = single_data1,
                        nn.name = "weighted.nn",
                        #assay = "RNA",
                        verbose = TRUE,
                        reduction.name="umap.biMod",
                        reduction.key = "UMAP.biMod_"
                        #dims = 1:50
                        )

############ tsne
single_data1 <- RunTSNE(object = single_data1,
                        nn.name = "weighted.nn",
                        #assay = "RNA",
                        verbose = TRUE,
                        reduction.name="tsne.biMod",
                        reduction.key = "TSNE.biMod_",
                        #dims = 1:50
                        )

################# 细胞注释 #########################
DefaultAssay(single_data1)<-"RNA"
pdf("wsnn_res.0.2.pdf")
DimPlot(single_data1,reduction = "umap.biMod",group.by = "wsnn_res.0.2",label = T)
dev.off()

marker_gene<-c("IL3RA","LILRA4",'PTGDS', ### pDC 
               "CLEC10A","FCER1A",  ### mDC
               "CD79A","MS4A1","IGHM","CD19","CD79B", ### B cells
               "CD14","S100A12",       ### CD14+ monocyte
               "FCGR3A","CX3CR1","MS4A7",      ### FCGR3A+ monocyte
               "CD4",   ### CD4+ T cells
               "CD8A",   ### CD8+ T cells
               "CD3D",'CD3E','CD2', ### T cells
               "GNLY","NKG7",'KLRF1','KLRD1' ### NK 
               )

pdf("marker for each cluster.pdf")
DotPlot(single_data1,features=unique(rev(marker_gene)),group.by="wsnn_res.0.2",assay="RNA")+
  scale_colour_gradient2(low="#3A71AA",mid="white",high="#B22028",midpoint=-0.4)+
  theme_bw()+
  coord_flip()+
  theme (axis.text.x = element_text (angle = 35, hjust = 1))
dev.off()

############ 注释
single_data1$cell_type<-single_data1$wsnn_res.0.2
single_data1$cell_type<-recode(single_data1$cell_type,
                               "0"="CD14+ Monocytes",
                               "1"="CD4+ T cells",
                               "2"="CD4+ T cells",
                               "3"="CD8+ T cells",
                               "4"="B cells",
                               "5"="NK cells",
                               "6"="NK cells",
                               "7"="FCGR3A+ Monocytes",
                               "8"="pDCs",
                               "9"="NK cells",
                               "10"="mDCs",
                               "11"="B cells")

single_data1$cell_type<-factor(single_data1$cell_type,levels = c("pDCs","mDCs","B cells","CD14+ Monocytes","FCGR3A+ Monocytes",
                                                                 "CD4+ T cells","CD8+ T cells","NK cells"))
colors<-c("#77A0EA","#ECB888","#F0BDD7","#A3C4A4","#B68ABC","#FC7540","#A032CB",
          "#CD8A75","#E9D7EA","#EED785","#6076B2","#F0BDD7","#F5D39C","#CE7460",
          "#7FC28D","#F7CBBD","#C76EA7","#FEF6AD","#EDE1E7","#E5B98D","#D8706F",
          "#C6D9DE","#E16EB8","#928EBC","#AF89BB","#DDDBEE","#F07392","#DEB9D4")

##############  umap图
pdf("cell_type_umap.pdf")
DimPlot(single_data1,reduction = "umap.biMod",group.by = "cell_type",
        cols = colors,label = T,repel = T,pt.size = 0.4,label.size = 3.5)+NoLegend()
dev.off()

############ dotplot for markers
pdf("marker for each celltype.pdf")
DotPlot(single_data1,features=unique(rev(marker_gene)),group.by="cell_type",assay="RNA")+
  scale_colour_gradient2(low="#3A71AA",mid="white",high="#B22028",midpoint=-0.4)+
  theme_bw()+
  coord_flip()+
  theme (axis.text.x = element_text (angle = 35, hjust = 1))
dev.off()

########### 保存单细胞数据
saveRDS(single_data1,file = "10X_PBMC.rds")

