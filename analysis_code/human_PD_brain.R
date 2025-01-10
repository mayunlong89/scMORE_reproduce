library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(cowplot)
library(readr)
library(tidyverse)
library(BSgenome.Hsapiens.UCSC.hg38)

#working location:
setwd('/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/')

human_PD_brain_10x <- readRDS("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/06-PD/PD_anno.rds")


###Load 10X format matrix
sample_list <- readr::read_csv('sample.csv')
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"

scmulti.list <- list()
for(i in c(1:length(sample_list$Samples))){
    sp <- sample_list$Samples[[i]]
    path1 <- glue::glue('./GSE193688/{sp}_filtered_feature_bc_matrix.h5')
    path2 <- glue::glue('./GSE193688/{sp}_atac_fragments.tsv.gz')
    sc_counts <- Read10X_h5(path1)
    sce <- CreateSeuratObject(counts = sc_counts$`Gene Expression`,
                                 assay = "RNA",min.cells = 0,
                                 min.features = 0,
                                 project = sp)
    sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
    sce[['peaks']] <- CreateChromatinAssay(counts = sc_counts$`Peaks`,
                                             annotation = annotations,
                                             fragments = path2,
                                             sep = c(":", "-"),
                                             genome = 'hg38')
    ############ 在ATAC矩阵下，进行质控 ################
DefaultAssay(sce)<-"peaks"
############ 核小体条带模式
sce <- NucleosomeSignal(object = sce)
############ 转录起始位点富集分数(TSS)
sce <- TSSEnrichment(object = sce, fast = FALSE)
############ 计算blacklist_fraction
sce$blacklist_fraction <- FractionCountsInRegion(object = sce, 
                                                         assay = 'peaks',
                                                         regions = blacklist_hg38)
sce<-subset(sce,subset = nCount_peaks < 60000 &
                       nCount_peaks > 1000 &
                       nFeature_peaks > 1000 &
                       nFeature_peaks < 25000 &
                       nCount_RNA < 20000 &
                       nCount_RNA > 500 &
                       nFeature_RNA > 100 &
                       nFeature_RNA < 7000 &
                       nucleosome_signal < 4 &
                       TSS.enrichment > 2 &
                       blacklist_fraction < 0.02 &
                       percent.mt < 5)
    print(sce)
DefaultAssay(sce) <- 'RNA'
sce <- NormalizeData(sce)
sce <- FindVariableFeatures(sce)
sce <- ScaleData(sce)
    scmulti.list <- append(scmulti.list,sce)
}

###Merge Seurat5 Object

combined <- merge(scmulti.list[[1]],scmulti.list[-1])
print(combined)

saveRDS(combined,'combined.rds')

###Unified peak list for all samples and then redo fragment counting based on the new peak list

peaks <- reduce(unlist(as(c(scmulti.list[[1]]@assays$peaks@ranges,
scmulti.list[[2]]@assays$peaks@ranges,
scmulti.list[[3]]@assays$peaks@ranges,
scmulti.list[[4]]@assays$peaks@ranges,
scmulti.list[[5]]@assays$peaks@ranges,
scmulti.list[[6]]@assays$peaks@ranges,
scmulti.list[[7]]@assays$peaks@ranges,
scmulti.list[[8]]@assays$peaks@ranges,
scmulti.list[[9]]@assays$peaks@ranges,
scmulti.list[[10]]@assays$peaks@ranges,
scmulti.list[[11]]@assays$peaks@ranges,
scmulti.list[[12]]@assays$peaks@ranges,
scmulti.list[[13]]@assays$peaks@ranges,
scmulti.list[[14]]@assays$peaks@ranges,
scmulti.list[[15]]@assays$peaks@ranges,
scmulti.list[[16]]@assays$peaks@ranges,
scmulti.list[[17]]@assays$peaks@ranges,
scmulti.list[[18]]@assays$peaks@ranges,
scmulti.list[[19]]@assays$peaks@ranges,
scmulti.list[[20]]@assays$peaks@ranges,
scmulti.list[[21]]@assays$peaks@ranges,
scmulti.list[[22]]@assays$peaks@ranges,
scmulti.list[[23]]@assays$peaks@ranges,
scmulti.list[[24]]@assays$peaks@ranges,
scmulti.list[[25]]@assays$peaks@ranges,
scmulti.list[[26]]@assays$peaks@ranges,
scmulti.list[[27]]@assays$peaks@ranges,
scmulti.list[[28]]@assays$peaks@ranges,
scmulti.list[[29]]@assays$peaks@ranges,
scmulti.list[[30]]@assays$peaks@ranges,
scmulti.list[[31]]@assays$peaks@ranges,
scmulti.list[[32]]@assays$peaks@ranges,
scmulti.list[[33]]@assays$peaks@ranges,
scmulti.list[[34]]@assays$peaks@ranges),
                          "GRangesList")))

peakwidths <- width(peaks)
peaks <- peaks[peakwidths < 10000 & peakwidths > 20]

main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
keep.peaks <- subset(peaks,seqnames(granges(peaks)) %in% main.chroms)

counts_atac_merged <- FeatureMatrix(seurat@assays$peaks@fragments,
                                    features = keep.peaks,
                                    cells = colnames(seurat))
seurat[['peaks']] <- CreateChromatinAssay(counts_atac_merged,
                                         fragments = seurat@assays$peaks@fragments,
                                         annotation = seurat@assays$peaks@annotation,
                                         sep = c(":","-"),
                                         genome = "hg38")

seurat <- JoinLayers(seurat)

###Save joinlayers for Pando
saveRDS(seurat,'Final_combined.rds')

###Split Matrix for CCA merge
obj <- readRDS('Final_combined.rds')
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$orig.ident)

obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- IntegrateLayers(
  object = obj, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = FALSE
)

###Calculate cca clusters
obj <- FindNeighbors(obj, dims = 1:20, reduction = "integrated.cca")
obj <- FindClusters(obj, resolution = 0.8, cluster.name = "cca_clusters")
obj <- RunUMAP(obj, reduction = "integrated.cca", dims = 1:20, reduction.name = "umap.cca")
p1 <- DimPlot(
  obj,
  reduction = "umap.cca",
  group.by = c("orig.ident",  "cca_clusters"),
  combine = FALSE, label.size = 2
)

levels(obj@meta.data$cca_clusters)

###Check markers from article

markers <- c(c('SLC17A7','SLC17A6','SLC32A1','TBR1','STMN2','MOBP','MOG','SOX2','OLIG1','OLIG2','SLC1A2','GFAP','AQP4','FLT1','VWF','PDGFRB','TREM2','GPR34','CD3E','CD8A'))
DotPlot(obj,features = markers)+coord_flip()

###Set new annotation
New_ids <- c('ODC','ODC','ODC','ODC','AS','MG','N','EC','ODC','OPC','T','AS','ODC','OPC','ODC','ODC','ODC','OPC','ODC','ODC','MG')

obj@meta.data$cell_type <- as.factor(obj@meta.data$cca_clusters)
levels(obj@meta.data$cell_type) <- New_ids
p2 <- DimPlot(
  obj,
  reduction = "umap.cca",
  group.by =   "cell_type",
  combine = FALSE, label.size = 2,label = T
)

###Check markers for new annotation
Idents(obj) <- obj@meta.data$cell_type
DotPlot(obj,features = markers)+coord_flip()

saveRDS(obj,'PD_anno.rds')

















