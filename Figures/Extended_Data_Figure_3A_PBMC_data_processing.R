library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)

# load the RNA and ATAC data
counts <- Read10X_h5("pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5")
fragpath <- "pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))

# create a Seurat object containing the RNA adata
pbmc <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# create ATAC assay and add it to the object
pbmc[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)

pbmc

DefaultAssay(pbmc) <- "ATAC"

pbmc <- NucleosomeSignal(pbmc)
pbmc <- TSSEnrichment(pbmc)

DensityScatter(pbmc, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)

VlnPlot(
  object = pbmc,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)

# filter out low quality cells
pbmc <- subset(
  x = pbmc,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 1800 &
    nCount_RNA > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1
)
pbmc

DefaultAssay(pbmc) <- "RNA"
pbmc <- SCTransform(pbmc)
pbmc <- RunPCA(pbmc)
var_genes <- VariableFeatures(pbmc)

DefaultAssay(pbmc) <- "ATAC"
pbmc <- FindTopFeatures(pbmc, min.cutoff = 5)
pbmc <- RunTFIDF(pbmc)
pbmc <- RunSVD(pbmc)
var_peaks <- VariableFeatures(pbmc)

library(SeuratDisk)

# load PBMC reference
reference <- LoadH5Seurat("pbmc_multimodal.h5seurat")
reference <- UpdateSeuratObject(reference)

DefaultAssay(pbmc) <- "SCT"

# transfer cell type labels from reference to query
transfer_anchors <- FindTransferAnchors(
  reference = reference,
  query = pbmc,
  normalization.method = "SCT",
  reference.reduction = "spca",
  recompute.residuals = FALSE,
  dims = 1:50
)

predictions <- TransferData(
  anchorset = transfer_anchors, 
  refdata = reference$celltype.l2,
  weight.reduction = pbmc[['pca']],
  dims = 1:50
)

pbmc <- AddMetaData(
  object = pbmc,
  metadata = predictions
)

# set the cell identities to the cell type predictions
Idents(pbmc) <- "predicted.id"

# remove low-quality predictions
pbmc <- pbmc[, pbmc$prediction.score.max > 0.5]

# build a joint neighbor graph using both assays
pbmc <- FindMultiModalNeighbors(
  object = pbmc,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 2:40),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

# build a joint UMAP visualization
pbmc <- RunUMAP(
  object = pbmc,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE
)

p <- DimPlot(pbmc, label = TRUE, repel = TRUE, reduction = "umap") + NoLegend()
ggsave('pbmc_predict_ct.pdf',p)
saveRDS(pbmc,'pbmc_scmulti_anno.rds')
write.table(var_genes,'var_genes.txt',quote=F,col.names = F,row.names = F)
write.table(var_peaks,'var_peaks.txt',quote=F,col.names = F,row.names = F)

#Select top 10 celltypes

cell_type_counts <- table(pbmc@meta.data$predicted.id)
sorted_cell_types <- sort(cell_type_counts, decreasing = TRUE)
top_10_cell_types <- names(sorted_cell_types)[1:10]

top_10_pbmc <- pbmc[, pbmc@meta.data$predicted.id %in% top_10_cell_types]

#Sample nCells to 5000
cell_type_counts <- table(top_10_pbmc@meta.data$predicted.id)
total_cells <- sum(cell_type_counts)
cells_to_sample <- round(cell_type_counts * 5000 / total_cells)

if (sum(cells_to_sample) < 5000) {
  while (sum(cells_to_sample) < 5000) {
    max_cell_type <- which.max(cells_to_sample)
    cells_to_sample[max_cell_type] <- cells_to_sample[max_cell_type] + 1
  }
}

sampled_cells <- NULL
for (ctype in names(cells_to_sample)) {
  cell_indices <- which(top_10_pbmc@meta.data$predicted.id == ctype)
  sampled_indices <- sample(cell_indices, cells_to_sample[ctype])
  if (is.null(sampled_cells)) {
    sampled_cells <- sampled_indices
  } else {
    sampled_cells <- c(sampled_cells, sampled_indices)
  }
}

sampled_top_10_pbmc <- top_10_pbmc[, sampled_cells]
saveRDS(sampled_top_10_pbmc,'pbmc_n5000_ct10.rds')








