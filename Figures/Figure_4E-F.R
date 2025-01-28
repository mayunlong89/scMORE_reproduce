###------Figure 4E

### Subset Monocytes

single_cell_mono <- subset(single_cell_8900, idents="Monocytes")

table(single_cell_mono$cell_type)

DimPlot(single_cell_mono,reduction = "umap.biMod",label = T, cols=pal, pt.size = 0.1, repel = T)


#re-clustering the UMAP
single_cell_mono <- NormalizeData(single_cell_mono)
single_cell_mono <- FindVariableFeatures(single_cell_mono)
single_cell_mono <- ScaleData(single_cell_mono)
single_cell_mono <- RunPCA(single_cell_mono,reduction.name = "pca")
single_cell_mono <- FindNeighbors(single_cell_mono,
                                  dims = 1:20,reduction = "pca")
single_cell_mono <- FindClusters(single_cell_mono,
                                 resolution = 1)
#Run UMAP
single_cell_mono <- RunUMAP(single_cell_mono, dims = 1:20,
                            reduction = "pca",reduction.name = "umap.new")

#Visualization
Idents(single_cell_mono) <- single_cell_mono$cell_type

#DimPlot(single_cell_mono, reduction = "umap.biMod")
#DimPlot(single_cell_mono, reduction = "umap.new" )
DimPlot(single_cell_mono,label = T, cols=c("#A13B46","#67ADB7"), pt.size = 0.1, repel = T)
DimPlot(single_cell_mono,label = T, cols=c("#67ADB7","#A13B46"), pt.size = 0.1, repel = T)
DimPlot(single_cell_mono,label = T,reduction = "umap.atac", cols=c("#67ADB7","#A13B46"), pt.size = 0.1, repel = T)
DimPlot(single_cell_mono,label = T,reduction = "umap.rna", cols=c("#67ADB7","#A13B46"), pt.size = 0.1, repel = T)

#Idents(single_cell_mono) <- single_cell_mono$RNA_snn_res.0.5
DimPlot(single_cell_mono,label = T, pt.size = 0.1, repel = T)

Idents(single_cell_mono) <- single_cell_mono$RNA_snn_res.1
table(single_cell_mono$RNA_snn_res.1)
DimPlot(single_cell_mono,label = T, pt.size = 0.1, repel = T)


#CD14+ Monocytes,FCGR3A+ Monocytes

##------------Figure---4F
FeaturePlot(single_cell_mono,features = "CD14") # classical monocytes
FeaturePlot(single_cell_mono,features = "FCGR3A") #non-classical monocytes

