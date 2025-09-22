
##-------Figure 4H-----------sup_Figure_S18AB--50cellstates----------------------------------------------------------------------------------
# Cell state inferences

###Slingshot for monocytes trajectory inference------------------------------------
suppressPackageStartupMessages({
  library(scran)
  library(scater)
  library(igraph)
  library(slingshot)
})
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
#BiocManager::install("kstreet13/slingshot")

#https://bioconductor.org/packages/devel/bioc/vignettes/slingshot/inst/doc/vignette.html

###using above prossed monocyte subset single-cell data
#single_cell_mono <- subset(single_cell, idents="Monocytes")
table(single_cell_mono$cell_type)

sce  <- as.SingleCellExperiment(single_cell_mono)

#data("slingshotExample")

library(slingshot, quietly = FALSE)    





#Gene filtering
# at least M (15) reads in at least N (15) cells
geneFilter <- apply(assays(sce)$counts,1,function(x){
  sum(x >= 3) >= 10
})
sce <- sce[geneFilter, ]

#Normalization
FQnorm <- function(counts){   
  rk <- apply(counts,2,rank,ties.method='min')  
  counts.sort <- apply(counts,2,sort)  
  refdist <- apply(counts.sort,1,median)   
  norm <- apply(rk,2,function(r){ refdist[r] })    
  rownames(norm) <- rownames(counts)   
  return(norm)
}
assays(sce)$norm <- FQnorm(assays(sce)$counts)

#Dimensionality Reduction

pca <- prcomp(t(log1p(assays(sce)$norm)), scale. = FALSE)
rd1 <- pca$x[,1:2]

plot(rd1, col = rgb(0,0,0,.5), pch=16, asp = 1)


library(uwot)
rd2 <- uwot::umap(t(log1p(assays(sce)$norm)))
colnames(rd2) <- c('UMAP1', 'UMAP2')

plot(rd2, col = rgb(0,0,0,.5), pch=16, asp = 1)

reducedDims(sce) <- SimpleList(PCA = rd1, UMAP = rd2)

#clustering Cells
library(mclust, quietly = TRUE)

#Gaussian mixture modeling
cl1 <- Mclust(rd1)$classification
colData(sce)$GMM <- cl1

library(RColorBrewer)
plot(rd1, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 1)

#K-means
cl2 <- kmeans(rd1, centers = 4)$cluster
colData(sce)$kmeans <- cl2

plot(rd1, col = brewer.pal(9,"Set1")[cl2], pch=16, asp = 1)

#using Slingshot
sce <- slingshot(sce, clusterLabels = 'GMM', reducedDim = 'PCA')
summary(sce$slingPseudotime_1)

library(grDevices)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]

plot(reducedDims(sce)$PCA, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, col='black')

plot(reducedDims(sce)$PCA, col = brewer.pal(9,'Set1')[sce$GMM], pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, type = 'lineages', col = 'black')

sce$slingPseudotime_1
colnames(sce)
colnames(single_cell_mono)


##adding pseudotime to seurat data
single_cell_mono$pseudotime <-  sce$slingPseudotime_1

##pseudotime feature plot
FeaturePlot(single_cell_mono,features = "pseudotime")

##@@@@@=====----Suppl_Figure_S14A
library(scCustomize)
FeaturePlot_scCustom(seurat_object=single_cell_mono,
                     reduction = "umap.rna",
                     features = "pseudotime",
                     colors_use=viridis_plasma_dark_high,
                     na_color = "lightgray")

#extracting pseudotimes, cell names
x1 <- as.matrix(single_cell_mono$pseudotime)
x2 <- as.matrix(colnames(single_cell_mono))
data_mono <- cbind(x1,x2)
data_mono <- as.data.frame(data_mono)
colnames(data_mono) <- c("pseudotimes","cell_ID")
head(data_mono)

#取bins ntile()默认升序
data_mono$pseudotimes <- as.numeric(data_mono$pseudotimes)
binning <- data_mono %>% mutate(rank=ntile(data_mono$pseudotimes,50))
#binning <- data_mono %>% arrange(desc(data_mono$pseudotimes)) %>% mutate(rank=ntile(data_mono$pseudotimes,50))



#xx <- aggregate(binning,by = list(binning$rank), mean)

binning <- binning[which(binning$rank %in% c(1:50)),]

cor(as.numeric(binning$pseudotimes), as.numeric(binning$rank))
plot(as.numeric(binning$rank),as.numeric(binning$pseudotimes))


table(binning$rank)

#Defining 100 cell states based on the pseudotimes
single_cell_mono$cell_state <- binning$rank


##--------------boxplots for bins and pseudotimes----Suppl_Figure_S14B
##ggplot
ggplot(binning,aes(x=as.factor(rank),y=pseudotimes))+
  geom_boxplot()+
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9","#999999", "#E69F00"))+
  #theme(legend.position="right")+
  # labs(title="Plot of length  per dose",x="Dose (mg)", y = "Length")+
  #geom_dotplot(binaxis='y', stackdir='center', stackratio=1.0, dotsize=0.5)+
  theme_classic()



#length(binning$rank)

Idents(single_cell_mono) <- single_cell_mono$cell_state
DimPlot(single_cell_mono,label = F,
        reduction = "umap.rna", 
        pt.size = 1, repel = F)

single_cell_mono1 <- subset(single_cell_mono,idents = c(1:50))
Idents(single_cell_mono1) <- single_cell_mono1$cell_state
DimPlot(single_cell_mono1,label = F,
        reduction = "umap.rna", 
        pt.size = 1, repel = F)

##pseudotime feature plot
FeaturePlot(single_cell_mono,features = "pseudotime")

FeaturePlot_scCustom(seurat_object=single_cell_mono,
                     reduction = "umap.rna",
                     features = "pseudotime",
                     colors_use=viridis_plasma_dark_high,
                     na_color = "lightgray")

FeaturePlot_scCustom(seurat_object=single_cell_mono,
                     reduction = "umap.rna",
                     features = "pseudotime",
                     colors_use=scCustomize_Palette(num_groups = 2),
                     na_color = "lightgray")

FeaturePlot_scCustom(seurat_object=single_cell_mono,
                     reduction = "umap.rna",
                     features = "pseudotime",
                     colors_use=c("navy","orange"),
                     na_color = "lightgray")

FeaturePlot_scCustom(seurat_object=single_cell_mono,
                     reduction = "umap.rna",
                     features = "pseudotime",
                     colors_use=viridis_inferno_light_high,
                     na_color = "lightgray")

FeaturePlot_scCustom(seurat_object=single_cell_mono,
                     reduction = "umap.rna",
                     features = "pseudotime",
                     colors_use=viridis_dark_high,
                     na_color = "lightgray")



