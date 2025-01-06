

# Load necessary libraries
library(Seurat)
library(ggplot2)
library(patchwork) # For combining plots
#---load pbmc_10x_real_data_downsampled_2000
pbmc_10x_real_data_downsampled_2000 <- readRDS("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/10X_PBMC/10X_PBMC_downsample_2000cells.rds")

DimPlot(pbmc_10x_real_data_downsampled_2000, reduction = "umap.new",cols = c("#67ADB7","#E77A77"))


single_cell <-pbmc_10x_real_data_downsampled_2000

