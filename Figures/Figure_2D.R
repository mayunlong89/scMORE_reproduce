
###----------significant TFs---------regulon feature plot---------Figure 2D and Sup Figure S2C

# Load necessary libraries
library(Seurat)
library(ggplot2)
library(patchwork) # For combining plots
#---load pbmc_10x_real_data_downsampled_2000
pbmc_10x_real_data_downsampled_2000 <- readRDS("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/10X_PBMC/10X_PBMC_downsample_2000cells.rds")

DimPlot(pbmc_10x_real_data_downsampled_2000, reduction = "umap.new",cols = c("#67ADB7","#E77A77"))


single_cell <-pbmc_10x_real_data_downsampled_2000


# Load necessary libraries
library(Seurat)
library(ggplot2)
library(patchwork) # For combining plots

#Default use 'Signac' method for analysis
grn_outputs <- createRegulon(single_cell,
                             n_targets = 5,
                             peak2gene_method="Signac",
                             infer_method = 'glm')


regulons <- grn_outputs$grn

#1
Module <- regulons[regulons$TF == "ZEB1", ]
ZEB1_regulon <- c(unique(Module$TF), unique(Module$Target))
# Add module score for ZEB1 regulon
single_cell <- AddModuleScore(
  object = single_cell,
  features = list(ZEB1_regulon),  # Use only the genes present in the object
  name = "ZEB1_regulon"
)


FeaturePlot(object = single_cell, reduction = "umap.new", features = "ZEB1_regulon1")
library(RColorBrewer)
FeaturePlot(object = single_cell, reduction = "umap.new", features = "ZEB1_regulon1") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu")))

# 使用 FeaturePlot 并自定义颜色梯度
FeaturePlot(
  object = single_cell,
  reduction = "umap.new",
  features = "ZEB1_regulon1"
) +
  scale_colour_gradientn(
    colours = c("gray", "white", "red")  # 自定义颜色梯度
  ) +
  theme_classic()  # 可选：更改图形主题

# 使用 FeaturePlot 并自定义颜色梯度
FeaturePlot(
  object = single_cell,
  reduction = "umap.new",
  features = "ZEB1"
) +
  scale_colour_gradientn(
    colours = c("gray", "white", "red")  # 自定义颜色梯度
  ) +
  theme_classic()  # 可选：更改图形主题


#2
Module <- regulons[regulons$TF == "FOXO1", ]
FOXO1_regulon <- c(unique(Module$TF), unique(Module$Target))
# Add module score for FOXO1 regulon
single_cell <- AddModuleScore(
  object = single_cell,
  features = list(FOXO1_regulon),  # Use only the genes present in the object
  name = "FOXO1_regulon"
)

library(RColorBrewer)
FeaturePlot(object = single_cell, reduction = "umap.new", features = "FOXO1_regulon1") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu")))

# 使用 FeaturePlot 并自定义颜色梯度
FeaturePlot(
  object = single_cell,
  reduction = "umap.new",
  features = "FOXO1_regulon1"
) +
  scale_colour_gradientn(
    colours = c("gray", "white", "red")  # 自定义颜色梯度
  ) +
  theme_classic()  # 可选：更改图形主题

# 使用 FeaturePlot 并自定义颜色梯度
FeaturePlot(
  object = single_cell,
  reduction = "umap.new",
  features = "FOXO1"
) +
  scale_colour_gradientn(
    colours = c("gray", "white", "red")  # 自定义颜色梯度
  ) +
  theme_classic()  # 可选：更改图形主题


#3
Module <- regulons[regulons$TF == "IKZF1", ]
IKZF1_regulon <- c(unique(Module$TF), unique(Module$Target))
# Add module score for IKZF1 regulon
single_cell <- AddModuleScore(
  object = single_cell,
  features = list(IKZF1_regulon),  # Use only the genes present in the object
  name = "IKZF1_regulon"
)

library(RColorBrewer)
FeaturePlot(object = single_cell, reduction = "umap.new", features = "IKZF1_regulon1") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu")))

FeaturePlot(object = single_cell, reduction = "umap.new", features = "IKZF1") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu")))

# 使用 FeaturePlot 并自定义颜色梯度
FeaturePlot(
  object = single_cell,
  reduction = "umap.new",
  features = "IKZF1"
) +
  scale_colour_gradientn(
    colours = c("gray", "white", "red")  # 自定义颜色梯度
  ) +
  theme_classic()  # 可选：更改图形主题

#4
Module <- regulons[regulons$TF == "ZNF721", ]
ZNF721_regulon <- c(unique(Module$TF), unique(Module$Target))
# Add module score for ZNF721 regulon
single_cell <- AddModuleScore(
  object = single_cell,
  features = list(ZNF721_regulon),  # Use only the genes present in the object
  name = "ZNF721_regulon"
)

library(RColorBrewer)
FeaturePlot(object = single_cell, reduction = "umap.new", features = "ZNF721_regulon1") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu")))

FeaturePlot(object = single_cell, reduction = "umap.new", features = "ZNF721") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu")))

# 使用 FeaturePlot 并自定义颜色梯度
FeaturePlot(
  object = single_cell,
  reduction = "umap.new",
  features = "ZNF721"
) +
  scale_colour_gradientn(
    colours = c("gray", "white", "red")  # 自定义颜色梯度
  ) +
  theme_classic()  # 可选：更改图形主题

# 使用 FeaturePlot 并自定义颜色梯度
FeaturePlot(
  object = single_cell,
  reduction = "umap.new",
  features = "ZNF721_regulon1"
) +
  scale_colour_gradientn(
    colours = c("gray", "white", "red")  # 自定义颜色梯度
  ) +
  theme_classic()  # 可选：更改图形主题


#5
Module <- regulons[regulons$TF == "ETS1", ]
ETS1_regulon <- c(unique(Module$TF), unique(Module$Target))
# Add module score for FOXO1 regulon
single_cell <- AddModuleScore(
  object = single_cell,
  features = list(ETS1_regulon),  # Use only the genes present in the object
  name = "ETS1_regulon"
)

library(RColorBrewer)
FeaturePlot(object = single_cell, reduction = "umap.new", features = "ETS1_regulon1") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu")))

# 使用 FeaturePlot 并自定义颜色梯度
FeaturePlot(
  object = single_cell,
  reduction = "umap.new",
  features = "ETS1_regulon1"
) +
  scale_colour_gradientn(
    colours = c("gray", "white", "red")  # 自定义颜色梯度
  ) +
  theme_classic()  # 可选：更改图形主题

# 使用 FeaturePlot 并自定义颜色梯度
FeaturePlot(
  object = single_cell,
  reduction = "umap.new",
  features = "ETS1"
) +
  scale_colour_gradientn(
    colours = c("gray", "white", "red")  # 自定义颜色梯度
  ) +
  theme_classic()  # 可选：更改图形主题


#6
Module <- regulons[regulons$TF == "FOXP1", ]
FOXP1_regulon <- c(unique(Module$TF), unique(Module$Target))
# Add module score for FOXO1 regulon
single_cell <- AddModuleScore(
  object = single_cell,
  features = list(FOXP1_regulon),  # Use only the genes present in the object
  name = "FOXP1_regulon"
)

library(RColorBrewer)
FeaturePlot(object = single_cell, reduction = "umap.new", features = "FOXP1_regulon1") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu")))

# 使用 FeaturePlot 并自定义颜色梯度
FeaturePlot(
  object = single_cell,
  reduction = "umap.new",
  features = "FOXP1_regulon1"
) +
  scale_colour_gradientn(
    colours = c("gray", "white", "red")  # 自定义颜色梯度
  ) +
  theme_classic()  # 可选：更改图形主题

# 使用 FeaturePlot 并自定义颜色梯度
FeaturePlot(
  object = single_cell,
  reduction = "umap.new",
  features = "FOXP1"
) +
  scale_colour_gradientn(
    colours = c("gray", "white", "red")  # 自定义颜色梯度
  ) +
  theme_classic()  # 可选：更改图形主题





#7
Module <- regulons[regulons$TF == "LEF1", ]
LEF1_regulon <- c(unique(Module$TF), unique(Module$Target))
# Add module score for FOXO1 regulon
single_cell <- AddModuleScore(
  object = single_cell,
  features = list(LEF1_regulon),  # Use only the genes present in the object
  name = "LEF1_regulon"
)

library(RColorBrewer)
FeaturePlot(object = single_cell, reduction = "umap.new", features = "LEF1_regulon1") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu")))

# 使用 FeaturePlot 并自定义颜色梯度
FeaturePlot(
  object = single_cell,
  reduction = "umap.new",
  features = "LEF1_regulon1"
) +
  scale_colour_gradientn(
    colours = c("gray", "white", "red")  # 自定义颜色梯度
  ) +
  theme_classic()  # 可选：更改图形主题

# 使用 FeaturePlot 并自定义颜色梯度
FeaturePlot(
  object = single_cell,
  reduction = "umap.new",
  features = "LEF1"
) +
  scale_colour_gradientn(
    colours = c("gray", "white", "red")  # 自定义颜色梯度
  ) +
  theme_classic()  # 可选：更改图形主题

