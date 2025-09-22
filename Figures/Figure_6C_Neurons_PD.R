##-----------Figure 7C 

setwd("~/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/06-PD/")

###---pando_results for brain organoids single-cell data

PD_seu0 <- readRDS("./PD_anno.rds")

PD_seu <- JoinLayers(PD_seu0)


##read metadata

meta_anno <- read.csv("01_metadata_PD.csv", header=T)
head(meta_anno)

PD_seu$group <- meta_anno$group[match(PD_seu$orig.ident, meta_anno$ID)]
PD_seu$age <- meta_anno$age[match(PD_seu$orig.ident, meta_anno$ID)]
PD_seu$sex <- meta_anno$sex[match(PD_seu$orig.ident, meta_anno$ID)]


grn_outputs <- readRDS("grn_outputs_PD.rds")

regulons <- grn_outputs$grn
head(regulons)

##------Figure 7C
library(dplyr)
# Aging-relevant regulons for Ns
Aging_regulons_Ns <- c(
  "DMBX1", "BCL11A", "STAT4", "CUX2", "DACH1", "MYT1L", 
  "NPAS2", "TSHZ3", "PKNOX2", "ZFPM2", "EBF1", "SOX11", 
  "TSHZ1", "TOX", "THRB", "TOX2", "ETV1", "ZNF423", 
  "LEF1", "ZEB1", "RORA", "TRERF1"
)

# Initialize an empty data frame for storing results with predefined columns
results <- data.frame(
  cell_type = character(),
  group = character(),
  sex = character(),
  Module_regulon = numeric(),
  Regulon = character(),
  stringsAsFactors = FALSE
)

# Loop through each regulon
for (i in seq_along(Aging_regulons_Ns)) {
  # Filter regulons by the current TF
  Module <- regulons[regulons$TF == Aging_regulons_Ns[i], ]
  
  # Combine TF and its target genes
  Module_regulon <- c(unique(Module$TF), unique(Module$Target))
  
  # Calculate addModuleScore for each regulon
  DefaultAssay(PD_seu) <- "RNA"
  PD_seu <- AddModuleScore(
    PD_seu,
    features = list(Module_regulon),
    name = "Module_regulon"  # Unique name for each regulon
  )
  
  # Extract metadata
  add_cell <- as.data.frame(PD_seu@meta.data)
  
  # Calculate the mean regulon activity by cell_type, group, and sex
  temp <- add_cell %>%
    group_by(cell_type, group, sex) %>%
    summarise(
      Module_regulon = mean(get(paste0("Module_regulon", 1)), na.rm = TRUE),
      .groups = "drop"
    )
  
  # Add a column to indicate the regulon name
  temp$Regulon <- Aging_regulons_Ns[i]
  
  # Append the results
  results <- rbind(results, temp)
}

# View the final results
print(results)

results$group_sex <- paste0(results$group,"_",results$sex)

results_Ns <- results[which(results$cell_type=="N"),]
  


# Save the results to a CSV file
write.csv(results, "02_Aging_regulon_activity_by_group_and_sex.csv", row.names = FALSE)


library(tibble)
library(pheatmap)
library(dplyr)

heatmap_data <- results_Ns %>%
  dplyr::select(group_sex, Regulon, Module_regulon) %>%
  pivot_wider(names_from = Regulon, values_from = Module_regulon)

heatmap_data <- as.data.frame(heatmap_data)

# 替代 column_to_rownames 方法
rownames(heatmap_data) <- heatmap_data[,1]
heatmap_data <- heatmap_data[, -1]  # 删除 group_sex 列
heatmap_data <- as.matrix(heatmap_data)

# 绘制热图

pearson_heatmap<-pheatmap(heatmap_data,cellwidth =12, cellheight = 12,fontsize=8,
                          border_color="black",
                          scale="column",
                          fontsize_row=8,
                          color=colorRampPalette(rev(c( "firebrick3","#E77A77", "#ECBA84","white","lightblue",'#336699')))(102),
                          cluster_rows = T,
                          cluster_cols = T,
                           display_numbers = F)
                      

heatmap_data1 <- t(heatmap_data)

pearson_heatmap<-pheatmap(heatmap_data1,cellwidth =12, cellheight = 12,fontsize=8,
                          border_color="black",
                          scale="row",
                          fontsize_row=8,
                          color=colorRampPalette(rev(c( "firebrick3","#E77A77", "#ECBA84","white","lightblue",'#336699')))(102),
                          cluster_rows = T,
                          cluster_cols = F,
                          display_numbers = F)

pearson_heatmap<-pheatmap(heatmap_data1,cellwidth =12, cellheight = 12,fontsize=8,
                          border_color="black",
                          scale="none",
                          fontsize_row=8,
                          color=colorRampPalette(rev(c( "firebrick3","#E77A77", "#ECBA84","white","lightblue")))(102),
                          cluster_rows = T,
                          cluster_cols = F,
                          display_numbers = F)

