
###------------------------------------------------------------------------------------
#2025-01-09
###------------------------------------------------------------------------------------

###required load R packages
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
})


#devtools::install_github("mayunlong89/scMORE")

library(scMORE)
library(GenomicRanges)
library(IRanges)
library(Seurat)
library(Signac)


#load data on 10x pbmc example data
#Load scMultiomic data
pbmc_10x <- readRDS("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/10X_PBMC/10X_PBMC.rds")

length(unique(pbmc_10x$cell_type))

single_cell_8900<-pbmc_10x


grn_outputs8900 <- createRegulon(single_cell_8900, n_targets=5,
                             peak2gene_method = 'Signac',
                             infer_method = 'glm')

grn_outputs8900_glmnet <- createRegulon(single_cell_8900, n_targets = 5,
                              peak2gene_method="Signac",
                              infer_method = "glmnet")

grn_outputs8900_xgb <- createRegulon(single_cell_8900, n_targets=5,
                              peak2gene_method = 'Signac',
                              infer_method = 'xgb')



#load data on 10x human brain data
#Load scMultiomic data
human_brain_10x <- readRDS("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/10X_human_brain/10X_Human_Brain.rds")

length(unique(human_brain_10x$cell_type))

single_cell_2464   <- human_brain_10x

grn_outputs2464_glm <- createRegulon(single_cell_2464, n_targets=5,
                                 peak2gene_method = 'Signac',
                                 infer_method = 'glm')

grn_outputs8900_glmnet <- createRegulon(single_cell_2464, n_targets = 5,
                                        peak2gene_method="Signac",
                                        infer_method = "glmnet")

grn_outputs2464_xgb <- createRegulon(single_cell_2464, n_targets=5,
                                     peak2gene_method = 'Signac',
                                     infer_method = 'xgb')



#load data on 10x human B Cell lymphoma
#Load scMultiomic data
human_Bcell_10x <- readRDS("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/10X_B_Cell_Lymphoma/10X_B_Cell_Lymphoma.rds")

length(unique(human_Bcell_10x$cell_type))

single_cell_10478 <- human_Bcell_10x


grn_outputs10478_glm <- createRegulon(single_cell_10478, n_targets=5,
                                     peak2gene_method = 'Signac',
                                     infer_method = 'glm')

grn_outputs10478_glmnet <- createRegulon(single_cell_10478, n_targets = 5,
                                        peak2gene_method="Signac",
                                        infer_method = "glmnet")

grn_outputs10478_xgb <- createRegulon(single_cell_10478, n_targets=5,
                                     peak2gene_method = 'Signac',
                                     infer_method = 'xgb')




#load data on 10x human PD Brain data
#Load scMultiomic data
human_PD_brain_10x <- readRDS("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/06-PD/PD_anno.rds")

human_PD_brain_10x <- JoinLayers(human_PD_brain_10x)

length(unique(human_PD_brain_10x$cell_type))

single_cell_58949  <- human_PD_brain_10x

grn_outputs58949_glm <- createRegulon(single_cell_58949, n_targets=5,
                                      peak2gene_method = 'Signac',
                                      infer_method = 'glm')

grn_outputs58949_glmnet <- createRegulon(single_cell_58949, n_targets = 5,
                                         peak2gene_method="Signac",
                                         infer_method = "glmnet")

grn_outputs58949_xgb <- createRegulon(single_cell_58949, n_targets=5,
                                      peak2gene_method = 'Signac',
                                      infer_method = 'xgb')




#load data on 10x human organoid eye cells
#Load scMultiomic data
human_organoid_eye_10x <- readRDS("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/01_five_eye_diseases_and_organoids/Final_anno_organoids.rds")
Idents(human_organoid_eye_10x) <- human_organoid_eye_10x$Final_anno

length(unique(human_organoid_eye_10x$Final_anno))

single_cell_2827 <- human_organoid_eye_10x


grn_outputs2827_glm <- createRegulon(single_cell_2827, n_targets=5,
                                      peak2gene_method = 'Signac',
                                      infer_method = 'glm')

grn_outputs2827_glmnet <- createRegulon(single_cell_2827, n_targets = 5,
                                         peak2gene_method="Signac",
                                         infer_method = "glmnet")

grn_outputs2827_xgb <- createRegulon(single_cell_2827, n_targets=5,
                                      peak2gene_method = 'Signac',
                                      infer_method = 'xgb')




#Visualization of plot------------------------------------------------------------------
# Load necessary libraries
library(dplyr)
library(ggplot2)


cell2regulon_Number <- read.table("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/05_glm_glmnet_xgb_tissues_cell_number_regulons_count.txt",header = T)


head(cell2regulon_Number)

ggplot(cell2regulon_Number, aes(x = method, y = Regulon_counts, fill = method)) +
  geom_violin(trim = FALSE, position = position_dodge(0.8)) + # 小提琴图
  geom_point(aes(color = Tissues), position = position_dodge(0.8), size = 1.5, alpha = 0.8) + # 数据点，用Tissues区分颜色
  scale_fill_manual(values = c( "#BFEFFF", "#EEDFCC","#FFDAB9")) + # 小提琴填充色
  scale_color_manual(values = c("orange", "blue", "#668B8B", "purple", "brown")) + # 数据点颜色映射
    labs(x = "", y = "Regulon counts", fill = "Method", color = "Tissues") + # 添加图例标签
  theme_minimal() +
  theme(legend.position = "right") # 图例放置在右侧



# Generate the violin plot
ggplot(cell2regulon_Number, aes(x = method, y = Regulon_counts, fill = method)) +
  geom_violin(trim = FALSE, position = position_dodge(0.8)) + # Adjust position for side-by-side violins
  geom_jitter(position = position_jitter(width = 0.2), size = 1.5, alpha = 0.8) + # 数据点
 # geom_line(aes(group = Tissues), color = "gray50", alpha = 0.6) + # 连线
  scale_fill_manual(values = c("#CDB38B", "orange","#668B8B")) + # Match colors
  labs(x = "theta parameter", y = "E score", fill = "method") +
  theme_minimal() +
  theme(legend.position = "right") # Place legend on the right


ggplot(cell2regulon_Number, aes(x = method, y = Regulon_counts, fill = method)) +
  geom_violin(trim = FALSE, position = position_dodge(0.8)) + # 小提琴图
  geom_point(position = position_dodge(0.8), size = 1.5, alpha = 0.8) + # 数据点
  #geom_line(aes(group = Tissues), position = position_dodge(0.8), color = "gray50", alpha = 0.6) + # 连线
  scale_fill_manual(values = c("#CDB38B", "orange", "#668B8B")) + # 配色
  labs(x = "theta parameter", y = "E score", fill = "method") +
  theme_minimal() +
  theme(legend.position = "right") # 图例放置在右侧



###------Significant regulons---Proportion

sig_regulon_proportion <- read.table("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/05_glm_glmnet_xgb_proportion.txt",header = T)


head(sig_regulon_proportion)

# Generate the violin plot
ggplot(sig_regulon_proportion, aes(x = Method, y = proportion, fill = group)) +
  geom_violin(trim = FALSE, position = position_dodge(0.8)) + # Adjust position for side-by-side violins
  geom_point(position = position_dodge(0.8), size = 1.5, alpha = 0.8) + # 数据点
  scale_fill_manual(values = c("blue", "orange")) + # Match colors
  labs(x = "theta parameter", y = "E score", fill = "method") +
  theme_minimal() +
  theme(legend.position = "right") # Place legend on the right



##------Significant regulons identified by GLM show in GLMNET---based on cosine

##-2
consine_average10traits_heatmap_glmnet <- read.table("01_10blood_cell_traits_sig_heatmap_glmnet.txt",header = TRUE,sep = "\t")

library(pheatmap)

data_h2 <- consine_average10traits_heatmap_glmnet[,c(-1)]
rownames(data_h2) <- consine_average10traits_heatmap_glmnet[,1]
pearson_heatmap<-pheatmap(data_h2,cellwidth =15, cellheight =15,fontsize=8,
                          border_color="black",scale="row",
                          fontsize_row=8,
                          color=colorRampPalette(rev(c("red","#E77A77", "#ECBA84","white","lightblue",'#336699')))(102),
                          cluster_rows = T,
                          cluster_cols = T)





##------Significant regulons identified by GLM show in XGBoost---based on cosine

##-3
consine_average10traits_heatmap_xgb <- read.table("01_10blood_cell_traits_sig_heatmap_xgb.txt",header = TRUE,sep = "\t")

library(pheatmap)

data_h3 <- consine_average10traits_heatmap_xgb[,c(-1)]
rownames(data_h3) <- consine_average10traits_heatmap_xgb[,1]
pearson_heatmap<-pheatmap(data_h3,cellwidth =15, cellheight =15,fontsize=8,
                          border_color="black",scale="row",
                          fontsize_row=8,
                          color=colorRampPalette(rev(c("red","#E77A77", "#ECBA84","white","lightblue",'#336699')))(102),
                          cluster_rows = T,
                          cluster_cols = T)





##-1
consine_average10traits_heatmap_glm <- read.table("01_10blood_cell_traits_sig_heatmap.txt",header = TRUE,sep = "\t")

library(pheatmap)

data_h1 <- consine_average10traits_heatmap_glm[,c(-1)]
rownames(data_h1) <- consine_average10traits_heatmap_glm[,1]
pearson_heatmap<-pheatmap(data_h1,cellwidth =15, cellheight =15,fontsize=8,
                          border_color="black",scale="row",
                          fontsize_row=8,
                          color=colorRampPalette(rev(c("red","#E77A77", "#ECBA84","white","lightblue",'#336699')))(102),
                          cluster_rows = T,
                          cluster_cols = T)

##-----------------------------------------------------------------------------------------------
##------glm_glmnet_xgb data transformed for correlation-----------------------

library(reshape2)

consine_average10traits_heatmap_glm

consine_average10traits_heatmap_xgb

consine_average10traits_heatmap_glmnet


md_glm <- melt(consine_average10traits_heatmap_glm, id=c("regulons"))
md_xgb <- melt(consine_average10traits_heatmap_xgb, id=c("regulons"))
md_glmnet <- melt(consine_average10traits_heatmap_glmnet, id=c("regulons"))

data <- md_glm

colnames(data) <- c("regulons","traits","glm_TRS")
head(data)
data$xgb_TRS <- md_xgb$value
data$glmnet_TRS <- md_glmnet$value


cor.test(data$glm_TRS,data$xgb_TRS)
cor.test(data$glm_TRS,data$glmnet_TRS)
cor.test(data$glmnet_TRS,data$xgb_TRS)



ggplot(data, aes(x = glm_TRS, y = xgb_TRS,color=traits)) +
  geom_point(size = 2, alpha = 0.7) +  # 散点
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # 添加拟合线
  labs(title = "Correlation by Cell Type",
       x = "glm_TRS",
       y = "xgb_TRS") +
  theme_classic() +  # 经典主题
  facet_wrap(~ traits)  # 按细胞类型分面

# Perform correlation test for each cell type
cor_results <- data %>%
  dplyr::group_by(traits) %>%
  dplyr::summarise(
    correlation = cor(glm_TRS, xgb_TRS, use = "complete.obs"),  # Pearson correlation
    p_value = cor.test(glm_TRS, xgb_TRS)$p.value,  # Extract p-value
    conf_low = cor.test(glm_TRS, xgb_TRS)$conf.int[1],  # Lower bound of confidence interval
    conf_high = cor.test(glm_TRS, xgb_TRS)$conf.int[2]   # Upper bound of confidence interval
  )

print(cor_results)


ggplot(data, aes(x = glm_TRS, y = xgb_TRS,color="red")) +
  geom_point(size = 2, alpha = 0.7) +  # 散点
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed",color="blue") +  # 添加拟合线
  labs(title = "Correlation by Cell Type",
       x = "glm_TRS",
       y = "xgb_TRS") +
  theme_classic()

cor_results2 <- data %>%
  dplyr::summarise(
    correlation = cor(glm_TRS, xgb_TRS, use = "complete.obs"),  # Pearson correlation
    p_value = cor.test(glm_TRS, xgb_TRS)$p.value,  # Extract p-value
    conf_low = cor.test(glm_TRS, xgb_TRS)$conf.int[1],  # Lower bound of confidence interval
    conf_high = cor.test(glm_TRS, xgb_TRS)$conf.int[2]   # Upper bound of confidence interval
  )
print(cor_results2)




ggplot(data, aes(x = glm_TRS, y = glmnet_TRS,color="red")) +
  geom_point(size = 2, alpha = 0.7) +  # 散点
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed",color="blue") +  # 添加拟合线
  labs(title = "Correlation by Cell Type",
       x = "glm_TRS",
       y = "glmnet_TRS") +
  theme_classic()

cor_results3 <- data %>%
  dplyr::summarise(
    correlation = cor(glm_TRS, glmnet_TRS, use = "complete.obs"),  # Pearson correlation
    p_value = cor.test(glm_TRS, glmnet_TRS)$p.value,  # Extract p-value
    conf_low = cor.test(glm_TRS, glmnet_TRS)$conf.int[1],  # Lower bound of confidence interval
    conf_high = cor.test(glm_TRS, glmnet_TRS)$conf.int[2]   # Upper bound of confidence interval
  )
print(cor_results3)



ggplot(data, aes(x = glm_TRS, y = glmnet_TRS,color=traits)) +
  geom_point(size = 2, alpha = 0.7) +  # 散点
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # 添加拟合线
  labs(title = "Correlation by Cell Type",
       x = "glm_TRS",
       y = "glmnet_TRS") +
  theme_classic() +  # 经典主题
  facet_wrap(~ traits)  # 按细胞类型分面

# Perform correlation test for each cell type
cor_results <- data %>%
  dplyr::group_by(traits) %>%
  dplyr::summarise(
    correlation = cor(glm_TRS, glmnet_TRS, use = "complete.obs"),  # Pearson correlation
    p_value = cor.test(glm_TRS, glmnet_TRS)$p.value,  # Extract p-value
    conf_low = cor.test(glm_TRS, glmnet_TRS)$conf.int[1],  # Lower bound of confidence interval
    conf_high = cor.test(glm_TRS, glmnet_TRS)$conf.int[2]   # Upper bound of confidence interval
  )

print(cor_results)



##----heatmap_for lymp_count and lymp_percent results in CD8+T cells across three GRN methods

##-1
consine_heatmap_glm_glmnet_xgb <- read.table("01_lymp_count_percent_glm_glmnet_xgb_heatmap.txt",header = TRUE,sep = "\t")

library(pheatmap)

data_h4 <- consine_heatmap_glm_glmnet_xgb[,c(-1)]
rownames(data_h4) <- consine_heatmap_glm_glmnet_xgb[,1]
pearson_heatmap<-pheatmap(data_h4,cellwidth =15, cellheight =15,fontsize=8,
                          border_color="black",scale="none",
                          fontsize_row=8,
                          color=colorRampPalette(rev(c("red","#E77A77", "#ECBA84","white",  "lightblue","lightblue")))(102),
                          cluster_rows = F,
                          cluster_cols = F)

# 将数据从宽格式转换为长格式
data_long <- melt(consine_heatmap_glm_glmnet_xgb, id.vars = "regulons",
                  variable.name = "traits", value.name = "TRS")

# 绘制 dot plot
ggplot(data_long, aes(x = regulons, y = traits, size = TRS, color = "red")) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(x = "Regulons", y = "Traits", size = "TRS", color = "Model") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



