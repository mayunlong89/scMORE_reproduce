
#---2025-06-26


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

library(Seurat)

###required load R packages

#devtools::install_github("mayunlong89/scMORE")

library(scMORE)
library(GenomicRanges)
library(IRanges)
library(Seurat)
library(Signac)




###---------PhastCons and exclude exons -----for ablation studies


#load data on 10x pbmc example data-----------------------------------
#load data on 10x pbmc example data-----------------------------------
#load data on 10x pbmc example data-----------------------------------
#load data on 10x pbmc example data-----------------------------------
#load data on 10x pbmc example data-----------------------------------
#load data on 10x pbmc example data-----------------------------------
#Load scMultiomic data
pbmc_10x <- readRDS("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/10X_PBMC/10X_PBMC.rds")

length(unique(pbmc_10x$cell_type))
single_cell_8900<-pbmc_10x
cell_type2<- as.character(single_cell_8900$cell_type)
cell_type2[which(cell_type2 %in% c("CD14+ Monocytes","FCGR3A+ Monocytes"))] <- "Monocytes"
single_cell_8900$cell_type2 <- cell_type2
Idents(single_cell_8900) <- single_cell_8900$cell_type2


grn_outputs8900_phastCons_withoutExons <- createRegulon(single_cell_8900, n_targets=5,
                                                        peak2gene_method = 'Signac',
                                                        infer_method = 'glm',
                                                        exclude_exon_regions = TRUE,
                                                        conserved_regions = phastConsElements20Mammals.UCSC.hg38)


saveRDS(grn_outputs8900_phastCons_withoutExons,file="~/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2025-05-20-scMORE_revision/05_parameter_assessment/grn_outputs8900_phastCons_withoutExons.rds")


##----Celltype-specificity for each genes
target_scores <- suppressWarnings(getSpecificity(single_cell_8900))



#Five autoimmune disease GWASs
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/08_five_immunediseases/")


##-----IBD------------------------------
##-----IBD------------------------------
##-----IBD------------------------------
##-----IBD------------------------------
##-----IBD------------------------------
gene_info_IBD <- read.table("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/08_five_immunediseases/IBD_processed_magma_results.genes.out",header = TRUE)
#snp_info
snp_info_IBD <- read.csv("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/08_five_immunediseases/IBD_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

gene_info<- gene_info_IBD
snp_info <- snp_info_IBD


message("Step 3: Retrieving gene-level association scores...")
geneRiskScores <- getGeneScore(gene_info)



#1.phastCons_withoutExons
regulon2disease_results_alpha1_IBD <- regulon2disease(
  grn_outputs = grn_outputs8900_phastCons_withoutExons,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  infer_method = 'glm',
  top_n = 5
)

write.csv(regulon2disease_results_alpha1_IBD, file="results_IBD_scMore_trait_results_phastCons_withoutExons_perm5_alpha=1.csv",quote=F, row.names = F)






##-----PBC------------------------------
##-----PBC------------------------------
##-----PBC------------------------------
##-----PBC------------------------------
gene_info_PBC <- read.table("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/08_five_immunediseases/PBC_count_GWAS_final.hg19_SNP_Gene_Analysis_results.genes.out",header = TRUE)
#snp_info
snp_info_PBC <- read.table("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/08_five_immunediseases/PBC_maf0.01.txt",header=T,stringsAsFactors = FALSE)

gene_info<- gene_info_PBC
snp_info <- snp_info_PBC


message("Step 3: Retrieving gene-level association scores...")
geneRiskScores <- getGeneScore(gene_info)


#1.phastCons_withoutExons
regulon2disease_results_alpha1_PBC <- regulon2disease(
  grn_outputs = grn_outputs8900_phastCons_withoutExons,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  infer_method = 'glm',
  top_n = 5
)

write.csv(regulon2disease_results_alpha1_PBC, file="results_PBC_scMore_trait_results_phastCons_withoutExons_perm1000_alpha=1.csv",quote=F, row.names = F)




##-----RA------------------------------
##-----RA------------------------------
##-----RA------------------------------
##-----RA------------------------------
##-----RA------------------------------
gene_info_RA <- read.table("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/08_five_immunediseases/RA_count_GWAS_final.hg19_SNP_Gene_Analysis_results.genes.out",header = TRUE)
#snp_info
snp_info_RA <- read.table("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/08_five_immunediseases/RA_maf0.01.txt",header=T,stringsAsFactors = FALSE)

gene_info<- gene_info_RA
snp_info <- snp_info_RA


message("Step 3: Retrieving gene-level association scores...")
geneRiskScores <- getGeneScore(gene_info)


#1.phastCons_withoutExons
regulon2disease_results_alpha1_RA <- regulon2disease(
  grn_outputs = grn_outputs8900_phastCons_withoutExons,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  infer_method = 'glm',
  top_n = 5
)

write.csv(regulon2disease_results_alpha1_RA, file="results_RA_scMore_trait_results_phastCons_withoutExons_perm1000_alpha=1.csv",quote=F, row.names = F)





##-----SLE------------------------------
##-----SLE------------------------------
##-----SLE------------------------------
##-----SLE------------------------------
##-----SLE------------------------------
gene_info_SLE <- read.table("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/08_five_immunediseases/SLE_count_GWAS_final.hg19_SNP_Gene_Analysis_results.genes.out",header = TRUE)
#snp_info
snp_info_SLE <- read.table("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/08_five_immunediseases/SLE_maf0.01.txt",header=T,stringsAsFactors = FALSE)

gene_info<- gene_info_SLE
snp_info <- snp_info_SLE


message("Step 3: Retrieving gene-level association scores...")
geneRiskScores <- getGeneScore(gene_info)



#1.phastCons_withoutExons
regulon2disease_results_alpha1_SLE <- regulon2disease(
  grn_outputs = grn_outputs8900_phastCons_withoutExons,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  infer_method = 'glm',
  top_n = 5
)

write.csv(regulon2disease_results_alpha1_SLE, file="results_SLE_scMore_trait_results_phastCons_withoutExons_perm1000_alpha=1.csv",quote=F, row.names = F)





##-----UC------------------------------
##-----UC------------------------------
##-----UC------------------------------
##-----UC------------------------------
##-----UC------------------------------
gene_info_UC <- read.table("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/08_five_immunediseases/UC_processed_magma_results.genes.out",header = TRUE)
#snp_info
snp_info_UC <- read.table("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/08_five_immunediseases/UC_maf0.01.txt",header=T,stringsAsFactors = FALSE)

gene_info<- gene_info_UC
snp_info <- snp_info_UC


message("Step 3: Retrieving gene-level association scores...")
geneRiskScores <- getGeneScore(gene_info)



#1.phastCons_withoutExons
regulon2disease_results_alpha1_UC <- regulon2disease(
  grn_outputs = grn_outputs8900_phastCons_withoutExons,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  infer_method = 'glm',
  top_n = 5
)

write.csv(regulon2disease_results_alpha1_UC, file="results_UC_scMore_trait_results_phastCons_withoutExons_perm1000_alpha=1.csv",quote=F, row.names = F)






##-----PD------------------------------
##-----PD------------------------------
##-----PD------------------------------
##-----PD------------------------------
##-----PD------------------------------
human_PD_brain_10x <- readRDS("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/06-PD/PD_anno.rds")

human_PD_brain_10x <- JoinLayers(human_PD_brain_10x)

length(unique(human_PD_brain_10x$cell_type))

single_cell_58949  <- human_PD_brain_10x
Idents(single_cell_58949) <- human_PD_brain_10x$cell_type


grn_outputsBrain58949_phastCons_withoutExons <- createRegulon(single_cell_58949, n_targets=5,
                                                              peak2gene_method = 'Signac',
                                                              infer_method = 'glm',
                                                              exclude_exon_regions = TRUE,
                                                              conserved_regions = phastConsElements20Mammals.UCSC.hg38)

gene_info_PD <- read.table("~/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/06-PD/PD_processed_magma_results.genes.out",header = TRUE)
#snp_info
snp_info_PD <- read.csv("~/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/06-PD/PD_maf0.01.txt",header=T,sep="\t",stringsAsFactors = FALSE)

gene_info<- gene_info_PD
snp_info <- snp_info_PD

message("Step 3: Retrieving gene-level association scores...")
geneRiskScores <- getGeneScore(gene_info)

##----Celltype-specificity for each genes
target_scores <- suppressWarnings(getSpecificity(single_cell_58949))


#1.phastCons_withoutExons
regulon2disease_results_alpha1_PD <- regulon2disease(
  grn_outputs = grn_outputsBrain58949_phastCons_withoutExons,
  target_scores = target_scores,
  geneRiskScores = geneRiskScores,
  snp_info = snp_info,
  perm_n = 1000,
  theta = 0.5,
  alpha = 1,
  buffer = 500,
  infer_method = 'glm',
  top_n = 5
)

write.csv(regulon2disease_results_alpha1_PD, file="results_PD_scMore_trait_results_phastCons_withoutExons_perm1000_alpha=1.csv",quote=F, row.names = F)



regulon2disease_results_alpha1_PD <- read.csv("results_PD_scMore_trait_results_phastCons_withoutExons_perm1000_alpha=1.csv")





###------PD-----TRS vs Penalty
###------PD-----TRS vs Penalty
###------PD-----TRS vs Penalty



library(ggtern)
library(dplyr)
library(readr)


##------scatter plot for alpha=0



regulon2disease_results_alpha1_IBD <-  read.csv("results_IBD_scMore_trait_results_phastCons_withoutExons_perm1000_alpha=1.csv")



# ──────────────────────────
# 1) 读入Data
# ──────────────────────────
df <- regulon2disease_results_alpha1_UC
df <- regulon2disease_results_alpha1_PBC
df <- regulon2disease_results_alpha1_SLE
df <- regulon2disease_results_alpha1_RA
df <- regulon2disease_results_alpha1_IBD
df <- regulon2disease_results_alpha1_PD



###------1---TRS with penalty VS TRS without penalty
df$TRS_raw2 <- df$CTS_raw + df$GRS_raw
df$TRS_penalty <- df$TRS_raw2 - df$Penalty

ggplot(df, aes(x = TRS_raw, y = TRS_raw2)) +
  geom_point(aes(color = Significance), alpha = 0.6, size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(
    x = "TRS with penalty",
    y = "TRS without penalty",
    title = "Impact of Penalty Term on TRS",
    color = "Significance"
  ) +
  scale_color_manual(values = c("Significant" = "red", "Nonsignificant" = "gray")) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "top"
  )


###------2---TRS with penalty VS Penalty
ggplot(df, aes(x = TRS_raw, y = Penalty)) +
  geom_point(aes(color = Significance), alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(
    x = "TRS with penalty",
    y = "Penalty (SD)",
    title = "Penalty Term vs TRS",
    color = "Significance"
  ) +
  scale_color_manual(values = c("Significant" = "red", "Nonsignificant" = "gray")) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "top"
  )



###------3---Boxplot penalty contribution
library(ggplot2)
library(tidyr)
library(dplyr)

df$TRS_sum <- df$CTS_raw + df$GRS_raw
df$CTS_contrib <- df$CTS_raw / df$TRS_sum
df$GRS_contrib <- df$GRS_raw / df$TRS_sum
df$Penalty_contrib <- df$Penalty / df$TRS_sum


# 将数据转换为长格式
df_long <- df %>%
  dplyr::select(CTS_contrib, GRS_contrib, Penalty_contrib) %>%
  pivot_longer(cols = everything(),
               names_to = "Component",
               values_to = "Contribution")

# 可选：优化标签显示
df_long$Component <- factor(df_long$Component,
                            levels = c("CTS_contrib", "GRS_contrib", "Penalty_contrib"),
                            labels = c("CTS", "GRS", "Penalty"))

# 画图
ggplot(df_long, aes(x = Component, y = Contribution, fill = Component)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.3, color = "black", size = 0.6) +
  labs(title = "Relative Contribution of CTS / GRS / Penalty to TRS",
       x = "Component",
       y = "Relative Contribution (scaled by TRS)") +
  scale_fill_manual(values = c("CTS" = "#66c2a5", "GRS" = "#fc8d62", "Penalty" = "#8da0cb")) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "none"
  )




