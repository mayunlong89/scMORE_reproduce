
#2025-01-03

##-----top significant four regulons MC permutation plots of empirical distributions

#for lymphocyte count (Figure S2E-F)

library(scMORE)
library(GenomicRanges)
library(IRanges)
library(Seurat)
library(Signac)

  #---load pbmc_10x_real_data_downsampled_2000
  pbmc_10x_real_data_downsampled_2000 <- readRDS("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/10X_PBMC/10X_PBMC_downsample_2000cells.rds")
  DimPlot(pbmc_10x_real_data_downsampled_2000, reduction = "umap.new",cols = c("#67ADB7","#E77A77"))
  single_cell <-pbmc_10x_real_data_downsampled_2000

  ##--scMORE for lymphocyte count-----------------------------------------------------------------------------------------
  setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
  gene_info_lymp_count <- read.table("lymp_count_processed_magma_results.genes.out",header = TRUE)
  snp_info_lymp_count <- read.csv("lymphocyte_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

  gene_info <- gene_info_lymp_count
  snp_info <-snp_info_lymp_count

  ##--scMORE for lymphocyte percent-----------------------------------------------------------------------------------------
  setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")
  gene_info_lymp_percent <- read.table("lymp_percent_processed_magma_results.genes.out",header = TRUE)
  snp_info_lymp_percent <- read.csv("lymphocyte_percent_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

  gene_info <- gene_info_lymp_percent
  snp_info <-snp_info_lymp_percent


  #Default use 'Signac' method for analysis
  grn_outputs <- createRegulon(single_cell,
                               n_targets = 5,
                               peak2gene_method="Signac",
                               infer_method = 'glm')


  #'--Parameters setting
  perm_n = 1000
  theta = 0.5
  alpha = 1
  top_n = 5
  buffer = 500
  p1 = 0.05
  p2 = 0.05
  p3 = 0.05

  # Step 2: Calculating specificity scores for target genes
  message("Step 2: Calculating specificity scores for target genes...")
  target_scores <- suppressWarnings(getSpecificity(single_cell,method = "cosine"))

  # Step 3: Get gene-level association scores
  message("Step 3: Retrieving gene-level association scores...")
  geneRiskScores <- getGeneScore(gene_info)
  # Step 4.1: Map SNPs to TF-peaks and target genes
  peak2gene_strength <- peak2gene(grn_outputs,
                                  infer_method = 'glm')  # Map peaks to genes
  #snp2peak_map <- snp2peak(snp_info_lymp, peak2gene_strength, buffer = buffer)  # SNP-to-peak mapping
  #snp2peak_map <- snp2peak(snp_info_mono, peak2gene_strength, buffer = buffer)  # SNP-to-peak mapping
  snp2peak_map <- snp2peak(snp_info, peak2gene_strength, buffer = buffer)  # SNP-to-peak mapping


  # Add gene risk scores to SNP-to-peak mapping
  snp2peak_map$geneScores <- geneRiskScores$logP[match(snp2peak_map$Target, geneRiskScores$SYMBOL)]
  snp2peak_map <- getPeakScore(snp2peak_map)  # Calculate peak importance scores

  # Extract regulons containing SNP, peak, TF, target gene, and importance score
  regulons <- snp2peak_map[, c("snp_id", "peak_ids", "TF", "Target", "Importance_weighted")]

  all_celltype_names <- unique(target_scores[, "celltypes"])  # Cell types

  # Step 4.3: Perform scMORE analysis for each cell type and regulon
  message("Step 4.3: Calculate Trait-Asssociated Regulon Score and Signifiance...")

      # Extract regulon genes for the current TF
      Module <- regulons[regulons$TF == "LEF1", ]
      Module_regulon <- c(unique(Module$TF), unique(Module$Target))

      # Calculate importance scores for the regulon
      eachModule_Importance_score <- getRiskScore(Module, top_n = top_n)

      # Subset target scores for the current cell type
      target_scores_sub <- target_scores[target_scores[, "celltypes"] == all_celltype_names[2], ]

      # Filter specificity scores for regulon genes
      each_module_score <- target_scores_sub[match(target_scores_sub[, "genes"], Module_regulon, nomatch = 0) > 0, ]
      each_module_score$anno <- ifelse(each_module_score$genes == Module_regulon[1], "TF", "Gene")

      # Add importance scores to regulon
      each_module_score$Importance_weighted <- eachModule_Importance_score$Importance_weighted[match(
        each_module_score$genes,
        eachModule_Importance_score$Target
      )]

      # Compute regulon scores (TRS)
      TRS <- getRegulonScore(each_module_score, theta = theta, alpha = alpha)

      # Run Monte Carlo (MC) permutation analysis
      len_of_regulon <- length(Module_regulon)
      target_scores_background <- target_scores_sub
      real_specificity <- target_scores_sub$scores[match(target_scores_sub[, "genes"], regulons$Target, nomatch = 0) > 0]
      real_importance <- regulons$Importance_weighted

      # MC analysis
      perm_results <- replicate(
        perm_n,
        getRandomScore(
          tf_list,
          target_scores_background,
          real_specificity,
          real_importance,
          len_of_regulon,
          theta = theta,
          alpha = alpha,
          top_n = top_n
        )
      )

      # Extract and clean permutation scores
      perm_specificity <- as.numeric(perm_results["SpecificityScore", ])
      perm_importance <- as.numeric(perm_results["ImportanceWeightScore", ])
      perm_regulon <- as.numeric(perm_results["RegulonScore", ])

      # Remove NA values
      perm_specificity <- perm_specificity[!is.na(perm_specificity)]
      perm_importance <- perm_importance[!is.na(perm_importance)]
      perm_regulon <- perm_regulon[!is.na(perm_regulon)]

      # Calculate p-values and z-scores
      p_specificity <- (1 + sum(perm_specificity >= TRS$SpecificityScore)) / (1 + length(perm_specificity))
      p_importance <- (1 + sum(perm_importance >= TRS$GeneRiskScore)) / (1 + length(perm_importance))
      p_regulon <- (1 + sum(perm_regulon >= TRS$RegulonScore)) / (1 + length(perm_regulon))

      # Compute z-scores for each regulon
      z_specificity <- if (sd(perm_specificity) == 0) NA else (TRS$SpecificityScore - mean(perm_specificity)) / sd(perm_specificity)
      z_importance <- if (sd(perm_importance) == 0) NA else (TRS$GeneRiskScore - mean(perm_importance)) / sd(perm_importance)
      z_regulon <- if (sd(perm_regulon) == 0) NA else (TRS$RegulonScore - mean(perm_regulon)) / sd(perm_regulon)



####--------Empirical distribution
      ## 确保 perm_regulon 是数值型且无缺失值
      perm_regulon <- as.numeric(perm_regulon)
      perm_regulon <- perm_regulon[!is.na(perm_regulon)]

      # 创建数据框
      data <- data.frame(value = perm_regulon)

      # 真实值
      TRS1 <- TRS$RegulonScore

      # 绘制经验分布
      library(ggplot2)

      ggplot(data, aes(x = value)) +
        # 绘制直方图
        geom_histogram(
          binwidth = 0.1,  # 设置柱宽
          fill = "#FFB233",  # 填充颜色
          color = "black",   # 边框颜色
          alpha = 0.7        # 透明度
        ) +
        # 添加竖线标注 TRS
        geom_vline(
          xintercept = TRS1,
          color = "red",
          linetype = "dashed",
          size = 1
        ) +
        # 设置标题和坐标轴标签
        labs(
          title = "Empirical Distribution of perm_regulon",
          x = "Empirical TRS values",
          y = "Frequency"
        ) +
        # 设置主题
        theme_classic() +
        theme(
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  # 标题样式
          axis.title = element_text(size = 14),  # 坐标轴标签样式
          axis.text = element_text(size = 12)    # 坐标轴刻度文字样式
        )




