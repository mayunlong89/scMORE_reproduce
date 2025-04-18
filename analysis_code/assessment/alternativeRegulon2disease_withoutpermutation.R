#' @title Identify Cell Type-Specific Regulons Relevant to Disease
#'
#' @description
#' Function to identify cell type-specific regulons (TRS) that are associated with disease by integrating
#' cell type specificity scores(CTS) and genetic relevance scores (GRS).
#'
#' @param grn_outputs GRN (Gene Regulatory Network) outputs containing TF-gene relationships.
#' @param target_scores A matrix of genes and their TF specificity scores across cell types.
#'                      Columns include:
#'                      - `genes`: Gene names
#'                      - `scores`: TF specificity scores
#'                      - `celltypes`: Cell type annotations
#' @param snp_info SNP information used for the analysis.
#' @param geneRiskScores Gene-based genetic association results from MAGMA. Must include the following columns:
#'                       - `SYMBOL`: Gene symbols
#'                       - `logP`: Log-transformed p-values
#'                       - `ZSTAT`: Z-scores from MAGMA
#' @param perm_n Number of permutations for Monte Carlo simulation. Default: 1000.
#' @param theta Weighting factor for integrating TF and gene specificity scores. Range: 0.1 to 1. Default: 0.5.
#' @param alpha Flexibility parameter for penalization in the scoring model. Default: 1.
#' @param top_n Number of top targets for each TF used to calculate TF importance. Default: 5.
#' @param buffer Numeric value specifying the flanking region size for genomic peaks.
#'               Extends the peak range upstream and downstream by the specified value.
#'               For example, `buffer = 500` adds 500 bp on both sides of a peak. Default: 500 bp.
#' @param infer_method GRN inference method:
#'   - 'glm' (default): Generalized linear model
#'   - 'cv.glmnet': Regularized GLMs
#'   - 'brms': Bayesian regression models
#'   - 'bagging_ridge': Bagging ridge and Bayesian ridge (CellOracle)
#'   - 'xgb': XGBoost gradient-boosted Random Forest (SCENIC)
#' @param alternative A string specifying the method to use for scoring. Options are:
#'                "AUCell", "RSS", "VAM", "Seurat", "UCell". Default is "AUCell".
#' @param p1 Threshold for statistical significance of the cell type-specificity score (CTS)
#'           for each regulon. Default: 0.05.
#' @param p2 Threshold for statistical significance of the genetic relevance score (GRS)
#'           for each regulon. Default: 0.05.
#' @param p3 Threshold for statistical significance of the trait-associated regulon score (TRS).
#'           Default: 0.05.
#'
#' @return A data frame containing the following:
#'         - `specificity`: Cell type-specificity scores for each regulon (CTS).
#'         - `genetic_risk_score`: Genetic relevance scores for each regulon (GRS).
#'         - `regulon_score`: Trait-associated regulon scores (TRS).
#'         - `p_values`: Statistical significance values for specificity, risk, and regulon scores.
#'
#' @export
#'
alternativeRegulon2disease2 <- function(grn_outputs,
                                       geneRiskScores,
                                       snp_info,
                                       perm_n = 1000,
                                       theta = 0.5,
                                       alpha = 1,
                                       top_n = 5,
                                       buffer = 500,
                                       infer_method = 'glm',
                                       alternative = "AUCell",
                                       p1 = 0.05,
                                       p2 = 0.05,
                                       p3 = 0.05) {

  # Step 4.1: Map SNPs to TF-peaks and target genes
  message("Step 4.1: Mapping SNPs to TF-peaks and target genes...")
  peak2gene_strength <- peak2gene(grn_outputs,infer_method = infer_method)  # Map peaks to genes
  snp2peak_map <- snp2peak(snp_info, peak2gene_strength, buffer = buffer)  # SNP-to-peak mapping

  # Add gene risk scores to SNP-to-peak mapping
  snp2peak_map$geneScores <- geneRiskScores$logP[match(snp2peak_map$Target, geneRiskScores$SYMBOL)]
  snp2peak_map <- getPeakScore(snp2peak_map)  # Calculate peak importance scores

  # Extract regulons containing SNP, peak, TF, target gene, and importance score (GRS score of each node)
  regulons <- snp2peak_map[, c("snp_id", "peak_ids", "TF", "Target", "Importance_weighted")]

  # Step 4.2: Define Regulon list and cell types
  message("Step 4.2: Defining Regulon list and cell types...")
  tf_list <- grn_outputs$tf_names  # Transcription factors

  # Initialize results data frame
  all_regulon_results_df <- data.frame(
    RegulonID = character(),
    RegulonName = character(),
    SpecificityScore = numeric(),
    SpecificityScore_p = numeric(),
    GeneRiskScore = numeric(),
    ImportanceWeightScore_p = numeric(),
    RegulonScore = numeric(),
    RegulonScore_p = numeric(),
    Celltype = character(),
    stringsAsFactors = FALSE
  )

  # Step 4.3: Perform scMORE analysis for each cell type and regulon
  message("Step 4.3: Calculate Trait-Asssociated Regulon Score and Signifiance...")


  # Calculate the cell type-specificity score (CTS) of each regulon
  alternative_regulon_score <- alternativeSpecificityScore(single_cell = single_cell,
                                                           tf_list = tf_list,
                                                           regulons = regulons,
                                                           alternative = alternative)

  #alternative_regulon_score$scores[alternative_regulon_score$scores < 0] <- 0
  #alternative_regulon_score$scores <- log10(alternative_regulon_score$scores + 1e-6)
  alternative_regulon_score$scores <- max_min_scale(alternative_regulon_score$scores)

  # extract cell types
  all_celltype_names <- unique(alternative_regulon_score[, "celltypes"])  # Cell types

  # Open progress bar
  pb <- txtProgressBar(style = 3)
  start_time <- Sys.time()
  total_run <- length(all_celltype_names) * length(tf_list)
  count <- 0

  # Runing regulon2disease_alternative()
  for (i in seq_along(all_celltype_names)) {
    for (j in seq_along(tf_list)) {

      # Extract regulon genes for the current TF
      Module <- regulons[regulons$TF == tf_list[j], ]
      Module_regulon <- c(unique(Module$TF), unique(Module$Target))

      # Calculate importance scores for the regulon
      each_module_score <- getRiskScore(Module, top_n = top_n)
      each_module_score$anno <- ifelse(each_module_score$Target == Module_regulon[1], "TF", "Gene")
      #Here, header of each_module_score: genes, scores, celltypes, anno, Importance_weighted
      #scores: cell type-specificity (CTS) score of each gene (node)
      #Importance_weighted: genetic relevance scores (GRS) of each node

      # Cell type specificity score (CTS) using alternative methods
      alternative_regulon_score_sub <- alternative_regulon_score$scores[which(alternative_regulon_score$celltypes==all_celltype_names[i] & alternative_regulon_score$regulons==tf_list[j])]

      # Compute regulon scores (TRS)
      TRS <- alternativeRegulonScore(each_module_score,
                                     alternative_regulon_score_sub,
                                     theta = theta, alpha = alpha)


      # Compute z-scores for each regulon
      z_specificity <- TRS$SpecificityScore
      #z_importance <- TRS$GeneRiskScore
      #z_regulon <- TRS$RegulonScore

      # Append results to the data frame
      all_regulon_results_df <- rbind(
        all_regulon_results_df,
        data.frame(
          RegulonID = paste0("Regulon_", j),
          RegulonName = Module_regulon[1],
          RegulonScore = z_specificity,
          #SpecificityScore_p = p_specificity,
          #GeneRiskScore = z_importance,
          #ImportanceWeightScore_p = p_importance,
          #RegulonScore = z_regulon,
          #RegulonScore_p = p_regulon,
          Celltype = all_celltype_names[i]
        )
      )

      # Update progress bar
      count <- count + 1
      setTxtProgressBar(pb, count / total_run)
    }
  }



  # Close progress bar and record running time
  close(pb)
  end_time <- Sys.time()
  cat("Running time:", round(difftime(end_time, start_time, units = "secs"), 2), "seconds\n")

  return(all_regulon_results_df)
}
