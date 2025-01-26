
# Step 1: 计算每个细胞类型的基因平均表达值
celltype_averages <- as.data.frame(AverageExpression(single_cell)$RNA)

# Step 2: 添加基因在所有细胞中的平均表达值
celltype_averages$Overall_Mean <- rowMeans(celltype_averages)

# Step 3: 计算每个细胞类型的特异性
specificity <- sweep(celltype_averages[, -ncol(celltype_averages)], 1, celltype_averages$Overall_Mean, "/")

# Step 4: 将 specificity 加入结果数据框
celltype_averages <- cbind(celltype_averages, specificity)

# 重命名 specificity 列
colnames(celltype_averages)[(ncol(celltype_averages) - ncol(specificity) + 1):ncol(celltype_averages)] <- paste0(colnames(specificity), "_Specificity")

# 查看结果
head(celltype_averages)


# Convert specificity matrix to long format
specificity_long <- specificity %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "genes") %>% # Add gene names as a column
  tidyr::pivot_longer(
    cols = -genes, # All columns except "genes"
    names_to = "celltypes", # Column name for cell types
    values_to = "scores"    # Column name for specificity scores
  )

# 查看结果
head(specificity_long)



celltype_averages <- as.data.frame(AverageExpression(single_cell)$RNA)

# Step 2: Add the average expression of each gene across all cell types
celltype_averages$Overall_Mean <- rowMeans(celltype_averages, na.rm = TRUE) # Handle NA values
celltype_averages$Overall_Mean[celltype_averages$Overall_Mean == 0 | is.na(celltype_averages$Overall_Mean)] <- 1e-6

# Step 3: Calculate specificity for each cell type
specificity <- sweep(celltype_averages[, -ncol(celltype_averages)], 1, celltype_averages$Overall_Mean, "/")

# Step 4: Convert specificity matrix to long format
target_scores <- specificity %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "genes") %>% # Add gene names as a column
  tidyr::pivot_longer(
    cols = -genes, # All columns except "genes"
    names_to = "celltypes", # Column name for cell types
    values_to = "scores"    # Column name for specificity scores
  )


# Step 3: Replace specificity scores of -1 with 0
target_scores$scores[target_scores$scores == -1] <- 0

# Step 4: Log10-transformation and max-min scale the specificity scores
target_scores$scores <- log10(target_scores$scores + 1e-6)
target_scores$scores <- (target_scores$scores - min(target_scores$scores, na.rm = TRUE)) /
  (max(target_scores$scores, na.rm = TRUE) - min(target_scores$scores, na.rm = TRUE))

