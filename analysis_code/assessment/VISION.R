


##----2
library(VISION)
# 检查和安装依赖
#install.packages(c("irlba", "Matrix", "methods"))

Module_regulon

# 获取 RNA counts 并归一化
counts <- single_cell@assays$RNA@counts
n.umi <- colSums(counts)

# 使用矩阵操作进行归一化

library(Matrix)
# 确保输入是稀疏矩阵
scaled_counts <- as(Matrix(scaled_counts, sparse = TRUE), "dgCMatrix")
scaled_counts <- t(t(counts) / n.umi) * median(n.umi)

# Function to convert to GMT format and store in memory
output_gmt <- function(geneset) {
  # Convert each gene set to GMT format as a string
  gmt_data <- sapply(names(geneset), function(i) {
    paste(c(i, "temp", geneset[[i]]), collapse = "\t")
  })

  # Return as a character vector
  return(gmt_data)
}

#generating .gmt file
# Use the function to get GMT format in memory
geneset <- list(Module_regulon)
names(geneset) <- Module_regulon[1]
#print(names(geneset))
Module_regulon_gmt <- output_gmt(geneset)
gmt_file <- tempfile(fileext = ".gmt")
writeLines(Module_regulon_gmt, con = gmt_file)

# 使用 Vision 运行分析
vis <- Vision(scaled_counts,       # 归一化的表达矩阵
  signatures = gmt_file       # 临时生成的签名文件
)

# 执行 Vision 分析
vis <- analyze(vis)
#Viewing the results
#viewResults(vis)
getSignatureAutocorrelation(vis)



