#2025-02-05

# Permutation analysis for 1000 genes: Assessing correlation with PubMed search results using RISmed across 100 iterations of random gene weight selection.

# Set CRAN mirror for installing packages
options(repos = structure(c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")))

# Set working directory
setwd("~/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/06-PD/")


# Load the required data---AS-specific TRSs for PD
gene <- read.table("01_PD_AS_regulons_RISmed_list.txt", header = FALSE)
colnames(gene) <- c("Gene", "Weight")  # Rename columns for clarity



##-----running RISmed method

# Load the RISmed package for PubMed querying
library(RISmed)

# Initialize variables
PMID <- NULL  # 用于存储 PubMed 查询结果

m<-length(gene$Weight)

# 选取前76个基因
traits <- gene[order(gene$Weight, decreasing = TRUE)[1:m], ]

# 循环处理每个基因
for (i in 1:nrow(traits)) {
  gene_name <- traits[i, "Gene"]
  cat("Processing gene ", i, ": ", gene_name, "\n", sep = "")
  
  # 在 tryCatch 外定义 search_term，构造检索字符串
  # 如果需要使用更多与 Parkinson disease 相关的关键词，可扩展 OR 子句
  search_term <- paste0(
    "(", gene_name, ") AND (PD OR \"Parkinson disease\" OR \"Parkinson's disease\" OR ",
    "\"Idiopathic Parkinson's disease\" OR Parkinsonism OR \"Paralysis agitans\" OR ",
    "\"Shaking palsy\")"
  )
  
  # 使用 tryCatch 进行错误处理
  tryCatch({
    res <- EUtilsSummary(search_term, type = 'esearch', db = 'pubmed')
    PMID <- rbind(PMID, c(gene_name, search_term, QueryCount(res)))
  }, error = function(e) {
    # 第一次出错后暂停5秒重试
    Sys.sleep(5)
    tryCatch({
      res <- EUtilsSummary(search_term, type = 'esearch', db = 'pubmed')
      PMID <- rbind(PMID, c(gene_name, search_term, QueryCount(res)))
    }, error = function(e) {
      # 第二次出错后暂停20秒重试
      Sys.sleep(20)
      res <- EUtilsSummary(search_term, type = 'esearch', db = 'pubmed')
      PMID <- rbind(PMID, c(gene_name, search_term, QueryCount(res)))
    })
  })
}

# 转换 PMID 结果为数据框以便处理
PMID <- as.data.frame(PMID, stringsAsFactors = FALSE)
colnames(PMID) <- c("Gene", "SearchTerm", "ResultCount")


# Convert PMID results to a data frame for easier handling
PMID <- as.data.frame(PMID, stringsAsFactors = FALSE)
colnames(PMID) <- c("Gene", "SearchTerm", "ResultCount")

# Add log2-transformed PubMed result counts
PMID$logCount <- log2(as.numeric(PMID$ResultCount) + 1)

# Annotate the weights (TRS) from the traits table
PMID$TRS <- traits$Weight[match(PMID$Gene, traits$Gene)]

# Perform correlation test between TRS and log-transformed PubMed counts
cor_TRS_log2count <- cor.test(as.numeric(PMID$TRS), as.numeric(PMID$logCount))
print(cor_TRS_log2count)


write.csv(PMID, file="02_PD_AS_RISmed_search_results_v2.csv",quote=F, row.names = F)

# Scatter plot with ggplot2
library(ggplot2)

# Scatter plot with regression line
scatter_plot <- ggplot(PMID, aes(x = logCount, y = TRS)) +
  geom_point(color = "blue", size = 2, alpha = 0.7) +  # Scatter points
  geom_smooth(method = "lm", color = "red", se = TRUE) +  # Regression line with confidence interval
  labs(
    title = "Correlation Between TRS and PubMed log2 Count",
    x = "log2(PubMed Result Count + 1)",
    y = "Trait Relevance Score (TRS)"
  ) +
  theme_minimal()

# Print the scatter plot
print(scatter_plot)

