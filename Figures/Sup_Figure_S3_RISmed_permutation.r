
#2025-01-05

# Permutation analysis for 1000 genes: Assessing correlation with PubMed search results using RISmed across 100 iterations of random gene weight selection.

# Set CRAN mirror for installing packages
options(repos = structure(c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")))

# Set working directory
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")


# Load the required data---lymphocyte count trait CD8+T-cells TRSs
gene <- read.table("02_lymp_count_2000cells_regulons_list_CD8T.txt", header = FALSE)
colnames(gene) <- c("Gene", "Weight")  # Rename columns for clarity


# Load the required data---lymphocyte count trait monocytes TRSs
#gene <- read.table("02_lymp_count_2000cells_regulons_list_monocytes.txt", header = FALSE)
#colnames(gene) <- c("Gene", "Weight")  # Rename columns for clarity


# Load the required data---lymphocyte percent trait CD8+T-cells TRSs
gene <- read.table("02_lymp_percent_2000cells_regulons_list_CD8T.txt", header = FALSE)
colnames(gene) <- c("Gene", "Weight")  # Rename columns for clarity




##-----running RISmed method

# Load the RISmed package for PubMed querying
library(RISmed)

# Initialize variables
PMID <- NULL  # To store PubMed query results

# Select the top 76 genes based on actual weights
traits <- gene[order(gene$Weight, decreasing = TRUE)[1:76], ]

# Loop through each gene in the selected traits
for (i in 1:nrow(traits)) {
  print(paste0("Processing gene ", i, ": ", traits[i, "Gene"]))

  # Attempt PubMed query with error handling
  tryCatch({
    # Use PubMed to search for the gene and "lymphocyte"
    res <- EUtilsSummary(paste0("(", traits[i, "Gene"], ") AND (lymphocyte)"),
                         type = 'esearch', db = 'pubmed')
    # Record the gene name, search term, and query result count
    PMID <- rbind(PMID, c(traits[i, "Gene"], "lymphocyte", QueryCount(res)))
  }, error = function(e) {
    # If an error occurs, retry with a delay
    Sys.sleep(5)
    tryCatch({
      res <- EUtilsSummary(paste0("(", traits[i, "Gene"], ") AND (lymphocyte)"),
                           type = 'esearch', db = 'pubmed')
      PMID <- rbind(PMID, c(traits[i, "Gene"], "lymphocyte", QueryCount(res)))
    }, error = function(e) {
      Sys.sleep(20)
      res <- EUtilsSummary(paste0("(", traits[i, "Gene"], ") AND (lymphocyte)"),
                           type = 'esearch', db = 'pubmed')
      PMID <- rbind(PMID, c(traits[i, "Gene"], "lymphocyte", QueryCount(res)))
    })
  })
}

# Convert PMID results to a data frame for easier handling
PMID <- as.data.frame(PMID, stringsAsFactors = FALSE)
colnames(PMID) <- c("Gene", "SearchTerm", "ResultCount")

# Add log2-transformed PubMed result counts
PMID$logCount <- log2(as.numeric(PMID$ResultCount) + 1)

# Annotate the weights (TRS) from the traits table
PMID$TRS <- traits$Weight[match(PMID$Gene, traits$Gene)]

# Perform correlation test between TRS and log-transformed PubMed counts
cor_TRS_log2count <- cor.test(PMID$TRS, PMID$logCount)
print(cor_TRS_log2count)


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



# Save the results
write.table(PMID, "02_lymp_percent_real_results_CD8T_cells.txt", quote = FALSE, row.names = FALSE, sep = "\t")

print("PubMed query completed and results saved.")

# Calculate correlation: gene weights vs log-transformed PubMed result counts



###-----------Permutation---analysis of random traits

# Load gene data (assumes no header)
# Gene file format: two columns (Gene Name, Weight/Score)
gene <- read.table("02_lymp_count_2000cells_regulons_list.txt", header = FALSE)


# Load the RISmed package for PubMed querying
library(RISmed)


res <- EUtilsSummary(paste0("(", traits[i, 1], ") AND (lymphocyte)"),
                     type = 'esearch', db = 'pubmed')
# Record the gene name, search term, and result count
PMID <- rbind(PMID, c(traits[i, 1], "lymphocyte", QueryCount(res)))



# Set random seed for reproducibility
set.seed(123456)

# Initialize variables
temp <- gene  # Temporary storage for shuffled gene weights
cor <- NULL   # To store correlation results

# Start permutation and PubMed search
for (j in 1:100) {  # Repeat the process 10 times
  # Shuffle gene weights
  temp[, 2] <- gene[sample(1:nrow(gene), nrow(gene)), 2]

  # Select the top 76 genes based on the shuffled weights
  traits <- temp[order(temp[, 2], decreasing = TRUE)[1:76], ]

  # Initialize storage for PubMed search results
  PMID <- NULL

  # Iterate over each selected gene
  for (i in 1:nrow(traits)) {
    print(paste0("Iteration j = ", j))  # Print current iteration
    print(paste0("Processing gene i = ", i))  # Print current gene

    # Attempt PubMed search with error handling
    tryCatch({
      # Search PubMed for articles related to the gene and "lymphocyte"
      res <- EUtilsSummary(paste0("(", traits[i, 1], ") AND (lymphocyte)"),
                           type = 'esearch', db = 'pubmed')
      # Record the gene name, search term, and result count
      PMID <- rbind(PMID, c(traits[i, 1], "lymphocyte", QueryCount(res)))
    }, error = function(e) {
      # Retry with a delay if an error occurs
      tryCatch({
        Sys.sleep(5)  # Wait 5 seconds before retrying
        res <- EUtilsSummary(paste0("(", traits[i, 1], ") AND (lymphocyte)"),
                             type = 'esearch', db = 'pubmed')
        PMID <- rbind(PMID, c(traits[i, 1], "lymphocyte", QueryCount(res)))
      }, error = function(e) {
        # Retry again with a longer delay if the second attempt fails
        Sys.sleep(20)  # Wait 20 seconds before final retry
        res <- EUtilsSummary(paste0("(", traits[i, 1], ") AND (lymphocyte)"),
                             type = 'esearch', db = 'pubmed')
        PMID <- rbind(PMID, c(traits[i, 1], "lymphocyte", QueryCount(res)))
      })
    })
  }

  # Calculate correlation: gene weights vs log-transformed PubMed result counts
  matched_indices <- na.omit(match(PMID[, 1], traits[, 1]))  # Match genes with PMID results
  cor <- c(cor, cor(traits[matched_indices, 2], log2(as.numeric(PMID[, 3]) + 1)))
}

# Save correlation results to a text file
write.table(cor, "02_lymp_count_cor_100_randoms.txt", quote = FALSE, row.names = FALSE)

# Note: This script aims to maximize the efficiency of tools used in human genetic analysis.

###lymp_count
hist(cor, col="#D2691E",xlab="Counts of overlapped genes",xlim = c(-0.5,0.5),main=NULL)
abline(v=0.2374394,col="blue",lty="longdash")
P_value = length(cor[cor>0.2374394])/length(cor)
P_value

###lymp_percent
hist(cor, col="#D2691E",xlab="Counts of overlapped genes",xlim = c(-0.5,0.5),main=NULL)
abline(v=0.2382936,col="blue",lty="longdash")
P_value = length(cor[cor>0.2382936])/length(cor)
P_value



#########Permutation_plot for random selections

results_permut <- read.table("cor_all.txt",header = T)

hist(results_permut[,1], col="#D2691E",xlab="Counts of overlapped genes",xlim = c(-0.2,0.2),main=NULL)
abline(v=0.21,col="blue",lty="longdash")
P_value= length(results_permut$x[results_permut$x>0.21])/length(results_permut$x)



