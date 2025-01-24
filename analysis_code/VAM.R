########-----------------VAM
library(VAM)  #install.packages("VAM")
#BiocManager::install("qusage")
library(qusage)
#devtools::install_github("cellgeni/sceasy")
#install.packages('reticulate')
library(sceasy)
library(reticulate)



## Convert seurat format to sce format
#sceasy::convertFormat(single_cell,from = "seurat",to="sce",outFile = "organoids_brain_multi_sce.rds")
library(Seurat)
library(SeuratDisk)
library(SingleCellExperiment)

# 将 Seurat 对象转换为 SingleCellExperiment 对象
sce_multi <- as.SingleCellExperiment(single_cell)


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


gs <- qusage::read.gmt(gmt_file)

sce_multi <- importGeneSetsFromList(inSCE = sce_multi,geneSetList = gs,
                                    by = "rownames")

sce_multi <- runVAM(inSCE = sce_multi,
                    geneSetCollectionName = "GeneSetCollection",
                    useAssay = "logcounts")

#VAM values
sce_multi@int_colData$reducedDims$VAM_GeneSetCollection_Distance
results_VAM <- sce_multi@int_colData$reducedDims$VAM_GeneSetCollection_CDF
VAM_cell<-c()
VAM_cell <- as.data.frame(Idents(single_cell))

VAM_cell$vam_score <- results_VAM[,1]

colnames(VAM_cell) <- c("celltypes","scores")
results_vam2 <- aggregate(scores~celltypes,VAM_cell,mean)

results_vam2$regluons <- Module_regulon[1]

head(results_vam2)
head(results)
