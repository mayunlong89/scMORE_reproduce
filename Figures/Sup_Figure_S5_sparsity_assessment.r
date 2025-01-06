
#@ 2024-05-09

#@ Here, we leverage several single-cell datasets with cell numbers ranging from 2000 to 100000
#@ to evaluate the sparsity of each gene.


#load data on 10x pbmc example data
#Load scMultiomic data
pbmc_10x <- readRDS("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/10X_PBMC/10X_PBMC.rds")
single_cell<- pbmc_10x
Idents(single_cell) <- single_cell$cell_type

#single_cell$cell_type2[which(single_cell$cell_type2 %in% c("CD14+ Monocytes","FCGR3A+ Monocytes"))] <- "Monocytes"

cell_type2<- as.character(single_cell$cell_type)
cell_type2[which(cell_type2 %in% c("CD14+ Monocytes","FCGR3A+ Monocytes"))] <- "Monocytes"
single_cell$cell_type2 <- cell_type2

Idents(single_cell) <- single_cell$cell_type2 

#subset cells of two cell types
pbmc_10x_real_data <- subset(single_cell,idents=c("Monocytes","NK cells"))
#downsample cell list
cell.list <- WhichCells(pbmc_10x_real_data,idents = c("Monocytes","NK cells"),downsample = 1000)
pbmc_10x_real_data_downsampled_2000 <- pbmc_10x_real_data[,cell.list]
table(pbmc_10x_real_data_downsampled_2000$cell_type2)
pbmc_10x_real_data_downsampled_2000 #this real dataset used for assessing ctDRTF performance

DimPlot(pbmc_10x_real_data_downsampled_2000)


data_ASD
data_SCZ

pbmc_10x_Bcells <- readRDS("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/10X_B_Cell_Lymphoma/10X_B_Cell_Lymphoma.rds")

pbmc_10x_brain <- readRDS("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/10X_human_brain/10X_Human_Brain.rds")


#pbmc_10x_brain 
## Calculate the sparsity
data_s1 <-gene_cosr(pbmc_10x_brain)
sparsity_cosr <- length(data_s1$scores[which(data_s1$scores==0)])/length(data_s1$scores)
sparsity_cosr

data_s1_e <-gene_expr(pbmc_10x_brain)
sparsity_expr_e <- length(data_s1_e$scores[which(data_s1_e$scores==0)])/length(data_s1_e$scores)
sparsity_expr_e

library(ggplot2)
ggplot(data_s1,aes(x=scores))+
  geom_density(fill="lightblue")+
  geom_vline(aes(xintercept=0),color="red",linetype="dashed")+
  theme_classic()

data_s1_e  <- data_s1_e [-which(data_s1_e$scores>10),]

ggplot(data_s1_e,aes(x=scores))+
  geom_density(fill="lightblue")+
  geom_vline(aes(xintercept=0),color="red",linetype="dashed")+
  theme_classic()


## Calculate the sparsity
data_s1_8900 <-gene_cosr(single_cell)
sparsity_cosr_8900 <- length(data_s1_8900$scores[which(data_s1_8900$scores==0)])/length(data_s1_8900$scores)
sparsity_cosr_8900

data_s1_8900_e <-gene_expr(single_cell)
sparsity_expr_8900_e <- length(data_s1_8900_e$scores[which(data_s1_8900_e$scores==0)])/length(data_s1_8900_e$scores)
sparsity_expr_8900_e


library(ggplot2)
ggplot(data_s1_8900,aes(x=scores))+
  geom_density(fill="lightblue")+
  geom_vline(aes(xintercept=0),color="red",linetype="dashed")+
  theme_classic()

data_s1_8900_e <- data_s1_8900_e[-which(data_s1_8900_e$scores>10),]

ggplot(data_s1_8900_e,aes(x=scores))+
  geom_density(fill="lightblue")+
  geom_vline(aes(xintercept=0),color="red",linetype="dashed")+
  theme_classic()

max(data_s1_8900_e$scores)
length(data_s1_8900_e$scores)
length(data_s1_8900_e$scores[which(data_s1_8900_e$scores>100)])

###
## Calculate the sparsity
 
data_s1_Bcells <-gene_cosr(pbmc_10x_Bcells)
sparsity_cosr_Bcells <- length(data_s1_Bcells$scores[which(data_s1_Bcells$scores==0)])/length(data_s1_Bcells$scores)
sparsity_cosr_Bcells

data_s1_gene_Bcells <-gene_expr(pbmc_10x_Bcells)
sparsity_expr_Bcells <- length(data_s1_gene_Bcells$scores[which(data_s1_gene_Bcells$scores==0)])/length(data_s1_gene_Bcells$scores)
sparsity_expr_Bcells


library(ggplot2)
ggplot(data_s1_Bcells,aes(x=scores))+
  geom_density(fill="lightblue")+
  geom_vline(aes(xintercept=0),color="red",linetype="dashed")+
  theme_classic()

data_s1_gene_Bcells <- data_s1_gene_Bcells[-which(data_s1_gene_Bcells$scores>10),]

ggplot(data_s1_gene_Bcells,aes(x=scores))+
  geom_density(fill="lightblue")+
  geom_vline(aes(xintercept=0),color="red",linetype="dashed")+
  theme_classic()


## Calculate the sparsity
data_s1_1600 <-gene_cosr(pbmc_10x_real_data_downsampled_1600)
sparsity_cosr_1600 <- length(data_s1_1600$scores[which(data_s1_1600$scores==0)])/length(data_s1_1600$scores)
sparsity_cosr_1600

data_s1_1600_e <-gene_expr(pbmc_10x_real_data_downsampled_1600)
sparsity_1600_e <- length(data_s1_1600_e$scores[which(data_s1_1600_e$scores==0)])/length(data_s1_1600_e$scores)
sparsity_1600_e


library(ggplot2)
ggplot(data_s1_1600,aes(x=scores))+
  geom_density(fill="lightblue")+
  geom_vline(aes(xintercept=0),color="red",linetype="dashed")+
  theme_classic()

data_s1_1600_e <- data_s1_1600_e[-which(data_s1_1600_e$scores>10),]

ggplot(data_s1_1600_e,aes(x=scores))+
  geom_density(fill="lightblue")+
  geom_vline(aes(xintercept=0),color="red",linetype="dashed")+
  theme_classic()

 

temp_test <- read.csv("/Users/mayunlong/Desktop/outputs_final_test32_2000cells_lymphocyte_single_model.csv",header = TRUE)

final_data <- temp_test
 
 if(is.na(final_data$Tval)) {
   
   final_data$Tval_scale <- max_min_scale(final_data$Tval)
   
 }
na.omit()



###--------------------------------------------------------------------------------------
#@ 2) Calculating specificity scores
#single_cell is the input single-cell data.
#num_genes <- length(rownames(single_cell))

#@ mode 1 use COSR-calculated cell type-specificity score of each gene for the linear regression model
gene_cosr <- function(single_cell = single_cell){
  
  #Calculating the specificity scores of TFs and targeting genes across all cell types
  
  #All genes in single-cell data 
  num_genes <- length(rownames(single_cell))
  
  #COSG::cosg()
  celltype_markers <- COSG::cosg(single_cell,
                                 groups = "all",
                                 assay = "RNA",
                                 slot = "data",
                                 mu=1,
                                 n_genes_user = num_genes) #number of genes used for analysis
  
  #@ num_genes2 <- length(rownames(single_cell@assays$peaks))
  #celltype_markers2 <- COSG::cosg(single_cell,
  #                               groups = "all",
  #                               assay = "peaks",
  #                               slot = "data",
  #                               mu=1,
  #                               n_genes_user = num_genes2) #number of peaks used for analysis
  
  #extracting cell type names
  all_celltype_names <- colnames(celltype_markers$names)
  
  ##-----Extracting the specificity score of each gene or TF in a given regulon
  #Empty data.frame
  data_s1 <- data.frame(matrix(ncol = 3, nrow = 0))
  names(data_s1) <- c("genes","scores","celltypes")
  for (i in 1:length(all_celltype_names)){
    
    data_s <- data.frame(celltype_markers$names[i],celltype_markers$scores[i])
    data_s[,3] <- rep(all_celltype_names[i],length(data_s[,1]))
    names(data_s) <- c("genes","scores","celltypes")
    
    #collecting all the cell types
    data_s1 <- rbind(data_s1,data_s)
  }
  
  #If specificity score = -1, transforming -1 to 0
  #data_s1[,2][which(data_s1[,2] == -1)] <- 0
  
  return(data_s1)
  
}

## Calculate the sparsity
data_s1 <-gene_cosr(pbmc_10x_real_data_downsampled)
sparsity_cosr <- length(data_s1$scores[which(data_s1$scores==0)])/length(data_s1$scores)
sparsity_cosr

## Calculate the sparsity
data_s1_2000 <-gene_cosr(pbmc_10x_real_data_downsampled_2000)
sparsity_cosr <- length(data_s1_2000$scores[which(data_s1_2000$scores==0)])/length(data_s1_2000$scores)
sparsity_cosr


DimPlot(single_cell, reduction = "umap.rna")
DimPlot(single_cell, reduction = "umap.atac")
DimPlot(single_cell, reduction = "umap.biMod")



## Calculate the sparsity
data_s1_8900 <-gene_cosr(single_cell)
sparsity_cosr <- length(data_s1_8900$scores[which(data_s1_8900$scores==0)])/length(data_s1_8900$scores)
sparsity_cosr

data_s1_4443 <-gene_cosr(pbmc_10x_real_data)
sparsity_cosr <- length(data_s1_8900$scores[which(data_s1_8900$scores==0)])/length(data_s1_8900$scores)
sparsity_cosr


data_s1_1600 <-gene_cosr(pbmc_10x_real_data_downsampled_1600)
sparsity_cosr <- length(data_s1_1600$scores[which(data_s1_1600$scores==0)])/length(data_s1_1600$scores)
sparsity_cosr





###--------------------------------------------------------------------------------------
#@ mode 2 use expression value of each gene for the linear regression model
gene_expr <- function(single_cell = single_cell){
  
  #Calculating the averarge expression value of TFs and targeting genes across all cell types
  
  #extracting cell type names
  all_celltype_names <- as.vector(unique(Idents(single_cell)))
  
  #Calculating the averarge expression value or average peak value
  celltype.averages <- AverageExpression(single_cell)
  #1) celltype.averages$RNA 
  celltype_rna <- celltype.averages[[1]]
  
  #2) celltype.averages$peaks
  #celltype.averages[[2]] 
  
  ##-----Extracting the average expression value of each gene or TF in a given regulon
  #Empty data.frame
  data_s1 <- data.frame(matrix(ncol = 3, nrow = 0))
  names(data_s1) <- c("genes","scores","celltypes")
  for (i in 1:length(all_celltype_names)){
    
    data_s <- data.frame(rownames(celltype_rna),celltype_rna[,i])
    data_s[,3] <- rep(all_celltype_names[i],length(celltype_rna[,1]))
    names(data_s) <- c("genes","scores","celltypes")
    
    #collecting all the cell types
    data_s1 <- rbind(data_s1,data_s)
  }
  
  #If specificity score = -1, transforming -1 to 0
  #data_s1[,2][which(data_s1[,2] == -1)] <- 0
  
  return(data_s1)
  
}

## Calculate the sparsity
data_s1 <-gene_expr(pbmc_10x_real_data_downsampled)
sparsity_cosr <- length(data_s1$scores[which(data_s1$scores==0)])/length(data_s1$scores)
sparsity_cosr


## Calculate the sparsity
data_s1_2000 <-gene_expr(pbmc_10x_real_data_downsampled_2000)
sparsity_cosr <- length(data_s1_2000$scores[which(data_s1_2000$scores==0)])/length(data_s1_2000$scores)
sparsity_cosr


## Calculate the sparsity
data_s1_8900 <-gene_expr(single_cell)
sparsity_cosr <- length(data_s1_8900$scores[which(data_s1_8900$scores==0)])/length(data_s1_8900$scores)
sparsity_cosr

data_s1_4443 <-gene_expr(pbmc_10x_real_data)
sparsity_cosr <- length(data_s1_8900$scores[which(data_s1_8900$scores==0)])/length(data_s1_8900$scores)
sparsity_cosr




