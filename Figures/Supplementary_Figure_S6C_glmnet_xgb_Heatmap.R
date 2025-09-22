
##------Significant regulons identified by GLM show in GLMNET---based on cosine--Figure S8C

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
