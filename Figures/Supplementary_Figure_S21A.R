##---------------------heatmap-----Figure 5F---- Supplementary Figure S21A
library(pheatmap)

setwd("~/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/06-PD/")

# top 5 for Figure 5F
heatmap_8diseases <- read.csv("01_8disease_heatmap_top5.csv",header = TRUE)

#--all for Supplementary_Figure S21
heatmap_8diseases <- read.csv("01_8disease_heatmap.csv",header = TRUE)

head(heatmap_8diseases)


data_h1 <- heatmap_8diseases[,-1]
rownames(data_h1) <- heatmap_8diseases[,1]


data_h1_score <- data_h1[,c(1,3,5,7,9,11,13)]

data_h1_p <- data_h1[,c(2,4,6,8,10,12,14)]

data_h1_p2 <- filter(data_h1_p, AS_p < 0.05 | OPC_p <0.05 | ODC_p < 0.05 | T_p < 0.05 | N_p < 0.05 |  MG_p < 0.05 | EC_p < 0.05)
data_h1_score2 <- data_h1_score[which(rownames(data_h1_score) %in% rownames(data_h1_p2)),]

p_vals<-(data_h1_p2)
heatmap_data <- (data_h1_score2)


pearson_heatmap<-pheatmap(heatmap_data,cellwidth =8, cellheight = 8,fontsize=8,
                          border_color="black",
                          scale="row",
                          fontsize_row=8,
                          color=colorRampPalette(rev(c( "firebrick3","#E77A77", "#ECBA84","white","lightblue",'#336699')))(102),
                          cluster_rows = F,
                          cluster_cols = T,
                          display_numbers = matrix(ifelse(p_vals <0.05 & heatmap_data > 3, "**", ifelse(p_vals<0.05 & heatmap_data>1.5,"*","")), nrow(p_vals)))



pearson_heatmap<-pheatmap(heatmap_data,cellwidth =8, cellheight = 1,fontsize=8,
                          border_color="black",
                          scale="row",
                          fontsize_row=3,
                          color=colorRampPalette(rev(c( "firebrick3","#E77A77", "#ECBA84","white","lightblue",'#336699')))(102),
                          cluster_rows = F,
                          cluster_cols = T,
                          display_numbers = matrix(ifelse(p_vals <0.05 & heatmap_data > 3, "**", ifelse(p_vals<0.05 & heatmap_data>1.5,"*","")), nrow(p_vals)))



pearson_heatmap<-pheatmap(heatmap_data,cellwidth =8, cellheight = 1,fontsize=3,
                          border_color="black",
                          scale="row",
                          fontsize_row=1.5,
                          color=colorRampPalette(rev(c( "firebrick3","#E77A77", "#ECBA84","white","lightblue",'#336699')))(102),
                          cluster_rows = F,
                          cluster_cols = T)



p_vals<-t(data_h1_p2)
heatmap_data <- t(data_h1_score2)


pearson_heatmap<-pheatmap(heatmap_data,cellwidth =10, cellheight = 10,fontsize=8,
                          border_color="black",
                          scale="column",
                          fontsize_row=8,
                          color=colorRampPalette(rev(c( "firebrick3","#E77A77", "#ECBA84","white","lightblue",'#336699')))(102),
                          cluster_rows = T,
                          cluster_cols = T,
                          display_numbers = matrix(ifelse(p_vals <0.05 & heatmap_data > 3, "**", ifelse(p_vals<0.05 & heatmap_data>1.5,"*","")), nrow(p_vals)))


pearson_heatmap<-pheatmap(heatmap_data,cellwidth =10, cellheight = 10,fontsize=8,
                          border_color="black",
                          scale="none",
                          fontsize_row=8,
                          color=colorRampPalette(rev(c( "firebrick3","#E77A77", "#ECBA84","white","lightblue",'#336699')))(102),
                          cluster_rows = T,
                          cluster_cols = T,
                          display_numbers = matrix(ifelse(p_vals <0.001 & heatmap_data > 5, "**", ifelse(p_vals<0.05 & heatmap_data>3,"*","")), nrow(p_vals)))





pearson_heatmap<-pheatmap(heatmap_data,cellwidth =1.551, cellheight = 8,fontsize=3,
                          border_color="black",
                          scale="column",
                          fontsize_row=8,
                          color=colorRampPalette(rev(c( "firebrick3","#E77A77", "#ECBA84","white","lightblue",'#336699')))(102),
                          cluster_rows = T,
                          cluster_cols = T,
                          display_numbers = matrix(ifelse(p_vals <0.001 & heatmap_data > 5, "**", ifelse(p_vals<0.05 & heatmap_data>5,"*","")), nrow(p_vals)))


pearson_heatmap<-pheatmap(heatmap_data,cellwidth =1.551, cellheight = 8,fontsize=3,
                          border_color="black",
                          scale="column",
                          fontsize_row=8,
                          color=colorRampPalette(rev(c( "firebrick3","#E77A77", "#ECBA84","white","lightblue",'#336699')))(102),
                          cluster_rows = T,
                          cluster_cols = T,
                          display_numbers = matrix(ifelse(p_vals <0.05 & heatmap_data > 3, "**", ifelse(p_vals<0.05 & heatmap_data>1.5,"*","")), nrow(p_vals)))


pearson_heatmap<-pheatmap(heatmap_data,cellwidth =1.551, cellheight = 8,fontsize=3,
                          border_color="black",
                          scale="column",
                          fontsize_row=8,
                          color=colorRampPalette(rev(c( "firebrick3","#E77A77", "#ECBA84","white","lightblue",'#336699')))(102),
                          cluster_rows = T,
                          cluster_cols = F)
