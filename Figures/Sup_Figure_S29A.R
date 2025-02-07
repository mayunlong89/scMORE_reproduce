

##---------------------heatmap---for 8 aging-related diseases---Supplemental Figure S29A
library(pheatmap)

setwd("~/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/06-PD/")

heatmap_8diseases <- read.csv("01_8disease_heatmap_PD_sig.csv",header = TRUE)

#heatmap_8diseases <- read.csv("01_8disease_heatmap_PD_neurons.csv",header = TRUE)

head(heatmap_8diseases)


data_h1 <- heatmap_8diseases[,-1]
rownames(data_h1) <- heatmap_8diseases[,1]


data_h1_score <- data_h1[,c(1,3,5,7,9,11,13,15)]

data_h1_p <- data_h1[,c(2,4,6,8,10,12,14,16)]

data_h1_p2 <- filter(data_h1_p, PD_P < 0.05 | EAA_P <0.05 | FI_P < 0.05 | Lifespan_P < 0.05 | Healthspan_P < 0.05 |  longevity_90_P < 0.05)
data_h1_score2 <- data_h1_score[which(rownames(data_h1_score) %in% rownames(data_h1_p2)),]

p_vals<-(data_h1_p2)
heatmap_data <- (data_h1_score2)


pearson_heatmap<-pheatmap(heatmap_data,cellwidth =8, cellheight = 2.8,fontsize=8,
                          border_color="black",
                          scale="row",
                          fontsize_row=3,
                          color=colorRampPalette(rev(c( "firebrick3","#E77A77", "#ECBA84","white","lightblue",'#336699')))(102),
                          cluster_rows = T,
                          cluster_cols = T,
                          display_numbers = matrix(ifelse(p_vals <0.05 & heatmap_data > 3, "**", ifelse(p_vals<0.05 & heatmap_data>1.5,"*","")), nrow(p_vals)))


