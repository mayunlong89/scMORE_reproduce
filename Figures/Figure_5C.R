


##---------------------heatmap---for Eight psychiatric diseases---Figure 5C
library(pheatmap)

setwd("~/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/10-8psychiatric_disorders//")

heatmap_8diseases <- read.csv("01_8disease_heatmap.csv",header = TRUE)
head(heatmap_8diseases)


data_h1 <- heatmap_8diseases[,-1]
rownames(data_h1) <- heatmap_8diseases[,1]


data_h1_score <- data_h1[,c(1,3,5,7,9,11,13,15)]

data_h1_p <- data_h1[,c(2,4,6,8,10,12,14,16)]

data_h1_p2 <- filter(data_h1_p, ADHD_P < 0.05 | AN_P <0.05 | ASD_P < 0.05 | BIP_P < 0.05 | MDD_P < 0.05 | OCD_P < 0.05 | SCZ_P < 0.05 | TS_P < 0.05)
data_h1_score2 <- data_h1_score[which(rownames(data_h1_score) %in% rownames(data_h1_p2)),]

p_vals<-(data_h1_p2)
heatmap_data <- (data_h1_score2)


pearson_heatmap<-pheatmap(heatmap_data,cellwidth =8, cellheight = 7.5,fontsize=8,
                          border_color="black",
                          scale="row",
                          fontsize_row=8,
                          color=colorRampPalette(rev(c( "firebrick3","#E77A77", "#ECBA84","white","lightblue",'#336699')))(102),
                          cluster_rows = T,
                          cluster_cols = T,
                          display_numbers = matrix(ifelse(p_vals <0.05 & heatmap_data > 3, "**", ifelse(p_vals<0.05 & heatmap_data>1.5,"*","")), nrow(p_vals)))


pearson_heatmap<-pheatmap(heatmap_data,cellwidth =8, cellheight =8,fontsize=8,
                          border_color="black",
                          scale="none",
                          fontsize_row=8,
                          color=colorRampPalette(rev(c( "firebrick3","#E77A77", "#ECBA84","white")))(102),
                          cluster_rows = F,
                          cluster_cols = T,
                          display_numbers = matrix(ifelse(p_vals <0.05 & heatmap_data > 3, "**", ifelse(p_vals<0.05 & heatmap_data>1.5,"*","")), nrow(p_vals)))

