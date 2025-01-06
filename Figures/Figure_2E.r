

##----------Heatmap plots------------------------------------Figure 2E
#@ Heatmap
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/02-10-blood_traits/scMORE_benchmark_10blood_trait_magma_vcf/")

library(pheatmap)


lym_heatmap <- read.table("lymp_count_percent_heatmap.txt",header = TRUE,sep = "\t")


data_h1 <- lym_heatmap[,-1]
rownames(data_h1) <- lym_heatmap[,1]


data_h1_score <- data_h1[,c(1,3,5,7)]

data_h1_p <- data_h1[,c(2,4,6,8)]

data_h1_p2 <- filter(data_h1_p,Lymp_count_CD8T_P < 0.05 | Lymp_percent_CD8T_P <0.05 | Lymp_count_monocyte_P < 0.05 | Lymp_percent_monocyte_P < 0.05 )
data_h1_score2 <- data_h1_score[which(rownames(data_h1_score) %in% rownames(data_h1_p2)),]

p_vals<-t(data_h1_p2)
heatmap_data <- t(data_h1_score2)

p_vals<- data_h1_p2
heatmap_data <- data_h1_score2

pearson_heatmap<-pheatmap(heatmap_data,cellwidth =15, cellheight =15,fontsize=8,
                          border_color="black",scale="none",
                          fontsize_row=8,
                          color=colorRampPalette(rev(c("#E77A77", "#ECBA84","white","lightblue",'#336699')))(102),
                          cluster_rows = F,
                          cluster_cols = F,
                          display_numbers = matrix(ifelse(p_vals <0.01 & heatmap_data > 1.96, "**",
                                                          ifelse(p_vals<0.05 & heatmap_data>1.96,"*","")), nrow(p_vals)))

