
##-2
consine_average10traits_heatmap <- read.table("01_10blood_cell_traits_sig_heatmap.txt",header = TRUE,sep = "\t")

library(pheatmap)

data_h1 <- consine_average10traits_heatmap[,c(-1)]
rownames(data_h1) <- consine_average10traits_heatmap[,1]
pearson_heatmap<-pheatmap(data_h1,cellwidth =15, cellheight =15,fontsize=8,
                          border_color="black",scale="row",
                          fontsize_row=8,
                          color=colorRampPalette(rev(c("red","#E77A77", "#ECBA84","white","lightblue",'#336699')))(102),
                          cluster_rows = T,
                          cluster_cols = T)
