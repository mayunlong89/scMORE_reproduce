
#2025-07-17


library(Seurat)

lymp_cout <- readRDS("scMore_pbmc_lymphocyte_count.rds")
write.csv(lymp_cout$scMore_trait_results,file="scMore_pbmc_lymphocyte_count.csv",quote = F, row.names = F)


lymp_percent <- readRDS("scMore_pbmc_lymphocyte_percent.rds")
write.csv(lymp_percent$scMore_trait_results,file="scMore_pbmc_lymp_percent.csv",quote = F, row.names = F)



lymp_cout2 <- readRDS("scMore_pbmc_sub_lymphocyte_count.rds")
write.csv(lymp_cout2$scMore_trait_results,file="scMore_pbmc_sub_lymphocyte_count.csv",quote = F, row.names = F)


lymp_percent2 <- readRDS("scMore_pbmc_sub_lymphocyte_percent.rds")
write.csv(lymp_percent2$scMore_trait_results,file="scMore_pbmc_sub_lymphocyte_percent.csv",quote = F, row.names = F)



####-visualization

##-1-----Figure B

library(pheatmap)
library(reshape2)

consine_heatmap_new_PBMC <- read.table("01_lymp_count_percent_New_PBMC_dataset.txt",header = TRUE,sep = "\t")

# 将数据从宽格式转换为长格式
data_long <- melt(consine_heatmap_new_PBMC, id.vars = "regulons",
                  variable.name = "traits", value.name = "TRS")

# 绘制 dot plot
ggplot(data_long, aes(x = regulons, y = traits, size = TRS, color = "red")) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(x = "Regulons", y = "Traits", size = "TRS", color = "Model") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



data_h4 <- consine_heatmap_new_PBMC[,c(-1)]
rownames(data_h4) <- consine_heatmap_new_PBMC[,1]
pearson_heatmap<-pheatmap(data_h4,cellwidth =15, cellheight =15,fontsize=8,
                          border_color="black",scale="none",
                          fontsize_row=8,
                          color=colorRampPalette(rev(c("red","#E77A77", "#ECBA84","white",  "lightblue","lightblue")))(102),
                          cluster_rows = F,
                          cluster_cols = F)

