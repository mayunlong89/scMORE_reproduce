
##----heatmap_for lymp_count and lymp_percent results in CD8+T cells across three GRN methods

##-1
consine_heatmap_glm_glmnet_xgb <- read.table("01_lymp_count_percent_glm_glmnet_xgb_heatmap.txt",header = TRUE,sep = "\t")

library(pheatmap)

data_h4 <- consine_heatmap_glm_glmnet_xgb[,c(-1)]
rownames(data_h4) <- consine_heatmap_glm_glmnet_xgb[,1]
pearson_heatmap<-pheatmap(data_h4,cellwidth =15, cellheight =15,fontsize=8,
                          border_color="black",scale="none",
                          fontsize_row=8,
                          color=colorRampPalette(rev(c("red","#E77A77", "#ECBA84","white",  "lightblue","lightblue")))(102),
                          cluster_rows = F,
                          cluster_cols = F)

# 将数据从宽格式转换为长格式
data_long <- melt(consine_heatmap_glm_glmnet_xgb, id.vars = "regulons",
                  variable.name = "traits", value.name = "TRS")

# 绘制 dot plot
ggplot(data_long, aes(x = regulons, y = traits, size = TRS, color = "red")) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(x = "Regulons", y = "Traits", size = "TRS", color = "Model") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


