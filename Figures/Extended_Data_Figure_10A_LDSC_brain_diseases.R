rm(list=ls())

library(readxl)
library(tidyr)
library(dplyr)
library(corrplot)
library(tibble)


# 读取数据
df <- read_excel("/share/pub/guiyy/Cerebral_comorbidities/LDSC/result/plot/ldsc_brain.xlsx")

corr_data <- df[,c("p1","p2","rg")]
corr_data$rg[corr_data$rg > 1] <- 1
corr_data <- unique(corr_data)
corr_wide <- spread(corr_data,p2,rg) 
rownames_vec <- corr_wide[[1]]  
corr_wide <- corr_wide[, -1]  
rownames(corr_wide) <- rownames_vec  


corr_p <- df[,c("p1","p2","p")]
corr_p <- unique(corr_p)
corr_wide_p <- spread(corr_p,p2,p)
rownames_p_vec <- corr_wide_p[[1]]  
corr_wide_p <- corr_wide_p[, -1]  
rownames(corr_wide_p) <- rownames_p_vec  

corr_wide_p <- as.data.frame(lapply(corr_wide_p, as.numeric), row.names = rownames(corr_wide_p))
corr_wide_p <- as.matrix(corr_wide_p)

col1=colorRampPalette(colors =c("#004E71","white","#B1182D"),space="Lab") 
pdf("/share/pub/guiyy/Cerebral_comorbidities/LDSC/result/plot/LDSC_brain.pdf",width=10,height=10)
corrplot(corr=as.matrix(corr_wide), p.mat = as.matrix(corr_wide_p),method = "circle",
         type = "upper", tl.pos="lt", tl.col="black", tl.srt = 45, insig="label_sig", 
         sig.level = c(.001, .01, .05),pch.cex = 0.8, col = col1(10))

corrplot(corr=as.matrix(corr_wide),type="lower",add=TRUE,method="number", col = col1(6),
         tl.pos = "n",cl.pos = "n",diag=FALSE ) 

dev.off()


