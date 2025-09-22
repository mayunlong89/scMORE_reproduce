

##------FOXP2---GO-term enrichments

library(ggplot2)

# Data
go_data <- data.frame(
  Term = c("cell fate commitment",
           "muscle tissue development",
           "cardiocyte differentiation",
           "ameboidal-type cell migration",
           "stem cell differentiation",
           "regulation of developmental growth",
           "cognition",
           "muscle cell proliferation",
           "forebrain development"),
  FDR = c(1.331623403, 1.343336529, 1.343336529, 1.343336529,
          1.369857002, 1.642503478, 1.642503478, 1.642503478, 8.57024772)
)

# 设置 GO term 顺序（按照 FDR 值从大到小排序）
go_data$Term <- factor(go_data$Term, levels = go_data$Term[order(go_data$FDR)])

# 画横向条形图
ggplot(go_data, aes(x = Term, y = FDR)) +
  geom_bar(stat = "identity", fill = "#B3CDE3", color = "black", width = 0.6) +
  coord_flip() +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold")) +
  ylab("-Log10(FDR)") +
  xlab("") +
  ggtitle("GO Terms Enrichment")
