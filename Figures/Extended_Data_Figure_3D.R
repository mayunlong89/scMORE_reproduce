

library(ggplot2)

# --OpenTarget-Lymphocyte count
df <- data.frame(
  Regulon = c("ETS1", "FOXO1", "FOXP1", "IKZF1","LEF1","ZEB1","ZNF721"),
  Overlap = c(6,5,5,10,4,4,1),
  Regulon_Size = c(23, 26, 13, 51,12 ,9,7),
  Pvalue = c(0.0002617551,0.0001821074,0.0001821074,8.407115e-05,0.001566536, 0.0004449588, 0.2779678)
)

# 添加新列
df$logP <- -log10(df$Pvalue)
df$Ratio <- df$Overlap / df$Regulon_Size

# 绘图
ggplot(df, aes(x = Ratio, y = logP, size = Regulon_Size, color = Regulon)) +
  geom_point(alpha = 0.8) +
  geom_text(aes(label = Regulon), vjust = -1.2, size = 4.5) +
  #geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", size = 1) +
  scale_size_continuous(name = "Regulon Size", range = c(5, 12)) +
  labs(x = "Overlap Ratio", y = "-log10(P-value)",
       title = "OpenTarget lymph_count Gene Enrichment in Core Regulons") +
  theme_classic() +
  theme(legend.position = "right")




library(ggplot2)

# --OpenTarget-Lymphocyte percent
df <- data.frame(
  Regulon = c("ETS1", "FOXO1", "FOXP1", "IKZF1","LEF1","ZEB1","ZNF721"),
  Overlap = c(5,5,6,8,3,4,0),
  Regulon_Size = c(23, 26, 13, 51,12 ,9,7),
  Pvalue = c(0.002125947,0.005694227,1.132249e-05,0.001983851,0.0151365,0.0004449588, 0.999)
)

# 添加新列
df$logP <- -log10(df$Pvalue)
df$Ratio <- df$Overlap / df$Regulon_Size

# 绘图
ggplot(df, aes(x = Ratio, y = logP, size = Regulon_Size, color = Regulon)) +
  geom_point(alpha = 0.8) +
  geom_text(aes(label = Regulon), vjust = -1.2, size = 4.5) +
  #geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", size = 1) +
  scale_size_continuous(name = "Regulon Size", range = c(5, 12)) +
  labs(x = "Overlap Ratio", y = "-log10(P-value)",
       title = "OpenTarget lymph_percent Gene Enrichment in Core Regulons") +
  theme_classic() +
  theme(legend.position = "right")


