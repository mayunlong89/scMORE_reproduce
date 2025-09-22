
library(ggplot2)

# Data
df <- data.frame(
  Regulon = c("ZFPM2", "ZEB1", "MEF2C", "NFATC2"),
  Overlap = c(26, 20, 12, 11),
  Regulon_Size = c(349, 230, 138, 145),
  Pvalue = c(0.009639885, 0.004302228, 0.023577914, 0.066764335)
)

# add new column
df$logP <- -log10(df$Pvalue)
df$Ratio <- df$Overlap / df$Regulon_Size

# plot
ggplot(df, aes(x = Ratio, y = logP, size = Regulon_Size, color = Regulon)) +
  geom_point(alpha = 0.8) +
  geom_text(aes(label = Regulon), vjust = -1.2, size = 4.5) +
  #geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", size = 1) +
  scale_size_continuous(name = "Regulon Size", range = c(5, 12)) +
  labs(x = "Overlap Ratio", y = "-log10(P-value)",
       title = "OpenTarget PD Gene Enrichment in Core Regulons") +
  theme_classic() +
  theme(legend.position = "right")



