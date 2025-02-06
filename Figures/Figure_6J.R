##-----------Figure 6J

setwd("~/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/06-PD/")

##read metadata

meta_anno <- read.csv("01_metadata_PD.csv", header=T)
head(meta_anno)

PD_seu$group <- meta_anno$group[match(PD_seu$orig.ident, meta_anno$ID)]
PD_seu$age <- meta_anno$age[match(PD_seu$orig.ident, meta_anno$ID)]
PD_seu$sex <- meta_anno$sex[match(PD_seu$orig.ident, meta_anno$ID)]


grn_outputs <- readRDS("grn_outputs_PD.rds")


regulons <- grn_outputs$grn
head(regulons)

Module <- regulons[regulons$TF == "PRDM16",]
Module <- regulons[regulons$TF == "TSHZ2",]
Module <- regulons[regulons$TF == "RFX4",]
Module <- regulons[regulons$TF == "ID4",]
Module <- regulons[regulons$TF == "AR",]
Module <- regulons[regulons$TF == "SOX6",]
Module <- regulons[regulons$TF == "NFATC4",]
Module <- regulons[regulons$TF == "CCDC88A",]

Module <- regulons[regulons$TF == "ASCL1",]


Module <- regulons[regulons$TF == "BCL11A",]
Module <- regulons[regulons$TF == "NR2F2",]
Module <- regulons[regulons$TF == "DMBX1",]

Module <- regulons[regulons$TF == "SPI1",]
Module <- regulons[regulons$TF == "RUNX1",]
Module <- regulons[regulons$TF == "IKZF1",]
Module <- regulons[regulons$TF == "PRDM1",]



Module <- regulons[regulons$TF == "HIVEP3",]




Module_regulon <- c(unique(Module$TF), unique(Module$Target))


#Calculating addModuleScore for each regulon
DefaultAssay(PD_seu) <- "RNA"
PD_seu <- AddModuleScore(PD_seu,
                              features = list(Module_regulon),
                              name="Module_regulon")


#add_cell <- as.data.frame(Idents(PD_seu))
#add_cell$add_score <- PD_seu$Module_regulon1

add_cell <- as.data.frame(PD_seu@meta.data)
head(add_cell)


##---boxplot
add_cell$group <- factor(add_cell$group, levels = c("Young", "Aged","PD"))

ggplot(add_cell, aes(x = group, y = Module_regulon1, fill = group)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "Distribution of Proportion by Cell Type and Sex",
    x = "",
    y = "Cell proportion",
    fill = "Group"
  ) +
  scale_fill_manual(values = c("PD" = "#e31a1c","Aged"="#85C89C","Young" = "#1f78b4")) +
  facet_wrap(~cell_type, scales = "free_y")


# 按细胞类型分组并执行 Wilcoxon 检验

add_cell2 <- add_cell[-which(add_cell$group=="Aged"),]
results <- add_cell2 %>%
  group_by(cell_type) %>%
  summarise(
    p_value = wilcox.test(Module_regulon1 ~ group)$p.value,
    statistic = wilcox.test(Module_regulon1 ~ group)$statistic
  )

# 查看结果
print(results)

add_cell2 <- add_cell[-which(add_cell$group=="Young"),]
results <- add_cell2 %>%
  group_by(cell_type) %>%
  summarise(
    p_value = wilcox.test(Module_regulon1 ~ group)$p.value,
    statistic = wilcox.test(Module_regulon1 ~ group)$statistic
  )

# 查看结果
print(results)


add_cell2 <- add_cell[-which(add_cell$group=="PD"),]
results <- add_cell2 %>%
  group_by(cell_type) %>%
  summarise(
    p_value = wilcox.test(Module_regulon1 ~ group)$p.value,
    statistic = wilcox.test(Module_regulon1 ~ group)$statistic
  )

# 查看结果
print(results)




