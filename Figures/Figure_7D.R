##------Figure 7D

Module <- regulons[regulons$TF == "DMBX1",]
Module <- regulons[regulons$TF == "BCL11A",]
Module <- regulons[regulons$TF == "STAT4",]
Module <- regulons[regulons$TF == "CUX2",]


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


##---boxplot----cell type---age x sex interaction
add_cell$group <- factor(add_cell$group, levels = c("Young", "Aged", "PD"))

ggplot(add_cell, aes(x = group, y = Module_regulon1, fill = sex)) + 
  geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.75)) + 
  theme_minimal() + 
  labs(
    title = "Distribution of regulon activity by Cell Type, Group, and Sex",
    x = "",
    y = "Regulon Activity",
    fill = "Sex"
  ) + 
  scale_fill_manual(values = c("M" = "#e31a1c", "F" = "#1f78b4")) + 
  facet_wrap(~cell_type, scales = "free_y") + 
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )





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

