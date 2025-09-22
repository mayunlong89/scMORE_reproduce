
##---Ns---Figure 7B
Module <- regulons[regulons$TF == "STAT4",]
Module <- regulons[regulons$TF == "BCL11A",]
Module <- regulons[regulons$TF == "DMBX1",]


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
    title = "Distribution of regulon activity by Cell Type and Sex",
    x = "",
    y = "regulon activity",
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



