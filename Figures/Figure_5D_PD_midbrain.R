
##----Figure 5D-----(These codes also contain old Sup_Figure S24 and S25, which have been deprecated)

cell_counts <- table(PD_seu$cell_type)


cell_proportions <- prop.table(cell_counts)

PD_seu@meta.data$cell_type <- factor(PD_seu@meta.data$cell_type,levels=c("ODC","AS","MG","N","EC","OPC","T"))

library(ggplot2)

# 转换为数据框格式
df <- as.data.frame(cell_proportions)
colnames(df) <- c("Cell_Type", "Proportion")

# 绘制条形图
ggplot(df, aes(x = Cell_Type, y = Proportion, fill = Cell_Type)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Cell Type Proportions", x = "Cell Type", y = "Proportion") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



library(dplyr)

# 按样本和细胞类型分组，计算细胞数量
cell_counts_by_sample <- PD_seu@meta.data %>%
  group_by(group, cell_type) %>%
  summarise(Count = n()) %>%
  mutate(Proportion = Count / sum(Count))

pal <- c("#A13B46","#67ADB7","#36600E","#6A8473","#C0BFDF","#7E6148FF",
         "#85C89C")

# 按特定顺序设置 group 的因子水平
cell_counts_by_sample$group <- factor(cell_counts_by_sample$group, levels = c("Young", "Aged", "PD"))

ggplot(cell_counts_by_sample, aes(x = group, y = Proportion, fill = cell_type)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_minimal() +
  labs(title = "Cell Type Proportions by Sample", x = "Sample", y = "Proportion") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_manual(values = pal)


library(dplyr)

# 按样本、性别和细胞类型分组，计算细胞数量
cell_counts_by_sample2 <- PD_seu@meta.data %>%
  group_by(orig.ident, group, sex, cell_type) %>%
  summarise(Count = n(), .groups = "drop") %>%
  mutate(Proportion_average = Count / sum(Count))



# Step 2: 将比例加入到元数据中
# 为每个细胞找到对应的比例值
PD_seu@meta.data <- PD_seu@meta.data %>%
  left_join(cell_counts_by_sample2, by = c("group", "cell_type"))

# Step 3: 检查是否成功加入
head(PD_seu@meta.data,20)




# 按特定顺序设置 group 的因子水平
cell_counts_by_sample2$group <- factor(cell_counts_by_sample2$group, levels = c("Young", "Aged", "PD"))

# 定义颜色调色板
pal <- c("#A13B46", "#67ADB7", "#36600E", "#6A8473", "#C0BFDF", "#7E6148FF", "#85C89C")

# 绘制堆叠条形图并按性别分面
ggplot(cell_counts_by_sample2, aes(x = group, y = Proportion_average, fill = cell_type)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_minimal() +
  labs(title = "Cell Type Proportions by Sample and Sex", x = "Sample Group", y = "Proportion") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = pal) +
  facet_wrap(~sex)




library(dplyr)
# 提取相关信息
meta_data <- as.data.frame(PD_seu@meta.data)

data_for_analysis <- meta_data[,c("group","age","Proportion_average","sex.y","cell_type")]

# 去重数据
data_for_analysis_unique <- data_for_analysis %>%
  distinct()

# 检查去重后的数据
head(data_for_analysis_unique)

data_for_analysis <- data_for_analysis_unique

write.csv(data_for_analysis,file = "data_for_analysis_cell_proportion.csv",quote=F, row.names = F)

# 检查提取的数据
head(data_for_analysis,36)
length(unique(data_for_analysis$cell_type))
length(unique(data_for_analysis$group))
length(unique(data_for_analysis$age))
length(data_for_analysis$group)


ggplot(data_for_analysis, aes(x = age, y = Proportion_average, color = sex.y)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal() +
  labs(
    title = "Correlation between Age and Proportion by Cell Type",
    x = "Age",
    y = "Proportion",
    color = "Sex"
  ) +
  scale_color_manual(values = c("M" = "#e31a1c", "F" = "#1f78b4")) +
  facet_wrap(~cell_type, scales = "free_y")



##---boxplot
data_for_analysis$group <- factor(data_for_analysis$group, levels = c("Young", "Aged","PD"))

ggplot(data_for_analysis, aes(x = group, y = Proportion_average, fill = group)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "Distribution of Proportion by Cell Type and Sex",
    x = "",
    y = "Cell proportion",
    fill = "Group"
  ) +
  #scale_fill_manual(values = c("M" = "#e31a1c", "F" = "#1f78b4")) +
  facet_wrap(~cell_type, scales = "free_y")


##--------去除---PD group
data_for_analysis <- data_for_analysis[-which(data_for_analysis$group =="PD"),]

ggplot(data_for_analysis, aes(x = age, y = Proportion_average, color = sex.y)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal() +
  labs(
    title = "Correlation between Age and Proportion by Cell Type",
    x = "Age",
    y = "Proportion",
    color = "Sex"
  ) +
  scale_color_manual(values = c("M" = "#e31a1c", "F" = "#1f78b4")) +
  facet_wrap(~cell_type, scales = "free_y")





ggplot(data_for_analysis, aes(x = group, y = Proportion_average, fill = sex.y)) +
  geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.75)) +
  theme_minimal() +
  labs(
    title = "Distribution of Proportion by Cell Type and Sex",
    x = "",
    y = "Cell proportion",
    fill = "Sex"
  ) +
  #scale_fill_manual(values = c("Aged" = "#e31a1c", "Young" = "#1f78b4")) +
  facet_wrap(~cell_type, scales = "free_y") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 10)
  )


library(dplyr)

# 按细胞类型分组并执行 Wilcoxon 检验

data_for_analysis2 <- data_for_analysis[-which(data_for_analysis$group=="Aged"),]
results <- data_for_analysis2 %>%
  group_by(cell_type) %>%
  summarise(
    p_value = wilcox.test(Proportion_average ~ group)$p.value,
    statistic = wilcox.test(Proportion_average ~ group)$statistic
  )

# 查看结果
print(results)

data_for_analysis2 <- data_for_analysis[-which(data_for_analysis$group=="Young"),]
results <- data_for_analysis2 %>%
  group_by(cell_type) %>%
  summarise(
    p_value = wilcox.test(Proportion_average ~ group)$p.value,
    statistic = wilcox.test(Proportion_average ~ group)$statistic
  )

# 查看结果
print(results)


data_for_analysis2 <- data_for_analysis[-which(data_for_analysis$group=="PD"),]
results <- data_for_analysis2 %>%
  group_by(sex.y,cell_type) %>%
  summarise(
    p_value = wilcox.test(Proportion_average ~ group)$p.value,
    statistic = wilcox.test(Proportion_average ~ group)$statistic
  )

# 查看结果
print(results)




library(dplyr)

# 按 cell_type 和 sex 分组，计算 Age 和 Proportion 的相关系数
correlation_results <- data_for_analysis %>%
  group_by(cell_type, sex.y) %>%
  summarise(
    correlation = cor(Proportion_average, age, method = "pearson", use = "complete.obs"),
    p_value = cor.test(Proportion_average, age, method = "pearson")$p.value,
    .groups = "drop"
  ) %>%
  mutate(FDR=p.adjust(p_value, method="fdr"))

# 查看相关系数结果
print(correlation_results)


library(ggplot2)

# 将相关性结果转为热图所需的格式
correlation_results$sex <- as.factor(correlation_results$sex.y)

ggplot(correlation_results, aes(x = cell_type, y = sex.y, fill = correlation)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#1f78b4", mid = "white", high = "#e31a1c", midpoint = 0) +
  theme_minimal() +
  labs(
    title = "Correlation between Age and Proportion by Cell Type and Sex",
    x = "Cell Type",
    y = "Sex",
    fill = "Correlation"
  )







ggplot(correlation_results, aes(x = cell_type, y = correlation, fill = sex.y)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(
    title = "Correlation between Age and Proportion by Cell Type and Sex",
    x = "Cell Type",
    y = "Correlation",
    fill = "Sex"
  )





