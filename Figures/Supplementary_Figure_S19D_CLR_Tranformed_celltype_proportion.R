
##=========修稿的时候--------加上CLR转换--------
##=========修稿的时候--------加上CLR转换--------
##=========修稿的时候--------加上CLR转换--------
##=========修稿的时候--------加上CLR转换--------
###_------Figure S D---Cell proportion CLR transformation
#🧩 Load R packages
library(dplyr)
library(tidyr)
library(compositions)
library(ggplot2)
library(ggpubr)


# Step 1: 从 meta.data 提取必要列
meta_df <- PD_seu@meta.data %>%
  dplyr::select(orig.ident, group, cell_type)

# Step 2: 每个 sample 内统计每种 cell_type 的数量
cell_counts <- meta_df %>%
  group_by(orig.ident, group, cell_type) %>%
  summarise(count = n(), .groups = "drop")

# Step 3: 计算每个 sample 的 cell type 组成比例
cell_prop <- cell_counts %>%
  group_by(orig.ident) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

# Step 4: pivot_wider，确保每个样本一行、每列是 cell type
prop_wide <- cell_prop %>%
  pivot_wider(
    id_cols = c(orig.ident, group),
    names_from = cell_type,
    values_from = proportion,
    values_fill = 0
  )

##---test sample total proportion = 1
prop_wide %>% count(orig.ident)


# === Step 5: CLR 转换 ===
meta_info <- prop_wide %>% dplyr::select(orig.ident, group)
prop_matrix <- prop_wide %>% dplyr::select(-orig.ident, -group)

rowSums(prop_matrix) %>% summary()
cell_prop %>% count(orig.ident) %>% summary()


clr_matrix <- clr(prop_matrix + 1e-10)
clr_df <- cbind(meta_info, clr_matrix)

# === Step 6: 转换为长格式，用于绘图 ===
clr_long <- clr_df %>%
  pivot_longer(cols = -c(orig.ident, group),
               names_to = "cell_type",
               values_to = "clr_value")

clr_long$group <- factor(clr_long$group, levels = c("Young", "Aged", "PD"))

# === Step 7: 画图：按 group 分组、cell_type 分面 ===
ggplot(clr_long, aes(x = group, y = clr_value, fill = group)) +
  geom_boxplot(outlier.size = 0.5) +
  stat_compare_means(comparisons = list(c("Young", "Aged"), c("Young", "PD"), c("Aged", "PD")),
                     method = "wilcox.test", label = "p.signif") +
  facet_wrap(~cell_type, scales = "free_y") +
  scale_fill_manual(values = c("Young" = "#8DA0CB", "Aged" = "#FC8D62", "PD" = "#66C2A5")) +
  labs(
    title = "CLR-transformed Cell Type Proportions by Group",
    x = "Group",
    y = "CLR value"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 12, face = "bold"))




###_-----因为图中有outliers, 我想去除outliers, 再计算wilcox test
##---单独分析---一个细胞类型内的两组:
remove_outliers <- function(x) {
  q1 <- quantile(x, 0.25, na.rm = TRUE)
  q3 <- quantile(x, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  lower <- q1 - 1.5 * iqr
  upper <- q3 + 1.5 * iqr
  x[x >= lower & x <= upper]
}


library(dplyr)

df_sub <- clr_long %>% filter(cell_type == "AS", group %in% c("Young", "PD"))

y <- remove_outliers(df_sub$clr_value[df_sub$group == "Young"])
a <- remove_outliers(df_sub$clr_value[df_sub$group == "PD"])

wilcox.test(y, a)

df_sub <- clr_long %>% filter(cell_type == "AS", group %in% c("Young", "Aged"))

y <- remove_outliers(df_sub$clr_value[df_sub$group == "Young"])
a <- remove_outliers(df_sub$clr_value[df_sub$group == "Aged"])

wilcox.test(y, a)

df_sub <- clr_long %>% filter(cell_type == "AS", group %in% c("Aged", "PD"))

y <- remove_outliers(df_sub$clr_value[df_sub$group == "Aged"])
a <- remove_outliers(df_sub$clr_value[df_sub$group == "PD"])

wilcox.test(y, a)


df_sub <- clr_long %>% filter(cell_type == "MG", group %in% c("Young", "PD"))

y <- remove_outliers(df_sub$clr_value[df_sub$group == "Young"])
a <- remove_outliers(df_sub$clr_value[df_sub$group == "PD"])

wilcox.test(y, a)

df_sub <- clr_long %>% filter(cell_type == "MG", group %in% c("Young", "Aged"))

y <- remove_outliers(df_sub$clr_value[df_sub$group == "Young"])
a <- remove_outliers(df_sub$clr_value[df_sub$group == "Aged"])

wilcox.test(y, a)

df_sub <- clr_long %>% filter(cell_type == "MG", group %in% c("Aged", "PD"))

y <- remove_outliers(df_sub$clr_value[df_sub$group == "Aged"])
a <- remove_outliers(df_sub$clr_value[df_sub$group == "PD"])

wilcox.test(y, a)



df_sub <- clr_long %>% filter(cell_type == "N", group %in% c("Young", "Aged"))

y <- remove_outliers(df_sub$clr_value[df_sub$group == "Young"])
a <- remove_outliers(df_sub$clr_value[df_sub$group == "Aged"])

wilcox.test(y, a)


df_sub <- clr_long %>% filter(cell_type == "N", group %in% c("Young", "PD"))

y <- remove_outliers(df_sub$clr_value[df_sub$group == "Young"])
a <- remove_outliers(df_sub$clr_value[df_sub$group == "PD"])

wilcox.test(y, a)

df_sub <- clr_long %>% filter(cell_type == "N", group %in% c("Aged", "PD"))

y <- remove_outliers(df_sub$clr_value[df_sub$group == "Aged"])
a <- remove_outliers(df_sub$clr_value[df_sub$group == "PD"])

wilcox.test(y, a)






df_sub <- clr_long %>% filter(cell_type == "OPC", group %in% c("Young", "Aged"))

y <- remove_outliers(df_sub$clr_value[df_sub$group == "Young"])
a <- remove_outliers(df_sub$clr_value[df_sub$group == "Aged"])

wilcox.test(y, a)


df_sub <- clr_long %>% filter(cell_type == "OPC", group %in% c("Young", "PD"))

y <- remove_outliers(df_sub$clr_value[df_sub$group == "Young"])
a <- remove_outliers(df_sub$clr_value[df_sub$group == "PD"])

wilcox.test(y, a)

df_sub <- clr_long %>% filter(cell_type == "OPC", group %in% c("Aged", "PD"))

y <- remove_outliers(df_sub$clr_value[df_sub$group == "Aged"])
a <- remove_outliers(df_sub$clr_value[df_sub$group == "PD"])

wilcox.test(y, a)



df_sub <- clr_long %>% filter(cell_type == "ODC", group %in% c("Young", "PD"))

y <- remove_outliers(df_sub$clr_value[df_sub$group == "Young"])
a <- remove_outliers(df_sub$clr_value[df_sub$group == "PD"])

wilcox.test(y, a)


#✅ 方法 1：Kruskal-Wallis 检验（三组非参数比较）
library(ggpubr)

ggboxplot(clr_long, x = "group", y = "clr_value", fill = "group") +
  facet_wrap(~cell_type, scales = "free_y") +
  stat_compare_means(method = "kruskal.test", label = "p.format") +
  theme_minimal()


