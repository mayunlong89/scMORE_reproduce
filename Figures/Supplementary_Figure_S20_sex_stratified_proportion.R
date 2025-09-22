
####-------按照sex stratified proportion
####-------按照sex stratified proportion
####-------按照sex stratified proportion
####-------按照sex stratified proportion
library(dplyr)
library(tidyr)
library(compositions)
library(ggplot2)
library(ggpubr)

# === Step 1: 从 meta.data 提取必要列，包括 sex ===
meta_df <- PD_seu@meta.data %>%
  dplyr::select(orig.ident, group, sex, cell_type)

# === Step 2: 每个 sample 内统计每种 cell_type 的数量 ===
cell_counts <- meta_df %>%
  group_by(orig.ident, group, sex, cell_type) %>%
  summarise(count = n(), .groups = "drop")

# === Step 3: 计算每个 sample 的 cell type 组成比例 ===
cell_prop <- cell_counts %>%
  group_by(orig.ident) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

# === Step 4: pivot_wider，每行一个样本 ===
prop_wide <- cell_prop %>%
  pivot_wider(
    id_cols = c(orig.ident, group, sex),
    names_from = cell_type,
    values_from = proportion,
    values_fill = 0
  )

# === Step 5: CLR 转换 ===
meta_info <- prop_wide %>% dplyr::select(orig.ident, group, sex)
prop_matrix <- prop_wide %>% dplyr::select(-orig.ident, -group, -sex)

clr_matrix <- clr(prop_matrix + 1e-10)
clr_df <- cbind(meta_info, clr_matrix)

# === Step 6: 转换为长格式 + 筛选去除PD组 ===
clr_long <- clr_df %>%
  pivot_longer(cols = -c(orig.ident, group, sex),
               names_to = "cell_type",
               values_to = "clr_value") %>%
  filter(group != "PD")  # 只保留 Young 和 Aged

# 设置顺序
clr_long$group <- factor(clr_long$group, levels = c("Young", "Aged"))
clr_long$sex <- factor(clr_long$sex, levels = c("M", "F"))

# === Step 7: 按 sex 分组画图 ===

# 筛选只有 Young 和 Aged
clr_sub <- clr_long %>% filter(group %in% c("Young", "Aged"))

# 分别计算每个 cell_type + sex 中 Young vs Aged 的 p 值
pvals_df <- clr_sub %>%
  group_by(cell_type, sex) %>%
  summarise(
    p = tryCatch(wilcox.test(clr_value ~ group)$p.value, error = function(e) NA),
    .groups = "drop"
  ) %>%
  mutate(
    y.position = max(clr_sub$clr_value, na.rm = TRUE) + 1,
    group1 = "Young",
    group2 = "Aged",
    p.signif = case_when(
      p <= 0.001 ~ "***",
      p <= 0.01 ~ "**",
      p <= 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

ggplot(clr_sub, aes(x = group, y = clr_value, fill = sex)) +
  geom_boxplot(outlier.size = 0.5) +
  stat_pvalue_manual(
    pvals_df,
    label = "p.signif",
    tip.length = 0.01,
    bracket.size = 0.3
  ) +
  facet_wrap(~cell_type, scales = "free_y") +
  scale_fill_manual(values = c("M" = "#E78AC3", "F" = "#66C2A5")) +
  labs(
    title = "CLR-transformed Cell Type Proportions by Group (Young vs Aged within each Sex)",
    x = "Group",
    y = "CLR value"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 12, face = "bold")
  )


# ==== 1. 定义去除异常值函数 ====
remove_outliers <- function(x) {
  q1 <- quantile(x, 0.25, na.rm = TRUE)
  q3 <- quantile(x, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  lower <- q1 - 1.5 * iqr
  upper <- q3 + 1.5 * iqr
  x[x >= lower & x <= upper]
}

# 手动设置 cell type 和 sex
cell <- "T"
sex_group <- "M"

# 手动提取子集：cell_type + sex + group（Young vs Aged）
df_sub <- clr_long %>% filter(cell_type == cell, sex == sex_group, group %in% c("Young", "Aged"))

# 分组并去除异常值
y <- remove_outliers(df_sub$clr_value[df_sub$group == "Young"])
a <- remove_outliers(df_sub$clr_value[df_sub$group == "Aged"])

# Wilcoxon 检验
wilcox.test(y, a)
