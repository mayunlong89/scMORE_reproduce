
###----alpha = 0
# 读取数据（确保为CSV格式）
#df <- read.csv("your_file.csv")  # ← 改为你的实际文件路径
df <- results_IBD_phastCons_withoutExons_alpha0$scMore_trait_results
# 计算三项指标的比例值
df <- df %>%
  mutate(sum_score = CTS + GRS + TRS,
         CTS_prop = CTS / sum_score,
         GRS_prop = GRS / sum_score,
         TRS_prop = TRS / sum_score)

# 绘制三元图，按Significance 上色
# 绘制三元图
ggtern(data = df, aes(x = CTS_prop, y = GRS_prop, z = TRS_prop)) +
  geom_point(aes(color = Significance), size = 2.5, alpha = 0.85) +
  scale_color_manual(values = c("Significant" = "red", "Nonsignificant" = "grey70")) +
  theme_bw() +
  labs(title = "Ternary Plot of CTS / GRS / TRS",
       x = "CTS", y = "GRS", z = "TRS", color = "Significance") +
  theme(legend.position = "right")


# 选择某个细胞类型，例如 Monocytes
df_mono <- df %>% filter(Celltype == "Monocytes")

# 计算三个指标的相对比例（归一化处理）
df_mono <- df_mono %>%
  mutate(sum_score = CTS + GRS + TRS,
         CTS_prop = CTS / sum_score,
         GRS_prop = GRS / sum_score,
         TRS_prop = TRS / sum_score)


ggtern(data = df_mono, aes(x = CTS_prop, y = GRS_prop, z = TRS_prop)) +
  geom_point(aes(color = Significance), size = 2.5, alpha = 0.85) +
  scale_color_manual(values = c("Significant" = "red", "Nonsignificant" = "grey70")) +
  theme_bw() +
  labs(title = "Ternary Plot of CTS / GRS / TRS",
       x = "CTS", y = "GRS", z = "TRS", color = "Significance") +
  theme(legend.position = "right")




###----alpha = 1

df <- read.csv("~/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2025-05-20-scMORE_revision/05_parameter_assessment/results_IBD_scMore_trait_results_phastCons_withoutExons.csv")
# 计算三项指标的比例值
df <- df %>%
  mutate(sum_score = CTS + GRS + TRS,
         CTS_prop = CTS / sum_score,
         GRS_prop = GRS / sum_score,
         TRS_prop = TRS / sum_score)

# 绘制三元图，按 Significance 上色
ggtern(data = df, aes(x = CTS_prop, y = GRS_prop, z = TRS_prop)) +
  geom_point(aes(color = Significance), size = 2.5, alpha = 0.85) +
  scale_color_manual(values = c("Significant" = "red", "Nonsignificant" = "grey70")) +
  theme_bw() +
  labs(title = "Ternary Plot of CTS / GRS / TRS",
       x = "CTS", y = "GRS", z = "TRS", color = "Significance") +
  theme(legend.position = "right")




# 选择某个细胞类型，例如 Monocytes
df_mono <- df %>% filter(Celltype == "Monocytes")

# 计算三个指标的相对比例（归一化处理）
df_mono <- df_mono %>%
  mutate(sum_score = CTS + GRS + TRS,
         CTS_prop = CTS / sum_score,
         GRS_prop = GRS / sum_score,
         TRS_prop = TRS / sum_score)


ggtern(data = df_mono, aes(x = CTS_prop, y = GRS_prop, z = TRS_prop)) +
  geom_point(aes(color = Significance), size = 2.5, alpha = 0.85) +
  scale_color_manual(values = c("Significant" = "red", "Nonsignificant" = "grey70")) +
  theme_bw() +
  labs(title = "Ternary Plot of CTS / GRS / TRS",
       x = "CTS", y = "GRS", z = "TRS", color = "Significance") +
  theme(legend.position = "right")

