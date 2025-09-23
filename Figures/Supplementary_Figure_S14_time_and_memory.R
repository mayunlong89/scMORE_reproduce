library(ggplot2)
library(dplyr)
library(stringr)


time_PD <- read.csv("/share/pub/guiyy/Cerebral_comorbidities/time_PD.csv", header = TRUE)

df_summary <- time_PD %>%
  group_by(cell, method) %>%
  summarise(
    mean_time = mean(Time),
    se_time = sd(Time)/sqrt(n()),
    .groups = "drop"
  )

stat.test <- time_PD %>%
  group_by(cell) %>%
  wilcox_test(Time ~ method) %>%
  add_xy_position(x = "cell", fun = "mean_sd", dodge = 0.9) 

p <- ggplot(df_summary, aes(x=factor(cell), y=mean_time, fill=method)) +
  geom_bar(stat="identity", position=position_dodge(width=0.9)) +
  geom_errorbar(aes(ymin=mean_time - se_time, ymax=mean_time + se_time),
                width=0.2, position=position_dodge(width=0.9)) +
  stat_pvalue_manual(stat.test, label = "p.adj", 
                     tip.length = 0.0, bracket.nudge.y = -2, 
                     inherit.aes = FALSE) +
  theme_classic(base_size=14) +
  labs(x="Number of Cells", y="Time (s)") +
  theme(
    legend.position="right",
    legend.key.size = unit(0.3, "cm"),
    axis.text.x = element_text(angle=45, hjust=1)
  )
p
