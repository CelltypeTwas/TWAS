

file1 <- "adjust1_1.extend1MB.part_all_ld.txt"
file2 <- "all.txt"

data1 <- scan(file1)
data2 <- scan(file2)

summary(data1)
summary(data2)

t_res <- t.test(data1, data2)
print(t_res)

library(ggplot2)
df <- data.frame(
  value = c(data1, data2),
  group = factor(c(rep("Simulation", length(data1)), rep("Observed", length(data2))))
)

p1 <- ggplot(df, aes(x = group, y = value, fill = group)) +
  geom_boxplot(alpha = 0.6, width = 0.4, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.6) +
  theme_bw(base_size = 14) +
  labs(title = "Comparison of R² Distributions", y = "R² Value", x = "") +
  scale_fill_manual(values = c("#619CFF", "#F8766D")) +
  theme(legend.position = "none")

p2 <- ggplot(df, aes(x = value, fill = group, color = group)) +
  geom_density(alpha = 0.4) +
  theme_bw(base_size = 14) +
  labs(title = "Density Plot of R² Values", x = "R² Value", y = "Density") +
  scale_fill_manual(values = c("#619CFF", "#F8766D")) +
  scale_color_manual(values = c("#619CFF", "#F8766D"))

print(p1)
print(p2)

ggsave("ld_comparison_boxplot.pdf", p1, width = 6, height = 5)
ggsave("ld_comparison_density.pdf", p2, width = 6, height = 5)