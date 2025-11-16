
library(ggplot2)

CellType <- c("BNaive","BMem","CD4Naive","CD4EM","CD4Treg","CD8Naive",
              "CD8GZMH","CD8GZMK","MAIT","γδT","NKDim","NKBright","cM","ncM")
A <- c(6357, 6535, 8004, 8004, 5386, 6535, 7074, 6650, 5409, 5572, 7410, 5409, 6430, 6133)

df <- data.frame(
  CellType = factor(CellType, levels = CellType),
  A = A
)


med <- median(df$A)

p <- ggplot(df, aes(x = CellType, y = A)) +
  geom_col(fill = "#377eb8", color = "black", width = 0.7) +
  geom_text(aes(label = A), vjust = -0.25, size = 3) +
  geom_hline(yintercept = med, color = "red", linetype = "dashed", linewidth = 0.8) +
  labs(x = "Cell Type", y = "A") +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    plot.margin = margin(10, 20, 10, 10)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.07))) +
  coord_cartesian(clip = "off")


print(p)

# ----- 导出 -----
ggsave("CellType_A_Bar_withMedian.pdf", p, width = 8, height = 5)
ggsave("CellType_A_Bar_withMedian.png", p, width = 8, height = 5, dpi = 300)

