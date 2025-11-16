


library(ggplot2)
library(ggpmisc)

cell_type_order <- c("BNaive", "BMem", "CD4Naive", "CD4EM", "CD4Treg", "CD8Naive",
                     "CD8GZMH", "CD8GZMK", "MAIT", "γδT", "NKDim", "NKBright", "cM", "ncM")

color_mapping <- c(
  "BNaive" = "#01579B",
  "BMem" = "#A5EDFF",
  "CD4EM" = "#1B5E20",
  "CD4Naive" = "#8BC34A",
  "CD4Treg" = "#D4E157",
  "CD8GZMH" = "#FF5722",
  "CD8GZMK" = "#FFCDD2",
  "CD8Naive" = "#F9A825",
  "MAIT" = "#FFEE58",
  "cM" = "#6A1B9A",
  "ncM" = "#E1BEE7",
  "NKBright" = "#82581F",
  "NKDim" = "#D7CCC8",
  "γδT" = "#4DD0E1"
)

B <- c(6357, 6535, 8004, 8004, 5386, 6535, 7074, 6650, 5409, 5572, 7410, 5409, 6430, 6133)
A <- c(873, 866, 972, 977, 748, 755, 921, 934, 296, 255, 974, 289, 569, 446)

df <- data.frame(
  CellType = cell_type_order,
  A = A,
  B = B
)

result <- cor.test(A, B, method = "pearson")
correlation <- result$estimate
correlation_label <- paste0("Pearson's r = ", round(correlation, 2), "\nP = ", signif(result$p.value, 3))

p1 <- ggplot(df, aes(x = A, y = B)) +
  geom_point(aes(color = CellType), size = 5, shape = 21, fill = "white", stroke = 2) +
  scale_color_manual(values = color_mapping[df$CellType]) +
  stat_poly_eq(formula = y ~ x, aes(label = ..eq.label..), size = 4, color = "darkred") +
  stat_poly_line(formula = y ~ x, color = "darkblue", size = 1.5) +
  annotate("text", x = Inf, y = Inf, label = correlation_label, hjust = 1.1, vjust = 1.5, size = 5, color = "darkgreen") +
  theme_classic() +
  xlab("A") +
  ylab("B") +
  guides(color = FALSE) +
  theme(
    plot.background = element_rect(fill = "white", color = NA), 
    panel.background = element_rect(fill = "white", color = NA),
    axis.title = element_text(face = "bold", size = 16, color = "black"),
    axis.text = element_text(size = 22, color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 18, color = "darkblue"),
    plot.subtitle = element_text(size = 22, hjust = 0.5, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    axis.line = element_line(color = "black", size = 0.6)
  ) +
  scale_y_continuous(limits = c(0, max(B)*1.1)) 
print(p1)
ggsave("prenumber_vs_cellnumber_correlation.pdf", p1, height = 7, width = 7)
