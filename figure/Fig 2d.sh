



library(ggplot2)
files <- list.files(pattern = "*.csv")
all_data <- data.frame(CellType = character(), Value = numeric())

for (f in files) {
  dat <- read.table(f, header = FALSE)
  celltype <- sub(".csv", "", f)
  temp_df <- data.frame(CellType = celltype, Value = dat$V1)
  all_data <- rbind(all_data, temp_df)
}

cell_type_order <- c("BNaive", "BMem", "CD4Naive", "CD4EM", "CD4Treg", 
                     "CD8Naive", "CD8GZMH", "CD8GZMK", "MAIT", "γδT", 
                     "NKDim", "NKBright", "cM", "ncM")

all_data$CellType <- factor(all_data$CellType, levels = cell_type_order)
overall_median <- median(all_data$Value, na.rm = TRUE)
new_celltype_pal <- c("#01579B","#A5EDFF", "#1B5E20","#8BC34A","#D4E157",
                      "#FF5722", "#FFCDD2","#F9A825","#FFEE58","#4DD0E1",
                      "#D7CCC8","#82581F","#6A1B9A","#E1BEE7")

p <- ggplot(all_data, aes(x = CellType, y = Value, color = CellType)) +
  geom_violin(
    fill = NA, trim = FALSE, linewidth = 1, scale = "width"
  ) +
  geom_boxplot(
    width = 0.15, fill = NA, outlier.shape = NA, linewidth = 1
  ) +
  geom_hline(
    yintercept = overall_median, color = "red", linetype = "solid", linewidth = 1
  ) +
  scale_color_manual(values = new_celltype_pal) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.position = "none"
  ) +
  labs(x = "Cell Type", y = "Value") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

print(p)
ggsave("violin_plot.pdf", plot = p, width = 10, height = 6, dpi = 300, units = "in", bg = "white")

