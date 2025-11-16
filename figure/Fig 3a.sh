



library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(tidyr)

data <- read.table("sig.with_hits.4cols.sx.tsv", header = FALSE, 
                   col.names = c("Trait", "Cell", "Value"))


trait_order <- c("BAS", "EOS", "Hb", "Ht", "LYM", "MCH", "MCHC", "MCV", 
                 "MON", "NEU", "PLT", "RBC", "WBC", "height", "BMI", 
                 "HLD", "HTN", "T2D", "MetS", "RA", "SLE", "T1D", 
                 "CD", "UC", "MS", "GD", "BC", "LC", "PCa", "BCC")


cell_order <- c("BNaive", "BMem", "CD4Naive", "CD4EM", "CD4Treg", 
                "CD8Naive", "CD8GZMH", "CD8GZMK", "MAIT", "γδT", 
                "NKDim", "NKBright", "cM", "ncM")

actual_traits <- intersect(trait_order, unique(data$Trait))
actual_cells <- intersect(cell_order, unique(data$Cell))

data$Trait <- factor(data$Trait, levels = actual_traits)
data$Cell <- factor(data$Cell, levels = actual_cells)

create_size_only_dotplot <- function(data, value_col = "Value") {

  size_range <- range(data[[value_col]], na.rm = TRUE)

  all_combinations <- expand.grid(
    Cell = levels(data$Cell),
    Trait = levels(data$Trait)
  )
  
  p <- ggplot(all_combinations, aes(x = Cell, y = Trait)) +
    geom_tile(fill = "white", color = "grey80", size = 0.5) +  
    geom_point(aes(size = get(value_col)), 
               data = data,  
               shape = 21, fill = "grey70", color = "black", stroke = 0.2) +
    scale_size_continuous(name = "Proportion", 
                          range = c(0.2, 6),
                          limits = c(0, max(size_range) * 1.1)) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white"),
      legend.position = "right"
    ) +
    labs(title = "Trait-Cell Association")
  
  return(p)
}


p_size_only <- create_size_only_dotplot(data)
print(p_size_only)

ggsave("trait_cell_heatmap_size_only.pdf", p_size_only, width = 14, height = 10)
ggsave("trait_cell_heatmap_size_only.png", p_size_only, width = 14, height = 10, dpi = 300)


create_simple_dotplot <- function(data, value_col = "Value") {
  
  size_range <- range(data[[value_col]], na.rm = TRUE)
  
  p <- ggplot(data, aes(x = Cell, y = Trait)) +
    geom_point(aes(size = get(value_col)), 
               shape = 21, fill = "grey70", color = "black", stroke = 0.2) +
    scale_size_continuous(name = "Proportion", 
                          range = c(0.2, 6),
                          limits = c(0, max(size_range) * 1.1)) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_blank(),
      panel.grid.major = element_line(color = "grey90", size = 0.5),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    ) +
    labs(title = "Trait-Cell Association (Simple)")
  
  return(p)
}


p_simple <- create_simple_dotplot(data)
print(p_simple)
ggsave("trait_cell_heatmap_simple.pdf", p_simple, width = 14, height = 10)


cat("\n数据摘要:\n")
cat("总行数:", nrow(data), "\n")
cat("Value 范围:", range(data$Value), "\n")
cat("实际Traits数:", length(actual_traits), "\n")
cat("实际Cells数:", length(actual_cells), "\n")