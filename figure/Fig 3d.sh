


library(reshape2)
library(ggplot2)
library(patchwork)
library(scales)

df <- read.table("all_3", header = FALSE,
                 col.names = c("Trait","col2","col3","col4","col5","total"))
df$part2 <- df$col4
df$part1 <- df$col2 - df$col4
df$sum   <- df$part1 + df$part2
df$part1_pct <- ifelse(df$sum > 0, df$part1 / df$sum * 100, 0)
df$part2_pct <- ifelse(df$sum > 0, df$part2 / df$sum * 100, 0)

bef <- read.table("bef_table_sort", header = FALSE,
                  col.names = c("Trait","Count"), stringsAsFactors = FALSE)

bef <- merge(df["Trait"], bef, by = "Trait", all.x = TRUE)
bef$Count[is.na(bef$Count)] <- 0

ord <- df$Trait 
bef$Trait <- factor(bef$Trait, levels = ord)
df$Trait  <- factor(df$Trait,  levels = ord)

col_orange <- "#ef8636"
col_blue   <- "#1f77b4"

p0 <- ggplot(bef, aes(Trait, Count)) +
  geom_col(fill = col_orange, color = "black", width = 0.7) +
  geom_text(aes(label = Count), vjust = -0.25, size = 3) +
  labs(y = "Count") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )

df_melt_count <- melt(df[, c("Trait","part1","part2")], id.vars = "Trait")
df_melt_count$variable <- factor(df_melt_count$variable, levels = c("part1","part2"))

p1 <- ggplot(df_melt_count, aes(Trait, value, fill = variable)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  scale_fill_manual(values = c(part1 = col_orange, part2 = col_blue),
                    labels = c("Part1","Part2"), name = "Component") +
  labs(y = "Count") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "top"
  )

df_melt_pct <- melt(df[, c("Trait","part1_pct","part2_pct")], id.vars = "Trait")
df_melt_pct$variable <- factor(df_melt_pct$variable, levels = c("part1_pct","part2_pct"))

p2 <- ggplot(df_melt_pct, aes(Trait, value, fill = variable)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  scale_y_continuous(labels = percent_format(scale = 1), limits = c(0, 100)) +
  scale_fill_manual(values = c(part1_pct = col_orange, part2_pct = col_blue),
                    labels = c("Part1","Part2"), name = "Component") +
  labs(x = "Trait", y = "Percentage") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    legend.position = "none"
  )

combined_plot <- p0 / p1 / p2 + plot_layout(heights = c(0.9, 1.1, 1.1))
combined_plot 
ggsave("combined_plot_original_order.pdf", combined_plot, width = 11, height = 10)
ggsave("combined_plot_original_order.png", combined_plot, width = 11, height = 10, dpi = 300)
