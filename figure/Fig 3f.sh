



library(ggplot2)
library(dplyr)

sharing_data <- data.frame(
  Category = c("Different lineage\nDifferent trait",
               "Different lineage\nSame trait",
               "Same lineage\nDifferent trait",
               "Same lineage\nSame trait"),
  Count = c(17, 23, 11, 36),
  Percentage = c(19.5, 26.4, 12.6, 41.4)
)

sharing_data <- sharing_data %>%
  mutate(
    Label = paste0(Count, "\n(", Percentage, "%)"),
    Category = factor(Category, levels = Category)
  )

colors <- c(
  "Different lineage\nDifferent trait" = "#5DBAA4",  
  "Different lineage\nSame trait" = "#F4976C",       
  "Same lineage\nDifferent trait" = "#7C9BC6",       
  "Same lineage\nSame trait" = "#E291B5"            
)


p1 <- ggplot(sharing_data, aes(x = "", y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 1, color = "white", size = 2) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = colors, name = "Sharing Pattern") +
  geom_text(aes(label = Label), 
            position = position_stack(vjust = 0.5),
            size = 5, fontface = "bold", color = "white") +
  theme_void() +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  ) +
  labs(title = "Sharing Pattern")


ggsave("sharing_pattern_ggplot2.png", p1, width = 10, height = 6, dpi = 300)
ggsave("sharing_pattern_ggplot2.pdf", p1, width = 10, height = 6)



png("sharing_pattern_base.png", width = 10, height = 6, units = "in", res = 300)

par(mar = c(2, 2, 3, 8), xpd = TRUE)


pie_data <- sharing_data$Count
pie_labels <- sharing_data$Label
pie_colors <- c("#5DBAA4", "#F4976C", "#7C9BC6", "#E291B5")

pie(pie_data, 
    labels = pie_labels,
    col = pie_colors,
    border = "white",
    lwd = 3,
    cex = 1.2,
    radius = 0.9,
    main = "Sharing Pattern",
    cex.main = 1.5,
    font.main = 2)

legend("right", 
       legend = c("Different lineage\nDifferent trait",
                  "Different lineage\nSame trait",
                  "Same lineage\nDifferent trait",
                  "Same lineage\nSame trait"),
       fill = pie_colors,
       border = "white",
       bty = "n",
       cex = 1.1,
       title = "Sharing Pattern",
       title.adj = 0)

dev.off()



  library(plotly)
  
  p3 <- plot_ly(sharing_data, 
                labels = ~Category, 
                values = ~Count,
                type = 'pie',
                textposition = 'inside',
                textinfo = 'label+percent',
                insidetextfont = list(color = '#FFFFFF', size = 14),
                marker = list(colors = pie_colors,
                             line = list(color = '#FFFFFF', width = 2)),
                showlegend = TRUE) %>%
    layout(title = list(text = "Sharing Pattern", 
                       font = list(size = 18, family = "Arial", color = "black")),
           xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
           yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
 


p4 <- ggplot(sharing_data, aes(x = reorder(Category, -Count), y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = colors) +
  geom_text(aes(label = paste0(Count, "\n(", Percentage, "%)")), 
            vjust = -0.5, size = 4, fontface = "bold") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "none",
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "Sharing Pattern - Bar Chart",
    x = "Category",
    y = "Count"
  ) +
  ylim(0, max(sharing_data$Count) * 1.15)

ggsave("sharing_pattern_barplot.png", p4, width = 10, height = 6, dpi = 300)

