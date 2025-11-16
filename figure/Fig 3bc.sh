
library(ggplot2)
library(gridExtra)


df <- read.table("cor.txt", header=FALSE)
colnames(df) <- c("col1", "col2", "trait", "col4", "col5")


df$sum <- df$col2

res1 <- cor.test(df$sum, df$col4, method="pearson")
res2 <- cor.test(df$sum, df$col5, method="pearson")

cat("Correlation between sum and col4:\n")
cat("Pearson R =", res1$estimate, "\n")
cat("P-value =", res1$p.value, "\n\n")

cat("Correlation between sum and col5:\n")
cat("Pearson R =", res2$estimate, "\n")
cat("P-value =", res2$p.value, "\n")

p1 <- ggplot(df, aes(x=sum, y=col4)) +
  geom_point(color="#1f77b4") +
  geom_smooth(method="lm", color="red") +
  labs(
    title="Sum vs Col4",
    subtitle=paste0("R = ", round(res1$estimate, 3), 
                    ", P = ", signif(res1$p.value, 3)),
    x="Sum (col1 + col2)", y="Col4"
  ) +
  theme_bw() +
  theme(panel.grid = element_blank())

p2 <- ggplot(df, aes(x=sum, y=col5)) +
  geom_point(color="#ff7f0e") +
  geom_smooth(method="lm", color="red") +
  labs(
    title="Sum vs Col5",
    subtitle=paste0("R = ", round(res2$estimate, 3), 
                    ", P = ", signif(res2$p.value, 3)),
    x="Sum (col1 + col2)", y="Col5"
  ) +
  theme_bw() +
  theme(panel.grid = element_blank())


combined_plot <- grid.arrange(p1, p2, ncol=1)
print(combined_plot) 

pdf("correlation_plot.pdf", width=8, height=10)
grid.arrange(p1, p2, ncol=1)
dev.off()
