# Load necessary libraries
rm(list = ls())
library(ggplot2)
library(cowplot)

# Define the range for x and k
x <- seq(0, 1, by = 0.01)
k <- seq(0, 50, by = 0.01)
kmax <- max(k)
xmax <- max(x)

# 1. Fix k = 0.5 and vary x
kfix <- 9
df1 <- data.frame(x = x, k = kfix)
df1$z <- with(df1, 1 / (1 - x * (1 - 1/(k + 1))))
ymax <- 1 / (1 - (1 - xmax * 1/(kfix + 1)))
df1$ratio <- df1$z/ymax

plot1 <- ggplot(df1, aes(x=x, y=ratio)) +
  geom_line(color="#D77186") +
  labs(title = "Given 90% unlabeled data") +
  geom_hline(yintercept = 1, linetype="dashed", color = "#E2C43F", size=1) +
  xlab(expression("Imputation " * r^2)) +
  ylab(expression(N[eff]/ N[total])) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.1, 0.25, 0.5, 0.75, 1), labels = scales::percent_format(accuracy = 1)) +
  theme_classic()

# 2. Fix x = 0.5 and vary k
xfix <- 0.5
df2 <- data.frame(x = xfix, k = k)
df2$z <- with(df2, 1 / (1 - x * (1 - 1/(k + 1))))
ymax <- 1/(1 - xfix)

plot2 <- ggplot(df2, aes(x=k, y=z)) +
  geom_line(color="#5c94c4") +
  labs(title=expression(paste("Given imputation " * r^2 == 0.5))) +
  geom_hline(yintercept = ymax, linetype="dashed", color = "#66A64F", size=1) +
  xlab(expression(N[unlab]/N[lab])) +
  ylab(expression(N[eff]/ N[lab])) +
  scale_y_continuous(limits = c(1, 2), breaks = c(1, 1.25, 1.5, 1.75, 2), labels = scales::percent_format(accuracy = 1)) +
  theme_classic() 

# Display the plots
print(plot1)
print(plot2)

p_all <- plot_grid(plot1, plot2, labels = c('a', 'b'), nrow = 1)

FigDir="./fig4/"
pdf_out_dir <-  paste0(FigDir, "/Fig4.pdf")
pdf(pdf_out_dir, width = 8, height = 4)
par(mar=c(5,5,4,4)+.1)
print(p_all)
dev.off()