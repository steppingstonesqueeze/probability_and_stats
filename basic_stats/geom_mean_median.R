#!/usr/bin/env Rscript
# geom_mean_median.R
# Plot mean vs. median of the geometric distribution

# Load required package
if(!require(ggplot2)) install.packages("ggplot2", repos="https://cloud.r-project.org")

# Sequence of p values
p_vals <- seq(0.01, 0.99, length.out = 200)

# Mean = 1/p
means <- 1 / p_vals

# Median = ceil(log(0.5) / log(1-p))
medians <- ceiling(log(0.5) / log(1 - p_vals))

df <- data.frame(p = p_vals, Mean = means, Median = medians)

# Plot
library(ggplot2)
ggplot(df, aes(x=p)) +
  geom_line(aes(y=Mean, color="Mean (1/p)"), size=1) +
  geom_line(aes(y=Median, color="Median"), linetype="dashed", size=1) +
  coord_cartesian(ylim=c(0,20), xlim=c(0,1)) +
  labs(title="Mean vs Median of Geometric Distribution",
       x="Success probability p", y="Value") +
  scale_color_manual("", values=c("Mean (1/p)"="blue", "Median"="red")) +
  theme_minimal(base_size=14)
