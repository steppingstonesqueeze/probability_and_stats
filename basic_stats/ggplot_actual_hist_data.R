# simple trick to get histogram data (actual histo data) from ggplot

library(ggplot2)

x <- rnorm(1000, 0, 1)
df <- data.frame(x = x)

gg <- ggplot(df, aes(x = x)) + geom_histogram(binwidth = 0.5, 
                                              fill = "red", 
                                              colour = "black",
                                              alpha = 0.4)
built <- ggplot_build(gg)
counts <- built$data[[1]]   # the histogram data
head(counts)

gg