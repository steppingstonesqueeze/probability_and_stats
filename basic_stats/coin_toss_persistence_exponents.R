# coin toss persistence exponents

library(ggplot2)
library(dplyr)

# Coin toss - histo of distance between successsive heads

sample_data <- c("H", "T")
N <- 100000
data <- sample(sample_data,
               N,
               replace = TRUE)

# Head locations
pos_H <- which(data == "H")
diff_pos_H <- diff(pos_H)

# Table
cat("Count table \n\n")
table(diff_pos_H)
Sys.sleep(2)

#Histo
h <- hist(diff_pos_H, col = "red", breaks = 20)

#Histo details

cat("Histogram counts \n")
h$counts
cat("Histogram breaks \n")
h$breaks
cat("Histogram midpoints \n")
h$mids
cat("Histogram density \n")
h$density

# lets now examine persistence

# heads fractions > 0.5 

data1 <- ifelse(data == "H", 1, 0)
heads_fraction <- cumsum(data1) / c(1:length(data1))
head(data1)
head(heads_fraction)

# simple plot
plot(heads_fraction, col = "blue")

heads_win <- ifelse(heads_fraction > 0.5, 1, 0)
head(heads_win)

plot(heads_win, col = "green")

# Use run length encoding to look into 0s and 1s runs

runs <- rle(data1)

results <- data.frame(
  group = runs$values,
  runs = runs$lengths
)

head(results)

# plot overlapping histos?

ggplot(
  data = results
) + geom_histogram(aes(x = runs, fill = as.factor(group)), alpha = 0.5)





