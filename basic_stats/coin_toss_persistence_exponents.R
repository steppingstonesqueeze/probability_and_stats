# coin toss persistence exponents

library(ggplot2)
library(dplyr)
set.seed(1234)
EPS <- 1.0e-8 # prevent infinities #
# Coin toss - histo of distance between successsive heads

sample_data <- c("H", "T")
N <- 100000
data <- sample(sample_data,
               N,
               replace = TRUE)

data1 <- ifelse(data == "H", 1, 0)

runs <- rle(data1)

results <- data.frame(
  group = runs$values,
  runs = runs$lengths
)

head(results)

#Histo 

ggplot(
  data = results
) + geom_histogram(aes(x = runs, fill = as.factor(group)), alpha = 0.5, bins = 100) + 
  ggtitle("Runs of 0s and 1s in coin tossing")

# lets now examine persistence

# heads fractions > 0.5 

heads_fraction <- cumsum(data1) / c(1:length(data1))
head(data1)
head(heads_fraction)

# heads_win defines when heads fraction > 0.5 ; then we look at ITS runs below
heads_win <- ifelse(heads_fraction > 0.5, 1, 0)
head(heads_win)


h_df <- data.frame(
  index = c(1:length(data1)),
  hf = heads_fraction,
  hw = heads_win
)

ggplot(
  data = h_df
) + geom_point(aes(x = index, y = hf), colour = "red") + 
  ggtitle("heads fraction for fair coin tossing")


# plot heads win
ggplot(
  data = h_df
) + geom_point(aes(x = index, y = hw), colour = "green") + 
  ggtitle("heads fraction > 0.5 for fair coin tossing is a 1, otherwise 0")

# Use run length encoding to look into 0s and 1s runs for the heads win data

runs <- rle(heads_win)

results <- data.frame(
  group = runs$values,
  runs = runs$lengths
)

head(results)

# plot overlapping histos?

ggplot(
  data = results
) + geom_histogram(aes(x = runs, fill = as.factor(group)), alpha = 0.5, bins = 100)

# Behooves us to do a log plot on this

ggplot(
  data = results
) + geom_histogram(aes(x = log(runs+EPS), fill = as.factor(group)), alpha = 0.5, bins = 100)


