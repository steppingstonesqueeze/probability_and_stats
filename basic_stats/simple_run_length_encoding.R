x <- c(1,1,1,1,0,0,0,1,1,1)
runs <- rle(x)

result <- data.frame(
  group = runs$values,
  runs = runs$lengths
)

result