

# Coin toss - histo of distance between successsive heads

sample_data <- c("H", "T")
N <- 10000
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
h <- hist(diff_pos_H, col = "red", breaks = 10)

#Histo details

cat("Histogram counts \n")
h$counts
cat("Histogram breaks \n")
h$breaks
cat("Histogram midpoints \n")
h$mids
cat("Histogram density \n")
h$density
