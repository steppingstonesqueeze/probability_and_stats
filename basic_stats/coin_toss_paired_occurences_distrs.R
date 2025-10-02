# coin toss - paired occurences distribution and a complication on the
# geometric distribution

library(ggplot2)
EPS <- 1.0e-10
p_H <- 0.7
p_T <- 1.0 - p_H
N <- 100000
N1 <- N - 1

options <- c("H", "T")

data <- sample(options,
               N,
               replace = TRUE,
               prob = c(p_H, p_T))

# code : 1 = "HH", 2 = "HT", 3 = "TH", 4 = "TT" ; 4 choices for the pairs

occurences <- numeric(length = 4)

paired_H_occurences <- NULL
paired_T_occurences <- NULL

for ( i in 1:N1) {
  j <- i + 1
  
  if (data[i] == "H" && data[j] == "H") {
    occurences[1] <- occurences[1] + 1
    # note position of first H down
    paired_H_occurences <- c(paired_H_occurences, i)
  } else if (data[i] == "H" && data[j] == "T") {
    occurences[2] <- occurences[2] + 1
  } else if (data[i] == "T" && data[j] == "H") {
    occurences[3] <- occurences[3] + 1
  } else {
    occurences[4] <- occurences[4] + 1
    # note position of first T down
    paired_T_occurences <- c(paired_T_occurences, i)
  }
}

occ_df <- data.frame(index = c(1:4), 
                     pairs = c("HH", "HT", "TH", "TT"),
                     occurences = occurences,
                     stringsAsFactors = FALSE)

ggplot(
  data = occ_df, aes(x = pairs, y = occurences)
) + geom_bar(stat = "identity", fill = factor(occ_df$index))

# distances between the starting point of successive pairs #

diff_paired_H_occurences <- diff(paired_H_occurences)
diff_paired_T_occurences <- diff(paired_T_occurences)

# plot their histos
par(mfrow=c(1,2))
hist(diff_paired_T_occurences, col = "red", main = "Histo of diffs paired T", xlab = "Diffs paired T")
hist(diff_paired_H_occurences, col = "blue", main = "Histo of diffs paired H", xlab = "Diffs paired H")


# log plots
par(mfrow=c(1,2))
hist(log(diff_paired_T_occurences+EPS), col = "red", main = "Histo of log diffs paired T", xlab = "Log Diffs paired T")
hist(log(diff_paired_H_occurences+EPS), col = "blue", main = "Histo of log diffs paired H", xlab = "Log Diffs paired H")



