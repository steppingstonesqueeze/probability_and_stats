# breiman inspired #
# Toss a coin N times - N even ; if exactly same number of heads and tails, victory 
# Plot number of times we get exact H,T against N

pi <- 3.141592654
N_seq <- c(100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000)
N_seq_2 <- N_seq / 2

len_N_seq <- length(N_seq)

N_trials <- 5000

sample_space <- c(0, 1)

same_ratio <- numeric(length = len_N_seq)

for (i in 1:len_N_seq) {
  same <- 0
  for (t in 1:N_trials) {
    data <- sample(sample_space, N_seq[i], replace = TRUE)
    if (sum(data) == N_seq_2[i]) {
      same <- same + 1
    }
  }
  
  same_ratio[i] <- same / N_trials
  
}

# This is trivial to compute using Stirling formula for n! 
# Assume 2n trials and then simply compute 2n Choose n - simplifying gives the below

theoretical_same_ratio <- sqrt(1.0 / (pi * N_seq_2))

comparison_df <- data.frame(N = N_seq, 
                            theoretical_same_ratio = theoretical_same_ratio,
                            same_ratio = same_ratio)



ggplot(
  data = comparison_df
) + geom_point(aes(x = N, y = theoretical_same_ratio), colour = "red") + 
  geom_smooth(aes(x = N, y = same_ratio), colour = "blue", se = FALSE, method = "loess") + 
  geom_line(aes(x = N, y = same_ratio), colour = "black")


