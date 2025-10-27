counts

ggplot(
  data = counts, aes(x = x, y = density)
) + geom_point(colour = "black") + geom_smooth(se = FALSE, colour = "green")


dice1 <- sample(c(1:6), 1000, replace = TRUE)
dice2 <- sample(c(1:6), 1000, replace = TRUE)

score <- sum(ifelse(dice1 + dice2 > 9, 1, 0))

score / 1000.0

x <- c(99, 77, 67, 70)
t.test(x, mu = 50, alternative = "greater")

(mean(x) - 50) / (sd(x) / sqrt(length(x)))

p_vec <- rep(1./6., 6)
p_vec

p_act <- c(0.12, 0.18, 0.2, 0.17, 0.19, 0.14)

sum(p_act)

t.test(x = p_vec,
       y = p_act)

p_act_2 <- c(0.85, 0.15, 0, 0, 0, 0)

t.test(x = p_vec,
       y = p_act_2)

t.test(p_vec - p_act)

pop <- rnorm(10000, 0, 3)

hist(pop)

#v sample from this and look at the samplong distribution of sample variances 

sampling_variance <- numeric(length = 1000)
for (i in 1:1000) {
  s <- sample(pop, 100, replace = TRUE)
  sampling_variance[i] <- sd(s)
}

hist(sampling_variance)



#v sampling distribution of Linear regression lines

# setup
# original data x_i, y_i i = 1,2,...N ; all x_is unique

# Expt 1 : M < N
# Sample points withour replacement from orig data - fit and get slope and intercept
# plot the histos


x <- seq(0, 10, by = 0.1)
x

num_pts <- length(x)

actual_slope <- 0.5
actual_intercept <- -0.3
noise_amplitude <- 0.85


y <- actual_slope * x + actual_intercept + noise_amplitude * runif(num_pts, -5, 5)

plot(x,y, col = "red")

cor(x,y)

# Sampling distribution of the slope via sampling from orig data pts

orig_lm <- lm(y ~ x)
attributes(orig_lm)

orig_intercept <- orig_lm$coefficients[1]
orig_slope <- orig_lm$coefficients[2]

num_sampled_pts <- 40

num_expts <- 1000

sample_slope <- numeric(length = num_expts)
sample_intercept <- numeric(length = num_expts)

# df[sample(nrow(df), 3), ]

df <- data.frame(x = x, y = y)
nrows_df <- nrow(df)

for (e in 1:num_expts) {
  
  data_sample <- df[sample(nrows_df, num_sampled_pts, replace = TRUE), ]
  
  x_s <- data_sample[,1]
  y_s <- data_sample[,2]
  
  # fit it
  lm1 <- lm(y_s ~ x_s)
  
  sample_intercept[e] <- lm1$coefficients[1]
  sample_slope[e] <- lm1$coefficients[2]
  
}

par(mfrow = c(2,1))

hist(sample_intercept, col = "blue", main = "sampling distribution of the intercept")
hist(sample_slope, col = "red", main = "sampling distribution of the slope")


### Noise dependence of good-ness



x <- seq(0, 10, by = 0.1)
x

num_pts <- length(x)

actual_slope <- 0.5
actual_intercept <- -0.3
noise_amplitude <- 2.3


y <- actual_slope * x + actual_intercept + noise_amplitude * runif(num_pts, -1, 1)

plot(x,y, col = "red")

cor(x,y)

# Sampling distribution of the slope via sampling from orig data pts

orig_lm <- lm(y ~ x)
attributes(orig_lm)

orig_intercept <- orig_lm$coefficients[1]
orig_slope <- orig_lm$coefficients[2]

num_sampled_pts <- 480

num_expts <- 1000

sample_slope <- numeric(length = num_expts)
sample_intercept <- numeric(length = num_expts)

# df[sample(nrow(df), 3), ]

df <- data.frame(x = x, y = y)
nrows_df <- nrow(df)

for (e in 1:num_expts) {
  
  data_sample <- df[sample(nrows_df, num_sampled_pts, replace = TRUE), ]
  
  x_s <- data_sample[,1]
  y_s <- data_sample[,2]
  
  # fit it
  lm1 <- lm(y_s ~ x_s)
  
  sample_intercept[e] <- lm1$coefficients[1]
  sample_slope[e] <- lm1$coefficients[2]
  
}

par(mfrow = c(2,1))

hist(sample_intercept, col = "blue", main = "sampling distribution of the intercept")
hist(sample_slope, col = "red", main = "sampling distribution of the slope")

mean(sample_intercept)
mean(sample_slope)
sd(sample_intercept)
sd(sample_slope)
