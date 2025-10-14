# type II error calculations in R #

set.seed(123)
n <- 20; mu0 <- 50; mu1 <- 55; sd <- 10
alpha <- 0.05
nsim <- 1e5

rejects <- replicate(nsim, {
  x <- rnorm(n, mean = mu1, sd = sd)
  t.test(x, mu = mu0)$p.value < alpha
})

power <- mean(rejects)
beta <- 1 - power
power; beta
