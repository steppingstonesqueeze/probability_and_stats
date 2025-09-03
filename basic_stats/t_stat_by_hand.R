# t-statistics by hand -- inpsiired by Albert


x <- rnorm(10,mean=50,sd=10)
y <- rnorm(10,mean=50,sd=10)

tstatistic <- function(x,y) {
  m <- length(x)
  n <- length(y)
  pooled_sd <- sqrt(((m-1)*sd(x)^2+(n-1)*sd(y)^2)/(m+n-2))
  t <- (mean(x)-mean(y))/(pooled_sd*sqrt(1/m+1/n))
  return(t)
}

# compute on some data 

alpha <- .1; m <- 100; n <- 100 # sets alpha, m, n
N <- 10000 # sets the number of simulations
n.reject <- 0 # counter of num. of rejections
for (i in 1:N)
{
  x <- rnorm(m,mean=0,sd=1) # simulates xs from population 1
  y <- rnorm(n,mean=0,sd=1) # simulates ys from population 2
  t <- tstatistic(x,y) # computes the t statistic
  if (abs(t)>qt(1-alpha/2,n+m-2))
    n.reject <- n.reject+1 # reject if |t| exceeds critical pt
}
true.sig.level <- n.reject/N #  

true.sig.level # should be at alpha ideally #

# To run a very large Monte Carlo sim - use replicate like so #
N <- 100000
sim_t <- replicate(N, {
  x <- rnorm(m, mean = 0, sd = 1)
  y <- rnorm(n, mean = 0, sd = 1)
  tstatistic(x, y)
})

true.sig.level <- mean(abs(sim_t) > qt(1 - alpha/2, m + n - 2))
true.sig.level






