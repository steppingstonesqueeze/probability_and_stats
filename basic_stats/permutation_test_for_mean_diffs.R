# robust permutation terst for mean differences

perm_test_once <- function(x, y, B = 1999, alternative = c("two.sided","less","greater")) {
  alternative <- match.arg(alternative)
  nx <- length(x); ny <- length(y)
  d_obs <- mean(x) - mean(y)
  z <- c(x, y)
  
  # generate permutation diffs
  diffs <- replicate(B, {
    idx <- sample.int(nx + ny, nx, replace = FALSE)
    mean(z[idx]) - mean(z[-idx])
  })
  
  if (alternative == "two.sided") {
    p <- mean(abs(diffs) >= abs(d_obs))
  } else if (alternative == "greater") {
    p <- mean(diffs >= d_obs)
  } else {
    p <- mean(diffs <= d_obs)
  }
  list(stat = d_obs, p.value = p)
}

simulate_perm_test <- function(
    rx, ry, m = 50, n = 50, N = 2000, alpha = 0.05,
    alternative = c("two.sided","less","greater"),
    B = 999, seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  alternative <- match.arg(alternative)
  rej <- logical(N)
  for (i in seq_len(N)) {
    x <- rx(m); y <- ry(n)
    p <- perm_test_once(x, y, B = B, alternative = alternative)$p.value
    rej[i] <- (p < alpha)
  }
  list(alpha = alpha, N = N, true_sig_level = mean(rej), reject = rej)
}

# quick demos
set.seed(1)
# Null normal
simulate_perm_test(function(n) rnorm(n), function(n) rnorm(n),
                   m=50,n=60,N=1000,alpha=0.05,alternative="two.sided",B=999)$true_sig_level

# Mean shift normal (greater)
simulate_perm_test(function(n) rnorm(n,0,1), function(n) rnorm(n,0.3,1),
                   m=40,n=40,N=1000,alpha=0.05,alternative="greater",B=999)$true_sig_level

# Exponential rates: rate_x=2, rate_y=1 ⇒ mean_x<mean_y ⇒ use "less"
simulate_perm_test(function(n) rexp(n, rate=2), function(n) rexp(n, rate=1),
                   m=60,n=60,N=1000,alpha=0.05,alternative="less",B=999)$true_sig_level
