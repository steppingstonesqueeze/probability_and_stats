#t-tests on streoids!


# --- Test statistics ---
tstat_pooled <- function(x, y) {
  m <- length(x); n <- length(y)
  sp <- sqrt(((m-1)*var(x) + (n-1)*var(y)) / (m+n-2))
  (mean(x) - mean(y)) / (sp * sqrt(1/m + 1/n))
}

tstat_welch <- function(x, y) {
  m <- length(x); n <- length(y)
  (mean(x) - mean(y)) / sqrt(var(x)/m + var(y)/n)
}

df_welch <- function(x, y) {
  m <- length(x); n <- length(y)
  vx <- var(x); vy <- var(y)
  num <- (vx/m + vy/n)^2
  den <- vx^2 / (m^2 * (m - 1)) + vy^2 / (n^2 * (n - 1))
  num / den
}

# --- General simulator ---
simulate_t_test <- function(
    rx, ry,                       # functions: n -> numeric vector
    m = 100, n = 100,             # sample sizes
    N = 10000,                    # simulations
    alpha = 0.10,
    alternative = c("two.sided", "less", "greater"),
    test = c("pooled", "welch")   # var.equal=TRUE vs Welch
) {
  alternative <- match.arg(alternative)
  test <- match.arg(test)
  
  sims <- replicate(N, {
    x <- rx(m); y <- ry(n)
    if (test == "pooled") {
      t <- tstat_pooled(x, y); df <- m + n - 2
    } else {
      t <- tstat_welch(x, y);  df <- df_welch(x, y)
    }
    c(t = t, df = df)
  })
  
  tvec <- sims["t", ]; dfvec <- sims["df", ]
  
  reject <- switch(
    alternative,
    "two.sided" = abs(tvec) > qt(1 - alpha/2, dfvec),
    "greater"   = tvec > qt(1 - alpha, dfvec),
    "less"      = tvec < qt(alpha, dfvec)
  )
  
  list(
    alpha = alpha, m = m, n = n, N = N,
    alternative = alternative, test = test,
    true_sig_level = mean(reject),
    t = tvec, df = dfvec
  )
}

# Usage

# Example 1
set.seed(1)
res <- simulate_t_test(
  rx = function(n) rnorm(n, 0, 1),
  ry = function(n) rnorm(n, 0, 1),
  m = 100, n = 100, N = 10000, alpha = 0.10,
  alternative = "two.sided", test = "pooled"
)
res$true_sig_level  # ~0.10 under the null


# Example 2 - non-normal but equal means

set.seed(1)
res_gamma_t <- simulate_t_test(
  rx = function(n) rgamma(n, shape = 4, scale = 2),     # mean = 8
  ry = function(n) rgamma(n, shape = 4, scale = 2),     # same mean
  m = 50, n = 70, N = 20000, alpha = 0.05,
  alternative = "two.sided", test = "welch"
)
res_gamma_t$true_sig_level


# Example 3 - unequal means
set.seed(1)
res_power <- simulate_t_test(
  rx = function(n) rnorm(n, 0, 1),
  ry = function(n) rnorm(n, 0.3, 1),   # mean shift of +0.3
  m = 40, n = 40, N = 20000, alpha = 0.05,
  alternative = "greater", test = "welch"
)
res_power$true_sig_level  # this is empirical power here








