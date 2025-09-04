# t tests on stroeids with power curves and viz p 2

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
    rx, ry,                    # functions: n -> numeric vector
    m = 100, n = 100,          # sample sizes
    N = 10000,                 # simulations
    alpha = 0.10,
    alternative = c("two.sided","less","greater"),
    test = c("pooled","welch"),
    return_p = TRUE,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
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
  
  crit_fun <- function(tv, dfv) {
    switch(alternative,
           "two.sided" = abs(tv) > qt(1 - alpha/2, dfv),
           "greater"   = tv > qt(1 - alpha, dfv),
           "less"      = tv < qt(alpha, dfv)
    )
  }
  reject <- crit_fun(tvec, dfvec)
  
  pvals <- if (return_p) {
    if (alternative == "two.sided") 2 * pt(-abs(tvec), dfvec)
    else if (alternative == "greater") 1 - pt(tvec, dfvec)
    else pt(tvec, dfvec) # "less"
  } else NULL
  
  list(
    alpha = alpha, m = m, n = n, N = N,
    alternative = alternative, test = test,
    true_sig_level = mean(reject),
    reject = reject, t = tvec, df = dfvec, p = pvals
  )
}

# --- Power curve helpers ---

# Power vs mean shift for Normals (Welch by default)
power_curve_normal <- function(deltas, m=50, n=50, N=5000, alpha=0.05,
                               alternative=c("two.sided","less","greater"),
                               test=c("welch","pooled"), seed=1) {
  alternative <- match.arg(alternative); test <- match.arg(test)
  out <- sapply(deltas, function(d) {
    res <- simulate_t_test(
      rx = function(nn) rnorm(nn, 0, 1),
      ry = function(nn) rnorm(nn, d, 1),
      m=m, n=n, N=N, alpha=alpha, alternative=alternative, test=test, seed=seed
    )
    res$true_sig_level
  })
  data.frame(delta=deltas, power=out)
}

# Power vs rate ratio for exponentials.
# rate_y is baseline; rate_x = rate_ratio * rate_y
power_curve_exp <- function(rate_y=1, rate_ratios=seq(0.25, 2, by=0.25),
                            m=50, n=50, N=5000, alpha=0.05,
                            alternative=c("two.sided","less","greater"),
                            test=c("welch","pooled"), seed=1) {
  alternative <- match.arg(alternative); test <- match.arg(test)
  out <- sapply(rate_ratios, function(rR) {
    res <- simulate_t_test(
      rx = function(nn) rexp(nn, rate = rR * rate_y),
      ry = function(nn) rexp(nn, rate = rate_y),
      m=m, n=n, N=N, alpha=alpha, alternative=alternative, test=test, seed=seed
    )
    res$true_sig_level
  })
  data.frame(rate_ratio = rate_ratios, power = out)
}


# Example 1: Null, normal, pooled
set.seed(1)
res <- simulate_t_test(
  rx=function(n) rnorm(n,0,1),
  ry=function(n) rnorm(n,0,1),
  m=100,n=100,N=10000,alpha=0.10,
  alternative="two.sided", test="pooled"
)
res$true_sig_level  # ~ 0.10

# Example 2: Null but non-normal (Gamma), Welch
set.seed(1)
res_gamma_t <- simulate_t_test(
  rx=function(n) rgamma(n, shape=4, scale=2),
  ry=function(n) rgamma(n, shape=4, scale=2),
  m=50,n=70,N=20000,alpha=0.05,
  alternative="two.sided", test="welch"
)
res_gamma_t$true_sig_level  # ~ 0.05

# Example 3: Power for normal mean shift (greater)
set.seed(1)
res_power <- simulate_t_test(
  rx=function(n) rnorm(n,0,1),
  ry=function(n) rnorm(n,0.3,1),
  m=40,n=40,N=20000,alpha=0.05,
  alternative="greater", test="welch"
)
res_power$true_sig_level  # empirical power

# Exponential, identical rates => ~ alpha
set.seed(1)
res_exp_null <- simulate_t_test(
  rx=function(n) rexp(n, 1),
  ry=function(n) rexp(n, 1),
  m=4000,n=4000,N=20000,alpha=0.05,
  alternative="greater", test="welch"
)
res_exp_null$true_sig_level  # ~ 0.05

# Exponential with rate_x = 1, rate_y = 0.01 (so mean_x = 1, mean_y = 100)
# Here mean(x) < mean(y), so use alternative = "less"
set.seed(1)
res_exp_diff <- simulate_t_test(
  rx=function(n) rexp(n, rate=1),
  ry=function(n) rexp(n, rate=0.01),
  m=4000,n=4000,N=20000,alpha=0.05,
  alternative="less", test="welch"
)
res_exp_diff$true_sig_level  # ~ 1 with such huge difference


# Normal mean-shift power curve
pc_norm <- power_curve_normal(seq(0, 0.5, by=0.05), m=40, n=40, N=4000,
                              alpha=0.05, alternative="greater", test="welch", seed=1)
head(pc_norm)

# Exponential, power vs rate ratio; if rate_ratio > 1 then mean_x < mean_y ⇒ use "less"
pc_exp <- power_curve_exp(rate_y=1, rate_ratios=seq(0.25, 2, by=0.25),
                          m=60, n=60, N=4000, alpha=0.05,
                          alternative="less", test="welch", seed=1)
head(pc_exp)
