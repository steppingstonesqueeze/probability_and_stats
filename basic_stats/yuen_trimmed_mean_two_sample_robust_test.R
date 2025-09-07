# Yuen's trimmed mean two sample test -- robust to outliers #

.yuen_core <- function(x, y, trim = 0.2, alternative = c("two.sided","less","greater")) {
  alternative <- match.arg(alternative)
  stopifnot(trim >= 0, trim < 0.5)
  tx <- sort(x); ty <- sort(y)
  nx <- length(x); ny <- length(y)
  gx <- floor(trim * nx); gy <- floor(trim * ny)
  hx <- nx - 2 * gx;       hy <- ny - 2 * gy
  if (hx < 2 || hy < 2) stop("Not enough data after trimming; reduce 'trim'.")
  
  # trimmed means
  x_t <- mean(tx[(gx+1):(nx-gx)])
  y_t <- mean(ty[(gy+1):(ny-gy)])
  
  # winsorize
  wx <- x; if (gx > 0) {
    lo <- tx[gx+1]; hi <- tx[nx-gx]
    wx[wx < lo] <- lo; wx[wx > hi] <- hi
  }
  wy <- y; if (gy > 0) {
    lo <- ty[gy+1]; hi <- ty[ny-gy]
    wy[wy < lo] <- lo; wy[wy > hi] <- hi
  }
  
  # winsorized variances
  swx2 <- var(wx); swy2 <- var(wy)
  
  # variance of trimmed means (Yuen’s formula)
  vx <- swx2 / (hx * (hx - 1))
  vy <- swy2 / (hy * (hy - 1))
  
  tval <- (x_t - y_t) / sqrt(vx + vy)
  
  # Welch-type df for Yuen
  df <- (vx + vy)^2 / (vx^2 / (hx - 1) + vy^2 / (hy - 1))
  
  p <- switch(alternative,
              "two.sided" = 2 * pt(-abs(tval), df),
              "greater"   = 1 - pt(tval, df),
              "less"      = pt(tval, df))
  
  list(statistic = tval, df = df, p.value = p,
       estimate = c(trimmed_mean_x = x_t, trimmed_mean_y = y_t),
       trim = trim, alternative = alternative)
}

yuen_test <- function(x, y, trim = 0.2, alternative = c("two.sided","less","greater")) {
  .yuen_core(x, y, trim = trim, alternative = alternative)
}

simulate_yuen <- function(
    rx, ry, m = 40, n = 40, N = 2000, alpha = 0.05,
    trim = 0.2, alternative = c("two.sided","less","greater"), seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  alternative <- match.arg(alternative)
  rej <- logical(N)
  for (i in seq_len(N)) {
    x <- rx(m); y <- ry(n)
    p <- yuen_test(x, y, trim = trim, alternative = alternative)$p.value
    rej[i] <- (p < alpha)
  }
  list(alpha = alpha, N = N, true_sig_level = mean(rej), reject = rej,
       trim = trim, alternative = alternative)
}

# demos
set.seed(1)
# Heavy-tailed Cauchy (null) — Yuen keeps size near alpha better than classic t
simulate_yuen(function(n) rcauchy(n), function(n) rcauchy(n),
              m=50,n=50,N=1000,alpha=0.05,trim=0.2,alternative="two.sided")$true_sig_level

# Outlier contamination — Yuen often has better power than Welch here
simulate_yuen(function(n) c(rnorm(n-2,0,1), 20, -20),
              function(n) rnorm(n, 0.3, 1),
              m=40,n=40,N=1000,alpha=0.05,trim=0.2,alternative="greater")$true_sig_level
