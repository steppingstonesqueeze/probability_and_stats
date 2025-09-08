# ================================================================
# Generalized CLT & LLN Failures: Pareto and Lévy (R, with viz)
# Safe numerics: no qstable(); robust density overlays; fewer warnings
# ================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("ggplot2",   quietly = TRUE)) install.packages("ggplot2")
  if (!requireNamespace("dplyr",     quietly = TRUE)) install.packages("dplyr")
  if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
  if (!requireNamespace("stabledist",quietly = TRUE)) install.packages("stabledist")
  if (!requireNamespace("tidyr",     quietly = TRUE)) install.packages("tidyr")
})

library(ggplot2)
library(dplyr)
library(patchwork)
library(stabledist)
library(tidyr)

# Optional: use libstableR if installed (more robust/fast pdf)
.have_libstable <- requireNamespace("libstableR", quietly = TRUE)

# -------------------------
# Pareto generators
# -------------------------
rpareto <- function(n, alpha, xm = 1) {
  # Type I Pareto on (xm, inf)
  xm / runif(n)^(1/alpha)
}
rpareto_sym <- function(n, alpha, xm = 1) {
  # Symmetric Pareto: +/- with equal prob, magnitude Pareto(alpha)
  sign <- ifelse(runif(n) < 0.5, -1, 1)
  sign * rpareto(n, alpha, xm)
}
pareto_mean <- function(alpha, xm = 1) if (alpha > 1) alpha * xm / (alpha - 1) else Inf

# -------------------------
# Lévy (α = 1/2) via stable sampler
# -------------------------
rlevy <- function(n, c = 1, mu = 0) {
  # Lévy increments are α-stable with α=1/2, β=1
  rstable(n, alpha = 0.5, beta = 1, gamma = c, delta = mu, pm = 1)
}

# -------------------------
# Stable helpers (safe numerics)
# -------------------------
safe_dstable <- function(x, alpha, beta, gamma, delta = 0) {
  if (.have_libstable) {
    # libstableR uses S0 parametrization (parametrization = 0)
    libstableR::dstable(x,
                        alpha = alpha, beta = beta,
                        sigma = gamma, mu = delta,
                        parametrization = 0L)
  } else {
    # stabledist density; suppress tail warnings
    suppressWarnings(stabledist::dstable(x,
                                         alpha = alpha, beta = beta,
                                         gamma = gamma, delta = delta, pm = 1))
  }
}

# Fit stable scale γ via IQR ratio, using simulation (no qstable)
fit_stable_scale_by_IQR_sim <- function(sample, alpha, beta, B = 20000L, seed = 123) {
  set.seed(seed)
  y <- stabledist::rstable(B, alpha = alpha, beta = beta, gamma = 1, delta = 0, pm = 1)
  iqr_samp <- diff(quantile(sample, c(0.25, 0.75), names = FALSE))
  iqr_theo <- diff(quantile(y,      c(0.25, 0.75), names = FALSE))
  as.numeric(iqr_samp / iqr_theo)
}

# Build a "theoretical" stable sample for Q–Q via simulation
stable_theoretical_sample <- function(n, alpha, beta, gamma, delta = 0, seed = 456) {
  set.seed(seed)
  stabledist::rstable(n, alpha = alpha, beta = beta, gamma = gamma, delta = delta, pm = 1)
}

# -------------------------
# Generalized CLT helpers
# -------------------------
# For α > 1 (finite mean): (S_n - nμ)/n^{1/α}
# For α ≤ 1 (infinite mean): S_n / n^{1/α}
normalized_sums <- function(rfun, alpha, n, reps, mean_fun = NULL) {
  Tn <- numeric(reps)
  mu <- if (is.null(mean_fun)) NA_real_ else mean_fun(alpha)
  for (j in seq_len(reps)) {
    s <- sum(rfun(n, alpha))
    if (!is.null(mean_fun) && is.finite(mu)) Tn[j] <- (s - n * mu) / (n^(1/alpha))
    else                                     Tn[j] <-  s / (n^(1/alpha))
  }
  Tn
}

# -------------------------
# Visualizations
# -------------------------

# 1) CLT failure under sqrt(n) scaling (Pareto α in (1,2))
plot_clt_failure <- function(alpha = 1.5, n = 2000, reps = 4000, one_sided = TRUE) {
  stopifnot(alpha > 1, alpha < 2)
  xm <- 1
  mu <- pareto_mean(alpha, xm = xm)
  Z <- replicate(reps, {
    x <- if (one_sided) rpareto(n, alpha, xm) else rpareto_sym(n, alpha, xm)
    m <- mean(x); s <- sd(x)
    sqrt(n) * (m - mu) / s
  })
  ggplot(data.frame(Z = Z), aes(Z)) +
    geom_histogram(aes(y = after_stat(density)), bins = 80, alpha = 0.7, fill = "#e67e22") +
    stat_function(fun = dnorm, linetype = 2, linewidth = 1) +
    labs(title = sprintf("CLT failure for Pareto(α=%.1f): √n scaling not Normal", alpha),
         subtitle = sprintf("n=%d, reps=%d — studentized means are heavy-tailed", n, reps),
         x = "√n (x̄ - μ) / s", y = "density") +
    theme_minimal(base_size = 12)
}

# 2) Generalized CLT: stable limit under n^{1/α} scaling
plot_generalized_clt <- function(alpha = 1.5, n = 5000, reps = 4000, one_sided = TRUE) {
  stopifnot(alpha > 1, alpha < 2)
  xm <- 1
  mu <- if (one_sided) pareto_mean(alpha, xm) else 0
  rfun <- if (one_sided) function(n, a) rpareto(n, a, xm)
  else            function(n, a) rpareto_sym(n, a, xm)
  
  Tn <- normalized_sums(rfun, alpha, n, reps, mean_fun = function(a) mu)
  beta <- if (one_sided) 1 else 0
  gamma_hat <- fit_stable_scale_by_IQR_sim(Tn, alpha = alpha, beta = beta)
  
  # Central overlay to avoid tail-integration nasties
  rng   <- quantile(Tn, c(0.02, 0.98), names = FALSE)
  xgrid <- seq(rng[1], rng[2], length.out = 400)
  dens  <- safe_dstable(xgrid, alpha = alpha, beta = beta, gamma = gamma_hat, delta = 0)
  
  p_hist <- ggplot(data.frame(Tn = Tn), aes(Tn)) +
    geom_histogram(aes(y = after_stat(density)), bins = 100, alpha = 0.75, fill = "#1f77b4") +
    geom_line(data = data.frame(x = xgrid, y = dens),
              aes(x = x, y = y), inherit.aes = FALSE,
              color = "#2ca02c", linewidth = 1.1) +
    labs(title = sprintf("Generalized CLT: (Sₙ - nμ)/n^{1/α} → α-stable (α=%.1f, β=%s)",
                         alpha, ifelse(beta==1,"+1","0")),
         subtitle = sprintf("Overlay via %s; γ ≈ %.3f",
                            if (.have_libstable) "libstableR" else "stabledist", gamma_hat),
         x = "stable-scaled sum", y = "density") +
    theme_minimal(base_size = 7)
  
  # Q–Q vs simulated α-stable (no qstable)
  y_theo <- stable_theoretical_sample(length(Tn), alpha, beta, gamma_hat, 0)
  dfqq <- data.frame(sample = sort(Tn), theo = sort(y_theo))
  p_qq <- ggplot(dfqq, aes(theo, sample)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    geom_point(alpha = 0.5, size = 0.8) +
    labs(title = "Q–Q vs simulated α-stable",
         x = "Theoretical α-stable quantiles", y = "Sample quantiles") +
    theme_minimal(base_size = 8)
  
  p_hist | p_qq
}

# 3) LLN failure for infinite-mean: running mean paths
plot_lln_failure_paths <- function(alpha = 0.8, n_path = 30000, n_paths = 8) {
  stopifnot(alpha <= 1)
  xm <- 1
  paths <- replicate(n_paths, {
    x <- rpareto(n_path, alpha, xm)
    cumsum(x) / seq_along(x)
  })
  df <- as.data.frame(paths)
  df$step <- seq_len(n_path)
  df_long <- tidyr::pivot_longer(df, -step, names_to = "path", values_to = "running_mean")
  ggplot(df_long, aes(step, running_mean, color = path)) +
    geom_line(alpha = 0.7) +
    scale_x_log10() +
    guides(color = "none") +
    labs(title = sprintf("LLN failure: running means for Pareto(α=%.1f) with infinite mean", alpha),
         subtitle = "Multiple paths, log-x axis; means do not stabilize",
         x = "n (log scale)", y = "running mean") +
    theme_minimal(base_size = 12)
}

# 4) Lévy sums: stable under n^{1/α} (α=1/2) and LLN failure
plot_levy_demo <- function(n = 2000, reps = 5000, c = 1) {
  alpha <- 0.5; beta <- 1
  # Simulate sums of Lévy increments, scale by n^{1/α}
  Tn <- replicate(reps, {
    x <- stabledist::rstable(n, alpha = alpha, beta = beta, gamma = c, delta = 0, pm = 1)
    sum(x) / (n^(1/alpha))
  })
  gamma_hat <- fit_stable_scale_by_IQR_sim(Tn, alpha = alpha, beta = beta)
  
  rng   <- quantile(Tn, c(0.02, 0.98), names = FALSE)
  xgrid <- seq(rng[1], rng[2], length.out = 400)
  dens  <- safe_dstable(xgrid, alpha = alpha, beta = beta, gamma = gamma_hat, delta = 0)
  
  p1 <- ggplot(data.frame(Tn = Tn), aes(Tn)) +
    geom_histogram(aes(y = after_stat(density)), bins = 100, alpha = 0.75, fill = "#d62728") +
    geom_line(data = data.frame(x = xgrid, y = dens),
              aes(x = x, y = y), inherit.aes = FALSE,
              color = "#2ca02c", linewidth = 1.1) +
    labs(title = "Lévy sums: Sₙ / n^{1/α} with α=1/2 → α-stable",
         subtitle = sprintf("Overlay via %s; γ ≈ %.3f",
                            if (.have_libstable) "libstableR" else "stabledist", gamma_hat),
         x = "stable-scaled sum", y = "density") +
    theme_minimal(base_size = 7)
  
  # Running mean path to show LLN failure (mean = ∞)
  x <- stabledist::rstable(30000, alpha = alpha, beta = beta, gamma = c, delta = 0, pm = 1)
  runmean <- cumsum(x) / seq_along(x)
  p2 <- ggplot(data.frame(step = seq_along(runmean), running_mean = runmean), aes(step, running_mean)) +
    geom_line(color = "#d62728") +
    scale_x_log10() +
    labs(title = "Running mean of Lévy(α=1/2): LLN fails (mean = ∞)",
         x = "n (log scale)", y = "running mean") +
    theme_minimal(base_size = 8)
  
  p1 | p2
}

# ------------------------------------------------
# EXAMPLES — run these blocks interactively
# ------------------------------------------------
  # (A) CLT failure (Pareto α=1.5) and generalized CLT overlay
  plot_clt_failure(alpha = 1.5, n = 2000, reps = 4000, one_sided = TRUE)
Sys.sleep(3)
  print(plot_generalized_clt(alpha = 1.5, n = 5000, reps = 4000, one_sided = TRUE))
  Sys.sleep(3)
  # (B) LLN failure: Pareto α=0.8 (infinite mean)
  print(plot_lln_failure_paths(alpha = 0.8, n_path = 30000, n_paths = 8))
  Sys.sleep(3)
  # (C) Lévy sums demo (α=1/2)
  print(plot_levy_demo(n = 2000, reps = 5000, c = 1))

