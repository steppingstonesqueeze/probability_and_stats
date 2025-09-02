# ---- Heavy-tail CLT demo with colourful viz (Pareto) ----
# Fixed version: avoids after_stat() inheritance issue in stat_function()

# Packages
pkgs <- c("ggplot2","dplyr","tidyr","purrr","tibble","patchwork","scales")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if(length(to_install)) install.packages(to_install, quiet = TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))

set.seed(42)

# Pareto generator (Type I), xm > 0, alpha > 0
rpareto <- function(n, alpha, xm = 1) xm / runif(n)^(1/alpha)

pareto_mean  <- function(alpha, xm=1) if (alpha > 1) alpha * xm / (alpha - 1) else NA_real_
pareto_var   <- function(alpha, xm=1) if (alpha > 2) alpha * xm^2 / ((alpha - 1)^2 * (alpha - 2)) else Inf

# Winsorize then trimmed-mean variance (for asymptotic SE of trimmed mean)
winsorize <- function(x, trim = 0.2){
  a <- quantile(x, probs = trim, type = 7)
  b <- quantile(x, probs = 1-trim, type = 7)
  pmin(pmax(x, a), b)
}
trimmed_mean_se <- function(x, trim = 0.2){
  # asymptotic Var(trimmed mean) ~ Var(Winsorized)/((1-2t)^2 * n)
  w <- winsorize(x, trim = trim)
  sw2 <- var(w)
  se <- sqrt(sw2) / ((1 - 2*trim) * sqrt(length(x)))
  se
}

# Simulation design
alphas <- c(2.5, 1.5)        # 2.5: CLT OK; 1.5: finite mean, infinite variance
Ns     <- c(20, 50, 100, 300, 1000, 3000)
reps   <- 4000
xm     <- 1

# ---------- Generate standardized sample mean draws ----------
# - For alpha > 2: Z_n = sqrt(n) * (xbar - mu) / sigma_true
# - For 1 < alpha <= 2: T_n = sqrt(n) * (xbar - mu) / s  (Studentized; still not N(0,1))

gen_stats <- function(alpha, n, reps){
  mu    <- pareto_mean(alpha, xm)
  s2    <- pareto_var(alpha, xm)
  Xbar  <- replicate(reps, mean(rpareto(n, alpha, xm)))
  s_hat <- replicate(reps, sd(  rpareto(n, alpha, xm)))
  out <- tibble(
    alpha = alpha, n = n,
    xbar = Xbar,
    s_hat = s_hat
  ) |>
    dplyr::mutate(
      mu = mu,
      sigma_true = sqrt(s2),
      stat_known = if (is.finite(s2)) sqrt(n) * (xbar - mu) / sigma_true else NA_real_,
      stat_t     = sqrt(n) * (xbar - mu) / s_hat
    )
  out
}

stats_df <- purrr::map_df(alphas, \(a) purrr::map_df(Ns, \(n) gen_stats(a, n, reps)))

# ---------- FIG A: Histograms of standardized means vs N(0,1) ----------
nice_levels <- c("alpha=2.5 (finite var)", "alpha=1.5 (∞ var)")
stats_plot <- stats_df |>
  dplyr::mutate(alpha_lab = factor(ifelse(alpha>2, nice_levels[1], nice_levels[2]),
                            levels = nice_levels),
         which_stat = ifelse(alpha>2, "known σ", "studentized s"),
         z = ifelse(alpha>2, stat_known, stat_t)) |>
  dplyr::filter(n %in% c(100, 2000)) # two sample sizes

# Use ..density.. for broad compatibility; add stat_function with inherit.aes = FALSE
p_hist <- ggplot2::ggplot(stats_plot, ggplot2::aes(x = z, y = ..density.., fill = factor(n))) +
  ggplot2::geom_histogram(bins = 60, alpha = 0.65, position = "identity") +
  ggplot2::stat_function(
    data = data.frame(x = c(-6, 6)),
    mapping = ggplot2::aes(x = x),
    fun = dnorm,
    linewidth = 1.1,
    linetype = 2,
    inherit.aes = FALSE
  ) +
  ggplot2::facet_grid(alpha_lab ~ which_stat, labeller = label_value) +
  ggplot2::scale_fill_brewer(palette = "Set2", name = "n") +
  ggplot2::labs(title = "CLT check via standardized sample means",
       subtitle = "Overlayed with standard normal density (dashed)",
       x = "standardized mean", y = "density") +
  ggplot2::theme_minimal(base_size = 12)

# ---------- FIG B: Q–Q plots vs Normal ----------
qq_samples <- stats_df |>
  dplyr::mutate(alpha_lab = factor(ifelse(alpha>2, nice_levels[1], nice_levels[2]),
                            levels = nice_levels),
         which_stat = ifelse(alpha>2, "known σ", "studentized s"),
         z = ifelse(alpha>2, stat_known, stat_t)) |>
  dplyr::filter(n %in% c(100, 2000)) |>
  dplyr::group_by(alpha_lab, which_stat, n) |>
  dplyr::mutate(q_theory = qnorm(ppoints(dplyr::n()), 0, 1),
         q_sample = sort(z)) |>
  dplyr::ungroup()

p_qq <- ggplot2::ggplot(qq_samples, ggplot2::aes(x = q_theory, y = q_sample, color = factor(n))) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2) +
  ggplot2::geom_point(alpha = 0.5, size = 0.7) +
  ggplot2::facet_grid(alpha_lab ~ which_stat) +
  ggplot2::scale_color_brewer(palette = "Set1", name = "n") +
  ggplot2::labs(title = "Q–Q against N(0,1): straight line = normality",
       x = "Theoretical normal quantiles", y = "Empirical quantiles") +
  ggplot2::theme_minimal(base_size = 12)

# ---------- FIG C: Coverage of nominal 95% CIs for the MEAN (plain t-interval) ----------
coverage_plain <- purrr::map_df(alphas, function(alpha){
  if (alpha <= 1) return(tibble())
  mu <- pareto_mean(alpha, xm)
  purrr::map_df(Ns, function(n){
    inside <- replicate(reps, {
      x <- rpareto(n, alpha, xm)
      m <- mean(x); s <- sd(x)
      se <- s / sqrt(n)
      tcrit <- qt(0.975, df = n - 1)
      (mu >= m - tcrit * se) && (mu <= m + tcrit * se)
    })
    tibble(alpha=alpha, n=n, coverage=mean(inside))
  })
})

p_cov_plain <- coverage_plain |>
  dplyr::mutate(alpha_lab = ifelse(alpha>2, nice_levels[1], nice_levels[2])) |>
  ggplot2::ggplot(ggplot2::aes(x = n, y = coverage, color = alpha_lab)) +
  ggplot2::geom_hline(yintercept = 0.95, linetype = 2) +
  ggplot2::geom_point(size = 2) + ggplot2::geom_line(linewidth = 1) +
  ggplot2::scale_x_log10(breaks = Ns) +
  ggplot2::scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
  ggplot2::scale_color_brewer(palette = "Dark2", name = "Tail regime") +
  ggplot2::labs(title = "Nominal 95% t-interval coverage for the mean",
       subtitle = "Heavy tails (α=1.5) under-cover even at large n",
       x = "n (log scale)", y = "empirical coverage") +
  ggplot2::theme_minimal(base_size = 12)

# ---------- FIG D: Robust 20% trimmed mean coverage vs plain mean ----------
coverage_trim <- purrr::map_df(alphas, function(alpha){
  if (alpha <= 1) return(tibble())
  mu <- pareto_mean(alpha, xm)
  purrr::map_df(Ns, function(n){
    # Plain t CI
    inside_plain <- replicate(reps, {
      x <- rpareto(n, alpha, xm)
      m <- mean(x); s <- sd(x)
      se <- s / sqrt(n)
      tcrit <- qt(0.975, df = n - 1)
      (mu >= m - tcrit * se) && (mu <= m + tcrit * se)
    })
    # 20% trimmed mean with Winsorized-variance SE (z CI)
    inside_trim <- replicate(reps, {
      x <- rpareto(n, alpha, xm)
      m <- mean(x, trim = 0.2)
      se <- trimmed_mean_se(x, trim = 0.2)
      zcrit <- qnorm(0.975)
      (mu >= m - zcrit * se) && (mu <= m + zcrit * se)
    })
    tibble(alpha=alpha, n=n,
           coverage_plain=mean(inside_plain),
           coverage_trim =mean(inside_trim))
  })
})

cov_long <- coverage_trim |>
  dplyr::mutate(alpha_lab = ifelse(alpha>2, nice_levels[1], nice_levels[2])) |>
  tidyr::pivot_longer(cols = dplyr::starts_with("coverage_"),
               names_to = "method", values_to = "coverage") |>
  dplyr::mutate(method = dplyr::recode(method,
                         coverage_plain = "Plain mean (t-CI)",
                         coverage_trim  = "20% trimmed (z-CI)"))

p_cov_comp <- ggplot2::ggplot(cov_long, ggplot2::aes(x = n, y = coverage, color = method))+
  ggplot2::geom_hline(yintercept = 0.95, linetype = 2) +
  ggplot2::geom_point(size = 2) + ggplot2::geom_line(linewidth = 1) +
  ggplot2::scale_x_log10(breaks = Ns) +
  ggplot2::scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
  ggplot2::scale_color_brewer(palette = "Set1", name = "") +
  ggplot2::facet_wrap(~ alpha_lab) +
  ggplot2::labs(title = "Coverage: robust trimmed mean vs plain mean",
       subtitle = "Trimming stabilizes inference under heavy tails (α=1.5)",
       x = "n (log scale)", y = "empirical coverage") +
  ggplot2::theme_minimal(base_size = 12)

# ---- Draw all figures ----
#(p_hist / p_qq) | (p_cov_plain / p_cov_comp)


(p_hist)
Sys.sleep(4)
(p_qq)
Sys.sleep(4)
(p_cov_plain)
Sys.sleep(4)
(p_cov_comp)
