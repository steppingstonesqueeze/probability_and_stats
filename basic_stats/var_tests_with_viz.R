# ================================================================
# Variance-Equality Tests from Scratch + Clean Visualizations (R)
#   - Two-sample F-test (normality-sensitive)
#   - Levene's test (center = mean)
#   - Brown–Forsythe test (center = median)
#   - Pretty ggplot diagnostics (violin + jitter, abs-dev boxplots, QQ)
#
# Usage examples are at the bottom of this file.
# ================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
  if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
  if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")
  if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
  if (!requireNamespace("purrr", quietly = TRUE)) install.packages("purrr")
})

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(purrr)

# ------------------------------------------------
# Helper: sample variance with (n-1) denominator
# ------------------------------------------------
sample_var <- function(x) {
  n <- length(x)
  if (n < 2) stop("Need at least 2 observations per group.")
  mean((x - mean(x))^2) * n/(n-1)
}

# ------------------------------------------------
# (1) Classical two-sample F-test for equality of variances
#     Assumes both groups are normally distributed.
# ------------------------------------------------
f_test_equal_var <- function(x, y, alternative = c("two.sided","greater","less"), conf.level = 0.95) {
  alternative <- match.arg(alternative)
  n1 <- length(x); n2 <- length(y)
  if (n1 < 2 || n2 < 2) stop("Each sample must have at least 2 obs.")
  s1 <- sample_var(x); s2 <- sample_var(y)
  Fstat <- s1/s2
  df1 <- n1 - 1; df2 <- n2 - 1

  if (alternative == "greater") {
    pval <- pf(Fstat, df1, df2, lower.tail = FALSE)
    ci <- c(qf(conf.level, df1, df2)^(-1) * Fstat, Inf)
  } else if (alternative == "less") {
    pval <- pf(Fstat, df1, df2, lower.tail = TRUE)
    ci <- c(0, qf(conf.level, df2, df1) * Fstat)
  } else {
    p_lower <- pf(Fstat, df1, df2, lower.tail = TRUE)
    p_upper <- pf(Fstat, df1, df2, lower.tail = FALSE)
    pval <- 2 * min(p_lower, p_upper)
    alpha <- 1 - conf.level
    ci <- c(Fstat / qf(1 - alpha/2, df1, df2), Fstat / qf(alpha/2, df1, df2))
  }

  list(
    method = "Two-sample F-test for equality of variances (from scratch)",
    estimate = c("s1^2" = s1, "s2^2" = s2, "ratio" = Fstat),
    parameter = c(df1 = df1, df2 = df2),
    statistic = c(F = Fstat),
    p.value = pval,
    conf.int = ci,
    alternative = alternative
  )
}

# ------------------------------------------------
# (2) Levene / Brown–Forsythe test (any number of groups)
#     center = "mean" (Levene) or "median" (Brown–Forsythe)
#     Performs one-way ANOVA on absolute deviations from group center.
# ------------------------------------------------
levene_like_test <- function(values, group, center = c("mean","median")) {
  center <- match.arg(center)
  df <- data.frame(values = values, group = factor(group))
  centers <- df %>% group_by(group) %>%
    summarize(center_val = if (center == "mean") mean(values) else stats::median(values), .groups = "drop")
  df <- df %>% left_join(centers, by = "group")
  df <- df %>% mutate(dev = abs(values - center_val))

  g <- nlevels(df$group)
  if (g < 2) stop("Need at least 2 groups.")
  N <- nrow(df)

  dev_summ <- df %>% group_by(group) %>% summarize(n_i = dplyr::n(), mean_i = mean(dev), .groups = "drop")
  grand_mean <- mean(df$dev)

  SS_between <- sum(dev_summ$n_i * (dev_summ$mean_i - grand_mean)^2)
  SS_within <- df %>% group_by(group) %>% summarize(ss = sum((dev - mean(dev))^2), .groups = "drop") %>% pull(ss) %>% sum()

  df1 <- g - 1
  df2 <- N - g
  MS_between <- SS_between / df1
  MS_within <- SS_within / df2
  Fstat <- MS_between / MS_within
  pval <- pf(Fstat, df1, df2, lower.tail = FALSE)

  list(
    method = if (center == "mean") "Levene's test (center = mean) from scratch" else "Brown–Forsythe test (center = median) from scratch",
    center = center,
    statistic = c(F = Fstat),
    parameter = c(df1 = df1, df2 = df2),
    p.value = pval,
    details = list(
      SS_between = SS_between, SS_within = SS_within,
      MS_between = MS_between, MS_within = MS_within,
      grand_mean = grand_mean
    )
  )
}

levene_test <- function(values, group) levene_like_test(values, group, center = "mean")
brown_forsythe_test <- function(values, group) levene_like_test(values, group, center = "median")

# ------------------------------------------------
# Visualization helpers
# ------------------------------------------------

# 1) Raw values: violin + boxplot + jitter
plot_values_by_group <- function(values, group, title = "Values by group") {
  df <- data.frame(values = values, group = factor(group))
  ggplot(df, aes(x = group, y = values, fill = group)) +
    geom_violin(alpha = 0.25, color = NA, trim = FALSE) +
    geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.55) +
    geom_jitter(width = 0.08, alpha = 0.6, size = 1.6) +
    guides(fill = "none") +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
}

# 2) Absolute deviations plot for a chosen center
plot_abs_deviations <- function(values, group, center = c("mean","median")) {
  center <- match.arg(center)
  df <- data.frame(values = values, group = factor(group))
  centers <- df %>% group_by(group) %>% summarize(center_val = if (center == "mean") mean(values) else stats::median(values), .groups = "drop")
  df <- df %>% left_join(centers, by = "group") %>% mutate(dev = abs(values - center_val))
  ggplot(df, aes(x = group, y = dev, fill = group)) +
    geom_boxplot(alpha = 0.6, width = 0.35, outlier.alpha = 0.35) +
    geom_jitter(width = 0.08, alpha = 0.5, size = 1.2) +
    guides(fill = "none") +
    labs(title = paste0("Abs. deviations (center = ", center, ")"), x = NULL, y = NULL) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
}

# 3) QQ plot by group (diagnose normality for F-test)
plot_qq_by_group <- function(values, group) {
  df <- data.frame(values = values, group = factor(group)) %>% arrange(group, values) %>% group_by(group)
  qq_df <- df %>% group_modify(~{
    x <- .x$values
    n <- length(x)
    tibble(
      sample = sort(x),
      theoretical = stats::qnorm(ppoints(n), mean = mean(x), sd = stats::sd(x))
    )
  })
  ggplot(qq_df, aes(x = theoretical, y = sample, color = group)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.6) +
    geom_point(alpha = 0.7, size = 1.5) +
    labs(title = "QQ plot vs fitted Normal (by group)", x = "Theoretical quantiles", y = "Sample quantiles") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
}

# 4) Assemble a diagnostic dashboard
`%||%` <- function(a, b) if (!is.null(a)) a else b

variance_dashboard <- function(values, group, title = NULL) {
  title <- title %||% "Variance diagnostics"
  p1 <- plot_values_by_group(values, group, title)
  p2 <- plot_abs_deviations(values, group, center = "mean")
  p3 <- plot_abs_deviations(values, group, center = "median")
  p4 <- plot_qq_by_group(values, group)
  (p1 | p4) / (p2 | p3)
}

# ------------------------------------------------
# Example generator
# ------------------------------------------------
make_example <- function(scenario = c("equal_normal","unequal_normal","non_normal"), n = 60, seed = 1) {
  set.seed(seed)
  scenario <- match.arg(scenario)
  n1 <- n2 <- as.integer(n/2)
  if (scenario == "equal_normal") {
    x <- rnorm(n1, mean = 0, sd = 1.0)
    y <- rnorm(n2, mean = 0.2, sd = 1.0)
  } else if (scenario == "unequal_normal") {
    x <- rnorm(n1, mean = 0, sd = 1.0)
    y <- rnorm(n2, mean = 0.2, sd = 2.0)
  } else {
    x <- rt(n1, df = 3) * 1.0 + 0
    y <- rt(n2, df = 3) * 1.0 + 0.2
  }
  tibble(values = c(x, y), group = factor(rep(c("A","B"), c(n1, n2))))
}

# ------------------------------------------------
# Pretty printer
# ------------------------------------------------
print_test <- function(res) {
  cat(res$method, "\n", sep = "")
  if (!is.null(res$center)) cat("center:", res$center, "\n")
  if (!is.null(res$estimate)) print(round(res$estimate, 6))
  cat(sprintf("F = %.4f, df1 = %d, df2 = %d, p = %.5f\n",
              res$statistic[["F"]], as.integer(res$parameter[["df1"]]), as.integer(res$parameter[["df2"]]), res$p.value))
  if (!is.null(res$conf.int)) cat("95% CI (variance ratio):", paste(round(res$conf.int, 4), collapse = " to "), "\n")
  cat("\n")
}

# ================================================================
# EXAMPLES — uncomment to run interactively
# ================================================================
if (interactive()) {
  # 1) Equal variances, Normal
  df1 <- make_example("equal_normal", n = 100, seed = 11)
  cat("--- Scenario 1: equal_normal ---\n")
  resF1 <- with(df1, f_test_equal_var(values[group=="A"], values[group=="B"], alternative = "two.sided"))
  resL1 <- with(df1, levene_test(values, group))
  resB1 <- with(df1, brown_forsythe_test(values, group))
  print_test(resF1); print_test(resL1); print_test(resB1)
  print(variance_dashboard(df1$values, df1$group, title = "Equal variances (Normal)"))

  # 2) Unequal variances, Normal
  df2 <- make_example("unequal_normal", n = 100, seed = 12)
  cat("--- Scenario 2: unequal_normal ---\n")
  resF2 <- with(df2, f_test_equal_var(values[group=="A"], values[group=="B"], alternative = "two.sided"))
  resL2 <- with(df2, levene_test(values, group))
  resB2 <- with(df2, brown_forsythe_test(values, group))
  print_test(resF2); print_test(resL2); print_test(resB2)
  print(variance_dashboard(df2$values, df2$group, title = "Unequal variances (Normal)"))

  # 3) Non-normal heavy tails (t3), equal variance
  df3 <- make_example("non_normal", n = 100, seed = 13)
  cat("--- Scenario 3: non_normal (heavy tails) ---\n")
  resF3 <- with(df3, f_test_equal_var(values[group=="A"], values[group=="B"], alternative = "two.sided"))
  resL3 <- with(df3, levene_test(values, group))
  resB3 <- with(df3, brown_forsythe_test(values, group))
  print_test(resF3); print_test(resL3); print_test(resB3)
  print(variance_dashboard(df3$values, df3$group, title = "Heavy tails (t3)"))
}
