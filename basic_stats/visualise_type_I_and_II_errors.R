plot_alpha_beta_ttest <- function(n, mu0, mu1, sd, alpha = 0.05,
                                  alternative = c("two.sided", "greater", "less")) {
  alternative <- match.arg(alternative)
  df  <- n - 1
  se  <- sd / sqrt(n)
  ncp <- (mu1 - mu0) / se  # noncentrality for t under H1
  
  # Critical values (central t under H0)
  crit <- switch(
    alternative,
    "two.sided" = qt(1 - alpha/2, df = df),
    "greater"   = qt(1 - alpha,   df = df),
    "less"      = qt(alpha,       df = df)
  )
  
  # Rejection and non-rejection regions in T-space
  rej_region <- switch(
    alternative,
    "two.sided" = list(left = -Inf, lcut = -crit, rcut = crit, right = Inf),
    "greater"   = list(left = crit, right = Inf),
    "less"      = list(left = -Inf, right = crit)
  )
  
  # Numeric alpha (under H0) and beta (under H1)
  if (alternative == "two.sided") {
    alpha_calc <- pt(-crit, df) + (1 - pt(crit, df))                     # should equal alpha
    beta_calc  <- pt(crit, df, ncp) - pt(-crit, df, ncp)                 # fail-to-reject under H1
  } else if (alternative == "greater") {
    alpha_calc <- 1 - pt(crit, df)
    beta_calc  <- pt(crit, df, ncp)
  } else { # "less"
    alpha_calc <- pt(crit, df)
    beta_calc  <- 1 - pt(crit, df, ncp)
  }
  
  # x-range to plot: cover both central and noncentral t supports
  left_q  <- min(qt(0.001, df), qt(0.001, df, ncp))
  right_q <- max(qt(0.999, df), qt(0.999, df, ncp))
  x <- seq(left_q, right_q, length.out = 4000)
  
  # Densities
  d0 <- dt(x, df)            # H0: central t
  d1 <- dt(x, df, ncp = ncp) # H1: noncentral t
  
  # Helper to polygon-shade regions
  shade_poly <- function(x, y, region) {
    polygon(x, y, border = NA, col = region$col)
  }
  
  # Start plot
  op <- par(no.readonly = TRUE); on.exit(par(op))
  par(mar = c(4.2, 4.2, 3.2, 1.2))
  plot(x, d0, type = "l", lwd = 2, col = "#1f77b4",
       xlab = "t statistic", ylab = "Density",
       main = sprintf("Type I (α, blue) and Type II (β, red) for one-sample t-test\nn=%d, df=%d, alpha=%.3f, alt=%s", 
                      n, df, alpha, alternative))
  lines(x, d1, lwd = 2, col = "#d62728")
  
  # Shade alpha (rejection under H0) in blue with transparency
  blue <- adjustcolor("#1f77b4", alpha.f = 0.35)
  red  <- adjustcolor("#d62728", alpha.f = 0.35)
  
  if (alternative == "two.sided") {
    # α regions: x <= -crit and x >= crit (under H0 density)
    idx_left  <- x <= -crit
    idx_right <- x >=  crit
    if (any(idx_left))  polygon(c(x[idx_left],  rev(x[idx_left])),
                                c(d0[idx_left], rep(0, sum(idx_left))), col = blue, border = NA)
    if (any(idx_right)) polygon(c(x[idx_right], rev(x[idx_right])),
                                c(d0[idx_right], rep(0, sum(idx_right))), col = blue, border = NA)
    
    # β region: -crit <= x <= crit (under H1 density)
    idx_beta <- (x >= -crit) & (x <= crit)
    if (any(idx_beta)) polygon(c(x[idx_beta],  rev(x[idx_beta])),
                               c(d1[idx_beta], rep(0, sum(idx_beta))), col = red, border = NA)
    
    abline(v = c(-crit, crit), lty = 2, lwd = 1.5)
    
  } else if (alternative == "greater") {
    # α region: x >= crit under H0
    idx_alpha <- x >= crit
    if (any(idx_alpha)) polygon(c(x[idx_alpha], rev(x[idx_alpha])),
                                c(d0[idx_alpha], rep(0, sum(idx_alpha))), col = blue, border = NA)
    # β region: x <= crit under H1
    idx_beta <- x <= crit
    if (any(idx_beta)) polygon(c(x[idx_beta], rev(x[idx_beta])),
                               c(d1[idx_beta], rep(0, sum(idx_beta))), col = red, border = NA)
    
    abline(v = crit, lty = 2, lwd = 1.5)
    
  } else { # "less"
    # α region: x <= crit under H0
    idx_alpha <- x <= crit
    if (any(idx_alpha)) polygon(c(x[idx_alpha], rev(x[idx_alpha])),
                                c(d0[idx_alpha], rep(0, sum(idx_alpha))), col = blue, border = NA)
    # β region: x >= crit under H1
    idx_beta <- x >= crit
    if (any(idx_beta)) polygon(c(x[idx_beta], rev(x[idx_beta])),
                               c(d1[idx_beta], rep(0, sum(idx_beta))), col = red, border = NA)
    
    abline(v = crit, lty = 2, lwd = 1.5)
  }
  
  legend("topright",
         legend = c(
           sprintf("H0: central t (blue)"),
           sprintf("H1: noncentral t (red), ncp = %.3f", ncp),
           sprintf("α ≈ %.4f", alpha_calc),
           sprintf("β ≈ %.4f  (power = %.4f)", beta_calc, 1 - beta_calc)
         ),
         lwd = c(2,2,NA,NA), col = c("#1f77b4", "#d62728", NA, NA),
         bty = "n", text.width = strwidth("H1: noncentral t (red), ncp = 0.000"),
         y.intersp = 1.1)
}

# -------------------------
# EXAMPLES
# -------------------------

# 1) Match your earlier data roughly: treat sample sd as population sd for illustration
x <- c(99, 77, 67, 70)
n  <- length(x)
mu0 <- 50
mu1 <- mean(x)      # "true" mean to visualize
sd_pop <- sd(x)     # using sample sd as proxy for pop sd for the picture

plot_alpha_beta_ttest(n, mu0, mu1, sd = sd_pop, alpha = 0.05, alternative = "greater")

# 2) Two-sided example with clearer separation
plot_alpha_beta_ttest(n = 25, mu0 = 0, mu1 = 0.4, sd = 1, alpha = 0.05, alternative = "two.sided")

# 3) Left-tailed example
plot_alpha_beta_ttest(n = 30, mu0 = 10, mu1 = 9.5, sd = 2, alpha = 0.01, alternative = "less")
