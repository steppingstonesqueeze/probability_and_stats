## ---------- LLN & CLT demo: nice vs nasty ----------
set.seed(123)

## Pareto generator (xm > 0, alpha > 0)
rpareto <- function(n, alpha = 1.5, xm = 1) {
  xm / (runif(n)^(1/alpha))
}

## True means/vars for the two dists we’ll compare
dist_specs <- list(
  exp = list(rfun = function(n) rexp(n, rate = 1),
             mean = 1, var = 1, has_finite_var = TRUE,
             label = "Exponential(rate=1)"),
  pareto = list(rfun = function(n) rpareto(n, alpha = 1.5, xm = 1),
                mean = 1.5/(0.5) * 1,  # alpha/(alpha-1) * xm with alpha=1.5, xm=1 => 3
                var = Inf, has_finite_var = FALSE,
                label = "Pareto(alpha=1.5, xm=1)")
)

## Helper: running mean path (LLN picture)
running_mean <- function(x) cumsum(x) / seq_along(x)

## CLT helper: standardized sample means
## - If finite var: Z_n = sqrt(n)*(xbar - mu)/sigma (known sigma)
## - If infinite var: T_n = sqrt(n)*(xbar - mu)/s (studentized; still fails nicely)
standardized_means <- function(rfun, mu, sigma2, n, reps, finite_var) {
  m <- replicate(reps, mean(rfun(n)))
  if (finite_var) {
    z <- sqrt(n) * (m - mu) / sqrt(sigma2)
  } else {
    s <- replicate(reps, sd(rfun(n)))
    z <- sqrt(n) * (m - mu) / s  # not N(0,1) when variance is infinite
  }
  z
}

## Coverage experiment:
## - If finite var: known-sigma z-CI for mean
## - If infinite var: t-CI using sample sd (will under-cover)
coverage_rate <- function(rfun, mu, sigma2, n, reps, finite_var, level = 0.95) {
  alpha <- 1 - level
  if (finite_var) {
    zcrit <- qnorm(1 - alpha/2)
    inside <- replicate(reps, {
      x <- rfun(n)
      m <- mean(x)
      se <- sqrt(sigma2 / n)
      (mu >= m - zcrit * se) && (mu <= m + zcrit * se)
    })
  } else {
    tcrit <- qt(1 - alpha/2, df = n - 1)
    inside <- replicate(reps, {
      x <- rfun(n)
      m <- mean(x)
      s <- sd(x)
      se <- s / sqrt(n)
      (mu >= m - tcrit * se) && (mu <= m + tcrit * se)
    })
  }
  mean(inside)
}

## ---------------- Run the demo ----------------
n_path <- 20000        # length of single LLN path
n_grid <- c(10, 30, 100, 500, 2000)
reps   <- 3000         # for CLT hist & coverage
par(mfrow = c(2, 3), mar = c(4.2, 4.2, 3, 1))

for (name in c("exp", "pareto")) {
  sp <- dist_specs[[name]]
  
  ## 1) LLN path
  x <- sp$rfun(n_path)
  mu <- sp$mean
  mu_hat <- running_mean(x)
  plot(mu_hat, type = "l", lwd = 2,
       main = paste0("LLN: ", sp$label, "\nRunning mean vs n"),
       xlab = "n", ylab = expression(bar(X)[n]))
  abline(h = mu, col = 2, lwd = 2, lty = 2)
  legend("topright", legend = c("running mean", "true mean"),
         lty = c(1, 2), lwd = c(2, 2), bty = "n")
  
  ## 2) CLT hist at moderate n
  n1 <- n_grid[3]  # 100
  z1 <- standardized_means(sp$rfun, mu, sp$var, n1, reps, sp$has_finite_var)
  hist(z1, breaks = 50, freq = FALSE,
       main = paste0("CLT @ n=", n1, "\n", sp$label),
       xlab = "standardized mean")
  curve(dnorm(x), add = TRUE, lwd = 2, lty = 2)
  
  ## 3) CLT hist at large n
  n2 <- tail(n_grid, 1)  # 2000
  z2 <- standardized_means(sp$rfun, mu, sp$var, n2, reps, sp$has_finite_var)
  hist(z2, breaks = 50, freq = FALSE,
       main = paste0("CLT @ n=", n2, "\n", sp$label),
       xlab = "standardized mean")
  curve(dnorm(x), add = TRUE, lwd = 2, lty = 2)
  
  ## Print coverage for a few n
  cat("\n", sp$label, "coverage (", if (sp$has_finite_var) "z-CI" else "t-CI",
      ", 95% target):\n", sep = "")
  for (n in n_grid) {
    cov <- coverage_rate(sp$rfun, mu, sp$var, n, reps, sp$has_finite_var, level = 0.95)
    cat(sprintf("  n=%4d : %.3f\n", n, cov))
  }
}

par(mfrow = c(1,1))
