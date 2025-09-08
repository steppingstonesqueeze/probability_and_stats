# ================================================================
# Heavy-tailed distributions: generators + quick tail plots (R)
# ================================================================

# ---- Optional packages (used only if available) ----
.have_actuar     <- requireNamespace("actuar",     quietly = TRUE)  # Benini, Benktander, Burr, Lomax, loglogis, etc.
.have_VGAM       <- requireNamespace("VGAM",       quietly = TRUE)  # discrete Zipf, Yule-Simon, Waring, etc.
.have_extraDistr <- requireNamespace("extraDistr", quietly = TRUE)  # Zipf (alt)
.have_poweRlaw   <- requireNamespace("poweRlaw",   quietly = TRUE)  # discrete power-law
.have_stabledist <- requireNamespace("stabledist", quietly = TRUE)  # α-stable (incl. Lévy)
.have_libstableR <- requireNamespace("libstableR", quietly = TRUE)  # robust α-stable pdf (not needed for sampling)

# ================================================================
# A. Power-law / Pareto-type (regularly varying)
# ================================================================

# Pareto Type I: xm > 0, alpha > 0
rpareto1 <- function(n, alpha, xm = 1) { xm / runif(n)^(1/alpha) }

# Lomax (Pareto II): scale > 0, alpha > 0 ; support x>=0
# F(x) = 1 - (1 + x/scale)^(-alpha)  -> inverse-CDF
rlomax <- function(n, alpha, scale = 1) { scale * (runif(n)^(-1/alpha) - 1) }

# Generalized Pareto (xi>0 only here; heavy tail): beta > 0
# F(x) = 1 - (1 + xi x / beta)^(-1/xi), x>=0
rgpd_pos <- function(n, xi, beta) {
  if (xi <= 0) stop("This helper is for xi>0 (heavy tails).")
  beta/xi * (runif(n)^(-xi) - 1)
}

# Burr XII (Singh–Maddala): c>0 (tail), k>0 (shape), scale>0
# F(x) = 1 - (1 + (x/scale)^c)^(-k) -> inverse-CDF
rburr12 <- function(n, c, k, scale = 1) {
  scale * ((runif(n)^(-1/k) - 1)^(1/c))
}

# Beta prime (a,b,scale) : X = scale * (G_a / G_b), G_* ~ Gamma(shape, rate=1)
rbetaprime <- function(n, a, b, scale = 1) {
  g1 <- rgamma(n, shape = a, rate = 1)
  g2 <- rgamma(n, shape = b, rate = 1)
  scale * (g1 / g2)
}

# F distribution (beta-prime special case)
rf_F <- function(n, df1, df2) stats::rf(n, df1, df2)

# Student t and Cauchy (built-in)
rt_   <- function(n, nu) stats::rt(n, df = nu)
rcauchy_ <- function(n, loc = 0, scale = 1) stats::rcauchy(n, location = loc, scale = scale)

# Lévy (α=1/2 one-sided): prefer stabledist; fallback via  μ + c / Z^2 with Z~N(0,1)
rlevy <- function(n, c = 1, mu = 0) {
  if (.have_stabledist) {
    stabledist::rstable(n, alpha = 0.5, beta = 1, gamma = c, delta = mu, pm = 1)
  } else {
    z <- rnorm(n)
    mu + c / (z^2)  # matches the classical Lévy(μ, c) parameterization
  }
}

# General α-stable (0<alpha<2): require stabledist (no robust CMS implementation here)
rstable_ <- function(n, alpha, beta = 0, gamma = 1, delta = 0) {
  if (!.have_stabledist) stop("Install 'stabledist' for general α-stable sampling.")
  stabledist::rstable(n, alpha = alpha, beta = beta, gamma = gamma, delta = delta, pm = 1)
}

# Log-logistic (Fisk): shape>0, scale>0 ; X = scale*(U/(1-U))^(1/shape)
rloglogis <- function(n, shape, scale = 1) { scale * ( (runif(n)/(1 - runif(n)))^(1/shape) ) }

# Fréchet: shape alpha>0, scale s>0 ; F(x)=exp(-(s/x)^alpha), x>0 -> inverse
rfrechet <- function(n, alpha, scale = 1) { scale / (-log(runif(n)))^(1/alpha) }

# Inverse-Gamma (shape a>0, rate b>0): X = 1 / Gamma(a, rate=b)
rinvgamma <- function(n, shape, rate) 1 / rgamma(n, shape = shape, rate = rate)

# Benini (if 'actuar' available)
rbenini_ <- function(n, scale = 1, shape1 = 1, shape2 = 1) {
  if (.have_actuar) actuar::rbenini(n, shape1 = shape1, shape2 = shape2, scale = scale)
  else stop("Install 'actuar' for Benini.")
}

# Benktander Type I & II (if 'actuar' available)
rbenktander1_ <- function(n, shape1, shape2) {
  if (.have_actuar) actuar::rbenktander1(n, shape1 = shape1, shape2 = shape2)
  else stop("Install 'actuar' for Benktander I.")
}
rbenktander2_ <- function(n, shape1, shape2) {
  if (.have_actuar) actuar::rbenktander2(n, shape1 = shape1, shape2 = shape2)
  else stop("Install 'actuar' for Benktander II.")
}

# Double-Pareto (simple symmetric |X| ~ Pareto(xm=1, alpha); sign ±)
rdoublepareto <- function(n, alpha, xm = 1) {
  s <- ifelse(runif(n) < 0.5, -1, 1)
  s * rpareto1(n, alpha = alpha, xm = xm)
}

# Double-Pareto Lognormal (dPLN) — requires a dedicated package; here: simple “DP over Lognormal scale” proxy
# NOTE: This is a pedagogical proxy, not the exact Reed–Jorgensen dPLN.
rdpln_proxy <- function(n, alpha = 1.5, mu = 0, sigma = 1) {
  # sample a lognormal scale then apply double-Pareto tail
  s <- ifelse(runif(n) < 0.5, -1, 1)
  L <- rlnorm(n, meanlog = mu, sdlog = sigma)
  s * (L * rpareto1(n, alpha, xm = 1))
}

# ================================================================
# B. Heavy but not power-law (subexponential)
# ================================================================
rlognormal_ <- function(n, meanlog = 0, sdlog = 1) stats::rlnorm(n, meanlog, sdlog)

# Weibull with shape<1 (stretched exponential tail)
rweibull_stretched <- function(n, shape = 0.7, scale = 1) stats::rweibull(n, shape = shape, scale = scale)

# Log-Gamma, Log-Weibull (log of light-tailed -> heavy-tailed on original scale)
rloggamma   <- function(n, shape, rate = 1) exp(rgamma(n, shape = shape, rate = rate))
rlogweibull <- function(n, shape, scale = 1) exp(rweibull(n, shape = shape, scale = scale))

# ================================================================
# C. Discrete heavy-tailed
# ================================================================

# Zipf/Zeta (discrete Pareto) — several options:
# 1) If extraDistr available: extraDistr::rzipf(n, N, s) (bounded support)
# 2) If VGAM available: VGAM::rzeta(n, s) (unbounded via zeta — careful)
# 3) Fallback: truncated Zipf with inverse CDF up to Nmax
# Discrete Zipf / Zeta sampler
# - If VGAM is available: unbounded rzeta()
# - Else: truncated inverse-CDF with controllable Nmax (fast & dependency-free)
rzipf_disc <- function(n, s, Nmax = 1e5) {
  stopifnot(s > 1)  # zeta needs s>1 for finite normalizing constant
  if (.have_VGAM && is.infinite(Nmax)) {
    return(VGAM::rzeta(n, shape = s))              # unbounded support
  }
  # Truncated inverse-CDF sampler (default Nmax = 1e5)
  Nmax <- as.integer(Nmax)
  k <- 1:Nmax
  w <- k^(-s)
  W <- cumsum(w); W <- W / W[length(W)]            # CDF on {1,...,Nmax}
  u <- runif(n)
  1L + findInterval(u, W)                          # returns in 1..Nmax
}

# Zipf–Mandelbrot (requires VGAM)
rzipf_mandelbrot_ <- function(n, s, b) {
  if (.have_VGAM) VGAM::rzipf.mandelbrot(n, N = Inf, shape = s, b = b)
  else stop("Install 'VGAM' for Zipf–Mandelbrot.")
}

# Yule–Simon (requires VGAM)
ryulesimon_ <- function(n, rho) { if (.have_VGAM) VGAM::ryules(n, rho) else stop("Install 'VGAM' for Yule–Simon.") }

# Waring (requires VGAM)
rwaring_ <- function(n, a, b) { if (.have_VGAM) VGAM::rwaring(n, a, b) else stop("Install 'VGAM' for Waring.") }

# Poisson–Lognormal (always easy)
rpois_lognormal <- function(n, meanlog = 0, sdlog = 1) {
  lambda <- rlnorm(n, meanlog, sdlog)
  rpois(n, lambda)
}

# Poisson–Pareto (heavy mixing; mean infinite if alpha<=1)
rpois_pareto <- function(n, alpha, xm = 1) {
  lambda <- rpareto1(n, alpha, xm)
  rpois(n, lambda)
}

# ================================================================
# D. Multivariate heavy-tailed
# ================================================================

# Multivariate t (p-dim): μ vector, Σ positive-definite, ν>0
rmvt_ <- function(n, mu, Sigma, df) {
  p <- length(mu)
  C <- chol(Sigma)
  Z <- matrix(rnorm(n * p), n, p) %*% C
  g <- rgamma(n, shape = df/2, rate = df/2)  # ~ χ^2_df / df
  sweep(Z / sqrt(g), 2, mu, FUN = "+")
}

# (Multivariate α-stable would need dedicated packages; omitted here.)

# ================================================================
# Tail visualization helpers
# ================================================================
# Empirical survival S(x)=P(X>x) on log–log; linear ~ power tail
tail_plot <- function(x, main = "Empirical tail (log–log)", add = FALSE, col = "#2c7fb8") {
  x <- sort(x)
  n <- length(x)
  S <- rev(seq_len(n)) / (n + 1)   # empirical survival
  if (!add) {
    plot(x, S, log = "xy", type = "p", pch = 16, col = col,
         xlab = "x (log)", ylab = "S(x)=P(X>x) (log)", main = main, cex = 0.6)
  } else {
    points(x, S, pch = 16, col = col, cex = 0.6)
  }
}

# Semilog tail (log S vs x) — straight line suggests exponential tail; curvature indicates heavier
tail_semilog <- function(x, main = "Empirical tail (semi-log)", add = FALSE, col = "#d95f02") {
  x <- sort(x)
  n <- length(x)
  S <- rev(seq_len(n)) / (n + 1)
  if (!add) {
    plot(x, S, log = "y", type = "p", pch = 16, col = col,
         xlab = "x", ylab = "S(x)=P(X>x) (log y)", main = main, cex = 0.6)
  } else {
    points(x, S, pch = 16, col = col, cex = 0.6)
  }
}

# ================================================================
# Mini demos (run by hand)
# ================================================================
# Example usage:
set.seed(1)
x1 <- rpareto1(5e4, alpha = 1.7, xm = 1)
x2 <- rlognormal_(5e4, 0, 1)
x3 <- rt_(5e4, nu = 2.5)
par(mfrow = c(1,2))
tail_plot(x1, main="Pareto vs Lognormal vs t"); tail_plot(x2, add=TRUE, col="#31a354"); tail_plot(x3, add=TRUE, col="#756bb1")
legend("bottomleft", bty="n", pch=16, col=c("#2c7fb8","#31a354","#756bb1"), legend=c("Pareto(1.7)","Lognormal","t(2.5)"))
tail_semilog(x1, main="Semilog tails"); tail_semilog(x2, add=TRUE, col="#31a354"); tail_semilog(x3, add=TRUE, col="#756bb1")

Sys.sleep(4)
# # Discrete Zipf (truncated) vs Poisson-lognormal:
z1 <- rzipf_disc(1e4, s = 1.4)          # truncated at 1e5 by default
z2 <- if (.have_VGAM) rzipf_disc(1e4, 1.4, Nmax = Inf) else NA

par(mfrow = c(1,2))
hist(log(z1), breaks = 200, col = "#fee08b", border = NA, main = "Zipf (truncated)", xlab = "k")
if (!is.na(z2[1])) hist(log(z2), breaks = 200, col = "#fc8d59", border = NA, main = "Zipf (VGAM rzeta)", xlab = "k")





