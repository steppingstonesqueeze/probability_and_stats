

# Coin toss - histo of distance between successsive heads

sample_data <- c("H", "T")
N <- 10000
data <- sample(sample_data,
               N,
               replace = TRUE)

# Head locations
pos_H <- which(data == "H")
diff_pos_H <- diff(pos_H)

# Table
cat("Count table \n\n")
table(diff_pos_H)
Sys.sleep(2)

#Histo
h <- hist(diff_pos_H, col = "red", breaks = 10)

#Histo details

cat("Histogram counts \n")
h$counts
cat("Histogram breaks \n")
h$breaks
cat("Histogram midpoints \n")
h$mids
cat("Histogram density \n")
h$density

# Fit to histogram densities - exponential distribution and a power law 

# Build a histogram without plotting and use densities as targets
mids <- h$mids
dens <- h$density
w <- h$counts  # weights = counts

# Exponential fit on binned density: dens ≈ λ exp(-λ m)
start_lam <- 1/mean(x)
fit_exp_nls <- nls(dens ~ lam * exp(-lam * mids),
                   start = list(lam = start_lam),
                   weights = w)
lam_nls <- coef(fit_exp_nls)[["lam"]]

# Power-law fit on binned density: dens ≈ (α-1) xmin^{α-1} m^{-α}, m ≥ xmin
start_alpha <- 2
start_xmin  <- min(mids[mids > 0])
fit_pl_nls <- nls(dens ~ (alpha - 1) * xmin^(alpha - 1) * mids^(-alpha),
                  start = list(alpha = start_alpha, xmin = start_xmin),
                  subset = mids >= start_xmin, weights = w)
alpha_nls <- coef(fit_pl_nls)[["alpha"]]
xmin_nls  <- coef(fit_pl_nls)[["xmin"]]

list(
  exponential_equation = sprintf("f(x) = %.4f * exp(-%.4f * x)", lam_nls, lam_nls),
  powerlaw_equation    = sprintf("f(x) = (%.4f-1) * %.4f^{%.4f-1} * x^{-%.4f},  x ≥ %.4f",
                                 alpha_nls, xmin_nls, alpha_nls, alpha_nls, xmin_nls)
)








