
# t-test Monte Carlo (R version) -------------------------------------------
# Pop1: N(0,1), Pop2: N(1,1), m=n=12, alpha=0.05, two-sided
# Methods: pooled t, Welch t, Yuen trimmed-mean t (trim=0.2)
# K=200 batches, T=1000 trials per batch

set.seed(123)
pop1 <- rnorm(100000, 0, 1)
pop2 <- rnorm(100000, 1, 1)

pooled_reject_rate <- function(xmat, ymat, alpha=0.05) {
  m <- ncol(xmat); n <- ncol(ymat)
  mx <- rowMeans(xmat); my <- rowMeans(ymat)
  vx <- apply(xmat, 1, var); vy <- apply(ymat, 1, var)
  sp2 <- ((m-1)*vx + (n-1)*vy) / (m + n - 2)
  se <- sqrt(sp2 * (1/m + 1/n))
  tvals <- (mx - my) / se
  df <- m + n - 2
  tcrit <- qt(1 - alpha/2, df)
  mean(abs(tvals) > tcrit)
}

welch_reject_rate <- function(xmat, ymat, alpha=0.05) {
  m <- ncol(xmat); n <- ncol(ymat)
  mx <- rowMeans(xmat); my <- rowMeans(ymat)
  vx <- apply(xmat, 1, var); vy <- apply(ymat, 1, var)
  se2 <- vx/m + vy/n
  tvals <- (mx - my) / sqrt(se2)
  df <- (se2^2) / ( (vx^2)/((m^2)*(m-1)) + (vy^2)/((n^2)*(n-1)) )
  tcrit <- qt(1 - alpha/2, df)
  mean(abs(tvals) > tcrit)
}

yuen_reject_rate <- function(xmat, ymat, trim=0.2, alpha=0.05) {
  m <- ncol(xmat); n <- ncol(ymat)
  gx <- floor(trim * m); gy <- floor(trim * n)
  hx <- m - 2*gx; hy <- n - 2*gy
  if (hx < 2 || hy < 2) stop("Not enough data after trimming; reduce trim.")
  sort_rows <- function(M) t(apply(M, 1, sort))
  sx <- sort_rows(xmat); sy <- sort_rows(ymat)
  tx <- sx[, (gx+1):(m-gx), drop=FALSE]
  ty <- sy[, (gy+1):(n-gy), drop=FALSE]
  xt <- rowMeans(tx); yt <- rowMeans(ty)
  wx <- sx; wy <- sy
  if (gx > 0) {
    wx[, 1:gx] <- sx[, gx+1, drop=FALSE]
    wx[, (m-gx+1):m] <- sx[, m-gx, drop=FALSE]
  }
  if (gy > 0) {
    wy[, 1:gy] <- sy[, gy+1, drop=FALSE]
    wy[, (n-gy+1):n] <- sy[, n-gy, drop=FALSE]
  }
  swx2 <- apply(wx, 1, var); swy2 <- apply(wy, 1, var)
  vx <- swx2 / (hx * (hx - 1)); vy <- swy2 / (hy * (hy - 1))
  tvals <- (xt - yt) / sqrt(vx + vy)
  df <- (vx + vy)^2 / ( (vx^2)/(hx - 1) + (vy^2)/(hy - 1) )
  tcrit <- qt(1 - alpha/2, df)
  mean(abs(tvals) > tcrit)
}

run_batched_mc <- function(pop1, pop2, m=12, n=12, alpha=0.05, K=200, T=1000, trim=0.2, seed=20250907) {
  set.seed(seed)
  N1 <- length(pop1); N2 <- length(pop2)
  pooled_rates <- numeric(K)
  welch_rates  <- numeric(K)
  yuen_rates   <- numeric(K)
  for (k in seq_len(K)) {
    idx_x <- matrix(sample.int(N1, T*m, replace=TRUE), nrow=T)
    idx_y <- matrix(sample.int(N2, T*n, replace=TRUE), nrow=T)
    xmat <- matrix(pop1[idx_x], nrow=T)
    ymat <- matrix(pop2[idx_y], nrow=T)
    pooled_rates[k] <- pooled_reject_rate(xmat, ymat, alpha)
    welch_rates[k]  <- welch_reject_rate(xmat, ymat, alpha)
    yuen_rates[k]   <- yuen_reject_rate(xmat, ymat, trim, alpha)
  }
  list(pooled=pooled_rates, welch=welch_rates, yuen=yuen_rates)
}

summarize <- function(v) {
  c(mean=mean(v), sd=sd(v), min=min(v),
    p2.5=quantile(v, 0.025), p25=quantile(v, 0.25),
    median=median(v), p75=quantile(v, 0.75), p97.5=quantile(v, 0.975),
    max=max(v))
}

# ---- Run ----
res <- run_batched_mc(pop1, pop2, m=12, n=12, alpha=0.05, K=200, T=1000, trim=0.2)

summary <- rbind(
  "Pooled t"       = summarize(res$pooled),
  "Welch t"        = summarize(res$welch),
  "Yuen (trim=0.2)"= summarize(res$yuen)
)
summary_df <- data.frame(method=rownames(summary), round(summary, 6), row.names = NULL)

outdir <- "ttest_mc_outputs_R"
dir.create(outdir, showWarnings = FALSE)
write.csv(summary_df, file.path(outdir, "reject_rate_summary.csv"), row.names = FALSE)

png(file.path(outdir, "hist_pooled.png"), width=1200, height=700, res=150)
hist(res$pooled, breaks=20, main="Pooled t: batched reject rates (m=n=12, alpha=0.05)",
     xlab="Batched reject rate")
dev.off()

png(file.path(outdir, "hist_welch.png"), width=1200, height=700, res=150)
hist(res$welch, breaks=20, main="Welch t: batched reject rates (m=n=12, alpha=0.05)",
     xlab="Batched reject rate")
dev.off()

png(file.path(outdir, "hist_yuen.png"), width=1200, height=700, res=150)
hist(res$yuen, breaks=20, main="Yuen (trim=0.2): batched reject rates (m=n=12, alpha=0.05)",
     xlab="Batched reject rate")
dev.off()
