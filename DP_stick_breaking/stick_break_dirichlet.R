#!/usr/bin/env Rscript
# stick_break_dirichlet.R
# Sample from Dirichlet(alpha) via finite stick-breaking.
#
# Usage examples:
#   Rscript stick_break_dirichlet.R --alphas 1,2,3 --n 5 --seed 42
#   Rscript stick_break_dirichlet.R --alphas 0.5,0.5,0.5,0.5 --n 10000 --summary --seed 1 --output samples.csv

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  # Simple flag parser
  out <- list(alphas=NULL, n=1L, seed=NULL, output=NULL, summary=FALSE)
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (key == "--alphas") {
      i <- i + 1L; out$alphas <- args[[i]]
    } else if (key == "--n") {
      i <- i + 1L; out$n <- as.integer(args[[i]])
    } else if (key == "--seed") {
      i <- i + 1L; out$seed <- as.integer(args[[i]])
    } else if (key == "--output") {
      i <- i + 1L; out$output <- args[[i]]
    } else if (key == "--summary") {
      out$summary <- TRUE
    } else {
      stop(paste0("Unknown argument: ", key))
    }
    i <- i + 1L
  }
  if (is.null(out$alphas)) stop("--alphas is required, e.g. --alphas 1,2,3")
  out
}

parse_alphas <- function(s) {
  vals <- as.numeric(strsplit(s, ",")[[1]])
  if (any(!is.finite(vals)) || any(vals <= 0)) stop("All alphas must be finite and > 0.")
  if (length(vals) < 2) stop("Need at least 2 alphas for a Dirichlet.")
  vals
}

stick_break_once <- function(alpha) {
  K <- length(alpha)
  pies <- numeric(K)
  remaining <- 1.0
  if (K > 1) {
    for (i in 1:(K-1)) {
      a_i <- alpha[i]
      b_i <- sum(alpha[(i+1):K])
      V_i <- stats::rbeta(1L, shape1 = a_i, shape2 = b_i)
      pies[i] <- remaining * V_i
      remaining <- remaining * (1 - V_i)
    }
  }
  pies[K] <- remaining
  pies
}

main <- function() {
  opts <- parse_args()
  alpha <- parse_alphas(opts$alphas)
  if (!is.null(opts$seed)) set.seed(opts$seed)
  K <- length(alpha)
  n <- opts$n
  samples <- matrix(0.0, nrow=n, ncol=K)
  for (i in 1:n) {
    samples[i,] <- stick_break_once(alpha)
  }
  # Check row sums
  row_sums <- rowSums(samples)
  if (!all(abs(row_sums - 1.0) < 1e-10)) {
    warning("Some rows do not sum to 1 within tolerance.")
  }

  if (!is.null(opts$output)) {
    colnames(samples) <- paste0("pi_", seq_len(K))
    # write.csv adds row names by default; disable:
    utils::write.csv(samples, file = opts$output, row.names = FALSE, quote = FALSE)
    cat(sprintf("Wrote %d samples to %s\n", n, opts$output))
  } else {
    max_show <- min(n, 10L)
    for (i in 1:max_show) {
      cat(paste(sprintf("%.6f", samples[i,]), collapse = ", "), "\n", sep="")
    }
    if (n > max_show) cat(sprintf("... (%d more rows)\n", n - max_show))
  }

  if (isTRUE(opts$summary)) {
    emp_mean <- colMeans(samples)
    th_mean <- alpha / sum(alpha)
    cat("\nSummary (empirical vs theoretical means):\n", sep="")
    for (j in 1:K) {
      cat(sprintf("pi_%d: emp=%.6f  th=%.6f\n", j, emp_mean[j], th_mean[j]))
    }
    cat(sprintf("Total alpha: %.6f\n", sum(alpha)))
  }
}

if (identical(environment(), globalenv())) {
  tryCatch(main(), error=function(e) { message("Error: ", e$message); quit(status=1) })
}
