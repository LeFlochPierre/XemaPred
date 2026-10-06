# -----------------------------------------------------------------------
# XemaPred particle-filter utilities
#
# Purpose:
# - Provide particle-filter-specific helper functions used by both
#   population-level and patient-specific XemaPred models.
# - Keep SMC utilities separate from the model-specific likelihoods.
# -----------------------------------------------------------------------

systematic_resample_idx <- function(w) {
  n <- length(w)

  if (n <= 0L) {
    return(integer(0))
  }

  if (n == 1L) {
    return(1L)
  }

  w <- as.numeric(w)
  w[!is.finite(w)] <- 0
  w[w < 0] <- 0

  s <- sum(w)

  if (!is.finite(s) || s <= 0) {
    w <- rep(1 / n, n)
  } else {
    w <- w / s
  }

  cdf <- cumsum(w)
  cdf[!is.finite(cdf)] <- 0
  cdf <- cummax(cdf)

  if (anyNA(cdf) || is.unsorted(cdf, strictly = FALSE) || cdf[n] <= 0) {
    return(sample.int(n, size = n, replace = TRUE, prob = w))
  }

  cdf[n] <- 1

  u0 <- runif(1) / n
  u <- u0 + (0:(n - 1)) / n

  # idx <- tryCatch(
  #   findInterval(
  #     u,
  #     cdf,
  #     rightmost.closed = TRUE,
  #     all.inside = TRUE
  #   ) + 1L,
  #   error = function(e) {
  #     sample.int(n, size = n, replace = TRUE, prob = w)
  #   }
  # )

  idx <- tryCatch(
    {
      j <- findInterval(u, cdf, rightmost.closed = TRUE) + 1L
      pmin(pmax(j, 1L), n)
    },
    error = function(e) {
      sample.int(n, size = n, replace = TRUE, prob = w)
    }
  )

  idx
}

liu_west_move <- function(theta_u, a = 0.98) {
  theta_u <- as.matrix(theta_u)

  np <- nrow(theta_u)
  P <- ncol(theta_u)

  cn <- colnames(theta_u)

  if (np < 2L || P < 1L) {
    return(theta_u)
  }

  m <- colMeans(theta_u)
  V <- stats::cov(theta_u)

  if (any(!is.finite(V))) {
    return(theta_u)
  }

  V <- V + diag(1e-10, P)

  h2 <- 1 - a^2
  R <- chol(h2 * V)

  Z <- matrix(rnorm(np * P), np, P)

  out <- a * theta_u +
    matrix((1 - a) * m, np, P, byrow = TRUE) +
    Z %*% R

  if (!is.null(cn)) {
    colnames(out) <- cn
  }

  out
}

liu_west_move_select <- function(theta_u, a = 0.98, cols_keep = NULL) {
  theta_u <- as.matrix(theta_u)

  cn <- colnames(theta_u)

  if (is.null(cn)) {
    stop("theta_u must have colnames")
  }

  if (is.null(cols_keep)) {
    cols_keep <- intersect(
      cn,
      c(
        "log_sigma_meas",
        "log_sigma_lat",
        "mu_y0",
        "log_sigma_y0"
      )
    )
  } else {
    cols_keep <- intersect(cols_keep, cn)
  }

  if (!length(cols_keep)) {
    return(theta_u)
  }

  out <- theta_u

  out[, cols_keep] <- liu_west_move(
    theta_u[, cols_keep, drop = FALSE],
    a = a
  )

  colnames(out) <- cn

  out
}

rtrunc_norm <- function(n, mean, sd, lo = -Inf, hi = Inf) {
  out <- numeric(n)
  i <- 1L

  while (i <= n) {
    z <- rnorm(n - i + 1L, mean, sd)
    z <- z[z >= lo & z <= hi]

    if (length(z)) {
      take <- min(length(z), n - i + 1L)
      out[i:(i + take - 1L)] <- z[1:take]
      i <- i + take
    }
  }

  out
}