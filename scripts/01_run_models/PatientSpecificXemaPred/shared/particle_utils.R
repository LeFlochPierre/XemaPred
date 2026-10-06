# -----------------------------------------------------------------------
# Particle filtering utilities for PatientSpecificXemaPred
# -----------------------------------------------------------------------

log_normalize <- function(logw) {
  logw <- as.numeric(logw)
  logw[!is.finite(logw)] <- -Inf

  m <- max(logw)
  if (!is.finite(m)) {
    n <- length(logw)
    return(rep(1 / n, n))
  }

  w <- exp(logw - m)
  s <- sum(w)

  if (!is.finite(s) || s <= 0) {
    n <- length(logw)
    return(rep(1 / n, n))
  }

  w / s
}

normalize_weights <- function(w) {
  w <- as.numeric(w)
  w[!is.finite(w)] <- 0
  w[w < 0] <- 0

  s <- sum(w)
  if (!is.finite(s) || s <= 0) {
    return(rep(1 / length(w), length(w)))
  }

  w / s
}

systematic_resample_idx <- function(w) {
  n <- length(w)
  if (n <= 0L) return(integer(0))
  if (n == 1L) return(1L)

  w <- normalize_weights(w)
  cdf <- cumsum(w)
  cdf[!is.finite(cdf)] <- 0
  cdf <- cummax(cdf)

  if (anyNA(cdf) || is.unsorted(cdf, strictly = FALSE) || cdf[n] <= 0) {
    return(sample.int(n, size = n, replace = TRUE, prob = w))
  }

  cdf[n] <- 1
  u0 <- runif(1) / n
  u <- u0 + (0:(n - 1)) / n

  # tryCatch(
  #   findInterval(u, cdf, rightmost.closed = TRUE, all.inside = TRUE) + 1L,
  #   error = function(e) sample.int(n, size = n, replace = TRUE, prob = w)
  # )
  tryCatch(
    {
      # findInterval returns 0 when u < cdf[1], and +1 maps that to particle 1.
      # all.inside = TRUE would clamp it to 1..n-1, making particle 1
      # unreachable and handing its weight to particle 2.
      j <- findInterval(u, cdf, rightmost.closed = TRUE) + 1L
      pmin(pmax(j, 1L), n)
    },
    error = function(e) sample.int(n, size = n, replace = TRUE, prob = w)
  )
  
}

liu_west_move <- function(theta_u, a = 0.98) {
  theta_u <- as.matrix(theta_u)
  np <- nrow(theta_u)
  P <- ncol(theta_u)
  cn <- colnames(theta_u)

  if (np < 2L || P < 1L) return(theta_u)

  m <- colMeans(theta_u)
  V <- stats::cov(theta_u)
  if (any(!is.finite(V))) return(theta_u)

  V <- V + diag(1e-10, P)
  h2 <- 1 - a^2
  R <- chol(h2 * V)
  Z <- matrix(rnorm(np * P), np, P)

  out <- a * theta_u +
    matrix((1 - a) * m, np, P, byrow = TRUE) +
    Z %*% R

  if (!is.null(cn)) colnames(out) <- cn
  out
}

liu_west_move_select <- function(theta_u, a = 0.98, cols_keep = NULL) {
  theta_u <- as.matrix(theta_u)
  cn <- colnames(theta_u)
  if (is.null(cn)) stop("theta_u must have colnames")

  if (is.null(cols_keep)) {
    cols_keep <- intersect(cn, c("log_sigma_meas", "log_sigma_lat", "mu_y0", "log_sigma_y0"))
  } else {
    cols_keep <- intersect(cols_keep, cn)
  }

  if (!length(cols_keep)) return(theta_u)

  out <- theta_u
  out[, cols_keep] <- liu_west_move(theta_u[, cols_keep, drop = FALSE], a = a)
  colnames(out) <- cn
  out
}
