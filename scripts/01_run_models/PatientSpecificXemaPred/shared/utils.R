# -----------------------------------------------------------------------
# Shared utilities for PatientSpecificXemaPred
# -----------------------------------------------------------------------

`%||%` <- function(a, b) if (!is.null(a)) a else b

inv_logit <- function(z) {
  z_dim <- dim(z)

  out <- numeric(length(z))
  pos <- z >= 0
  out[pos] <- 1 / (1 + exp(-z[pos]))
  ez <- exp(z[!pos])
  out[!pos] <- ez / (1 + ez)

  if (!is.null(z_dim)) dim(out) <- z_dim
  out
}

logit <- function(p) {
  p <- pmin(pmax(p, 1e-12), 1 - 1e-12)
  log(p / (1 - p))
}

clip <- function(x, lo, hi) pmax(lo, pmin(hi, x))

as_matrix <- function(x) {
  if (is.matrix(x)) return(x)
  x <- as.matrix(x)
  if (is.null(nrow(x))) x <- matrix(x, nrow = 1)
  x
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

log_sum_exp <- function(v) {
  v <- v[is.finite(v)]
  if (!length(v)) return(-Inf)
  m <- max(v)
  m + log(sum(exp(v - m)))
}

log_mean_exp <- function(logv) {
  logv <- logv[is.finite(logv)]
  if (!length(logv)) return(-Inf)
  m <- max(logv)
  m + log(mean(exp(logv - m)))
}

safe_read_rds <- function(f) {
  tryCatch(
    readRDS(f),
    error = function(e) {
      message("[WARN] readRDS failed: ", f, " | ", conditionMessage(e))
      NULL
    }
  )
}

elapsed_seconds <- function(a, b) {
  if (any(is.na(c(a, b)))) return(NA_real_)
  as.numeric(difftime(b, a, units = "secs"))
}
