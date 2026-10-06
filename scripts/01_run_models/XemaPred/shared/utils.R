# -----------------------------------------------------------------------
# XemaPred utility functions
#
# Purpose:
# - Provide small helper functions used across XemaPred cohort and
#   patient-specific workflows.
# - Keep these functions independent from the specific intensity,
#   subjective, and extent model implementations.
# -----------------------------------------------------------------------

`%||%` <- function(a, b) {
  if (!is.null(a)) a else b
}

as_matrix <- function(x) {
  if (is.matrix(x)) {
    return(x)
  }

  x <- as.matrix(x)

  if (is.null(nrow(x))) {
    x <- matrix(x, nrow = 1)
  }

  x
}

inv_logit <- function(z) {
  z_dim <- dim(z)

  out <- numeric(length(z))

  pos <- z >= 0
  out[pos] <- 1 / (1 + exp(-z[pos]))

  ez <- exp(z[!pos])
  out[!pos] <- ez / (1 + ez)

  if (!is.null(z_dim)) {
    dim(out) <- z_dim
  }

  out
}

logit <- function(p) {
  p <- pmin(pmax(p, 1e-12), 1 - 1e-12)
  log(p / (1 - p))
}

clip <- function(x, lo, hi) {
  pmax(lo, pmin(hi, x))
}

rhalfnorm <- function(n, sd) {
  abs(rnorm(n, mean = 0, sd = sd))
}

log_sum_exp <- function(v) {
  v <- v[is.finite(v)]

  if (!length(v)) {
    return(-Inf)
  }

  m <- max(v)
  m + log(sum(exp(v - m)))
}

log_mean_exp <- function(logv) {
  logv <- logv[is.finite(logv)]

  if (!length(logv)) {
    return(-Inf)
  }

  m <- max(logv)
  m + log(mean(exp(logv - m)))
}

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

compute_rps <- function(y_true, samples, M_max) {
  y_true <- as.integer(y_true)
  samples <- as.integer(samples)

  y_true <- max(0L, min(M_max, y_true))
  samples <- pmin(M_max, pmax(0L, samples))

  n <- length(samples)

  if (n == 0L) {
    return(NA_real_)
  }

  counts <- tabulate(samples + 1L, nbins = M_max + 1L)
  pmf <- counts / n

  cdf_pred <- cumsum(pmf)

  k <- 0:M_max
  cdf_true <- as.numeric(k >= y_true)

  sum((cdf_pred - cdf_true)^2)
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

wait_for_expected_patient_finals <- function(
    final_dir,
    patients,
    filename = "predictions.rds",
    tries = 30L,
    sleep_sec = 2
) {
  expected <- length(patients)

  for (k in seq_len(tries)) {
    files <- list.files(
      final_dir,
      pattern = paste0("^", filename, "$"),
      recursive = TRUE,
      full.names = FALSE
    )

    pats_present <- unique(sub(paste0("/", filename, "$"), "", files))
    n_present <- sum(grepl("^patient_[0-9]+$", pats_present))

    message(
      "[INFO] ",
      filename,
      " present: ",
      n_present,
      "/",
      expected,
      " (attempt ",
      k,
      "/",
      tries,
      ")"
    )

    if (n_present >= expected) {
      return(invisible(TRUE))
    }

    Sys.sleep(sleep_sec)
  }

  message(
    "[WARN] Timeout waiting for all ",
    filename,
    "; proceeding with whatever is visible."
  )

  invisible(FALSE)
}

collect_all_patients <- function(
    final_dir,
    filename,
    out_name,
    tries = 30L,
    sleep_sec = 2
) {
  for (k in seq_len(tries)) {
    files <- list.files(
      final_dir,
      pattern = paste0("^", filename, "$"),
      recursive = TRUE,
      full.names = TRUE
    )

    if (length(files) > 0) {
      objs <- purrr::map(files, safe_read_rds)
      objs <- objs[!vapply(objs, is.null, logical(1))]

      if (length(objs) > 0) {
        df <- dplyr::bind_rows(objs)
        out <- file.path(final_dir, out_name)

        saveRDS(df, out)

        message(
          "[INFO] Saved combined ",
          filename,
          " -> ",
          out,
          " (n_files=",
          length(objs),
          ")"
        )

        return(invisible(out))
      }
    }

    message(
      "[INFO] Collector waiting for ",
      filename,
      " (attempt ",
      k,
      "/",
      tries,
      ") ..."
    )

    Sys.sleep(sleep_sec)
  }

  message("[WARN] Could not create ", out_name, " (no readable inputs).")

  invisible(NULL)
}