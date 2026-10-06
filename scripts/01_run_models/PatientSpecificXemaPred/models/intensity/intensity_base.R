# ======================================================================
# PatientSpecificXemaPred
#
# OrderedRW mathematical definition matched to population XemaPred.
# Patient-specific inference architecture is retained.
# ======================================================================

suppressPackageStartupMessages({
  library(stats)
})


# ----------------------------------------------------------------------
# Cutpoint utilities
# ----------------------------------------------------------------------

make_ct_mat <- function(delta, M_max) {
  stopifnot(is.matrix(delta), ncol(delta) == (M_max - 1L))

  np <- nrow(delta)
  ct <- matrix(0.0, nrow = np, ncol = M_max)

  if (M_max >= 2L) {
    ct[, 2:M_max] <- t(apply(delta, 1, cumsum))
  }

  # Matches EczemaPred's make_ct() in stan/include/functions_OrderedRW.stan:
  #   ct = cumulative_sum(append_row(0, delta)) * (M - 1) + 0.5
  ct * (M_max - 1) + 0.5
}


rdirichlet_mat <- function(n, alpha) {
  K <- length(alpha)

  x <- matrix(
    rgamma(
      n * K,
      shape = rep(alpha, each = n),
      rate = 1
    ),
    nrow = n,
    ncol = K
  )

  x / rowSums(x)
}


as_matrix <- function(x) {
  if (is.matrix(x)) {
    return(x)
  }

  x <- as.matrix(x)

  if (is.null(nrow(x))) {
    x <- matrix(
      x,
      nrow = 1L
    )
  }

  x
}


# ----------------------------------------------------------------------
# Priors
# ----------------------------------------------------------------------

get_priors_intensity <- function(M_max) {
  list(
    # Relative cutpoint increments lie on a simplex.
    prior_delta = rep(
      1,
      M_max - 1L
    ),

    mu0_mean = 0.5,
    mu0_sd = 0.25,
    sigma0_sd = 0.125,

    # sigma_lat / M ~ LogNormal(...)
    prior_sigma_lat = c(
      -log(10),
      0.5 * log(4)
    ),

    # sigma_meas / M ~ LogNormal(...)
    prior_sigma_meas = c(
      -log(10),
      0.5 * log(4)
    )
  )
}


theta_colnames_int <- function(M_max) {
  c(
    "log_sigma_meas",
    "log_sigma_lat",
    "mu_y0",
    "log_sigma_y0",
    paste0(
      "log_delta_",
      seq_len(M_max - 1L)
    )
  )
}


# ----------------------------------------------------------------------
# Initialise parameter particles
# ----------------------------------------------------------------------

init_theta_particles_int <- function(
    np,
    priors,
    M_max
) {
  # Same Dirichlet cutpoint prior as population XemaPred.
  delta <- rdirichlet_mat(
    np,
    priors$prior_delta
  )

  sigma_meas <- M_max * rlnorm(
    np,
    meanlog = priors$prior_sigma_meas[1],
    sdlog = priors$prior_sigma_meas[2]
  )

  sigma_meas <- pmax(
    sigma_meas,
    1e-8
  )

  sigma_lat <- M_max * rlnorm(
    np,
    meanlog = priors$prior_sigma_lat[1],
    sdlog = priors$prior_sigma_lat[2]
  )

  sigma_lat <- pmax(
    sigma_lat,
    1e-8
  )

  mu_y0 <- M_max * rnorm(
    np,
    mean = priors$mu0_mean,
    sd = priors$mu0_sd
  )

  sigma_y0 <- abs(
    M_max * rnorm(
      np,
      mean = 0,
      sd = priors$sigma0_sd
    )
  )

  sigma_y0 <- pmax(
    sigma_y0,
    1e-8
  )

  list(
    delta = delta,
    sigma_meas = sigma_meas,
    sigma_lat = sigma_lat,
    mu_y0 = mu_y0,
    sigma_y0 = sigma_y0
  )
}


# ----------------------------------------------------------------------
# Parameter transformations
# ----------------------------------------------------------------------

pack_theta_int <- function(
    delta,
    sigma_meas,
    sigma_lat,
    mu_y0,
    sigma_y0,
    M_max
) {
  # Centred log coordinates for the simplex.
  log_delta <- log(
    pmax(
      delta,
      1e-12
    )
  )

  log_delta <- log_delta -
    rowMeans(log_delta)

  out <- cbind(
    log_sigma_meas = log(
      pmax(
        sigma_meas,
        1e-12
      )
    ),
    log_sigma_lat = log(
      pmax(
        sigma_lat,
        1e-12
      )
    ),
    mu_y0 = mu_y0,
    log_sigma_y0 = log(
      pmax(
        sigma_y0,
        1e-12
      )
    ),
    log_delta
  )

  colnames(out) <- theta_colnames_int(
    M_max
  )

  out
}


unpack_theta_int <- function(
    theta_u,
    M_max
) {
  theta_u <- as_matrix(theta_u)

  delta_cols <- paste0(
    "log_delta_",
    seq_len(M_max - 1L)
  )

  log_delta <- theta_u[
    ,
    delta_cols,
    drop = FALSE
  ]

  log_delta <- log_delta -
    rowMeans(log_delta)

  log_delta <- log_delta -
    apply(
      log_delta,
      1,
      max
    )

  delta_raw <- exp(log_delta)

  delta <- delta_raw /
    rowSums(delta_raw)

  list(
    sigma_meas = exp(
      theta_u[, "log_sigma_meas"]
    ),
    sigma_lat = exp(
      theta_u[, "log_sigma_lat"]
    ),
    mu_y0 = theta_u[, "mu_y0"],
    sigma_y0 = exp(
      theta_u[, "log_sigma_y0"]
    ),
    delta = delta
  )
}


# ----------------------------------------------------------------------
# Ordered-logistic likelihood
# ----------------------------------------------------------------------

lik_intensity_ordlogit <- function(
    y_cat,
    x,
    sigma_meas,
    delta,
    M_max
) {
  y_cat <- as.integer(y_cat)

  K <- M_max + 1L

  if (y_cat < 1L || y_cat > K) {
    stop(
      "y_cat out of range: must be 1..K"
    )
  }

  np <- length(x)

  if (!is.matrix(delta)) {
    delta <- as.matrix(delta)
  }

  if (ncol(delta) != M_max - 1L) {
    stop(
      "delta must have M_max - 1 columns"
    )
  }

  if (nrow(delta) == 1L && np > 1L) {
    delta <- delta[
      rep(1L, np),
      ,
      drop = FALSE
    ]
  } else if (nrow(delta) != np) {
    stop(
      "delta must be 1 x (M_max-1) ",
      "or np x (M_max-1)"
    )
  }

  if (length(sigma_meas) == 1L) {
    s <- rep(
      pmax(
        as.numeric(sigma_meas),
        1e-12
      ) * sqrt(3) / pi,
      np
    )
  } else if (length(sigma_meas) == np) {
    s <- pmax(
      as.numeric(sigma_meas),
      1e-12
    ) * sqrt(3) / pi
  } else {
    stop(
      "sigma_meas must be scalar or length np"
    )
  }

  ct <- make_ct_mat(
    delta,
    M_max
  )

  z_lat <- x / s
  z_ct <- ct / s

  Fk <- inv_logit(
    z_ct - z_lat
  )

  if (y_cat == 1L) {
    p <- Fk[, 1L]
  } else if (y_cat == K) {
    p <- 1 - Fk[, M_max]
  } else {
    p <-
      Fk[, y_cat] -
      Fk[, y_cat - 1L]
  }

  pmax(
    p,
    1e-12
  )
}

# ----------------------------------------------------------------------
# Forecast from the current patient-specific OrderedRW state
# ----------------------------------------------------------------------

pf_forecast_intensity <- function(
    x_last,
    w_last,
    theta_last_u,
    time_last,
    time_test,
    np_sim = 2000L,
    M_max
) {
  time_test <- as.integer(time_test)
  np_sim <- as.integer(np_sim)

  H <- length(time_test)
  K <- M_max + 1L

  theta_last_u <- as_matrix(theta_last_u)

  if (length(x_last) != length(w_last)) {
    stop(
      "x_last and w_last must have the same length."
    )
  }

  if (nrow(theta_last_u) != length(x_last)) {
    stop(
      "theta_last_u must have one row per latent particle."
    )
  }

  idx <- sample.int(
    n = length(x_last),
    size = np_sim,
    replace = TRUE,
    prob = normalize_weights(w_last)
  )

  x <- as.numeric(
    x_last[idx]
  )

  theta_u <- theta_last_u[
    idx,
    ,
    drop = FALSE
  ]

  x_fc <- matrix(
    NA_real_,
    nrow = np_sim,
    ncol = H
  )

  y_fc <- matrix(
    NA_integer_,
    nrow = np_sim,
    ncol = H
  )

  t_prev <- time_last

  if (
    is.null(t_prev) ||
    is.na(t_prev)
  ) {
    t_prev <- -1L
  }

  t_prev <- as.integer(t_prev)

  for (h in seq_len(H)) {
    t_now <- as.integer(
      time_test[h]
    )

    dt <- if (t_prev < 0L) {
      0L
    } else {
      max(
        0L,
        t_now - t_prev
      )
    }

    th <- unpack_theta_int(
      theta_u,
      M_max
    )

    # Same OrderedRW propagation as the live filter.
    if (dt > 0L) {
      x <- x +
        sqrt(dt) *
        th$sigma_lat *
        rnorm(np_sim)
    }

    x_fc[, h] <- x

    # Calculate all K ordered-logistic category probabilities through
    # the same likelihood helper used during filtering.
    probs <- matrix(
      0,
      nrow = np_sim,
      ncol = K
    )

    for (yc in seq_len(K)) {
      probs[, yc] <- lik_intensity_ordlogit(
        y_cat = yc,
        x = x,
        sigma_meas = th$sigma_meas,
        delta = th$delta,
        M_max = M_max
      )
    }

    probs <- pmax(
      probs,
      1e-12
    )

    probs <- probs /
      rowSums(probs)

    # Sample one categorical observation for each particle.
    cdf <- t(
      apply(
        probs,
        1,
        cumsum
      )
    )

    u <- runif(np_sim)

    y_fc[, h] <- pmin(
      K,
      rowSums(
        cdf < u
      ) + 1L
    )

    t_prev <- t_now
  }

  list(
    x_forecast = x_fc,
    y_forecast = y_fc,
    theta_forecast_u = theta_u
  )
}


# ----------------------------------------------------------------------
# Log predictive density
# ----------------------------------------------------------------------

compute_lpd_intensity <- function(
    yc_true,
    x_particles,
    theta_u_mat,
    M_max
) {
  theta_u_mat <- as_matrix(
    theta_u_mat
  )

  if (nrow(theta_u_mat) != length(x_particles)) {
    stop(
      "theta_u_mat must have one row per x particle."
    )
  }

  th <- unpack_theta_int(
    theta_u_mat,
    M_max
  )

  pmf <- lik_intensity_ordlogit(
    y_cat = yc_true,
    x = x_particles,
    sigma_meas = th$sigma_meas,
    delta = th$delta,
    M_max = M_max
  )

  log_mean_exp(
    log(
      pmax(
        pmf,
        1e-12
      )
    )
  )
}