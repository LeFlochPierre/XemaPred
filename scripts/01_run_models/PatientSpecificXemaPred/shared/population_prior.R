# -----------------------------------------------------------------------
# Population-informed prior helpers for PatientSpecificXemaPred
#
# prior_mode = "population", prior_power = 1.0:
# use a Gaussian approximation to the FINAL external population-XemaPred
# posterior as the prior for the target patient's parameter vector.
#
# Family mapping:
# - intensity: population and patient-specific parameterisations match,
#   except simplex cutpoints are transformed to independent ALR coordinates.
# - subjective: population and patient-specific theta_u match directly.
# - extent: population theta_u = (log_sigma, mu10, log_sigma10), whereas
#   patient-specific theta_u = (logit_p10, log_sigma).  We therefore use
#   the population posterior predictive prior for a new patient's logit_p10:
#       logit_p10 | mu10,sigma10 ~ N(mu10, sigma10^2)
#   and moment-match the induced joint distribution with log_sigma.
# -----------------------------------------------------------------------

normalize_prior_weights <- function(w) {
  w <- as.numeric(w)
  w[!is.finite(w) | w < 0] <- 0
  s <- sum(w)
  if (!is.finite(s) || s <= 0) {
    return(rep(1 / length(w), length(w)))
  }
  w / s
}

weighted_mean_cov <- function(X, w) {
  X <- as.matrix(X)
  w <- normalize_prior_weights(w)

  if (nrow(X) != length(w)) {
    stop("[ERROR] Population-prior theta/weight dimension mismatch.")
  }

  mu <- colSums(sweep(X, 1, w, "*"))
  Xc <- sweep(X, 2, mu, "-")
  Sigma <- crossprod(Xc, Xc * w)

  list(mean = mu, cov = Sigma)
}

regularize_cov <- function(Sigma) {
  Sigma <- as.matrix(Sigma)

  if (nrow(Sigma) != ncol(Sigma)) {
    stop("[ERROR] Prior covariance is not square.")
  }

  Sigma <- (Sigma + t(Sigma)) / 2

  d <- diag(Sigma)
  positive_scale <- d[is.finite(d) & d > 0]
  scale <- if (length(positive_scale)) mean(positive_scale) else 1

  ridge <- max(1e-8, 1e-6 * scale)
  Sigma + diag(ridge, nrow(Sigma))
}

rmvn_prior <- function(n, mu, Sigma) {
  mu <- as.numeric(mu)
  names(mu) <- names(mu)

  Sigma <- regularize_cov(Sigma)
  R <- chol(Sigma)

  Z <- matrix(
    rnorm(as.integer(n) * length(mu)),
    nrow = as.integer(n),
    ncol = length(mu)
  )

  out <- sweep(Z %*% R, 2, mu, "+")
  colnames(out) <- names(mu)
  out
}

dmvn_prior_log <- function(X, mu, Sigma) {
  X <- as.matrix(X)
  mu <- as.numeric(mu)

  Sigma <- regularize_cov(Sigma)
  R <- chol(Sigma)

  Xc <- sweep(X, 2, mu, "-")
  whitened <- t(forwardsolve(t(R), t(Xc)))

  mahal <- rowSums(whitened^2)
  logdet <- 2 * sum(log(diag(R)))
  p <- length(mu)

  -0.5 * (p * log(2 * pi) + logdet + mahal)
}

# -----------------------------------------------------------------------
# Intensity coordinates
# -----------------------------------------------------------------------

theta_to_prior_z_intensity <- function(theta_u, M_max) {
  theta_u <- as_matrix(theta_u)

  delta_cols <- paste0(
    "log_delta_",
    seq_len(M_max - 1L)
  )

  missing_cols <- setdiff(
    c(
      "log_sigma_meas",
      "log_sigma_lat",
      "mu_y0",
      "log_sigma_y0",
      delta_cols
    ),
    colnames(theta_u)
  )

  if (length(missing_cols)) {
    stop(
      "[ERROR] Intensity population prior is missing theta columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  log_delta <- theta_u[, delta_cols, drop = FALSE]
  log_delta <- log_delta - rowMeans(log_delta)

  # Additive log-ratio coordinates using final simplex element as reference.
  n_delta <- ncol(log_delta)

  if (n_delta > 1L) {
    alr <- sweep(
      log_delta[, seq_len(n_delta - 1L), drop = FALSE],
      1,
      log_delta[, n_delta],
      "-"
    )
    colnames(alr) <- paste0("delta_alr_", seq_len(n_delta - 1L))
  } else {
    alr <- matrix(
      numeric(0),
      nrow = nrow(theta_u),
      ncol = 0L
    )
  }

  cbind(
    log_sigma_meas = theta_u[, "log_sigma_meas"],
    log_sigma_lat = theta_u[, "log_sigma_lat"],
    mu_y0 = theta_u[, "mu_y0"],
    log_sigma_y0 = theta_u[, "log_sigma_y0"],
    alr
  )
}

prior_z_to_theta_intensity <- function(z, M_max) {
  z <- as.matrix(z)

  n <- nrow(z)
  n_delta <- M_max - 1L
  n_alr <- n_delta - 1L

  if (n_alr > 0L) {
    alr_cols <- paste0("delta_alr_", seq_len(n_alr))

    eta <- cbind(
      z[, alr_cols, drop = FALSE],
      delta_reference = 0
    )

    eta <- eta - apply(eta, 1, max)
    delta_raw <- exp(eta)
    delta <- delta_raw / rowSums(delta_raw)
  } else {
    delta <- matrix(1, nrow = n, ncol = 1L)
  }

  log_delta <- log(pmax(delta, 1e-12))
  log_delta <- log_delta - rowMeans(log_delta)
  colnames(log_delta) <- paste0(
    "log_delta_",
    seq_len(n_delta)
  )

  out <- cbind(
    log_sigma_meas = z[, "log_sigma_meas"],
    log_sigma_lat = z[, "log_sigma_lat"],
    mu_y0 = z[, "mu_y0"],
    log_sigma_y0 = z[, "log_sigma_y0"],
    log_delta
  )

  as_matrix(out)
}

# -----------------------------------------------------------------------
# Convert saved population-XemaPred posterior bank object into the
# effective prior required by the patient-specific model.
# -----------------------------------------------------------------------

make_effective_population_prior <- function(
    raw_prior,
    model_family,
    M_max = NULL
) {
  theta <- as.matrix(raw_prior$theta_u)
  w <- normalize_prior_weights(raw_prior$weights)

  if (nrow(theta) != length(w)) {
    stop("[ERROR] Population-prior theta_u/weights mismatch.")
  }

  if (identical(model_family, "intensity")) {
    if (is.null(M_max)) {
      stop("[ERROR] M_max is required for intensity population prior.")
    }

    z <- theta_to_prior_z_intensity(theta, M_max)
    mom <- weighted_mean_cov(z, w)

    return(list(
      source_dataset = raw_prior$source_dataset,
      score = raw_prior$score,
      model_family = model_family,
      source_iteration = raw_prior$source_iteration,
      source_training_day = raw_prior$source_training_day,
      prior_power = 1.0,
      coordinate_names = colnames(z),
      mean = mom$mean,
      cov = regularize_cov(mom$cov)
    ))
  }

  if (identical(model_family, "subjective")) {
    required <- c(
      "log_sigma",
      "mu0",
      "log_sigma0"
    )

    missing_cols <- setdiff(required, colnames(theta))
    if (length(missing_cols)) {
      stop(
        "[ERROR] Subjective population prior missing columns: ",
        paste(missing_cols, collapse = ", ")
      )
    }

    z <- theta[, required, drop = FALSE]
    mom <- weighted_mean_cov(z, w)

    return(list(
      source_dataset = raw_prior$source_dataset,
      score = raw_prior$score,
      model_family = model_family,
      source_iteration = raw_prior$source_iteration,
      source_training_day = raw_prior$source_training_day,
      prior_power = 1.0,
      coordinate_names = required,
      mean = mom$mean,
      cov = regularize_cov(mom$cov)
    ))
  }

  if (identical(model_family, "extent")) {
    # Population extent parameterisation:
    #   log_sigma, mu10, log_sigma10
    #
    # Patient-specific parameterisation:
    #   logit_p10, log_sigma
    #
    # New-patient predictive prior:
    #   logit_p10 | mu10,sigma10 ~ N(mu10, sigma10^2)

    required <- c(
      "log_sigma",
      "mu10",
      "log_sigma10"
    )

    missing_cols <- setdiff(required, colnames(theta))
    if (length(missing_cols)) {
      stop(
        "[ERROR] Extent population prior missing columns: ",
        paste(missing_cols, collapse = ", ")
      )
    }

    log_sigma <- theta[, "log_sigma"]
    mu10 <- theta[, "mu10"]
    sigma10 <- exp(theta[, "log_sigma10"])

    mean_p10 <- sum(w * mu10)
    mean_log_sigma <- sum(w * log_sigma)

    var_p10 <- sum(
      w * (sigma10^2 + mu10^2)
    ) - mean_p10^2

    var_log_sigma <- sum(
      w * (log_sigma - mean_log_sigma)^2
    )

    cov_p10_log_sigma <- sum(
      w *
        (mu10 - mean_p10) *
        (log_sigma - mean_log_sigma)
    )

    mu <- c(
      logit_p10 = mean_p10,
      log_sigma = mean_log_sigma
    )

    Sigma <- matrix(
      c(
        var_p10,
        cov_p10_log_sigma,
        cov_p10_log_sigma,
        var_log_sigma
      ),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(names(mu), names(mu))
    )

    return(list(
      source_dataset = raw_prior$source_dataset,
      score = raw_prior$score,
      model_family = model_family,
      source_iteration = raw_prior$source_iteration,
      source_training_day = raw_prior$source_training_day,
      prior_power = 1.0,
      coordinate_names = names(mu),
      mean = mu,
      cov = regularize_cov(Sigma)
    ))
  }

  stop("[ERROR] Unknown model family: ", model_family)
}

# -----------------------------------------------------------------------
# Load + validate prior bank object for one target run.
# -----------------------------------------------------------------------

load_population_prior_for_run <- function(run_info) {

  # ---------------------------------------------------------------------
  # Base prior: nothing external to load
  # ---------------------------------------------------------------------

  if (identical(run_info$prior_mode, "base")) {
    return(NULL)
  }

  if (!identical(run_info$prior_mode, "population")) {
    stop(
      "[ERROR] Unknown prior_mode: ",
      run_info$prior_mode
    )
  }


  # ---------------------------------------------------------------------
  # Load population posterior used as the informative prior
  # ---------------------------------------------------------------------

  raw <- readRDS(
    run_info$population_prior_file
  )


  # ---------------------------------------------------------------------
  # Validate population-prior metadata
  # ---------------------------------------------------------------------

  if (
    !identical(
      as.character(raw$source_dataset),
      as.character(run_info$prior_source_dataset)
    )
  ) {

    stop(
      "[ERROR] Population-prior source dataset mismatch."
    )
  }


  if (
    !identical(
      as.character(raw$score),
      as.character(run_info$score)
    )
  ) {

    stop(
      "[ERROR] Population-prior score mismatch."
    )
  }


  if (
    !identical(
      as.character(raw$model_family),
      as.character(run_info$model_family)
    )
  ) {

    stop(
      "[ERROR] Population-prior model-family mismatch."
    )
  }


  # ---------------------------------------------------------------------
  # Convert population posterior into the effective patient-specific prior
  # ---------------------------------------------------------------------

  effective <- make_effective_population_prior(
    raw_prior = raw,
    model_family = run_info$model_family,
    M_max = run_info$M_max
  )

  effective$prior_combination <-
    run_info$prior_combination

  # ---------------------------------------------------------------------
  # Record provenance
  # ---------------------------------------------------------------------

  effective$source_fold <-
    run_info$prior_source_fold

  effective$source_file <-
    normalizePath(
      run_info$population_prior_file,
      mustWork = FALSE
    )


  # ---------------------------------------------------------------------
  # Unique prior ID
  #
  # For ordinary external-dataset priors:
  #
  #   Derexyl::dryness::iter20
  #
  # For cross-fold priors:
  #
  #   PFDC::dryness::fold1::iter20
  #   PFDC::dryness::fold2::iter20
  #
  # This prevents cached patient states fitted with Fold 1 from being
  # reused accidentally for a Fold 2 prior.
  # ---------------------------------------------------------------------

  prior_id_parts <- c(
    effective$source_dataset,
    effective$score
  )

  if (
    !is.null(run_info$prior_source_fold) &&
    !is.na(run_info$prior_source_fold)
  ) {

    prior_id_parts <- c(
      prior_id_parts,
      paste0(
        "fold",
        run_info$prior_source_fold
      )
    )
  }

  prior_id_parts <- c(
    prior_id_parts,
    paste0(
      "iter",
      effective$source_iteration
    )
  )

  effective$prior_id <-
    paste(
      prior_id_parts,
      collapse = "::"
    )


  # ---------------------------------------------------------------------
  # Logging
  # ---------------------------------------------------------------------

  cat(
    "[INFO] Population prior loaded:\n",
    "  Source dataset: ",
    effective$source_dataset,
    "\n",
    "  Source fold: ",
    if (
      !is.null(effective$source_fold) &&
      !is.na(effective$source_fold)
    ) {
      effective$source_fold
    } else {
      "NA"
    },
    "\n",
    "  Source file: ",
    effective$source_file,
    "\n",
    "  Prior ID: ",
    effective$prior_id,
    "\n",
    sep = ""
  )


  effective
}

# -----------------------------------------------------------------------
# Generic sampling / density in patient-specific theta_u coordinates.
# -----------------------------------------------------------------------

sample_population_prior_theta <- function(
    n,
    population_prior,
    model_family,
    M_max = NULL
) {
  z <- rmvn_prior(
    n = n,
    mu = population_prior$mean,
    Sigma = population_prior$cov
  )

  colnames(z) <- population_prior$coordinate_names

  if (identical(model_family, "intensity")) {
    return(
      prior_z_to_theta_intensity(
        z,
        M_max = M_max
      )
    )
  }

  if (identical(model_family, "subjective")) {
    return(
      as_matrix(
        z[, c(
          "log_sigma",
          "mu0",
          "log_sigma0"
        ), drop = FALSE]
      )
    )
  }

  if (identical(model_family, "extent")) {
    return(
      as_matrix(
        z[, c(
          "logit_p10",
          "log_sigma"
        ), drop = FALSE]
      )
    )
  }

  stop("[ERROR] Unknown model family: ", model_family)
}

log_population_prior_theta <- function(
    theta_u,
    population_prior,
    model_family,
    M_max = NULL
) {
  theta_u <- as_matrix(theta_u)

  if (identical(model_family, "intensity")) {
    z <- theta_to_prior_z_intensity(
      theta_u,
      M_max = M_max
    )
  } else if (identical(model_family, "subjective")) {
    z <- theta_u[, c(
      "log_sigma",
      "mu0",
      "log_sigma0"
    ), drop = FALSE]
  } else if (identical(model_family, "extent")) {
    z <- theta_u[, c(
      "logit_p10",
      "log_sigma"
    ), drop = FALSE]
  } else {
    stop("[ERROR] Unknown model family: ", model_family)
  }

  z <- z[
    ,
    population_prior$coordinate_names,
    drop = FALSE
  ]

  dmvn_prior_log(
    X = z,
    mu = population_prior$mean,
    Sigma = population_prior$cov
  )
}

# =======================================================================
# Base-prior sampling in patient-specific theta_u coordinates
# =======================================================================

sample_base_prior_theta <- function(
    n,
    model_family,
    priors,
    M_max = NULL
) {

  n <- as.integer(n)

  if (identical(model_family, "intensity")) {

    if (is.null(M_max)) {
      stop("[ERROR] M_max required for intensity base prior.")
    }

    th0 <- init_theta_particles_int(
      n,
      priors,
      M_max
    )

    return(
      pack_theta_int(
        delta = th0$delta,
        sigma_meas = th0$sigma_meas,
        sigma_lat = th0$sigma_lat,
        mu_y0 = th0$mu_y0,
        sigma_y0 = th0$sigma_y0,
        M_max = M_max
      )
    )
  }


  if (identical(model_family, "subjective")) {

    th0 <- init_theta_particles_subj(
      n,
      priors
    )

    return(
      pack_theta_subj(
        th0$sigma,
        th0$mu0,
        th0$sigma0
      )
    )
  }


  if (identical(model_family, "extent")) {

    th0 <- init_theta_particles_extent(
      n,
      priors
    )

    return(
      pack_theta_extent(
        th0$logit_p10,
        th0$sigma
      )
    )
  }


  stop(
    "[ERROR] Unknown model family: ",
    model_family
  )
}


# =======================================================================
# Robust mixture-prior sampler
#
# p_mix(theta) =
#   (1 - w) p_base(theta) +
#   w       p_population(theta)
#
# Draws EXACTLY from that mixture, therefore initial particle weights
# remain uniform.
# =======================================================================

sample_mixture_prior_theta <- function(
    n,
    population_weight,
    population_prior,
    model_family,
    priors,
    M_max = NULL
) {

  n <- as.integer(n)
  w <- as.numeric(population_weight)

  if (
    !is.finite(w) ||
    w < 0 ||
    w > 1
  ) {
    stop(
      "[ERROR] population_weight must satisfy 0 <= w <= 1."
    )
  }


  # Start with base-prior particles.
  theta <- sample_base_prior_theta(
    n = n,
    model_family = model_family,
    priors = priors,
    M_max = M_max
  )


  # Decide independently which particles come from population prior.
  use_population <- runif(n) < w

  n_population <- sum(use_population)


  if (n_population > 0L) {

    theta_population <- sample_population_prior_theta(
      n = n_population,
      population_prior = population_prior,
      model_family = model_family,
      M_max = M_max
    )


    if (
      !setequal(
        colnames(theta),
        colnames(theta_population)
      )
    ) {
      stop(
        "[ERROR] Base and population theta coordinates do not match."
      )
    }


    theta_population <- theta_population[
      ,
      colnames(theta),
      drop = FALSE
    ]


    theta[use_population, ] <- theta_population
  }


  theta
}


# =======================================================================
# Stable log-density of robust mixture prior
#
# log[
#   (1-w) p_base(theta)
#   + w p_population(theta)
# ]
# =======================================================================

log_mixture_prior <- function(
    lp_base,
    lp_population,
    population_weight
) {

  w <- as.numeric(population_weight)

  if (w <= 0) {
    return(lp_base)
  }

  if (w >= 1) {
    return(lp_population)
  }


  a <- log1p(-w) + lp_base
  b <- log(w) + lp_population

  m <- pmax(a, b)

  m +
    log(
      exp(a - m) +
        exp(b - m)
    )
}