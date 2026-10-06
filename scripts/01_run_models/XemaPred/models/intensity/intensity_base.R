# intensity_functions.R  (EczemaPred OrderedRW — PF approximation)

suppressPackageStartupMessages({ library(stats) })


# Writeup cutpoints: c0=0, deltas > 0, cutpoints are cumulative sums
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
  x <- matrix(rgamma(n * K, shape = rep(alpha, each = n), rate = 1),
              nrow = n, ncol = K)
  x / rowSums(x)
}

as_matrix <- function(x) {
  if (is.matrix(x)) return(x)
  x <- as.matrix(x)
  if (is.null(nrow(x))) x <- matrix(x, nrow = 1)
  x
}

# -------------------------------------------------------------------
# Priors: EXACT default_prior.OrderedRW
# -------------------------------------------------------------------
get_priors_intensity <- function(M_max) {
  list(
    prior_delta = rep(1, M_max - 1L),

    mu0_mean = 0.5,
    mu0_sd   = 0.25,
    sigma0_sd = 0.125,

    # Stan: sigma_lat / M ~ lognormal(prior_sigma_lat[1], prior_sigma_lat[2])
    prior_sigma_lat  = c(-log(10), 0.5 * log(4)),

    # Stan: sigma_meas / M ~ lognormal(prior_sigma_meas[1], prior_sigma_meas[2])
    prior_sigma_meas = c(-log(10), 0.5 * log(4))
  )
}



theta_colnames_int <- function(M_max) {
  c("log_sigma_meas", "log_sigma_lat", "mu_y0", "log_sigma_y0",
    paste0("log_delta_", seq_len(M_max - 1L)))
}

init_theta_particles_int <- function(np, priors, M_max) {

  # NEW: delta simplex like Stan
  delta <- rdirichlet_mat(np, priors$prior_delta)   # np x (M_max-1), rows sum to 1

  # sigma_meas (EczemaPred prior) — THIS GOES HERE
  sigma_meas <- M_max * rlnorm(np,
    meanlog = priors$prior_sigma_meas[1],
    sdlog   = priors$prior_sigma_meas[2]
  )
  sigma_meas <- pmax(sigma_meas, 1e-8)

  # sigma_lat — scale to latent space (see note below)
  sigma_lat <- M_max * rlnorm(np,
    meanlog = priors$prior_sigma_lat[1],
    sdlog   = priors$prior_sigma_lat[2]
  )
  sigma_lat <- pmax(sigma_lat, 1e-8)


  # mu_y0 — scale to latent space
  mu_y0 <- M_max * rnorm(np, priors$mu0_mean, priors$mu0_sd)

  # sigma_y0 — scale to latent space
  sigma_y0 <- abs(M_max * rnorm(np, 0, priors$sigma0_sd))
  sigma_y0 <- pmax(sigma_y0, 1e-8)

  list(delta = delta, sigma_meas = sigma_meas,
       sigma_lat = sigma_lat, mu_y0 = mu_y0, sigma_y0 = sigma_y0)
}


pack_theta_int <- function(delta, sigma_meas, sigma_lat, mu_y0, sigma_y0, M_max) {
  # Store delta on a centred unconstrained log scale.
  # unpack_theta_int() maps this back to the simplex.
  log_delta <- log(pmax(delta, 1e-12))
  log_delta <- log_delta - rowMeans(log_delta)

  out <- cbind(
    log_sigma_meas = log(pmax(sigma_meas, 1e-12)),
    log_sigma_lat  = log(pmax(sigma_lat,  1e-12)),
    mu_y0          = mu_y0,
    log_sigma_y0   = log(pmax(sigma_y0,  1e-12)),
    log_delta
  )

  colnames(out) <- theta_colnames_int(M_max)
  out
}

unpack_theta_int <- function(theta_u, M_max) {
  theta_u <- as_matrix(theta_u)
  dcols <- paste0("log_delta_", seq_len(M_max - 1L))

  # Numerically stable softmax back to simplex.
  log_delta <- theta_u[, dcols, drop = FALSE]
  log_delta <- log_delta - rowMeans(log_delta)
  log_delta <- log_delta - apply(log_delta, 1, max)

  delta_raw <- exp(log_delta)
  delta <- delta_raw / rowSums(delta_raw)

  list(
    sigma_meas = exp(theta_u[, "log_sigma_meas"]),
    sigma_lat  = exp(theta_u[, "log_sigma_lat"]),
    mu_y0      = theta_u[, "mu_y0"],
    sigma_y0   = exp(theta_u[, "log_sigma_y0"]),
    delta      = delta
  )
}

# -------------------------------------------------------------------
# Likelihood: ordered logistic like Stan
# y_cat is 1..(M_max+1)
# -------------------------------------------------------------------
lik_intensity_ordlogit <- function(y_cat, x, sigma_meas, delta, M_max) {
  y_cat <- as.integer(y_cat)
  K <- M_max + 1L
  if (y_cat < 1L || y_cat > K) stop("y_cat out of range: must be 1..K")

  np <- length(x)

  # --- coerce delta to matrix ---
  if (!is.matrix(delta)) delta <- as.matrix(delta)

  # Allow either:
  #   - delta is 1 x (M_max-1): global cutpoints -> replicate to np
  #   - delta is np x (M_max-1): particle-specific cutpoints (rare)
  if (ncol(delta) != (M_max - 1L)) stop("delta must have (M_max-1) columns")
  if (nrow(delta) == 1L && np > 1L) {
    delta <- delta[rep(1L, np), , drop = FALSE]
  } else if (nrow(delta) != np) {
    stop("delta must be 1 x (M_max-1) or np x (M_max-1)")
  }

  # sigma_meas: allow scalar or length np
  if (length(sigma_meas) == 1L) {
    s <- rep(pmax(as.numeric(sigma_meas), 1e-12) * sqrt(3) / pi, np)
  } else if (length(sigma_meas) == np) {
    s <- pmax(as.numeric(sigma_meas), 1e-12) * sqrt(3) / pi  # convert to logistic scale
  } else {
    stop("sigma_meas must be scalar or length np")
  }

  ct <- make_ct_mat(delta, M_max)  # now np x M_max

  z_lat <- x / s
  z_ct  <- ct / s
  Fk <- inv_logit(z_ct - z_lat)    # np x M_max

  if (y_cat == 1L) {
    p <- Fk[, 1L]
  } else if (y_cat == K) {
    p <- 1 - Fk[, M_max]
  } else {
    p <- Fk[, y_cat] - Fk[, y_cat - 1L]
  }

  pmax(p, 1e-12)
}


compute_lpd_intensity <- function(yc_true, x_particles, theta_u_mat, M_max) {
  theta_u_mat <- as_matrix(theta_u_mat)
  # theta_u_mat <- assert_theta_u_int(theta_u_mat, M_max, where = "compute_lpd_intensity")
  th <- unpack_theta_int(theta_u_mat, M_max)
  pmf <- lik_intensity_ordlogit(yc_true, x_particles, th$sigma_meas, th$delta, M_max)
  log_mean_exp(log(pmf))
}

