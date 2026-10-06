# Refactored from legacy SoloFast/patient-specific XemaPred code.
# Model mathematics intentionally preserved.

# -------------------------------------------------------------------
# Extent (binomial Markov chain) — MATCH EczemaPred Stan per patient
# -------------------------------------------------------------------

get_priors_extent <- function() {
  list(
    # Stan priors:
    # sigma ~ N+(0, (0.25 log(5))^2)
    sigma_loc = 0.0,
    sigma_scale = 0.25 * log(5),

    # SOLO (no pooling): logit(p10) ~ N(0, 1.5^2)
    logit_p10_loc = 0.0,
    logit_p10_scale = 1.5,

    # logit_tss1_0 ~ N(-1, 1)
    logit_tss1_0_loc = -1.0,
    logit_tss1_0_scale = 1.0
  )
}

# --- theta_u contains: logit_p10, log_sigma
init_theta_particles_extent <- function(np, priors) {
  logit_p10 <- rnorm(np, priors$logit_p10_loc, priors$logit_p10_scale)

  sigma <- abs(rnorm(np, priors$sigma_loc, priors$sigma_scale))
  sigma <- pmax(sigma, 1e-12)

  list(
    logit_p10 = logit_p10,
    sigma = sigma
  )
}

pack_theta_extent <- function(logit_p10, sigma) {
  cbind(
    logit_p10 = logit_p10,
    log_sigma = log(sigma)
  )
}


unpack_theta_extent <- function(theta_u) {
  list(
    logit_p10 = theta_u[, "logit_p10"],
    sigma = exp(theta_u[, "log_sigma"])
  )
}

# EczemaPred mapping: ss1 = inv_logit(logit_tss1)/(1+p10)
extent_derived <- function(logit_tss1, logit_p10) {
  p10 <- inv_logit(logit_p10)
  p11 <- 1 - p10
  ss1 <- inv_logit(logit_tss1) / (1 + p10)
  ss1 <- clip(ss1, 1e-12, 1 - 1e-12)

  p01 <- p10 * ss1 / (1 - ss1)
  p01 <- clip(p01, 1e-12, 1 - 1e-12)

  list(p10 = p10, p11 = p11, ss1 = ss1, p01 = p01)
}

pf_extent <- function(y_obs, day_obs, t_day, priors, np = 2000L, ess_frac = 0.5, a = 0.98) {
  # Particle parameters
  th0 <- init_theta_particles_extent(np, priors)
  theta_u <- pack_theta_extent(th0$logit_p10, th0$sigma)
  stopifnot(all(colnames(theta_u) == c("logit_p10", "log_sigma")))

  # Initial latent: logit_tss1_0 ~ N(-1, 1)  (Stan prior_logit_tss1_0)
  logit_tss1 <- rnorm(np, priors$logit_tss1_0_loc, priors$logit_tss1_0_scale)

  # logit_p10 per particle (SOLO / no pooling): logit(p10) ~ N(0, 1.5^2)
  th <- unpack_theta_extent(theta_u)
  logit_p10 <- th$logit_p10

  der <- extent_derived(logit_tss1, logit_p10)

  # Stan: y_lat[t0] = ss1[t0]
  y_lat <- der$ss1

  w <- rep(1 / np, np)

  obs_map <- tibble(day = as.integer(day_obs), y = as.integer(y_obs)) %>% distinct(day, .keep_all = TRUE)

  resample_move <- function(y_lat, logit_tss1, w, theta_u) {
    ess <- 1 / sum(w^2)
    if (ess < ess_frac * length(w)) {
      idx <- systematic_resample_idx(w)
      y_lat <- y_lat[idx]
      logit_tss1 <- logit_tss1[idx]
      theta_u <- theta_u[idx, , drop = FALSE]
      w <- rep(1 / length(w), length(w))
      theta_u <- liu_west_move(theta_u, a = a)
    }
    list(y_lat = y_lat, logit_tss1 = logit_tss1, w = w, theta_u = theta_u)
  }

  T_steps <- t_day + 1L
  y_store     <- vector("list", T_steps)
  tss1_store  <- vector("list", T_steps)
  w_store     <- vector("list", T_steps)
  theta_store <- vector("list", T_steps)

  for (d in 0:t_day) {
    if (d > 0) {
      th <- unpack_theta_extent(theta_u)

      # recompute logit_p10 after Liu–West move
      logit_p10 <- th$logit_p10

      # Stan RW: logit_tss1[t] = logit_tss1[t-1] + sigma * eta_t
      logit_tss1 <- logit_tss1 + th$sigma * rnorm(np)

      der <- extent_derived(logit_tss1, logit_p10)

      # Markov chain recursion (Stan y_lat):
      # y_lat[t] = p01[t]*(1-y_lat[t-1]) + p11*y_lat[t-1]
      y_lat <- der$p01 * (1 - y_lat) + der$p11 * y_lat
      y_lat <- clip(y_lat, 1e-12, 1 - 1e-12)
    }

    row <- obs_map %>% filter(day == d)
    if (nrow(row) == 1) {
      y <- as.integer(pmin(100L, pmax(0L, row$y[1])))
      ll   <- dbinom(y, size = 100L, prob = clip(y_lat, 1e-12, 1 - 1e-12), log = TRUE)
      logw <- log(pmax(w, 1e-12)) + ll
      w    <- log_normalize(logw)

      rr <- resample_move(y_lat, logit_tss1, w, theta_u)
      y_lat <- rr$y_lat
      logit_tss1 <- rr$logit_tss1
      w <- rr$w
      theta_u <- rr$theta_u
    }

    idx <- d + 1L
    y_store[[idx]]     <- y_lat
    tss1_store[[idx]]  <- logit_tss1
    w_store[[idx]]     <- w
    theta_store[[idx]] <- theta_u
  }

  # Return "x" as y_lat for compatibility with your runner naming
  list(x = y_store, logit_tss1 = tss1_store, w = w_store, theta_u = theta_store)
}

pf_forecast_extent <- function(x_last, logit_tss1_last, w_last, theta_last_u, day_last, day_test_vec, np_sim = 2000L) {
  H <- length(day_test_vec)
  idx <- sample.int(length(x_last), size = np_sim, replace = TRUE, prob = normalize_weights(w_last))


  y_lat <- x_last[idx]
  logit_tss1 <- logit_tss1_last[idx]
  theta_u <- theta_last_u[idx, , drop = FALSE]

  x_mat <- matrix(NA_real_, nrow = np_sim, ncol = H)
  y_mat <- matrix(NA_integer_, nrow = np_sim, ncol = H)

  day_prev <- as.integer(day_last)
  for (h in 1:H) {
    target_day <- as.integer(day_test_vec[h])
    dt <- as.integer(target_day - day_prev)
    if (dt > 0L) {
      for (k in 1:dt) {
        th <- unpack_theta_extent(theta_u)
        logit_p10 <- th$logit_p10

        logit_tss1 <- logit_tss1 + th$sigma * rnorm(np_sim)
        der <- extent_derived(logit_tss1, logit_p10)

        y_lat <- der$p01 * (1 - y_lat) + der$p11 * y_lat
        y_lat <- clip(y_lat, 1e-12, 1 - 1e-12)
      }
    }

    x_mat[, h] <- y_lat
    y_mat[, h] <- rbinom(np_sim, size = 100L, prob = y_lat)
    day_prev <- target_day
  }
  list(x_forecast = x_mat, y_forecast = y_mat, theta_u = theta_u, logit_tss1 = logit_tss1)
}

compute_lpd_extent <- function(y_true, x_particles) {
  y_true <- as.integer(pmin(100L, pmax(0L, y_true)))
  logpmf <- dbinom(y_true, size = 100L, prob = clip(x_particles, 1e-12, 1 - 1e-12), log = TRUE)
  log_mean_exp(logpmf)
}
