# ============================================================
# subjective_functions.R
# (EczemaPred BinRW — PF approximation)
#
# Matches:
#   y(t) ~ Binomial(M, p(t))
#   x(t) = logit(p(t))
#   x(t+1) ~ Normal(x(t), sigma^2)
#
# Priors:
#   sigma  ~ N+(0, (0.25*log(5))^2)
#   mu0    ~ N(0, 1)
#   sigma0 ~ N+(0, 1.5^2)
#   x(t0)  ~ N(mu0, sigma0^2)
#
# Notes:
# - We carry parameters (sigma, mu0, sigma0) as particles with Liu–West.
# - Latent state is x(t) = logit(p(t)).
# - Supports irregular time gaps via sqrt(dt) scaling.
#
# Depends on shared utilities:
#   inv_logit, log_normalize, log_mean_exp, normalize_weights,
#   systematic_resample_idx, liu_west_move
# ============================================================

suppressPackageStartupMessages({ library(stats) })

# ------------------------------------------------------------
# Priors (EXACT from writeup)
# ------------------------------------------------------------
get_priors_subjective <- function() {
  list(
    mu0_mean = 0,
    mu0_sd = 1,
    sigma0_sd = 1.5,
    sigma_sd = 0.25 * log(5)
  )
}

# half-normal sampler
rhalfnorm <- function(n, sd) abs(rnorm(n, mean = 0, sd = sd))

# ------------------------------------------------------------
# Theta: sigma, mu0, sigma0
# pack into unconstrained space:
#   log_sigma, mu0, log_sigma0
# ------------------------------------------------------------
init_theta_particles_subj <- function(np, priors) {
  sigma  <- pmax(rhalfnorm(np, priors$sigma_sd),  1e-8)
  mu0    <- rnorm(np, mean = priors$mu0_mean, sd = priors$mu0_sd)
  sigma0 <- pmax(rhalfnorm(np, priors$sigma0_sd), 1e-8)
  list(sigma = sigma, mu0 = mu0, sigma0 = sigma0)
}

pack_theta_subj <- function(sigma, mu0, sigma0) {
  out <- cbind(
    log_sigma  = log(pmax(sigma,  1e-12)),
    mu0        = mu0,
    log_sigma0 = log(pmax(sigma0, 1e-12))
  )
  colnames(out) <- c("log_sigma", "mu0", "log_sigma0")
  out
}

unpack_theta_subj <- function(theta_u) {
  theta_u <- as.matrix(theta_u)
  list(
    sigma  = exp(theta_u[, "log_sigma"]),
    mu0    = theta_u[, "mu0"],
    sigma0 = exp(theta_u[, "log_sigma0"])
  )
}

# ------------------------------------------------------------
# Likelihood: y ~ Bin(M, inv_logit(x))
# ------------------------------------------------------------
lik_subj_binom <- function(y, x, M = 100L) {
  p <- inv_logit(x)
  lik <- dbinom(y, size = M, prob = p)
  pmax(lik, 1e-12)
}

loglik_subj_binom <- function(y, x, M = 100L) {
  p <- inv_logit(x)
  dbinom(y, size = M, prob = p, log = TRUE)
}

# ------------------------------------------------------------
# PF (BinRW)
# ------------------------------------------------------------
pf_subjective <- function(y_train, time_train, priors,
                          np = 2000, ess_frac = 0.5, a = 0.98, M = 100L) {
  T <- length(y_train)
  stopifnot(length(time_train) == T)

  # parameter particles
  th0 <- init_theta_particles_subj(np, priors)
  theta_u <- pack_theta_subj(th0$sigma, th0$mu0, th0$sigma0)

  # latent initial: x(t0) ~ N(mu0, sigma0^2)
  th <- unpack_theta_subj(theta_u)
  x <- rnorm(np, mean = th$mu0, sd = th$sigma0)

  x_store <- vector("list", T)
  w_store <- vector("list", T)
  theta_store <- vector("list", T)

  w <- rep(1 / np, np)

  resample_move <- function(x, w, theta_u) {
    ess <- 1 / sum(w^2)
    if (ess < ess_frac * length(w)) {
      idx <- systematic_resample_idx(w)
      x <- x[idx]
      theta_u <- theta_u[idx, , drop = FALSE]
      w <- rep(1 / length(w), length(w))
      theta_u <- liu_west_move(theta_u, a = a)
    }
    list(x = x, w = w, theta_u = theta_u)
  }

  for (t in seq_len(T)) {
    if (t > 1L) {
      dt <- as.integer(time_train[t] - time_train[t - 1L])
      if (dt > 0L) {
        th <- unpack_theta_subj(theta_u)
        x <- rnorm(np, mean = x, sd = th$sigma * sqrt(dt))
      }
    }

    ll <- loglik_subj_binom(as.integer(y_train[t]), x, M = M)
    logw <- log(pmax(w, 1e-12)) + ll
    w <- log_normalize(logw)

    rr <- resample_move(x, w, theta_u)
    x <- rr$x; w <- rr$w; theta_u <- rr$theta_u

    x_store[[t]] <- x
    w_store[[t]] <- w
    theta_store[[t]] <- theta_u
  }

  list(x = x_store, w = w_store, theta_u = theta_store)
}

# ------------------------------------------------------------
# Forecast
# ------------------------------------------------------------
pf_forecast_subjective <- function(x_last, w_last, theta_last_u,
                                   time_last, time_test,
                                   np_sim = 2000, M = 100L) {
  H <- length(time_test)

  idx <- sample.int(length(x_last), size = np_sim, replace = TRUE,
                    prob = normalize_weights(w_last))
  x <- x_last[idx]
  theta_u <- theta_last_u[idx, , drop = FALSE]

  x_fc <- matrix(NA_real_, nrow = np_sim, ncol = H)
  y_fc <- matrix(NA_integer_, nrow = np_sim, ncol = H)

  t_prev <- as.integer(time_last)

  for (h in seq_len(H)) {
    t_now <- as.integer(time_test[h])
    dt <- as.integer(t_now - t_prev)

    if (dt > 0L) {
      th <- unpack_theta_subj(theta_u)
      x <- rnorm(np_sim, mean = x, sd = th$sigma * sqrt(dt))
    }

    x_fc[, h] <- x
    p <- inv_logit(x)
    y_fc[, h] <- rbinom(np_sim, size = M, prob = p)

    t_prev <- t_now
  }

  list(x_forecast = x_fc, y_forecast = y_fc, theta_forecast_u = theta_u)
}

# ------------------------------------------------------------
# LPD
# ------------------------------------------------------------
compute_lpd_subjective <- function(y_true, x_particles, M = 100L) {
  ll <- loglik_subj_binom(as.integer(y_true), x_particles, M = M)
  log_mean_exp(ll)
}
