# -------------------------------------------------------------------
# Extent (binomial Markov chain)
# -------------------------------------------------------------------

get_priors_extent <- function() {
  list(
    # Stan priors:
    # sigma ~ N+(0, (0.25 log(5))^2)
    sigma_loc = 0.0,
    sigma_scale = 0.25 * log(5),

    # logit(p10) ~ N(0, 1.5^2)
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

compute_lpd_extent <- function(y_true, x_particles) {
  y_true <- as.integer(pmin(100L, pmax(0L, y_true)))
  logpmf <- dbinom(y_true, size = 100L, prob = clip(x_particles, 1e-12, 1 - 1e-12), log = TRUE)
  log_mean_exp(logpmf)
}
