# Refactored from legacy SoloFast/patient-specific XemaPred code.
# Model mathematics intentionally preserved.

# subjective_cache_functions.R
suppressPackageStartupMessages({ library(stats) })

# State contains:
#   x        : latent particles (logit p)
#   w        : weights
#   theta_u  : parameter particles (log_sigma, mu0, log_sigma0)
#   time_last: last assimilated time (integer)
#   np, ess_frac, M, priors
#   train_cache, pmmh_nx, proposal_sd, accept, move_attempt

pf_subjective_init <- function(
    np,
    priors,
    ess_frac = 0.5,
    M = 100L,
    pmmh_nx = 128L,
    pmmh_window_n = Inf,
    proposal_sd = c(0.03, 0.10, 0.03),
    prior_mode = "base",
    prior_combination = "power",
    prior_power = 0,
    population_prior = NULL,
    pmmh_prior_mode = "matched"
) {

  # ------------------------------------------------------------
  # Check prior configuration
  # ------------------------------------------------------------

  if (!pmmh_prior_mode %in% c("matched", "base")) {
    stop(
      "[ERROR] pmmh_prior_mode must be ",
      "'matched' or 'base'."
    )
  }

  if (!prior_mode %in% c("base", "population")) {
    stop(
      "[ERROR] prior_mode must be 'base' or 'population'."
    )
  }

  if (identical(prior_mode, "base")) {

    prior_combination <- "base"

  } else {

    if (is.null(population_prior)) {
      stop(
        "[ERROR] population_prior is required ",
        "when prior_mode = 'population'."
      )
    }

    if (!prior_combination %in% c("power", "mixture")) {
      stop(
        "[ERROR] prior_combination must be ",
        "'power' or 'mixture' when prior_mode = 'population'."
      )
    }

    if (identical(prior_combination, "power")) {

      if (
        !is.finite(prior_power) ||
        prior_power <= 0 ||
        prior_power > 1
      ) {
        stop(
          "[ERROR] Power population prior requires ",
          "0 < prior_power <= 1."
        )
      }

    } else if (identical(prior_combination, "mixture")) {

      if (
        !is.finite(prior_power) ||
        prior_power < 0 ||
        prior_power > 1
      ) {
        stop(
          "[ERROR] Mixture population prior requires ",
          "0 <= prior_power <= 1."
        )
      }
    }
  }

  # ------------------------------------------------------------
  # Initialise parameter particles
  # ------------------------------------------------------------

  w_init <- rep(1 / np, np)

  if (identical(prior_mode, "population")) {

    # ----------------------------------------------------------
    # Robust arithmetic mixture
    # ----------------------------------------------------------

    if (identical(prior_combination, "mixture")) {

      theta_u <- sample_mixture_prior_theta(
        n = np,
        population_weight = prior_power,
        population_prior = population_prior,
        model_family = "subjective",
        priors = priors
      )

      # Exact mixture sampling -> uniform initial weights.
      w_init <- rep(1 / np, np)

    } else {

      # --------------------------------------------------------
      # Original geometric / power prior
      # --------------------------------------------------------

      theta_u <- sample_population_prior_theta(
        n = np,
        population_prior = population_prior,
        model_family = "subjective"
      )

      if (prior_power < 1.0) {

        lp_base <- vapply(
          seq_len(np),
          function(i) {
            log_prior_theta_u_subj(
              theta_u[i, , drop = FALSE],
              priors
            )
          },
          numeric(1)
        )

        lp_pop <- vapply(
          seq_len(np),
          function(i) {
            log_population_prior_theta(
              theta_u = theta_u[i, , drop = FALSE],
              population_prior = population_prior,
              model_family = "subjective"
            )
          },
          numeric(1)
        )

        logw_init <- (1 - prior_power) * (
          lp_base - lp_pop
        )

        w_init <- log_normalize(logw_init)
      }
    }

  } else {

    # ----------------------------------------------------------
    # Original base prior
    # ----------------------------------------------------------

    th0 <- init_theta_particles_subj(
      np,
      priors
    )

    theta_u <- pack_theta_subj(
      th0$sigma,
      th0$mu0,
      th0$sigma0
    )
  }

  # ------------------------------------------------------------
  # Initialise latent state
  # ------------------------------------------------------------

  th <- unpack_theta_subj(theta_u)

  x <- rnorm(
    np,
    mean = th$mu0,
    sd = th$sigma0
  )

  list(
    x = x,
    w = w_init,
    theta_u = as_matrix(theta_u),

    time_last = -1L,
    np = np,
    ess_frac = ess_frac,
    M = as.integer(M),

    priors = priors,

    prior_mode = prior_mode,
    prior_combination = prior_combination,
    prior_power = prior_power,
    population_prior = population_prior,
    prior_id = population_prior$prior_id %||% "base",
    pmmh_prior_mode = pmmh_prior_mode,

    train_cache = data.frame(
      TimeRel = integer(0),
      y = integer(0)
    ),

    pmmh_nx = as.integer(pmmh_nx),

    pmmh_window_n = if (
      is.finite(pmmh_window_n)
    ) {
      as.integer(pmmh_window_n)
    } else {
      Inf
    },

    proposal_sd = proposal_sd,
    accept = integer(np),
    move_attempt = integer(np)
  )
}

log_prior_theta_u_subj <- function(theta_u, priors) {
  th <- unpack_theta_subj(theta_u)

  sigma  <- as.numeric(th$sigma[1])
  mu0    <- as.numeric(th$mu0[1])
  sigma0 <- as.numeric(th$sigma0[1])

  lp_sigma  <- dnorm(sigma, 0, priors$sigma_sd, log = TRUE) + log(2) + as.numeric(theta_u[1, "log_sigma"])
  lp_sigma0 <- dnorm(sigma0, 0, priors$sigma0_sd, log = TRUE) + log(2) + as.numeric(theta_u[1, "log_sigma0"])
  lp_mu0    <- dnorm(mu0, priors$mu0_mean, priors$mu0_sd, log = TRUE)

  as.numeric(lp_sigma + lp_sigma0 + lp_mu0)
}

log_effective_prior_subjective <- function(
    theta_u,
    state
) {

  # ------------------------------------------------------------
  # Original patient-specific/base prior
  # ------------------------------------------------------------

  lp_base <- log_prior_theta_u_subj(
    theta_u,
    state$priors
  )

  # ------------------------------------------------------------
  # Warm-start diagnostic:
  # population initialization + base-prior PMMH
  # ------------------------------------------------------------

  if (identical(
    state$pmmh_prior_mode,
    "base"
  )) {
    return(lp_base)
  }

  # ------------------------------------------------------------
  # Matched population-informed prior
  # ------------------------------------------------------------

  if (identical(
    state$prior_mode,
    "population"
  )) {

    lp_pop <- log_population_prior_theta(
      theta_u = theta_u,
      population_prior = state$population_prior,
      model_family = "subjective"
    )

    alpha <- state$prior_power

    # Robust arithmetic mixture.
    if (identical(
      state$prior_combination %||% "power",
      "mixture"
    )) {

      return(
        log_mixture_prior(
          lp_base = lp_base,
          lp_population = lp_pop,
          population_weight = alpha
        )
      )
    }

    # Original geometric / power prior.
    return(
      (1 - alpha) * lp_base +
        alpha * lp_pop
    )
  }

  lp_base
}

propose_theta_u_subj <- function(theta_u_curr, proposal_sd) {
  th <- as.matrix(theta_u_curr)
  th[1, "log_sigma"]  <- th[1, "log_sigma"]  + rnorm(1, 0, proposal_sd[1])
  th[1, "mu0"]        <- th[1, "mu0"]        + rnorm(1, 0, proposal_sd[2])
  th[1, "log_sigma0"] <- th[1, "log_sigma0"] + rnorm(1, 0, proposal_sd[3])
  th
}

pf_subjective_aux_propagate <- function(st, theta_j, dt) {
  dt <- as.integer(dt)
  if (dt <= 0L) return(st)

  th <- unpack_theta_subj(theta_j)
  st$x <- rnorm(st$Nx, mean = st$x, sd = th$sigma[1] * sqrt(dt))
  st
}

pf_subjective_aux_init <- function(Nx, theta_j, M) {
  th <- unpack_theta_subj(theta_j)

  list(
    x = rnorm(Nx, mean = th$mu0[1], sd = th$sigma0[1]),
    w = rep(1 / Nx, Nx),
    time_last = -1L,
    Nx = as.integer(Nx),
    M = as.integer(M)
  )
}

pf_subjective_aux_update_one <- function(st, theta_j, y_t, t_now) {
  t_now <- as.integer(t_now)
  t_prev <- as.integer(st$time_last)

  dt <- if (t_prev < 0L) 0L else max(0L, t_now - t_prev)
  st <- pf_subjective_aux_propagate(st, theta_j, dt)

  y_t <- as.integer(pmin(st$M, pmax(0L, y_t)))
  ll <- loglik_subj_binom(y_t, st$x, M = st$M)

  loglik_inc <- log_sum_exp(log(pmax(st$w, 1e-12)) + ll)
  st$w <- log_normalize(log(pmax(st$w, 1e-12)) + ll)

  ess <- 1 / sum(st$w^2)
  if (ess < 0.5 * st$Nx) {
    idx <- systematic_resample_idx(st$w)
    st$x <- st$x[idx]
    st$w <- rep(1 / st$Nx, st$Nx)
  }

  st$time_last <- t_now

  list(st = st, loglik_inc = loglik_inc)
}

pf_subjective_resample_move <- function(state) {
  w <- state$w
  np <- length(w)
  ess <- 1 / sum(w^2)

  state$n_check  <- (state$n_check %||% 0L) + 1L          
  state$ess_last <- ess                                   
  state$ess_min  <- min(state$ess_min %||% Inf, ess)       
  state$ess_sum  <- (state$ess_sum %||% 0) + ess           

  if (ess < state$ess_frac * np) {
    state$n_resample <- (state$n_resample %||% 0L) + 1L    
    idx <- systematic_resample_idx(w)
    state$x <- state$x[idx]
    state$theta_u <- state$theta_u[idx, , drop = FALSE]
    state$w <- rep(1 / np, np)

    if (!is.null(state$train_cache) && nrow(state$train_cache) > 0) {
      move_idx <- sample.int(np, size = min(50L, np), replace = FALSE)
      state$n_rejuvenated <- length(move_idx)                              
      state$n_moves_cum   <- (state$n_moves_cum %||% 0L) + length(move_idx) 
      for (j in move_idx) {
        state <- pmmh_move_one_particle_subjective(state, j)
      }
    }
  }

  state
}

pf_subjective_aux_update_many <- function(st, theta_j, y_vec, time_vec) {
  stopifnot(length(y_vec) == length(time_vec))
  if (!length(y_vec)) return(list(st = st, loglik_sum = 0.0))

  s <- 0.0
  for (i in seq_along(y_vec)) {
    rr <- pf_subjective_aux_update_one(st, theta_j, y_t = y_vec[i], t_now = time_vec[i])
    st <- rr$st
    s <- s + rr$loglik_inc
  }

  list(st = st, loglik_sum = s)
}

get_pmmh_window_subjective <- function(df, window_n = Inf) {
  if (is.null(df) || nrow(df) == 0) return(df)

  df <- df[order(df$TimeRel), , drop = FALSE]

  if (is.finite(window_n)) {
    n <- as.integer(window_n)
    if (nrow(df) > n) {
      df <- df[(nrow(df) - n + 1L):nrow(df), , drop = FALSE]
    }
  }

  df
}

rebuild_pf_and_loglik_subjective <- function(state, theta_u_one) {
  theta_u_one <- as_matrix(theta_u_one)

  st_aux <- pf_subjective_aux_init(
    Nx = state$pmmh_nx,
    theta_j = theta_u_one,
    M = state$M
  )

  df <- get_pmmh_window_subjective(
    state$train_cache,
    window_n = state$pmmh_window_n %||% Inf
  )
  if (!is.null(df) && nrow(df) > 0) {
    rr <- pf_subjective_aux_update_many(
      st = st_aux,
      theta_j = theta_u_one,
      y_vec = df$y,
      time_vec = df$TimeRel
    )
    st_aux <- rr$st
    loglik_hat <- rr$loglik_sum
  } else {
    loglik_hat <- 0.0
  }

  idx <- sample.int(st_aux$Nx, size = 1L, prob = normalize_weights(st_aux$w))

  list(
    loglik_hat = loglik_hat,
    x_last = st_aux$x[idx]
  )
}

pmmh_move_one_particle_subjective <- function(state, j) {
  state$move_attempt[j] <- state$move_attempt[j] + 1L

  theta_curr <- state$theta_u[j, , drop = FALSE]

  curr <- rebuild_pf_and_loglik_subjective(state, theta_curr)
  ll_curr <- curr$loglik_hat
  lp_curr <- log_effective_prior_subjective(theta_curr, state)

  theta_prop <- propose_theta_u_subj(theta_curr, state$proposal_sd)

  prop <- rebuild_pf_and_loglik_subjective(state, theta_prop)
  ll_prop <- prop$loglik_hat
  lp_prop <- log_effective_prior_subjective(theta_prop, state)

  log_alpha <- (ll_prop + lp_prop) - (ll_curr + lp_curr)

  if (is.finite(log_alpha) && log(runif(1)) < log_alpha) {
    state$theta_u[j, ] <- theta_prop
    state$x[j] <- prop$x_last
    state$accept[j] <- state$accept[j] + 1L
  } else {
    state$x[j] <- curr$x_last
  }

  state
}


pf_subjective_propagate <- function(state, dt) {
  dt <- as.integer(dt)
  if (dt <= 0L) return(state)

  th <- unpack_theta_subj(state$theta_u)
  # RW on logit-scale, scale with sqrt(dt) like your original
  state$x <- rnorm(length(state$x), mean = state$x, sd = th$sigma * sqrt(dt))
  state
}

pf_subjective_update_one <- function(state, y_t, dt) {
  state <- pf_subjective_propagate(state, dt)

  y_t <- as.integer(pmin(state$M, pmax(0L, y_t)))
  ll <- loglik_subj_binom(y_t, state$x, M = state$M)
  state$w <- log_normalize(log(pmax(state$w, 1e-12)) + ll)

  state <- pf_subjective_resample_move(state)
  state
}

pf_subjective_update_many <- function(state, y_vec, time_vec) {
  stopifnot(length(y_vec) == length(time_vec))
  if (!length(y_vec)) return(state)

  time_vec <- as.integer(time_vec)
  y_vec    <- as.integer(y_vec)
  t_prev   <- as.integer(state$time_last)

  new_df <- data.frame(TimeRel = time_vec, y = y_vec)
  if (is.null(state$train_cache) || nrow(state$train_cache) == 0) {
    state$train_cache <- new_df
  } else {
    state$train_cache <- dplyr::bind_rows(state$train_cache, new_df) %>%
      dplyr::distinct(TimeRel, .keep_all = TRUE) %>%
      dplyr::arrange(TimeRel)
  }

  for (i in seq_along(y_vec)) {
    t_now <- time_vec[i]
    dt <- max(0L, t_now - t_prev)
    state <- pf_subjective_update_one(state, y_t = y_vec[i], dt = dt)
    state$time_last <- t_now
    t_prev <- t_now
  }

  state
}

pf_forecast_subjective_from_state <- function(state, time_test, np_sim = 2000L) {
  H <- length(time_test)

  idx <- sample.int(length(state$x), size = np_sim, replace = TRUE,
                    prob = normalize_weights(state$w))

  x <- state$x[idx]
  theta_u <- state$theta_u[idx, , drop = FALSE]

  x_fc <- matrix(NA_real_, nrow = np_sim, ncol = H)
  y_fc <- matrix(NA_integer_, nrow = np_sim, ncol = H)

  t_prev <- as.integer(state$time_last)

  for (h in seq_len(H)) {
    t_now <- as.integer(time_test[h])
    dt <- max(0L, t_now - t_prev)

    if (dt > 0L) {
      th <- unpack_theta_subj(theta_u)
      x <- rnorm(np_sim, mean = x, sd = th$sigma * sqrt(dt))
    }

    x_fc[, h] <- x
    p <- inv_logit(x)
    y_fc[, h] <- rbinom(np_sim, size = state$M, prob = p)

    t_prev <- t_now
  }

  list(x_forecast = x_fc, y_forecast = y_fc, theta_forecast_u = theta_u)
}
