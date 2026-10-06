# Refactored from legacy SoloFast/patient-specific XemaPred code.
# Model mathematics intentionally preserved.

# extent_cache_functions.R
suppressPackageStartupMessages({ library(stats) })

pf_extent_init <- function(
    np,
    priors,
    ess_frac = 0.5,
    pmmh_nx = 128L,
    pmmh_window_n = Inf,
    proposal_sd = c(0.05, 0.10),
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
    #
    # p(theta) =
    #   (1-w) p_base(theta) +
    #   w     p_population(theta)
    # ----------------------------------------------------------

    if (identical(prior_combination, "mixture")) {

      theta_u <- sample_mixture_prior_theta(
        n = np,
        population_weight = prior_power,
        population_prior = population_prior,
        model_family = "extent",
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
        model_family = "extent"
      )

      if (prior_power < 1.0) {

        lp_base <- vapply(
          seq_len(np),
          function(i) {
            log_prior_theta_u_extent(
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
              model_family = "extent"
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

    th0 <- init_theta_particles_extent(
      np,
      priors
    )

    theta_u <- pack_theta_extent(
      th0$logit_p10,
      th0$sigma
    )
  }

  # ------------------------------------------------------------
  # Initialise latent state
  # ------------------------------------------------------------

  logit_tss1 <- rnorm(
    np,
    priors$logit_tss1_0_loc,
    priors$logit_tss1_0_scale
  )

  th <- unpack_theta_extent(theta_u)

  der <- extent_derived(
    logit_tss1,
    th$logit_p10
  )

  y_lat <- der$ss1

  list(
    y_lat = y_lat,
    logit_tss1 = logit_tss1,

    w = w_init,
    theta_u = as_matrix(theta_u),

    time_last = -1L,
    np = np,
    ess_frac = ess_frac,

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

log_prior_theta_u_extent <- function(theta_u, priors) {
  th <- unpack_theta_extent(theta_u)

  logit_p10 <- as.numeric(th$logit_p10[1])
  sigma     <- as.numeric(th$sigma[1])

  lp_logit_p10 <- dnorm(logit_p10, 0, 1.5, log = TRUE)
  lp_sigma     <- dnorm(sigma, 0, 0.25 * log(5), log = TRUE) + log(2) + as.numeric(theta_u[1, "log_sigma"])

  as.numeric(lp_logit_p10 + lp_sigma)
}

log_effective_prior_extent <- function(
    theta_u,
    state
) {

  # ------------------------------------------------------------
  # Original patient-specific/base prior
  # ------------------------------------------------------------

  lp_base <- log_prior_theta_u_extent(
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
      model_family = "extent"
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

propose_theta_u_extent <- function(theta_u_curr, proposal_sd) {
  th <- as.matrix(theta_u_curr)
  th[1, "logit_p10"] <- th[1, "logit_p10"] + rnorm(1, 0, proposal_sd[1])
  th[1, "log_sigma"] <- th[1, "log_sigma"] + rnorm(1, 0, proposal_sd[2])
  th
}

pf_extent_aux_init <- function(Nx, theta_j, priors) {
  th <- unpack_theta_extent(theta_j)

  logit_tss1 <- rnorm(Nx, priors$logit_tss1_0_loc, priors$logit_tss1_0_scale)
  der <- extent_derived(logit_tss1, th$logit_p10[1])
  y_lat <- clip(der$ss1, 1e-12, 1 - 1e-12)

  list(
    y_lat = y_lat,
    logit_tss1 = logit_tss1,
    w = rep(1 / Nx, Nx),
    time_last = -1L,
    Nx = as.integer(Nx)
  )
}

pf_extent_aux_step_one_day <- function(st, theta_j) {
  th <- unpack_theta_extent(theta_j)

  st$logit_tss1 <- st$logit_tss1 + th$sigma[1] * rnorm(st$Nx)
  der <- extent_derived(st$logit_tss1, th$logit_p10[1])

  st$y_lat <- der$p01 * (1 - st$y_lat) + der$p11 * st$y_lat
  st$y_lat <- clip(st$y_lat, 1e-12, 1 - 1e-12)

  st
}

pf_extent_aux_propagate <- function(st, theta_j, dt) {
  dt <- as.integer(dt)
  if (dt <= 0L) return(st)

  for (k in seq_len(dt)) {
    st <- pf_extent_aux_step_one_day(st, theta_j)
  }
  st
}

pf_extent_aux_update_one <- function(st, theta_j, y_t, t_now) {
  t_now <- as.integer(t_now)
  t_prev <- as.integer(st$time_last)

  dt <- if (t_prev < 0L) 0L else max(0L, t_now - t_prev)
  st <- pf_extent_aux_propagate(st, theta_j, dt)

  y_t <- as.integer(pmin(100L, pmax(0L, y_t)))
  ll <- dbinom(y_t, size = 100L, prob = clip(st$y_lat, 1e-12, 1 - 1e-12), log = TRUE)

  loglik_inc <- log_sum_exp(log(pmax(st$w, 1e-12)) + ll)
  st$w <- log_normalize(log(pmax(st$w, 1e-12)) + ll)

  ess <- 1 / sum(st$w^2)
  if (ess < 0.5 * st$Nx) {
    idx <- systematic_resample_idx(st$w)
    st$y_lat <- st$y_lat[idx]
    st$logit_tss1 <- st$logit_tss1[idx]
    st$w <- rep(1 / st$Nx, st$Nx)
  }

  st$time_last <- t_now

  list(st = st, loglik_inc = loglik_inc)
}

pf_extent_aux_update_many <- function(st, theta_j, y_vec, time_vec) {
  stopifnot(length(y_vec) == length(time_vec))
  if (!length(y_vec)) return(list(st = st, loglik_sum = 0.0))

  s <- 0.0
  for (i in seq_along(y_vec)) {
    rr <- pf_extent_aux_update_one(st, theta_j, y_t = y_vec[i], t_now = time_vec[i])
    st <- rr$st
    s <- s + rr$loglik_inc
  }

  list(st = st, loglik_sum = s)
}

get_pmmh_window_extent <- function(df, window_n = Inf) {
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

rebuild_pf_and_loglik_extent <- function(state, theta_u_one) {
  theta_u_one <- as_matrix(theta_u_one)

  st_aux <- pf_extent_aux_init(
    Nx = state$pmmh_nx,
    theta_j = theta_u_one,
    priors = state$priors
  )

  df <- get_pmmh_window_extent(
    state$train_cache,
    window_n = state$pmmh_window_n %||% Inf
  )
  if (!is.null(df) && nrow(df) > 0) {
    rr <- pf_extent_aux_update_many(
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
    y_lat_last = st_aux$y_lat[idx],
    logit_tss1_last = st_aux$logit_tss1[idx]
  )
}

pmmh_move_one_particle_extent <- function(state, j) {
  state$move_attempt[j] <- state$move_attempt[j] + 1L

  theta_curr <- state$theta_u[j, , drop = FALSE]

  curr <- rebuild_pf_and_loglik_extent(state, theta_curr)
  ll_curr <- curr$loglik_hat
  lp_curr <- log_effective_prior_extent(theta_curr, state)

  theta_prop <- propose_theta_u_extent(theta_curr, state$proposal_sd)

  prop <- rebuild_pf_and_loglik_extent(state, theta_prop)
  ll_prop <- prop$loglik_hat
  lp_prop <- log_effective_prior_extent(theta_prop, state)

  log_alpha <- (ll_prop + lp_prop) - (ll_curr + lp_curr)

  if (is.finite(log_alpha) && log(runif(1)) < log_alpha) {
    state$theta_u[j, ] <- theta_prop
    state$y_lat[j] <- prop$y_lat_last
    state$logit_tss1[j] <- prop$logit_tss1_last
    state$accept[j] <- state$accept[j] + 1L
  } else {
    state$y_lat[j] <- curr$y_lat_last
    state$logit_tss1[j] <- curr$logit_tss1_last
  }

  state
}

pf_extent_resample_move <- function(state) {
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
    state$y_lat <- state$y_lat[idx]
    state$logit_tss1 <- state$logit_tss1[idx]
    state$theta_u <- state$theta_u[idx, , drop = FALSE]
    state$w <- rep(1 / np, np)

    if (!is.null(state$train_cache) && nrow(state$train_cache) > 0) {
      move_idx <- sample.int(np, size = min(50L, np), replace = FALSE)
      state$n_rejuvenated <- length(move_idx)                              
      state$n_moves_cum   <- (state$n_moves_cum %||% 0L) + length(move_idx) 
      for (j in move_idx) {
        state <- pmmh_move_one_particle_extent(state, j)
      }
    }
  }

  state
}

pf_extent_step_one_day <- function(state) {
  th <- unpack_theta_extent(state$theta_u)
  logit_p10 <- th$logit_p10

  state$logit_tss1 <- state$logit_tss1 + th$sigma * rnorm(length(state$y_lat))
  der <- extent_derived(state$logit_tss1, logit_p10)

  state$y_lat <- der$p01 * (1 - state$y_lat) + der$p11 * state$y_lat
  state$y_lat <- clip(state$y_lat, 1e-12, 1 - 1e-12)

  state
}

pf_extent_propagate <- function(state, dt) {
  dt <- as.integer(dt)
  if (dt <= 0L) return(state)
  for (k in seq_len(dt)) state <- pf_extent_step_one_day(state)
  state
}

pf_extent_update_one <- function(state, y_t, dt) {
  # propagate forward dt days
  state <- pf_extent_propagate(state, dt)

  # weight update on observation y_t (0..100)
  y_t <- as.integer(pmin(100L, pmax(0L, y_t)))
  ll <- dbinom(y_t, size = 100L, prob = clip(state$y_lat, 1e-12, 1 - 1e-12), log = TRUE)
  state$w <- log_normalize(log(pmax(state$w, 1e-12)) + ll)

  # resample + rejuvenate if needed
  state <- pf_extent_resample_move(state)

  state
}

pf_extent_update_many <- function(state, y_vec, day_vec) {
  stopifnot(length(y_vec) == length(day_vec))
  if (!length(y_vec)) return(state)

  day_vec <- as.integer(day_vec)
  t_prev <- state$time_last

  new_df <- data.frame(TimeRel = day_vec, y = y_vec)
  if (is.null(state$train_cache) || nrow(state$train_cache) == 0) {
    state$train_cache <- new_df
  } else {
    state$train_cache <- dplyr::bind_rows(state$train_cache, new_df) %>%
      dplyr::distinct(TimeRel, .keep_all = TRUE) %>%
      dplyr::arrange(TimeRel)
  }

  for (i in seq_along(y_vec)) {
    t_now <- day_vec[i]
    dt <- max(0L, t_now - t_prev)
    state <- pf_extent_update_one(state, y_t = y_vec[i], dt = dt)
    state$time_last <- t_now  # IMPORTANT
    t_prev <- t_now
  }
  state
}

pf_forecast_extent_from_state <- function(state, day_test_vec, np_sim = 2000L) {
  H <- length(day_test_vec)
  idx <- sample.int(length(state$y_lat), size = np_sim, replace = TRUE,
                    prob = normalize_weights(state$w))

  y_lat     <- state$y_lat[idx]
  logit_tss1 <- state$logit_tss1[idx]
  theta_u   <- state$theta_u[idx, , drop = FALSE]

  x_mat <- matrix(NA_real_, nrow = np_sim, ncol = H)
  y_mat <- matrix(NA_integer_, nrow = np_sim, ncol = H)

  day_prev <- as.integer(state$time_last)

  for (h in seq_len(H)) {
    target_day <- as.integer(day_test_vec[h])
    dt <- max(0L, target_day - day_prev)

    if (dt > 0L) {
      for (k in seq_len(dt)) {
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

  list(x_forecast = x_mat, y_forecast = y_mat)
}
