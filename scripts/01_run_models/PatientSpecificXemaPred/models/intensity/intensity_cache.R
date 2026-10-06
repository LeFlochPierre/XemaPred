# Refactored from legacy SoloFast/patient-specific XemaPred code.
# Model mathematics intentionally preserved.

# intensity_cache_functions.R
suppressPackageStartupMessages({ library(stats) })

# assumes these already exist (from solofast_functions.R + intensity_functions.R):
# - as_matrix(), unpack_theta_int(), init_theta_particles_int(), pack_theta_int()
# - lik_intensity_ordlogit(), log_normalize(), systematic_resample_idx()
# - normalize_weights(), log_sum_exp()
# - get_priors_intensity()

pf_intensity_init <- function(
    np,
    priors,
    M_max,
    ess_frac = 0.5,
    pmmh_nx = 128L,
    pmmh_window_n = Inf,
    proposal_sd = c(0.074, 0.074, 0.25, 0.074, 0.074),
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
    # Robust mixture prior:
    #
    # p(theta) =
    #   (1 - w) p_base(theta) +
    #   w       p_population(theta)
    #
    # Here prior_power is interpreted as w.
    # ----------------------------------------------------------

    if (identical(prior_combination, "mixture")) {

      theta_u <- sample_mixture_prior_theta(
        n = np,
        population_weight = prior_power,
        population_prior = population_prior,
        model_family = "intensity",
        priors = priors,
        M_max = M_max
      )

      # Exact sampling from the mixture -> uniform initial weights.
      w_init <- rep(1 / np, np)

    } else {

      # --------------------------------------------------------
      # Original geometric / power prior:
      #
      # p(theta) proportional to
      # p_base(theta)^(1-alpha) *
      # p_population(theta)^alpha
      # --------------------------------------------------------

      theta_u <- sample_population_prior_theta(
        n = np,
        population_prior = population_prior,
        model_family = "intensity",
        M_max = M_max
      )

      # Importance correction from population proposal to
      # the tempered power prior.
      if (prior_power < 1.0) {

        lp_base <- vapply(
          seq_len(np),
          function(i) {
            log_prior_theta_u_intensity(
              theta_u[i, , drop = FALSE],
              priors,
              M_max
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
              model_family = "intensity",
              M_max = M_max
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

    th0 <- init_theta_particles_int(
      np,
      priors,
      M_max
    )

    theta_u <- pack_theta_int(
      delta = th0$delta,
      sigma_meas = th0$sigma_meas,
      sigma_lat = th0$sigma_lat,
      mu_y0 = th0$mu_y0,
      sigma_y0 = th0$sigma_y0,
      M_max = M_max
    )
  }

  # ------------------------------------------------------------
  # Initialise latent state
  # ------------------------------------------------------------

  th <- unpack_theta_int(theta_u, M_max)

  x <- rnorm(
    np,
    mean = th$mu_y0,
    sd = th$sigma_y0
  )

  list(
    x = x,
    w = w_init,
    theta_u = as_matrix(theta_u),
    time_last = -1L,
    np = np,
    ess_frac = ess_frac,
    M_max = M_max,
    priors = priors,

    prior_mode = prior_mode,
    prior_combination = prior_combination,
    prior_power = prior_power,
    population_prior = population_prior,
    pmmh_prior_mode = pmmh_prior_mode,
    prior_id = population_prior$prior_id %||% "base",

    train_cache = data.frame(
      TimeRel = integer(0),
      yc = integer(0)
    ),

    pmmh_nx = as.integer(pmmh_nx),
    pmmh_window_n = as.integer(pmmh_window_n),
    proposal_sd = proposal_sd,
    accept = integer(np),
    move_attempt = integer(np)
  )
}

log_effective_prior_intensity <- function(
    theta_u,
    state
) {

  # ------------------------------------------------------------
  # Original patient-specific/base prior
  # ------------------------------------------------------------

  lp_base <- log_prior_theta_u_intensity(
    theta_u,
    state$priors,
    state$M_max
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
      model_family = "intensity",
      M_max = state$M_max
    )

    alpha <- state$prior_power

    # ----------------------------------------------------------
    # Robust arithmetic mixture
    #
    # p(theta) =
    #   (1-alpha) p_base(theta) +
    #   alpha     p_population(theta)
    # ----------------------------------------------------------

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

    # ----------------------------------------------------------
    # Original geometric / power prior
    #
    # p(theta) proportional to
    # p_base(theta)^(1-alpha) *
    # p_population(theta)^alpha
    # ----------------------------------------------------------

    return(
      (1 - alpha) * lp_base +
        alpha * lp_pop
    )
  }

  lp_base
}

log_prior_theta_u_intensity <- function(
    theta_u,
    priors,
    M_max
) {
  theta_u <- as.matrix(theta_u)

  log_sigma_meas <-
    theta_u[, "log_sigma_meas"]

  log_sigma_lat <-
    theta_u[, "log_sigma_lat"]

  mu_y0 <-
    theta_u[, "mu_y0"]

  log_sigma_y0 <-
    theta_u[, "log_sigma_y0"]

  sigma_meas <- exp(
    log_sigma_meas
  )

  sigma_lat <- exp(
    log_sigma_lat
  )

  sigma_y0 <- exp(
    log_sigma_y0
  )

  # sigma_meas / M ~ LogNormal(...)
  lp_sigma_meas <- stats::dlnorm(
    sigma_meas / M_max,
    meanlog = priors$prior_sigma_meas[1],
    sdlog = priors$prior_sigma_meas[2],
    log = TRUE
  ) -
    log(M_max) +
    log_sigma_meas

  # sigma_lat / M ~ LogNormal(...)
  lp_sigma_lat <- stats::dlnorm(
    sigma_lat / M_max,
    meanlog = priors$prior_sigma_lat[1],
    sdlog = priors$prior_sigma_lat[2],
    log = TRUE
  ) -
    log(M_max) +
    log_sigma_lat

  # mu_y0 ~ Normal(M * mu0_mean, M * mu0_sd)
  lp_mu_y0 <- stats::dnorm(
    mu_y0,
    mean = M_max * priors$mu0_mean,
    sd = M_max * priors$mu0_sd,
    log = TRUE
  )

  # sigma_y0 ~ HalfNormal(0, M * sigma0_sd)
  lp_sigma_y0 <- stats::dnorm(
    sigma_y0,
    mean = 0,
    sd = M_max * priors$sigma0_sd,
    log = TRUE
  ) +
    log(2) +
    log_sigma_y0

  # delta ~ Dirichlet(prior_delta)
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

  alpha <- priors$prior_delta

  lp_delta <- apply(
    delta,
    1,
    function(d) {
      d <- pmax(
        d,
        1e-12
      )

      lp_dirichlet <-
        lgamma(sum(alpha)) -
        sum(lgamma(alpha)) +
        sum(
          (alpha - 1) *
            log(d)
        )

      # Jacobian for the transformed simplex coordinates.
      lp_jacobian <- sum(
        log(d)
      )

      lp_dirichlet +
        lp_jacobian
    }
  )

  as.numeric(
    lp_sigma_meas +
      lp_sigma_lat +
      lp_mu_y0 +
      lp_sigma_y0 +
      lp_delta
  )
}

propose_theta_u_intensity <- function(
    theta_u_curr,
    M_max,
    proposal_sd
) {
  th <- as.matrix(
    theta_u_curr
  )

  th[1, "log_sigma_meas"] <-
    th[1, "log_sigma_meas"] +
    rnorm(
      1,
      0,
      proposal_sd[1]
    )

  th[1, "log_sigma_lat"] <-
    th[1, "log_sigma_lat"] +
    rnorm(
      1,
      0,
      proposal_sd[2]
    )

  th[1, "mu_y0"] <-
    th[1, "mu_y0"] +
    rnorm(
      1,
      0,
      proposal_sd[3]
    )

  th[1, "log_sigma_y0"] <-
    th[1, "log_sigma_y0"] +
    rnorm(
      1,
      0,
      proposal_sd[4]
    )

  delta_cols <- paste0(
    "log_delta_",
    seq_len(M_max - 1L)
  )

  # Work on centred log-simplex coordinates.
  th[1, delta_cols] <-
    th[1, delta_cols] -
    mean(
      th[1, delta_cols]
    )

  th[1, delta_cols] <-
    th[1, delta_cols] +
    rnorm(
      M_max - 1L,
      0,
      proposal_sd[5]
    )

  th[1, delta_cols] <-
    th[1, delta_cols] -
    mean(
      th[1, delta_cols]
    )

  th
}

pf_intensity_aux_init <- function(Nx, theta_j, M_max) {
  th <- unpack_theta_int(theta_j, M_max)

  list(
    x = rnorm(Nx, mean = th$mu_y0[1], sd = th$sigma_y0[1]),
    w = rep(1 / Nx, Nx),
    time_last = -1L,
    Nx = as.integer(Nx),
    M_max = as.integer(M_max)
  )
}

pf_intensity_aux_propagate <- function(
    st,
    theta_j,
    dt
) {
  dt <- as.integer(dt)

  if (dt <= 0L) {
    return(st)
  }

  th <- unpack_theta_int(
    theta_j,
    st$M_max
  )

  st$x <- st$x +
    sqrt(dt) *
    as.numeric(th$sigma_lat[1]) *
    rnorm(st$Nx)

  st
}

pf_intensity_aux_update_one <- function(st, theta_j, yc_t, t_now) {
  t_now <- as.integer(t_now)
  t_prev <- as.integer(st$time_last)

  dt <- if (t_prev < 0L) 0L else max(0L, t_now - t_prev)
  st <- pf_intensity_aux_propagate(st, theta_j, dt)

  th <- unpack_theta_int(theta_j, st$M_max)

  sigma_rep <- rep(th$sigma_meas[1], st$Nx)
  delta_rep <- th$delta[rep(1L, st$Nx), , drop = FALSE]

  pmf <- lik_intensity_ordlogit(
    y_cat = yc_t,
    x = st$x,
    sigma_meas = sigma_rep,
    delta = delta_rep,
    M_max = st$M_max
  )
  pmf <- pmax(pmf, 1e-12)

  loglik_inc <- log_sum_exp(log(pmax(st$w, 1e-12)) + log(pmf))
  st$w <- log_normalize(log(pmax(st$w, 1e-12)) + log(pmf))

  ess <- 1 / sum(st$w^2)
  if (ess < 0.5 * st$Nx) {
    idx <- systematic_resample_idx(st$w)
    st$x <- st$x[idx]
    st$w <- rep(1 / st$Nx, st$Nx)
  }

  st$time_last <- t_now

  list(st = st, loglik_inc = loglik_inc)
}

pf_intensity_aux_update_many <- function(st, theta_j, yc_vec, time_vec) {
  stopifnot(length(yc_vec) == length(time_vec))
  if (!length(yc_vec)) return(list(st = st, loglik_sum = 0.0))

  s <- 0.0
  for (i in seq_along(yc_vec)) {
    rr <- pf_intensity_aux_update_one(st, theta_j, yc_t = yc_vec[i], t_now = time_vec[i])
    st <- rr$st
    s <- s + rr$loglik_inc
  }

  list(st = st, loglik_sum = s)
}

pf_intensity_resample_move <- function(state) {
  w <- state$w
  np <- length(w)
  ess <- 1 / sum(w^2)

  # --- diagnostics: record every ESS check, not just the ones that fire ---
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

    # PMMH rejuvenation instead of Liu–West
    if (!is.null(state$train_cache) && nrow(state$train_cache) > 0) {
      move_idx <- sample.int(np, size = min(50L, np), replace = FALSE)
      state$n_rejuvenated <- length(move_idx)               
      state$n_moves_cum   <- (state$n_moves_cum %||% 0L) + length(move_idx)
      for (j in move_idx) {
        state <- pmmh_move_one_particle_intensity(state, j)
      }
    }
  }

  state
}

get_pmmh_window_intensity <- function(df, window_n = Inf) {
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

rebuild_pf_and_loglik_intensity <- function(state, theta_u_one) {
  theta_u_one <- as_matrix(theta_u_one)

  st_aux <- pf_intensity_aux_init(
    Nx = state$pmmh_nx,
    theta_j = theta_u_one,
    M_max = state$M_max
  )

  df <- get_pmmh_window_intensity(
    state$train_cache,
    window_n = state$pmmh_window_n %||% Inf
  )

  if (!is.null(df) && nrow(df) > 0) {
    rr <- pf_intensity_aux_update_many(
      st = st_aux,
      theta_j = theta_u_one,
      yc_vec = df$yc,
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
    x_last = st_aux$x[idx],
    aux_state = st_aux
  )
}

pmmh_move_one_particle_intensity <- function(state, j) {
  state$move_attempt[j] <- state$move_attempt[j] + 1L

  theta_curr <- state$theta_u[j, , drop = FALSE]

  curr <- rebuild_pf_and_loglik_intensity(state, theta_curr)
  ll_curr <- curr$loglik_hat
  lp_curr <- log_effective_prior_intensity(theta_curr, state)

  theta_prop <- propose_theta_u_intensity(
    theta_u_curr = theta_curr,
    M_max = state$M_max,
    proposal_sd = state$proposal_sd
  )

  prop <- rebuild_pf_and_loglik_intensity(state, theta_prop)
  ll_prop <- prop$loglik_hat
  lp_prop <- log_effective_prior_intensity(theta_prop, state)

  log_alpha <- (ll_prop + lp_prop) - (ll_curr + lp_curr)

  if (is.finite(log_alpha) && log(runif(1)) < log_alpha) {
    state$theta_u[j, ] <- theta_prop
    state$x[j] <- prop$x_last
    state$accept[j] <- state$accept[j] + 1L
  } else {
    # keep current theta, but refresh x to be coherent with cached data
    state$x[j] <- curr$x_last
  }

  state
}

pf_intensity_propagate <- function(
    state,
    dt
) {
  dt <- as.integer(dt)

  if (dt <= 0L) {
    return(state)
  }

  np <- length(state$x)

  th <- unpack_theta_int(
    state$theta_u,
    state$M_max
  )

  state$x <- state$x +
    sqrt(dt) *
    th$sigma_lat *
    rnorm(np)

  state
}

pf_intensity_update_one <- function(state, yc_t, dt) {
  # 1) propagate latent by dt (if gaps)
  state <- pf_intensity_propagate(state, dt)

  # 2) weight update using likelihood at time t
  th <- unpack_theta_int(state$theta_u, state$M_max)
  pmf <- lik_intensity_ordlogit(yc_t, state$x, th$sigma_meas, th$delta, state$M_max)

  logw <- log(pmax(state$w, 1e-12)) + log(pmax(pmf, 1e-12))
  state$w <- log_normalize(logw)

  # 3) resample + rejuvenate if ESS low
  state <- pf_intensity_resample_move(state)

  # Return the updated state
  state
}

pf_intensity_update_many <- function(
    state,
    yc_vec,
    time_vec
) {
  stopifnot(
    length(yc_vec) == length(time_vec)
  )

  if (length(yc_vec) == 0L) {
    return(state)
  }

  time_vec <- as.integer(time_vec)
  yc_vec <- as.integer(yc_vec)

  t_prev <- state$time_last

  if (
    is.null(t_prev) ||
    is.na(t_prev)
  ) {
    t_prev <- -1L
  }

  t_prev <- as.integer(t_prev)

  for (i in seq_along(yc_vec)) {
    t_now <- as.integer(time_vec[i])
    yc_now <- as.integer(yc_vec[i])

    # Add only the observation currently being assimilated.
    #
    # This ensures that a PMMH move triggered during this update
    # can use observations up to and including t_now, but cannot
    # use future observations from the same input batch.
    current_df <- data.frame(
      TimeRel = t_now,
      yc = yc_now
    )

    if (
      is.null(state$train_cache) ||
      nrow(state$train_cache) == 0L
    ) {
      state$train_cache <- current_df
    } else {
      state$train_cache <-
        dplyr::bind_rows(
          state$train_cache,
          current_df
        ) %>%
        dplyr::distinct(
          .data$TimeRel,
          .keep_all = TRUE
        ) %>%
        dplyr::arrange(
          .data$TimeRel
        )
    }

    # No propagation before the first observed measurement.
    dt <- if (t_prev < 0L) {
      0L
    } else {
      max(
        0L,
        t_now - t_prev
      )
    }

    state <- pf_intensity_update_one(
      state = state,
      yc_t = yc_now,
      dt = dt
    )

    state$time_last <- t_now
    t_prev <- t_now
  }

  state
}