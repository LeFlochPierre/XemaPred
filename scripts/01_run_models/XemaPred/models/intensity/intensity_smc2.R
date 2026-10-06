# intensity_smc2_cache_functions.R
suppressPackageStartupMessages({ library(stats) })

# ============================================================
# INTENSITY SMC² (PF inside PF) — OrderedRW (EczemaPred-style)
#
# Outer particles (j = 1..J): global parameters
#   theta_u[j, ] encodes:
#     delta (cutpoint increments), sigma_meas, sigma_lat, mu_y0, sigma_y0
#
# Inner PF (per outer particle j, per patient pid):
#   latent x(t) particles (length Nx), weights, time_last
#
# IMPORTANT SHAPE RULE:
#   lik_intensity_ordlogit() expects:
#     x length = Nx
#     sigma_meas length = Nx
#     delta = Nx x (M_max-1)
#   Therefore we must broadcast outer delta/sigma_meas to Nx inside update/forecast.
# ============================================================


# -----------------------------
# Inner PF: init
# -----------------------------
pf_intensity_inner_init <- function(Nx, theta_j, priors, M_max,
                                    ess_x = 0.5, a_x = 0.98) {
  th <- unpack_theta_int(theta_j, M_max)

  # x0 ~ N(mu_y0, sigma_y0^2)
  x <- rnorm(Nx, mean = as.numeric(th$mu_y0[1]), sd = as.numeric(th$sigma_y0[1]))

  list(
    Nx = as.integer(Nx),
    ess_x = as.numeric(ess_x),
    a_x = as.numeric(a_x),
    M_max = as.integer(M_max),

    x = x,
    w = rep(1 / Nx, Nx),
    time_last = -1L
  )
}

log_prior_theta_u_intensity <- function(theta_u, priors, M_max) {
  theta_u <- as.matrix(theta_u)

  log_sigma_meas <- theta_u[, "log_sigma_meas"]
  log_sigma_lat  <- theta_u[, "log_sigma_lat"]
  mu_y0          <- theta_u[, "mu_y0"]
  log_sigma_y0   <- theta_u[, "log_sigma_y0"]

  sigma_meas <- exp(log_sigma_meas)
  sigma_lat  <- exp(log_sigma_lat)
  sigma_y0   <- exp(log_sigma_y0)

  # sigma_meas / M ~ LogNormal(prior_sigma_meas[1], prior_sigma_meas[2])
  lp_sigma_meas <- stats::dlnorm(
    sigma_meas / M_max,
    meanlog = priors$prior_sigma_meas[1],
    sdlog = priors$prior_sigma_meas[2],
    log = TRUE
  ) - log(M_max) + log_sigma_meas

  # sigma_lat / M ~ LogNormal(prior_sigma_lat[1], prior_sigma_lat[2])
  lp_sigma_lat <- stats::dlnorm(
    sigma_lat / M_max,
    meanlog = priors$prior_sigma_lat[1],
    sdlog = priors$prior_sigma_lat[2],
    log = TRUE
  ) - log(M_max) + log_sigma_lat

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
  ) + log(2) + log_sigma_y0

  # delta ~ Dirichlet(prior_delta)
  delta_cols <- paste0("log_delta_", seq_len(M_max - 1L))
  log_delta <- theta_u[, delta_cols, drop = FALSE]

  log_delta <- log_delta - rowMeans(log_delta)
  log_delta <- log_delta - apply(log_delta, 1, max)

  delta_raw <- exp(log_delta)
  delta <- delta_raw / rowSums(delta_raw)

  alpha <- priors$prior_delta

  lp_delta <- apply(delta, 1, function(d) {
    d <- pmax(d, 1e-12)

    lp_dirichlet <- lgamma(sum(alpha)) -
      sum(lgamma(alpha)) +
      sum((alpha - 1) * log(d))

    # Jacobian term for simplex/log-delta representation.
    lp_jacobian <- sum(log(d))

    lp_dirichlet + lp_jacobian
  })

  as.numeric(
    lp_sigma_meas +
      lp_sigma_lat +
      lp_mu_y0 +
      lp_sigma_y0 +
      lp_delta
  )
}

propose_theta_u_intensity <- function(theta_u_curr, M_max, proposal_sd) {
  th <- as.matrix(theta_u_curr)

  th[1, "log_sigma_meas"] <- th[1, "log_sigma_meas"] + rnorm(1, 0, proposal_sd[1])
  th[1, "log_sigma_lat"]  <- th[1, "log_sigma_lat"]  + rnorm(1, 0, proposal_sd[2])
  th[1, "mu_y0"]          <- th[1, "mu_y0"]          + rnorm(1, 0, proposal_sd[3])
  th[1, "log_sigma_y0"]   <- th[1, "log_sigma_y0"]   + rnorm(1, 0, proposal_sd[4])

  delta_cols <- paste0("log_delta_", seq_len(M_max - 1L))

  # Work on centred log-delta coordinates.
  # This removes the non-identifiable common shift.
  th[1, delta_cols] <- th[1, delta_cols] - mean(th[1, delta_cols])
  th[1, delta_cols] <- th[1, delta_cols] + rnorm(M_max - 1L, 0, proposal_sd[5])
  th[1, delta_cols] <- th[1, delta_cols] - mean(th[1, delta_cols])

  th
}



pf_intensity_inner_resample <- function(st) {
  w <- st$w
  Nx <- st$Nx
  ess <- 1 / sum(w^2)

  if (ess < st$ess_x * Nx) {
    idx <- systematic_resample_idx(w)
    st$x <- st$x[idx]
    st$w <- rep(1 / Nx, Nx)
  }
  st
}

pf_intensity_inner_propagate <- function(st, dt, sigma_lat) {
  dt <- as.integer(dt)
  if (dt <= 0L) return(st)

  # RW aggregation across dt steps (sqrt(dt) version)
  st$x <- st$x + sqrt(dt) * sigma_lat * rnorm(st$Nx)
  st
}

# -----------------------------
# Inner PF: update one obs
# -----------------------------
pf_intensity_inner_update_one <- function(st, theta_j, yc_t, t_now) {
  t_now <- as.integer(t_now)
  t_prev <- st$time_last
  if (is.null(t_prev) || is.na(t_prev)) t_prev <- -1L
  t_prev <- as.integer(t_prev)

  th <- unpack_theta_int(theta_j, st$M_max)

  dt <- if (t_prev < 0L) 0L else max(0L, t_now - t_prev)

  st <- pf_intensity_inner_propagate(
    st,
    dt = dt,
    sigma_lat = as.numeric(th$sigma_lat[1])
  )

  # ---- broadcast outer params to inner particle dimension (Nx) ----
  sigma_meas_j <- as.numeric(th$sigma_meas[1])
  sigma_rep <- rep(sigma_meas_j, st$Nx)

  delta_1row <- th$delta[1, , drop = FALSE]                  # 1 x (M_max-1)
  delta_rep  <- delta_1row[rep(1L, st$Nx), , drop = FALSE]   # Nx x (M_max-1)

  pmf <- lik_intensity_ordlogit(
    y_cat = as.integer(yc_t),
    x = st$x,                     # length Nx
    sigma_meas = sigma_rep,       # length Nx
    delta = delta_rep,            # Nx x (M_max-1)
    M_max = st$M_max
  )
  pmf <- pmax(pmf, 1e-12)

  loglik_inc <- log_sum_exp(log(pmax(st$w, 1e-12)) + log(pmf))
  st$w <- log_normalize(log(pmax(st$w, 1e-12)) + log(pmf))

  st <- pf_intensity_inner_resample(st)
  st$time_last <- t_now

  list(st = st, loglik_inc = loglik_inc)
}


pf_intensity_inner_update_many <- function(st, theta_j, yc_vec, t_vec) {
  stopifnot(length(yc_vec) == length(t_vec))
  if (!length(yc_vec)) return(list(st = st, loglik_sum = 0.0))

  s <- 0.0
  for (i in seq_along(yc_vec)) {
    rr <- pf_intensity_inner_update_one(st, theta_j, yc_t = yc_vec[i], t_now = t_vec[i])
    st <- rr$st
    s <- s + rr$loglik_inc
  }
  list(st = st, loglik_sum = s)
}

# -----------------------------
# Inner PF: forecast
# -----------------------------
pf_intensity_inner_forecast <- function(st, theta_j, time_test, np_sim = 2000L) {
  time_test <- as.integer(time_test)
  H <- length(time_test)

  th <- unpack_theta_int(theta_j, st$M_max)

  w_use <- normalize_weights(st$w)
  idx <- sample.int(st$Nx, size = np_sim, replace = TRUE, prob = w_use)
  x <- st$x[idx]

  t_prev <- st$time_last
  if (is.null(t_prev) || is.na(t_prev)) t_prev <- -1L
  t_prev <- as.integer(t_prev)

  x_fc <- matrix(NA_real_, nrow = np_sim, ncol = H)
  y_fc <- matrix(NA_integer_, nrow = np_sim, ncol = H)

  # broadcast delta, sigma_meas once for all forecast sims
  sigma_meas_j <- as.numeric(th$sigma_meas[1])
  sigma_rep <- rep(sigma_meas_j, np_sim)

  delta_1row <- th$delta[1, , drop = FALSE]
  delta_rep  <- delta_1row[rep(1L, np_sim), , drop = FALSE]

  K <- st$M_max + 1L

  for (h in seq_len(H)) {
    t_now <- as.integer(time_test[h])
    dt <- if (t_prev < 0L) 0L else max(0L, t_now - t_prev)

    if (dt > 0L) {
      x <- x + sqrt(dt) * as.numeric(th$sigma_lat[1]) * rnorm(np_sim)
    }

    x_fc[, h] <- x

    # sample categorical y using pmf across categories (cheap for M_max small)
    pmf_mat <- matrix(0, nrow = np_sim, ncol = K)
    for (yc in seq_len(K)) {
      pmf_mat[, yc] <- lik_intensity_ordlogit(
        y_cat = yc,
        x = x,
        sigma_meas = sigma_rep,
        delta = delta_rep,
        M_max = st$M_max
      )
    }
    pmf_mat <- pmax(pmf_mat, 1e-12)
    pmf_mat <- pmf_mat / rowSums(pmf_mat)

    u <- runif(np_sim)
    cdf <- t(apply(pmf_mat, 1, cumsum))
    y_fc[, h] <- vapply(seq_len(np_sim), function(i) which(u[i] <= cdf[i, ])[1], integer(1))

    t_prev <- t_now
  }

  list(x_forecast = x_fc, y_forecast = y_fc)
}

get_pmmh_window_intensity <- function(dfp, window_n = Inf) {
  if (is.null(dfp) || nrow(dfp) == 0) return(dfp)

  dfp <- dfp[order(dfp$TimeRel), , drop = FALSE]

  if (is.finite(window_n)) {
    n <- as.integer(window_n)
    if (nrow(dfp) > n) {
      dfp <- dfp[(nrow(dfp) - n + 1L):nrow(dfp), , drop = FALSE]
    }
  }

  dfp
}

rebuild_inner_and_loglik_intensity <- function(st, theta_u_one, window_n = Inf) {
  theta_u_one <- as.matrix(theta_u_one)

  inner_new <- setNames(vector("list", length(st$patient_ids)), st$patient_ids)
  ll_sum <- 0.0

  for (pid in st$patient_ids) {
    inn <- pf_intensity_inner_init(
      Nx = st$Nx,
      theta_j = theta_u_one,
      priors = st$priors,
      M_max = st$M_max,
      ess_x = st$ess_x,
      a_x = st$a_x
    )

    dfp <- get_pmmh_window_intensity(
      st$train_cache[[pid]],
      window_n = window_n
    )

    if (!is.null(dfp) && nrow(dfp) > 0) {
      rr <- pf_intensity_inner_update_many(
        st = inn,
        theta_j = theta_u_one,
        yc_vec = dfp$yc,
        t_vec  = dfp$TimeRel
      )
      inn <- rr$st
      ll_sum <- ll_sum + rr$loglik_sum
    }

    inner_new[[pid]] <- inn
  }

  list(inner = inner_new, loglik_hat = ll_sum)
}

pmmh_move_one_outer_intensity <- function(st, j, proposal_sd = NULL) {
  window_n <- st$pmmh_window_n %||% Inf

  st$move_attempt[j] <- st$move_attempt[j] + 1L
  if (is.null(proposal_sd)) proposal_sd <- st$proposal_sd

  theta_curr <- st$theta_u[j, , drop = FALSE]

  curr <- rebuild_inner_and_loglik_intensity(
    st,
    theta_u_one = theta_curr,
    window_n = window_n
  )

  ll_curr <- curr$loglik_hat
  lp_curr <- log_prior_theta_u_intensity(theta_curr, st$priors, st$M_max)

  theta_prop <- propose_theta_u_intensity(theta_curr, st$M_max, proposal_sd)

  prop <- rebuild_inner_and_loglik_intensity(
    st,
    theta_u_one = theta_prop,
    window_n = window_n
  )

  ll_prop <- prop$loglik_hat
  lp_prop <- log_prior_theta_u_intensity(theta_prop, st$priors, st$M_max)

  log_alpha <- (ll_prop + lp_prop) - (ll_curr + lp_curr)

  if (is.finite(log_alpha) && log(runif(1)) < log_alpha) {
    st$theta_u[j, ] <- theta_prop
    st$inner[[j]] <- prop$inner
    st$loglik_hat[j] <- ll_prop
    st$accept[j] <- st$accept[j] + 1L
  } else {
    st$inner[[j]] <- curr$inner
    st$loglik_hat[j] <- ll_curr
  }

  st
}

smc2_intensity_pmmh_move <- function(st, n_steps = 1L, n_workers = 1L, seed = NULL) {
  n_steps <- as.integer(n_steps)
  if (n_steps <= 0L) return(st)

  J <- st$J
  n_workers <- as.integer(n_workers %||% 1L)

  for (m in seq_len(n_steps)) {

    if (n_workers > 1L && foreach::getDoParWorkers() > 1L) {
      if (!is.null(seed)) doRNG::registerDoRNG(seed + m)

      res <- foreach::foreach(
        j = seq_len(J),
        .packages = c("stats"),
        .export = c(
          "pmmh_move_one_outer_intensity",
          "propose_theta_u_intensity",
          "rebuild_inner_and_loglik_intensity",
          "get_pmmh_window_intensity",
          "log_prior_theta_u_intensity",
          "pf_intensity_inner_init",
          "pf_intensity_inner_update_many",
          "pf_intensity_inner_update_one",
          "pf_intensity_inner_propagate",
          "pf_intensity_inner_resample",
          "unpack_theta_int",
          "lik_intensity_ordlogit",
          "inv_logit",
          "as_matrix",
          "make_ct_mat",
          "clip",
          "log_sum_exp",
          "log_normalize",
          "normalize_weights",
          "systematic_resample_idx"
        )
      ) %dorng% {
        st_j <- st
        st_j <- pmmh_move_one_outer_intensity(st_j, j, proposal_sd = st$proposal_sd)
        list(
          j = j,
          theta = st_j$theta_u[j, ],
          inner = st_j$inner[[j]],
          ll = st_j$loglik_hat[j],
          acc = st_j$accept[j],
          att = st_j$move_attempt[j]
        )
      }

      for (rr in res) {
        st$theta_u[rr$j, ] <- rr$theta
        st$inner[[rr$j]] <- rr$inner
        st$loglik_hat[rr$j] <- rr$ll
        st$accept[rr$j] <- rr$acc
        st$move_attempt[rr$j] <- rr$att
      }

    } else {
      for (j in seq_len(J)) {
        st <- pmmh_move_one_outer_intensity(st, j, proposal_sd = st$proposal_sd)
      }
    }
  }

  st
}


# ============================================================
# Outer SMC² state
# ============================================================

smc2_intensity_init <- function(J, Nx, priors, M_max, patient_ids,
                                ess_theta = 0.5, a_theta = 0.98,
                                ess_x = 0.5, a_x = 0.98) {
  patient_ids <- as.character(patient_ids)

  th0 <- init_theta_particles_int(J, priors, M_max)
  theta_u <- pack_theta_int(th0$delta, th0$sigma_meas, th0$sigma_lat, th0$mu_y0, th0$sigma_y0, M_max)
  theta_u <- as.matrix(theta_u)

  cat("theta_u colnames:\n")
  print(colnames(theta_u))

  W <- rep(1 / J, J)

  inner <- vector("list", J)
  for (j in seq_len(J)) {
    inner_j <- setNames(vector("list", length(patient_ids)), patient_ids)
    theta_j <- theta_u[j, , drop = FALSE]
    for (pid in patient_ids) {
      inner_j[[pid]] <- pf_intensity_inner_init(
        Nx = Nx, theta_j = theta_j, priors = priors, M_max = M_max,
        ess_x = ess_x, a_x = a_x
      )
    }
    inner[[j]] <- inner_j
  }

  list(
    J = as.integer(J),
    Nx = as.integer(Nx),
    M_max = as.integer(M_max),

    ess_theta = as.numeric(ess_theta),
    a_theta = as.numeric(a_theta),

    ess_x = as.numeric(ess_x),
    a_x = as.numeric(a_x),

    theta_u = theta_u,  # J x P
    W = W,              # length J
    inner = inner,       # list[[j]][[pid]]

    pmmh_window_n = Inf,
    priors = priors,
    patient_ids = patient_ids,
    train_cache = setNames(vector("list", length(patient_ids)), patient_ids), # for storing train obs if needed#

    loglik_hat = rep(0.0, J), # for storing latest loglik estimate per particle
    accept = integer(J),
    move_attempt = integer(J),
    proposal_sd =  c(0.054, 0.054, 0.18, 0.054, 0.054)
  )
}

smc2_intensity_resample_outer_only <- function(st) {
  W <- st$W
  J <- st$J
  ess <- 1 / sum(W^2)

  did_resample <- FALSE
  if (ess < st$ess_theta * J) {
    did_resample <- TRUE
    idx <- systematic_resample_idx(W)

    st$theta_u <- st$theta_u[idx, , drop = FALSE]
    st$inner <- st$inner[idx]
    st$loglik_hat <- st$loglik_hat[idx]

    st$accept <- st$accept[idx]
    st$move_attempt <- st$move_attempt[idx]

    st$W <- rep(1 / J, J)
  }

  attr(st, "did_resample") <- did_resample
  st
}


# obs_df columns: Patient, TimeRel, yc (yc in 1..K)
smc2_intensity_update_many <- function(st, obs_df, n_workers = 1L, seed = NULL) {
  if (is.null(obs_df) || nrow(obs_df) == 0) {
    return(list(st = st, loglik_vec = rep(0.0, st$J)))
  }

  t_update0 <- Sys.time()
  pmmh_time <- 0
  did_rs <- FALSE

  stopifnot(all(c("Patient", "TimeRel", "yc") %in% colnames(obs_df)))

  obs_df <- obs_df %>%
    dplyr::mutate(
      Patient = as.character(.data$Patient),
      TimeRel = as.integer(.data$TimeRel),
      yc      = as.integer(.data$yc)
    ) %>%
    dplyr::arrange(.data$Patient, .data$TimeRel)

  obs_split <- split(obs_df, obs_df$Patient)

  # cache the obs we assimilate (needed for rebuild after LW move)
  if (!is.null(st$train_cache)) {
    for (pid in names(obs_split)) {
      newp <- obs_split[[pid]] %>% dplyr::arrange(.data$TimeRel)
      oldp <- st$train_cache[[pid]]
      if (is.null(oldp) || nrow(oldp) == 0) {
        st$train_cache[[pid]] <- newp
      } else {
        st$train_cache[[pid]] <- dplyr::bind_rows(oldp, newp) %>%
          dplyr::distinct(.data$TimeRel, .keep_all = TRUE) %>%
          dplyr::arrange(.data$TimeRel)
      }
    }
  }

  # -------------------------
  # Outer update (parallel)
  # -------------------------
  if (!is.null(seed)) {
    RNGkind("L'Ecuyer-CMRG")
    set.seed(as.integer(seed))
  }

  if (n_workers <= 1L) {
    foreach::registerDoSEQ()
  }

  if (!is.null(seed)) doRNG::registerDoRNG(as.integer(seed))

  res <- foreach::foreach(
    j = seq_len(st$J),
    .inorder = TRUE,
    .packages = c("stats"),
    .export = c(
      "as_matrix",
      "pf_intensity_inner_update_many",
      "pf_intensity_inner_update_one",
      "pf_intensity_inner_propagate",
      "pf_intensity_inner_resample",
      "unpack_theta_int",
      "lik_intensity_ordlogit",
      "make_ct_mat",
      "inv_logit",
      "log_sum_exp",
      "log_normalize",
      "systematic_resample_idx",
      "normalize_weights"
    )
  ) %dorng% {
    theta_j <- st$theta_u[j, , drop = FALSE]
    inner_j <- st$inner[[j]]

    s <- 0.0
    for (pid in names(obs_split)) {
      dfp <- obs_split[[pid]]
      if (nrow(dfp) == 0) next

      rr <- pf_intensity_inner_update_many(
        st = inner_j[[pid]],
        theta_j = theta_j,
        yc_vec = dfp$yc,
        t_vec  = dfp$TimeRel
      )
      inner_j[[pid]] <- rr$st
      s <- s + rr$loglik_sum
    }

    list(inner_j = inner_j, loglik = s)
  }

  # stitch back
  loglik_j <- numeric(st$J)
  for (j in seq_len(st$J)) {
    st$inner[[j]] <- res[[j]]$inner_j
    loglik_j[j] <- res[[j]]$loglik
  }

  # 1) accumulate running loglik estimate (p_hat)
  st$loglik_hat <- st$loglik_hat + loglik_j

  # 2) standard SMC² outer weight update (NO tempering)
  logW_raw <- log(pmax(st$W, 1e-12)) + loglik_j
  st$W <- log_normalize(logW_raw)

  ess_before <- 1 / sum(st$W^2)
  maxw_before <- max(st$W)
  cat(sprintf("[DBG][outer] ESS(before resample)=%.1f/%d thr=%.1f maxW=%.3f\n",
              ess_before, st$J, st$ess_theta * st$J, maxw_before))

  # 3) resample outer if ESS low
  st <- smc2_intensity_resample_outer_only(st)

  did_rs <- isTRUE(attr(st, "did_resample"))
  ess_after <- 1 / sum(st$W^2)
  maxw_after <- max(st$W)
  cat(sprintf("[DBG][outer] resample=%s | ESS(after)=%.1f/%d maxW=%.3f\n",
              did_rs, ess_after, st$J, maxw_after))


  # 4) PMMH move if resampling happened (same policy as extent)
  did_rs <- isTRUE(attr(st, "did_resample"))

  if (did_rs) {
    t_pmmh0 <- Sys.time()

    st <- smc2_intensity_pmmh_move(
      st,
      n_steps = 1L,
      n_workers = n_workers,
      seed = seed
    )

    t_pmmh1 <- Sys.time()
    pmmh_time <- as.numeric(difftime(t_pmmh1, t_pmmh0, units = "secs"))
  }

  t_update1 <- Sys.time()
  update_time <- as.numeric(difftime(t_update1, t_update0, units = "secs"))

  list(
    st = st,
    loglik_vec = loglik_j,
    update_time = update_time,
    pmmh_time = pmmh_time,
    did_resample = did_rs,
    ess_before = ess_before
  )

}


# Mixture forecast for one patient:
# sample outer index j ~ W, then forecast from that inner PF state
# Returns: x_forecast, y_forecast, theta_forecast_u (so your scoring code still works)
smc2_intensity_forecast_pid <- function(st, pid, time_test, np_sim = 2000L) {
  pid <- as.character(pid)
  time_test <- as.integer(time_test)
  H <- length(time_test)

  J <- length(st$W)
  idx_theta <- sample.int(J, size = np_sim, replace = TRUE, prob = normalize_weights(st$W))

  x_out <- matrix(NA_real_, nrow = np_sim, ncol = H)
  y_out <- matrix(NA_integer_, nrow = np_sim, ncol = H)

  theta_out <- matrix(NA_real_, nrow = np_sim, ncol = ncol(st$theta_u))
  colnames(theta_out) <- colnames(st$theta_u)

  ujs <- sort(unique(idx_theta))
  for (j in ujs) {
    rows <- which(idx_theta == j)
    m <- length(rows)
    if (m == 0L) next

    theta_j <- st$theta_u[j, , drop = FALSE]
    inner_st <- st$inner[[j]][[pid]]

    fc <- pf_intensity_inner_forecast(inner_st, theta_j, time_test = time_test, np_sim = m)

    x_out[rows, ] <- fc$x_forecast
    y_out[rows, ] <- fc$y_forecast
    theta_out[rows, ] <- st$theta_u[rep(j, m), , drop = FALSE]
  }

  list(x_forecast = x_out, y_forecast = y_out, theta_forecast_u = theta_out)
}
