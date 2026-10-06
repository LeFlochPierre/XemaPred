# extent_smc2_cache_functions.R
suppressPackageStartupMessages({ library(stats) })

# ============================================================
# EXTENT SMC² (PF inside PF)
#
# Outer particles: hyperparameters
#   Theta = (sigma, mu10, sigma10)
#   represented as theta_u = [log_sigma, mu10, log_sigma10]
#
# Inner PF (per outer particle j, per patient pid):
#   - logit_tss1 (RW with sigma from outer)
#   - logit_p10  (static, drawn from N(mu10, sigma10^2); rejuvenated by PMMH)
#   - y_lat      (derived via extent_derived() and Markov recursion)
#   - w          (inner weights)
#   - time_last  (patient last assimilated time)
#
# This matches the EczemaPred hierarchical extent model,
# but uses SMC²-style sequential inference.
# ============================================================

# -----------------------------
# Hyper-priors
# -----------------------------
get_priors_extent_hyper <- function() {
  list(
    # sigma ~ N+(0, (0.25 log(5))^2)
    sigma_loc = 0.0,
    sigma_scale = 0.25 * log(5),

    # mu10 ~ N(0, 1)
    mu10_loc = 0.0,
    mu10_scale = 1.0,

    # sigma10 ~ N+(0, 1.5^2)
    sigma10_loc = 0.0,
    sigma10_scale = 1.5,

    # initial latent logit_tss1_0 ~ N(-1, 1)
    logit_tss1_0_loc = -1.0,
    logit_tss1_0_scale = 1.0
  )
}

log_prior_theta_u_extent <- function(theta_u, priors) {
  theta_u <- as.matrix(theta_u)
  th <- unpack_theta_extent_hyper(theta_u)

  sigma   <- as.numeric(th$sigma)
  mu10    <- as.numeric(th$mu10)
  sigma10 <- as.numeric(th$sigma10)

  # priors are half-normal via abs(N(0, sd)), loc=0 in your writeup
  lp_sigma   <- dnorm(sigma,   mean = priors$sigma_loc,   sd = priors$sigma_scale,   log = TRUE) + log(2)
  lp_mu10    <- dnorm(mu10,    mean = priors$mu10_loc,    sd = priors$mu10_scale,    log = TRUE)
  lp_sigma10 <- dnorm(sigma10, mean = priors$sigma10_loc, sd = priors$sigma10_scale, log = TRUE) + log(2)

  # Jacobian for exp transforms: sigma = exp(log_sigma), sigma10 = exp(log_sigma10)
  lp_jac <- as.numeric(theta_u[, "log_sigma"] + theta_u[, "log_sigma10"])

  as.numeric(lp_sigma + lp_mu10 + lp_sigma10 + lp_jac)
}


pack_theta_extent_hyper <- function(log_sigma, mu10, log_sigma10) {
  out <- cbind(log_sigma = log_sigma, mu10 = mu10, log_sigma10 = log_sigma10)
  out
}

unpack_theta_extent_hyper <- function(theta_u) {
  theta_u <- as.matrix(theta_u)
  list(
    sigma = exp(theta_u[, "log_sigma"]),
    mu10 = theta_u[, "mu10"],
    sigma10 = exp(theta_u[, "log_sigma10"])
  )
}

init_theta_particles_extent_hyper <- function(J, priors) {
  mu10 <- rnorm(J, priors$mu10_loc, priors$mu10_scale)

  sigma <- abs(rnorm(J, priors$sigma_loc, priors$sigma_scale))
  sigma <- pmax(sigma, 1e-12)

  sigma10 <- abs(rnorm(J, priors$sigma10_loc, priors$sigma10_scale))
  sigma10 <- pmax(sigma10, 1e-12)

  pack_theta_extent_hyper(
    log_sigma = log(sigma),
    mu10 = mu10,
    log_sigma10 = log(sigma10)
  )
}

propose_theta_u_extent <- function(theta_u_curr, proposal_sd = c(0.05, 0.10, 0.05)) {
  th <- as.matrix(theta_u_curr)
  th[1, "log_sigma"]   <- th[1, "log_sigma"]   + rnorm(1, 0, proposal_sd[1])
  th[1, "mu10"]        <- th[1, "mu10"]        + rnorm(1, 0, proposal_sd[2])
  th[1, "log_sigma10"] <- th[1, "log_sigma10"] + rnorm(1, 0, proposal_sd[3])
  th
}

pmmh_move_one_outer_extent <- function(st, j, proposal_sd = c(0.05, 0.10, 0.05)) {
  window_n <- st$pmmh_window_n %||% Inf
  st$move_attempt[j] <- st$move_attempt[j] + 1L

  theta_curr <- st$theta_u[j, , drop = FALSE]
  curr <- rebuild_inner_and_loglik_extent(
    st,
    theta_u_one = theta_curr,
    window_n = window_n
  )

  ll_curr <- curr$loglik_hat
  lp_curr <- log_prior_theta_u_extent(theta_curr, st$priors_h)

  theta_prop <- propose_theta_u_extent(theta_curr, proposal_sd = proposal_sd)

  prop <- rebuild_inner_and_loglik_extent(
    st,
    theta_u_one = theta_prop,
    window_n = window_n
  )
  ll_prop <- prop$loglik_hat
  lp_prop <- log_prior_theta_u_extent(theta_prop, st$priors_h)

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

# -----------------------------
# Inner PF kernels (conditional on outer theta)
# -----------------------------
pf_extent_inner_init <- function(Nx, theta_j, priors, ess_x = 0.5, a_x = 0.98) {
  th <- unpack_theta_extent_hyper(theta_j)
  mu10 <- as.numeric(th$mu10)
  sigma10 <- as.numeric(th$sigma10)

  logit_p10 <- rnorm(Nx, mean = mu10, sd = sigma10)
  logit_tss1 <- rnorm(Nx, priors$logit_tss1_0_loc, priors$logit_tss1_0_scale)

  der <- extent_derived(logit_tss1, logit_p10)
  y_lat <- clip(der$ss1, 1e-12, 1 - 1e-12)

  list(
    Nx = as.integer(Nx),
    ess_x = as.numeric(ess_x),
    a_x = as.numeric(a_x),

    logit_p10 = logit_p10,
    logit_tss1 = logit_tss1,
    y_lat = y_lat,
    w = rep(1 / Nx, Nx),

    time_last = -1L
  )
}


pf_extent_inner_resample_move <- function(st) {
  w <- st$w
  Nx <- st$Nx
  ess <- 1 / sum(w^2)

  if (ess < st$ess_x * Nx) {
    idx <- systematic_resample_idx(w)

    st$logit_p10 <- st$logit_p10[idx]
    st$logit_tss1 <- st$logit_tss1[idx]
    st$y_lat <- st$y_lat[idx]
    st$w <- rep(1 / Nx, Nx)
  }

  st
}

pf_extent_inner_step_one_day <- function(st, sigma) {
  st$logit_tss1 <- st$logit_tss1 + sigma * rnorm(st$Nx)

  der <- extent_derived(st$logit_tss1, st$logit_p10)
  st$y_lat <- clip(der$p01 * (1 - st$y_lat) + der$p11 * st$y_lat, 1e-12, 1 - 1e-12)

  st
}

pf_extent_inner_propagate <- function(st, dt, sigma) {
  dt <- as.integer(dt)
  if (dt <= 0L) return(st)
  for (k in seq_len(dt)) st <- pf_extent_inner_step_one_day(st, sigma = sigma)
  st
}

# Update one observation; returns list(st=..., loglik_inc=...)
pf_extent_inner_update_one <- function(st, y_t, t_now, sigma) {
  t_now <- as.integer(t_now)

  t_prev <- st$time_last
  if (is.null(t_prev) || is.na(t_prev)) t_prev <- -1L
  t_prev <- as.integer(t_prev)

  dt <- if (t_prev < 0L) 0L else max(0L, t_now - t_prev)
  st <- pf_extent_inner_propagate(st, dt, sigma = sigma)

  y_t <- as.integer(pmin(100L, pmax(0L, y_t)))
  ll <- dbinom(y_t, size = 100L, prob = clip(st$y_lat, 1e-12, 1 - 1e-12), log = TRUE)

  loglik_inc <- log_sum_exp(log(pmax(st$w, 1e-12)) + ll)
  st$w <- log_normalize(log(pmax(st$w, 1e-12)) + ll)

  st <- pf_extent_inner_resample_move(st)
  st$time_last <- t_now

  list(st = st, loglik_inc = loglik_inc)
}


pf_extent_inner_update_many <- function(st, y_vec, t_vec, sigma) {
  stopifnot(length(y_vec) == length(t_vec))
  if (!length(y_vec)) return(list(st = st, loglik_sum = 0.0))

  s <- 0.0
  for (i in seq_along(y_vec)) {
    rr <- pf_extent_inner_update_one(st, y_t = y_vec[i], t_now = t_vec[i], sigma = sigma)
    st <- rr$st
    s <- s + rr$loglik_inc
  }
  list(st = st, loglik_sum = s)
}

reset_or_rebuild_inner_filters_extent <- function(st, mode = c("rebuild","reset")) {
  mode <- match.arg(mode)
  has_cache <- !is.null(st$train_cache) &&
    any(vapply(st$train_cache, function(x) !is.null(x) && nrow(x) > 0, logical(1)))

  if (mode == "rebuild" && !has_cache) mode <- "reset"

  for (j in seq_len(st$J)) {
    theta_j <- st$theta_u[j, , drop = FALSE]
    thj <- unpack_theta_extent_hyper(theta_j)
    sigma_j <- as.numeric(thj$sigma)

    pids_j <- names(st$inner[[j]])
    if (mode == "rebuild") pids_j <- intersect(pids_j, names(st$train_cache))

    for (pid in pids_j) {
      inner_new <- pf_extent_inner_init(
        Nx = st$Nx, theta_j = theta_j, priors = st$priors_h,
        ess_x = st$ess_x, a_x = st$a_x
      )

      if (mode == "rebuild") {
        dfp <- st$train_cache[[pid]]
        if (!is.null(dfp) && nrow(dfp) > 0) {
          rr <- pf_extent_inner_update_many(
            inner_new,
            y_vec = dfp$y,
            t_vec = dfp$TimeRel,
            sigma = sigma_j
          )
          inner_new <- rr$st
        }
      }

      st$inner[[j]][[pid]] <- inner_new
    }
  }

  st
}

get_pmmh_window_extent <- function(dfp, window_n = Inf) {
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

rebuild_inner_and_loglik_extent <- function(st, theta_u_one, window_n = Inf) {
  theta_u_one <- as.matrix(theta_u_one)
  th <- unpack_theta_extent_hyper(theta_u_one)
  sigma <- as.numeric(th$sigma)

  inner_new <- setNames(vector("list", length(st$patient_ids)), st$patient_ids)
  ll_sum <- 0.0

  for (pid in st$patient_ids) {
    inn <- pf_extent_inner_init(
      Nx = st$Nx, theta_j = theta_u_one, priors = st$priors_h,
      ess_x = st$ess_x, a_x = st$a_x
    )

    dfp <- get_pmmh_window_extent(
      st$train_cache[[pid]],
      window_n = window_n
    )
    if (!is.null(dfp) && nrow(dfp) > 0) {
      rr <- pf_extent_inner_update_many(
        inn,
        y_vec = dfp$y,
        t_vec = dfp$TimeRel,
        sigma = sigma
      )
      inn <- rr$st
      ll_sum <- ll_sum + rr$loglik_sum
    }

    inner_new[[pid]] <- inn
  }

  list(inner = inner_new, loglik_hat = ll_sum)
}


pf_extent_inner_forecast <- function(st, time_test, sigma, np_sim = 2000L) {
  time_test <- as.integer(time_test)
  H <- length(time_test)

  w_use <- normalize_weights(st$w)
  idx <- sample.int(st$Nx, size = np_sim, replace = TRUE, prob = w_use)

  y_lat <- st$y_lat[idx]
  logit_tss1 <- st$logit_tss1[idx]
  logit_p10 <- st$logit_p10[idx]

  x_mat <- matrix(NA_real_, nrow = np_sim, ncol = H)
  y_mat <- matrix(NA_integer_, nrow = np_sim, ncol = H)

  t_prev <- st$time_last
  if (is.null(t_prev) || is.na(t_prev)) t_prev <- -1L
  t_prev <- as.integer(t_prev)

  for (h in seq_len(H)) {
    t_now <- as.integer(time_test[h])
    dt <- if (t_prev < 0L) 0L else max(0L, t_now - t_prev)

    if (dt > 0L) {
      for (k in seq_len(dt)) {
        logit_tss1 <- logit_tss1 + sigma * rnorm(np_sim)
        der <- extent_derived(logit_tss1, logit_p10)
        y_lat <- clip(der$p01 * (1 - y_lat) + der$p11 * y_lat, 1e-12, 1 - 1e-12)
      }
    }

    x_mat[, h] <- y_lat
    y_mat[, h] <- rbinom(np_sim, size = 100L, prob = y_lat)

    t_prev <- t_now
  }

  list(x_forecast = x_mat, y_forecast = y_mat)
}


# -----------------------------
# Outer SMC² state
# -----------------------------
smc2_extent_init <- function(J, Nx, priors_h, patient_ids,
                             ess_theta = 0.5, a_theta = 0.98,
                             ess_x = 0.5, a_x = 0.98) {

  patient_ids <- as.character(patient_ids)

  theta_u <- init_theta_particles_extent_hyper(J, priors_h)
  W <- rep(1 / J, J)

  # inner[[j]][[pid]] = inner PF state
  inner <- vector("list", J)
  for (j in seq_len(J)) {
    inner_j <- setNames(vector("list", length(patient_ids)), patient_ids)
    theta_j <- theta_u[j, , drop = FALSE]
    for (pid in patient_ids) {
      inner_j[[pid]] <- pf_extent_inner_init(
        Nx = Nx, theta_j = theta_j, priors = priors_h, ess_x = ess_x, a_x = a_x
      )
    }
    inner[[j]] <- inner_j
  }

  list(
    J = as.integer(J),
    Nx = as.integer(Nx),

    ess_theta = as.numeric(ess_theta),
    a_theta = as.numeric(a_theta),

    ess_x = as.numeric(ess_x),
    a_x = as.numeric(a_x),

    theta_u = theta_u,   # J x 3
    W = rep(1 / J, J),   # length J
    inner = inner,        # nested lists

    pmmh_window_n = Inf,   # PMMH window size (Inf for all data)
    priors_h = priors_h,
    patient_ids = patient_ids,
    train_cache = setNames(vector("list", length(patient_ids)), patient_ids),

    loglik_hat = rep(0.0, J), # running log p_hat(y_1:t | theta_j)
    accept = integer(J),        # PMMH accept count
    move_attempt = integer(J),   # PMMH attempt count
    proposal_sd = c(0.073, 0.145, 0.073)
  )
}

smc2_extent_resample_outer_only <- function(st) {
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


smc2_extent_pmmh_move <- function(st, n_steps = 1L, proposal_sd = c(0.05, 0.10, 0.05),
                                 n_workers = 1L, seed = NULL) {

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
          "pmmh_move_one_outer_extent",
          "propose_theta_u_extent",
          "rebuild_inner_and_loglik_extent",
          "get_pmmh_window_extent",
          "log_prior_theta_u_extent",
          "inv_logit",
          "pf_extent_inner_init",
          "pf_extent_inner_update_many",
          "pf_extent_inner_update_one",
          "pf_extent_inner_propagate",
          "pf_extent_inner_step_one_day",
          "pf_extent_inner_resample_move",
          "unpack_theta_extent_hyper",
          "extent_derived", "clip",
          "systematic_resample_idx",
          "log_sum_exp", "log_normalize", "normalize_weights"
        )
      ) %dorng% {
        # Work on ONE particle j only
        st_j <- st
        st_j <- pmmh_move_one_outer_extent(st_j, j, proposal_sd = proposal_sd)
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
        st <- pmmh_move_one_outer_extent(st, j, proposal_sd = proposal_sd)
      }
    }
  }

  st
}

# obs_df must contain Patient, TimeRel, y
# Updates all outer particles and their inner PFs, and updates outer weights.
# obs_df must contain Patient, TimeRel, y
smc2_extent_update_many <- function(st, obs_df, n_workers = 1L, seed = NULL) {
  if (is.null(obs_df) || nrow(obs_df) == 0) {
    return(list(
      st = st,
      loglik_vec = rep(0.0, st$J),
      update_time = 0,
      pmmh_time = 0,
      did_resample = FALSE
    ))
  }

  t_update0 <- Sys.time()
  pmmh_time <- 0
  did_rs <- FALSE

  stopifnot(all(c("Patient", "TimeRel", "y") %in% colnames(obs_df)))

  obs_df <- obs_df %>%
    dplyr::mutate(
      Patient = as.character(.data$Patient),
      TimeRel = as.integer(.data$TimeRel),
      y = as.integer(.data$y)
    ) %>%
    dplyr::arrange(.data$Patient, .data$TimeRel)

  obs_split <- split(obs_df, obs_df$Patient)

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


  J <- st$J
  loglik_j <- numeric(J)

  # --- worker function for a single outer particle j ---
  one_j <- function(j) {
    s <- 0.0

    thj <- unpack_theta_extent_hyper(st$theta_u[j, , drop = FALSE])
    sigma_j <- as.numeric(thj$sigma)

    inner_j <- st$inner[[j]]  # list[[pid]] for this j

    for (pid in names(obs_split)) {
      dfp <- obs_split[[pid]]
      if (nrow(dfp) == 0) next

      inner_st <- inner_j[[pid]]

      keep <- which(dfp$TimeRel > inner_st$time_last)
      if (!length(keep)) next

      rr <- pf_extent_inner_update_many(
        inner_st,
        y_vec = dfp$y[keep],
        t_vec = dfp$TimeRel[keep],
        sigma = sigma_j
      )

      inner_j[[pid]] <- rr$st
      s <- s + rr$loglik_sum
    }

    list(j = j, inner_j = inner_j, loglik = s)
  }

  n_workers <- as.integer(n_workers %||% 1L)

  if (n_workers <= 1L || foreach::getDoParWorkers() <= 1L) {
    foreach::registerDoSEQ()
  }

  if (n_workers > 1L && foreach::getDoParWorkers() > 1L) {
    if (!is.null(seed)) doRNG::registerDoRNG(seed)

    res <- foreach::foreach(
      j = seq_len(J),
      .packages = c("stats"),
      .export = c(
        # inner update chain
        "pf_extent_inner_update_many",
        "pf_extent_inner_update_one",
        "pf_extent_inner_propagate",
        "pf_extent_inner_step_one_day",
        "pf_extent_inner_resample_move",
        "unpack_theta_extent_hyper",
        "inv_logit",
        "pf_extent_inner_update_many",
        "log_sum_exp", "log_normalize", "normalize_weights",
        "systematic_resample_idx",
        "extent_derived", "clip"
      )
    ) %dorng% {
      one_j(j)
    }

    # stitch results back
    for (rr in res) {
      st$inner[[rr$j]] <- rr$inner_j
      loglik_j[rr$j] <- rr$loglik
    }

  } else {
    for (j in seq_len(J)) {
      rr <- one_j(j)
      st$inner[[j]] <- rr$inner_j
      loglik_j[j] <- rr$loglik
    }
  }

  # 1) accumulate running loglik estimate
  st$loglik_hat <- st$loglik_hat + loglik_j

  # 2) standard SMC² outer weight update (NO tempering)
  logW_raw <- log(pmax(st$W, 1e-12)) + loglik_j
  st$W <- log_normalize(logW_raw)

  ess_before <- 1 / sum(st$W^2)
  maxw_before <- max(st$W)
  cat(sprintf("[DBG][outer] ESS(before resample)=%.1f/%d thr=%.1f maxW=%.3f\n",
              ess_before, st$J, st$ess_theta * st$J, maxw_before))

  # 3) resample if ESS low
  st <- smc2_extent_resample_outer_only(st)

  did_rs <- isTRUE(attr(st, "did_resample"))
  ess_after <- 1 / sum(st$W^2)
  maxw_after <- max(st$W)
  cat(sprintf("[DBG][outer] resample=%s | ESS(after)=%.1f/%d maxW=%.3f\n",
              did_rs, ess_after, st$J, maxw_after))

  # 4) PMMH move (only if resampling happened OR always when ESS low)
  # simplest: always attempt 1 PMMH step after resampling decision
  if (isTRUE(did_rs)) {
    t_pmmh0 <- Sys.time()

    st <- smc2_extent_pmmh_move(
      st,
      n_steps = 1L,
      proposal_sd = st$proposal_sd %||% c(0.05, 0.10, 0.05),
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
# sample outer index j ~ W, then forecast from that inner PF state.
smc2_extent_forecast_pid <- function(st, pid, time_test, np_sim = 2000L) {
  pid <- as.character(pid)
  time_test <- as.integer(time_test)
  H <- length(time_test)

  stopifnot(!is.null(st$W))
  J <- length(st$W)
  stopifnot(J >= 1L)

  idx_theta <- sample.int(
    J, size = np_sim, replace = TRUE,
    prob = normalize_weights(st$W)
  )

  x_out <- matrix(NA_real_, nrow = np_sim, ncol = H)
  y_out <- matrix(NA_integer_, nrow = np_sim, ncol = H)

  ujs <- sort(unique(idx_theta))
  for (j in ujs) {
    rows <- which(idx_theta == j)
    m <- length(rows)
    if (m == 0L) next

    inner_st <- st$inner[[j]][[pid]]
    thj <- unpack_theta_extent_hyper(st$theta_u[j, , drop = FALSE])
    sigma_j <- as.numeric(thj$sigma)

    fc <- pf_extent_inner_forecast(inner_st, time_test = time_test, sigma = sigma_j, np_sim = m)

    x_out[rows, ] <- fc$x_forecast
    y_out[rows, ] <- fc$y_forecast
  }

  list(x_forecast = x_out, y_forecast = y_out)
}

