# ============================================================
# subjective_smc2_cache_functions.R
#
# SUBJECTIVE SMC² (PF inside PF)
#
# Model (per patient):
#   y(t) ~ Bin(M, p(t))
#   z(t) = logit(p(t))
#   z(t+1) ~ N(z(t), sigma^2)
#
# Global (outer) parameters:
#   sigma  ~ half-normal(0, 0.25*log(5))
#   mu0    ~ N(0, 1)
#   sigma0 ~ half-normal(0, 1.5)
#
# Inner PF state per patient:
#   z particles (Nx), weights, time_last
#
# Depends on shared utilities:
#   inv_logit, log_normalize, log_sum_exp, normalize_weights,
#   systematic_resample_idx, liu_west_move_select
# ============================================================

suppressPackageStartupMessages({
  library(stats)
  library(foreach)
  library(doRNG)
})
# -------------------------
# priors + packing
# -------------------------
get_priors_subjective_smc2 <- function() {
  list(
    sigma_sd = 0.25 * log(5),
    mu0_mean = 0.0,
    mu0_sd = 1.0,
    sigma0_sd = 1.5
  )
}

log_prior_theta_u_subj <- function(theta_u, priors) {
  th <- unpack_theta_subj_smc2(theta_u)

  sigma  <- as.numeric(th$sigma[1])
  mu0    <- as.numeric(th$mu0[1])
  sigma0 <- as.numeric(th$sigma0[1])

  # half-normal priors with Jacobian for log_sigma, log_sigma0
  # sigma ~ HalfNormal(0, priors$sigma_sd)
  lp_sigma  <- dnorm(sigma, 0, priors$sigma_sd, log = TRUE) + log(2) + as.numeric(theta_u[1, "log_sigma"])
  # sigma0 ~ HalfNormal(0, priors$sigma0_sd)
  lp_sigma0 <- dnorm(sigma0, 0, priors$sigma0_sd, log = TRUE) + log(2) + as.numeric(theta_u[1, "log_sigma0"])
  # mu0 ~ Normal(mu0_mean, mu0_sd)
  lp_mu0 <- dnorm(mu0, priors$mu0_mean, priors$mu0_sd, log = TRUE)

  as.numeric(lp_sigma + lp_mu0 + lp_sigma0)
}

propose_theta_u_subj <- function(theta_u_curr, proposal_sd) {
  th <- as.matrix(theta_u_curr)
  th[1, "log_sigma"]  <- th[1, "log_sigma"]  + rnorm(1, 0, proposal_sd[1])
  th[1, "mu0"]        <- th[1, "mu0"]        + rnorm(1, 0, proposal_sd[2])
  th[1, "log_sigma0"] <- th[1, "log_sigma0"] + rnorm(1, 0, proposal_sd[3])
  th
}

get_pmmh_window_subj <- function(dfp, window_n = Inf) {
  if (is.null(dfp) || nrow(dfp) == 0) return(dfp)

  dfp <- dfp[order(dfp$TimeRel), , drop = FALSE]

  if (is.finite(window_n)) {
    n_keep <- as.integer(window_n)
    if (nrow(dfp) > n_keep) {
      dfp <- dfp[(nrow(dfp) - n_keep + 1L):nrow(dfp), , drop = FALSE]
    }
  }

  dfp
}

rebuild_inner_and_loglik_subj <- function(st, theta_u_one, window_n = Inf) {
  theta_u_one <- as.matrix(theta_u_one)

  inner_new <- setNames(vector("list", length(st$patient_ids)), st$patient_ids)
  ll_sum <- 0.0

  for (pid in st$patient_ids) {
    inn <- pf_subj_inner_init(
      Nx = st$Nx,
      theta_j = theta_u_one,
      priors = st$priors,
      M = st$M,
      ess_x = st$ess_x,
      a_x = st$a_x
    )

    dfp <- get_pmmh_window_subj(
      st$train_cache[[pid]],
      window_n = window_n
    )

    if (!is.null(dfp) && nrow(dfp) > 0) {
      rr <- pf_subj_inner_update_many(
        st = inn,
        theta_j = theta_u_one,
        y_vec = dfp$y,
        t_vec = dfp$TimeRel
      )

      inn <- rr$st
      ll_sum <- ll_sum + rr$loglik_sum
    }

    inner_new[[pid]] <- inn
  }

  list(inner = inner_new, loglik_hat = ll_sum)
}

pmmh_move_one_outer_subj <- function(st, j, proposal_sd = NULL) {
  window_n <- st$pmmh_window_n %||% Inf

  st$move_attempt[j] <- st$move_attempt[j] + 1L
  if (is.null(proposal_sd)) proposal_sd <- st$proposal_sd

  theta_curr <- st$theta_u[j, , drop = FALSE]

  curr <- rebuild_inner_and_loglik_subj(
    st,
    theta_u_one = theta_curr,
    window_n = window_n
  )

  ll_curr <- curr$loglik_hat
  lp_curr <- log_prior_theta_u_subj(theta_curr, st$priors)

  theta_prop <- propose_theta_u_subj(theta_curr, proposal_sd)

  prop <- rebuild_inner_and_loglik_subj(
    st,
    theta_u_one = theta_prop,
    window_n = window_n
  )

  ll_prop <- prop$loglik_hat
  lp_prop <- log_prior_theta_u_subj(theta_prop, st$priors)

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


smc2_subjective_pmmh_move <- function(st, n_steps = 1L, n_workers = 1L, seed = NULL) {
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
          "pmmh_move_one_outer_subj",
          "propose_theta_u_subj",
          "rebuild_inner_and_loglik_subj",
          "get_pmmh_window_subj",
          "log_prior_theta_u_subj",
          "pf_subj_inner_init",
          "pf_subj_inner_update_many",
          "pf_subj_inner_update_one",
          "pf_subj_inner_propagate",
          "pf_subj_inner_resample",
          "unpack_theta_subj_smc2",
          "as_matrix",
          "inv_logit",
          "log_sum_exp",
          "log_normalize",
          "normalize_weights",
          "systematic_resample_idx"
        )
      ) %dorng% {
        st_j <- st

        st_j <- pmmh_move_one_outer_subj(
          st_j,
          j,
          proposal_sd = st$proposal_sd
        )

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
        st <- pmmh_move_one_outer_subj(
          st,
          j,
          proposal_sd = st$proposal_sd
        )
      }
    }
  }

  st
}


theta_colnames_subj <- function() c("log_sigma", "mu0", "log_sigma0")

rhalfnorm <- function(n, sd) abs(rnorm(n, 0, sd))

init_theta_particles_subj_smc2 <- function(J, priors) {
  sigma  <- pmax(rhalfnorm(J, priors$sigma_sd),  1e-12)
  mu0    <- rnorm(J, priors$mu0_mean, priors$mu0_sd)
  sigma0 <- pmax(rhalfnorm(J, priors$sigma0_sd), 1e-12)
  list(sigma = sigma, mu0 = mu0, sigma0 = sigma0)
}

pack_theta_subj_smc2 <- function(sigma, mu0, sigma0) {
  out <- cbind(
    log_sigma  = log(pmax(sigma,  1e-12)),
    mu0        = mu0,
    log_sigma0 = log(pmax(sigma0, 1e-12))
  )
  colnames(out) <- theta_colnames_subj()
  out
}

as_matrix <- function(x) {
  if (is.matrix(x)) return(x)
  x <- as.matrix(x)
  if (is.null(nrow(x))) x <- matrix(x, nrow = 1)
  x
}

unpack_theta_subj_smc2 <- function(theta_u) {
  theta_u <- as_matrix(theta_u)
  list(
    sigma  = exp(theta_u[, "log_sigma"]),
    mu0    = theta_u[, "mu0"],
    sigma0 = exp(theta_u[, "log_sigma0"])
  )
}

# -------------------------
# inner PF
# -------------------------
pf_subj_inner_init <- function(Nx, theta_j, priors, M, ess_x = 0.5, a_x = 0.98) {
  th <- unpack_theta_subj_smc2(theta_j)
  z <- rnorm(Nx, mean = th$mu0[1], sd = th$sigma0[1])

  list(
    Nx = as.integer(Nx),
    ess_x = as.numeric(ess_x),
    a_x = as.numeric(a_x),
    M = as.integer(M),
    z = z,
    w = rep(1 / Nx, Nx),
    time_last = -1L
  )
}

pf_subj_inner_resample <- function(st) {
  w <- st$w
  Nx <- st$Nx
  ess <- 1 / sum(w^2)

  if (ess < st$ess_x * Nx) {
    idx <- systematic_resample_idx(w)
    st$z <- st$z[idx]
    st$w <- rep(1 / Nx, Nx)
  }
  st
}

pf_subj_inner_propagate <- function(st, dt, sigma) {
  dt <- as.integer(dt)
  if (dt <= 0L) return(st)

  st$z <- st$z + sqrt(dt) * sigma * rnorm(st$Nx)
  st
}

pf_subj_inner_update_one <- function(st, theta_j, y_t, t_now) {
  t_now <- as.integer(t_now)
  t_prev <- st$time_last
  if (is.null(t_prev) || is.na(t_prev)) t_prev <- -1L
  t_prev <- as.integer(t_prev)

  th <- unpack_theta_subj_smc2(theta_j)
  sigma_j <- as.numeric(th$sigma[1])

  dt <- if (t_prev < 0L) 0L else max(0L, t_now - t_prev)
  st <- pf_subj_inner_propagate(st, dt, sigma = sigma_j)

  p <- inv_logit(st$z)
  pmf <- dbinom(x = as.integer(y_t), size = st$M, prob = p)
  pmf <- pmax(pmf, 1e-12)

  loglik_inc <- log_sum_exp(log(pmax(st$w, 1e-12)) + log(pmf))
  st$w <- log_normalize(log(pmax(st$w, 1e-12)) + log(pmf))

  st <- pf_subj_inner_resample(st)
  st$time_last <- t_now

  list(st = st, loglik_inc = loglik_inc)
}

pf_subj_inner_update_many <- function(st, theta_j, y_vec, t_vec) {
  stopifnot(length(y_vec) == length(t_vec))
  if (!length(y_vec)) return(list(st = st, loglik_sum = 0.0))

  s <- 0.0
  for (i in seq_along(y_vec)) {
    rr <- pf_subj_inner_update_one(st, theta_j, y_t = y_vec[i], t_now = t_vec[i])
    st <- rr$st
    s <- s + rr$loglik_inc
  }
  list(st = st, loglik_sum = s)
}

reset_or_rebuild_inner_filters_subj <- function(st, mode = c("rebuild","reset")) {
  mode <- match.arg(mode)
  has_cache <- !is.null(st$train_cache) &&
    any(vapply(st$train_cache, function(x) !is.null(x) && nrow(x) > 0, logical(1)))
  if (mode == "rebuild" && !has_cache) mode <- "reset"

  for (j in seq_len(st$J)) {
    theta_j <- st$theta_u[j, , drop = FALSE]
    pids_j <- names(st$inner[[j]])
    if (mode == "rebuild") pids_j <- intersect(pids_j, names(st$train_cache))

    for (pid in pids_j) {
      inner_new <- pf_subj_inner_init(
        Nx = st$Nx, theta_j = theta_j, priors = st$priors, M = st$M,
        ess_x = st$ess_x, a_x = st$a_x
      )

      if (mode == "rebuild") {
        dfp <- st$train_cache[[pid]]
        if (!is.null(dfp) && nrow(dfp) > 0) {
          rr <- pf_subj_inner_update_many(
            st = inner_new, theta_j = theta_j,
            y_vec = dfp$y, t_vec = dfp$TimeRel
          )
          inner_new <- rr$st
        }
      }

      st$inner[[j]][[pid]] <- inner_new
    }
  }
  st
}


pf_subj_inner_forecast <- function(st, theta_j, time_test, np_sim = 2000L) {
  time_test <- as.integer(time_test)
  H <- length(time_test)

  th <- unpack_theta_subj_smc2(theta_j)
  sigma_j <- as.numeric(th$sigma[1])

  w_use <- normalize_weights(st$w)
  idx <- sample.int(st$Nx, size = np_sim, replace = TRUE, prob = w_use)
  z <- st$z[idx]

  z_fc <- matrix(NA_real_, nrow = np_sim, ncol = H)
  y_fc <- matrix(NA_integer_, nrow = np_sim, ncol = H)

  t_prev <- st$time_last
  if (is.null(t_prev) || is.na(t_prev)) t_prev <- -1L
  t_prev <- as.integer(t_prev)

  for (h in seq_len(H)) {
    t_now <- as.integer(time_test[h])
    dt <- if (t_prev < 0L) 0L else max(0L, t_now - t_prev)

    if (dt > 0L) {
      z <- z + sqrt(dt) * sigma_j * rnorm(np_sim)
    }

    z_fc[, h] <- z
    p <- inv_logit(z)
    y_fc[, h] <- rbinom(np_sim, size = st$M, prob = p)

    t_prev <- t_now
  }

  list(z_forecast = z_fc, y_forecast = y_fc)
}

# -------------------------
# outer SMC² state
# -------------------------
smc2_subjective_init <- function(J, Nx, priors, M, patient_ids,
                                 ess_theta = 0.5, a_theta = 0.98,
                                 ess_x = 0.5, a_x = 0.98) {
  patient_ids <- as.character(patient_ids)

  th0 <- init_theta_particles_subj_smc2(J, priors)
  theta_u <- pack_theta_subj_smc2(th0$sigma, th0$mu0, th0$sigma0)
  theta_u <- as.matrix(theta_u)

  W <- rep(1 / J, J)

  inner <- vector("list", J)
  for (j in seq_len(J)) {
    inner_j <- setNames(vector("list", length(patient_ids)), patient_ids)
    theta_j <- theta_u[j, , drop = FALSE]
    for (pid in patient_ids) {
      inner_j[[pid]] <- pf_subj_inner_init(
        Nx = Nx, theta_j = theta_j, priors = priors, M = M,
        ess_x = ess_x, a_x = a_x
      )
    }
    inner[[j]] <- inner_j
  }

  list(
    J = as.integer(J),
    Nx = as.integer(Nx),
    M = as.integer(M),

    ess_theta = as.numeric(ess_theta),
    a_theta = as.numeric(a_theta),

    ess_x = as.numeric(ess_x),
    a_x = as.numeric(a_x),

    theta_u = theta_u,  # J x 3
    W = W,              # length J
    inner = inner,       # list[[j]][[pid]]

    pmmh_window_n = Inf,  # PMMH window size (Inf for all data)
    priors = priors,
    patient_ids = patient_ids,
    train_cache = setNames(vector("list", length(patient_ids)), patient_ids),

    loglik_hat = rep(0.0, J),  # loglik estimate per outer particle (updated in PMMH moves)
    accept = integer(J),       # count of accepted PMMH moves per outer particle
    move_attempt = integer(J), # count of attempted PMMH moves per outer particle
    proposal_sd = c(0.043, 0.145, 0.043) # PMMH proposal sd for each parameter
  )
}

smc2_subjective_resample_outer_only <- function(st) {
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


# obs_df: Patient, TimeRel, y (0..M)
smc2_subjective_update_many <- function(st, obs_df, n_workers = 1L, seed = NULL) {
  if (is.null(obs_df) || nrow(obs_df) == 0)
    return(list(st = st, loglik_vec = rep(0.0, st$J)))
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

  one_j <- function(j) {
    theta_j <- st$theta_u[j, , drop = FALSE]
    inner_j <- st$inner[[j]]
    s <- 0.0

    for (pid in names(obs_split)) {
      dfp <- obs_split[[pid]]
      if (nrow(dfp) == 0) next

      inner_st <- inner_j[[pid]]

      keep <- which(dfp$TimeRel > inner_st$time_last)
      if (!length(keep)) next

      rr <- pf_subj_inner_update_many(
        inner_st, theta_j,
        y_vec = dfp$y[keep],
        t_vec = dfp$TimeRel[keep]
      )

      inner_j[[pid]] <- rr$st
      s <- s + rr$loglik_sum
    }

    list(j = j, inner_j = inner_j, loglik = s)
  }

  n_workers <- as.integer(n_workers %||% 1L)

  if (n_workers > 1L && foreach::getDoParWorkers() > 1L) {
    if (!is.null(seed)) doRNG::registerDoRNG(seed)

    res <- foreach::foreach(
      j = seq_len(J),
      .packages = c("stats"),
      .export = c(
        "pf_subj_inner_update_many",
        "pf_subj_inner_update_one",
        "pf_subj_inner_propagate",
        "pf_subj_inner_resample",
        "unpack_theta_subj_smc2",
        "as_matrix",
        "inv_logit", "log_sum_exp", "log_normalize",
        "systematic_resample_idx", "normalize_weights"
      )
    ) %dorng% {
      one_j(j)
    }

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

  # 1) accumulate likelihood estimate per outer particle
  st$loglik_hat <- st$loglik_hat + loglik_j

  # 2) standard outer weight update
  logW_raw <- log(pmax(st$W, 1e-12)) + loglik_j
  st$W <- log_normalize(logW_raw)

  ess_before <- 1 / sum(st$W^2)
  maxw_before <- max(st$W)
  cat(sprintf("[DBG][outer] ESS(before resample)=%.1f/%d thr=%.1f maxW=%.3f\n",
              ess_before, st$J, st$ess_theta * st$J, maxw_before))

  # 3) resample outer only if ESS low
  st <- smc2_subjective_resample_outer_only(st)

  did_rs <- isTRUE(attr(st, "did_resample"))
  ess_after <- 1 / sum(st$W^2)
  maxw_after <- max(st$W)
  cat(sprintf("[DBG][outer] resample=%s | ESS(after)=%.1f/%d maxW=%.3f\n",
              did_rs, ess_after, st$J, maxw_after))

  # 4) PMMH rejuvenation only if resampled
  did_rs <- isTRUE(attr(st, "did_resample"))

  if (did_rs) {
    t_pmmh0 <- Sys.time()

    st <- smc2_subjective_pmmh_move(
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


# mixture forecast per patient:
# sample outer index j ~ W, then forecast from that inner PF
# returns theta_forecast_u so scoring code can reuse it
smc2_subjective_forecast_pid <- function(st, pid, time_test, np_sim = 2000L) {
  pid <- as.character(pid)
  time_test <- as.integer(time_test)
  H <- length(time_test)

  J <- length(st$W)
  idx_theta <- sample.int(J, size = np_sim, replace = TRUE,
                          prob = normalize_weights(st$W))

  y_out <- matrix(NA_integer_, nrow = np_sim, ncol = H)
  z_out <- matrix(NA_real_, nrow = np_sim, ncol = H)

  theta_out <- matrix(NA_real_, nrow = np_sim, ncol = ncol(st$theta_u))
  colnames(theta_out) <- colnames(st$theta_u)

  ujs <- sort(unique(idx_theta))
  for (j in ujs) {
    rows <- which(idx_theta == j)
    m <- length(rows)
    if (m == 0L) next

    theta_j <- st$theta_u[j, , drop = FALSE]
    inner_st <- st$inner[[j]][[pid]]

    fc <- pf_subj_inner_forecast(inner_st, theta_j, time_test = time_test, np_sim = m)

    z_out[rows, ] <- fc$z_forecast
    y_out[rows, ] <- fc$y_forecast
    theta_out[rows, ] <- st$theta_u[rep(j, m), , drop = FALSE]
  }

  list(z_forecast = z_out, y_forecast = y_out, theta_forecast_u = theta_out)
}
