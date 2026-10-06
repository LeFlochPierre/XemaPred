# -----------------------------------------------------------------------
# XemaPred forecast and scoring utilities
#
# Purpose:
# - Forecast one patient's test observations from the current SMC2 state.
# - Compute LPD, RPS, median prediction, and posterior predictive samples.
# - Return patient-level prediction rows to the cohort runner.
#
# Notes:
# - This file no longer writes per-patient files.
# - The cohort runner saves one iter-XXXX.rds and one diag-XXXX.rds per
#   forward-chaining iteration, matching the EczemaPred output layout.
# -----------------------------------------------------------------------

forecast_xemapred_patient <- function(state, pid, time_test, run_info) {
  family <- run_info$model_family

  if (family == "intensity") {
    return(
      smc2_intensity_forecast_pid(
        st = state,
        pid = pid,
        time_test = time_test,
        np_sim = run_info$n_particles
      )
    )
  }

  if (family == "subjective") {
    return(
      smc2_subjective_forecast_pid(
        st = state,
        pid = pid,
        time_test = time_test,
        np_sim = run_info$n_particles
      )
    )
  }

  if (family == "extent") {
    return(
      smc2_extent_forecast_pid(
        st = state,
        pid = pid,
        time_test = time_test,
        np_sim = run_info$n_particles
      )
    )
  }

  stop("[ERROR] Unknown XemaPred model family: ", family)
}

score_xemapred_forecast <- function(fc, test, run_info) {
  family <- run_info$model_family
  nrow_test <- nrow(test)

  lpd_vec <- numeric(nrow_test)
  rps_vec <- numeric(nrow_test)
  samples_list <- vector("list", nrow_test)
  y_pred_vec <- numeric(nrow_test)

  if (family == "intensity") {
    test_vals <- as.integer(round(test$Score / run_info$reso))
    test_vals <- pmin(run_info$M_max, pmax(0L, test_vals))
    yc_test <- test_vals + 1L

    for (i in seq_len(nrow_test)) {
      draws0 <- fc$y_forecast[, i] - 1L

      samples_list[[i]] <- draws0

      lpd_vec[i] <- compute_lpd_intensity(
        yc_true = yc_test[i],
        x_particles = fc$x_forecast[, i],
        theta_u_mat = fc$theta_forecast_u,
        M_max = run_info$M_max
      )

      rps_vec[i] <- compute_rps(
        y_true = test_vals[i],
        samples = draws0,
        M_max = run_info$M_max
      )

      y_pred_vec[i] <- median(draws0)
    }
  } else if (family == "subjective") {
    M_binom <- as.integer(run_info$M_max)

    y_test_int <- as.integer(round(test$Score / run_info$reso))
    y_test_int <- pmin(M_binom, pmax(0L, y_test_int))

    for (i in seq_len(nrow_test)) {
      draws0 <- fc$y_forecast[, i]

      samples_list[[i]] <- draws0

      lpd_vec[i] <- compute_lpd_subjective(
        y_true = y_test_int[i],
        x_particles = fc$z_forecast[, i],
        M = M_binom
      )

      rps_vec[i] <- compute_rps(
        y_true = y_test_int[i],
        samples = draws0,
        M_max = M_binom
      )

      y_pred_vec[i] <- median(draws0)
    }
  } else if (family == "extent") {
    y_test_int <- as.integer(round(test$Score))
    y_test_int <- pmin(100L, pmax(0L, y_test_int))

    for (i in seq_len(nrow_test)) {
      draws0 <- fc$y_forecast[, i]

      samples_list[[i]] <- draws0

      lpd_vec[i] <- compute_lpd_extent(
        y_true = y_test_int[i],
        x_particles = fc$x_forecast[, i]
      )

      rps_vec[i] <- compute_rps(
        y_true = y_test_int[i],
        samples = draws0,
        M_max = 100L
      )

      y_pred_vec[i] <- median(draws0)
    }
  } else {
    stop("[ERROR] Unknown XemaPred model family: ", family)
  }

  if (family %in% c("intensity", "subjective")) {
    samples_list <- lapply(samples_list, function(v) v * run_info$reso)
    y_pred_vec <- y_pred_vec * run_info$reso
  }

  list(
    lpd = lpd_vec,
    RPS = rps_vec,
    y_pred = y_pred_vec,
    Samples = samples_list
  )
}

make_xemapred_prediction_rows <- function(pid, it, test, score_obj) {
  test |>
    dplyr::mutate(
      lpd = score_obj$lpd,
      RPS = score_obj$RPS,
      y_pred = score_obj$y_pred,
      Samples = score_obj$Samples,
      xemapred_iter = as.integer(it)
    ) |>
    dplyr::select(
      -dplyr::any_of(c("LastTime", "LastScore"))
    )
}

forecast_score_xemapred_patient_iteration <- function(
    pid,
    it,
    state,
    data_obj,
    update_result,
    run_info,
    paths = NULL
) {
  set.seed(run_info$seed + 2000000L + 100000L * as.integer(pid) + as.integer(it))

  split <- split_xemapred_patient_iteration(
    data_obj = data_obj,
    pid = pid,
    it = it
  )

  test <- split$test

  if (nrow(test) == 0) {
    return(NULL)
  }

  time_test <- as.integer(test$TimeRel)

  t_fc0 <- Sys.time()

  fc <- forecast_xemapred_patient(
    state = state,
    pid = pid,
    time_test = time_test,
    run_info = run_info
  )

  t_fc1 <- Sys.time()

  t_sc0 <- Sys.time()

  score_obj <- score_xemapred_forecast(
    fc = fc,
    test = test,
    run_info = run_info
  )

  t_sc1 <- Sys.time()

  perf <- make_xemapred_prediction_rows(
    pid = pid,
    it = it,
    test = test,
    score_obj = score_obj
  )

  list(
    Patient = pid,
    perf = perf,
    forecast_time = as.numeric(difftime(t_fc1, t_fc0, units = "secs")),
    score_time = as.numeric(difftime(t_sc1, t_sc0, units = "secs")),
    n_test_rows = nrow(perf)
  )
}

make_xemapred_iteration_diag <- function(run_info, it, update_result, patient_outputs) {
  patient_outputs <- purrr::compact(patient_outputs)

  compute_time <- as.numeric(
    difftime(update_result$t_fit1, update_result$t_fit0, units = "secs")
  )

  forecast_time <- sum(
    vapply(
      patient_outputs,
      function(x) x$forecast_time %||% 0,
      numeric(1)
    ),
    na.rm = TRUE
  )

  score_time <- sum(
    vapply(
      patient_outputs,
      function(x) x$score_time %||% 0,
      numeric(1)
    ),
    na.rm = TRUE
  )

  n_test_rows <- sum(
    vapply(
      patient_outputs,
      function(x) x$n_test_rows %||% 0L,
      integer(1)
    ),
    na.rm = TRUE
  )

  st  <- update_result$state
    W   <- st$W %||% st$w
    np  <- st$J %||% NROW(st$theta_u)

    ess_after <- if (is.null(W) || !sum(W, na.rm = TRUE)) NA_real_
                else { wn <- W / sum(W); 1 / sum(wn^2) }

    n_uniq <- if (is.null(st$theta_u)) NA_integer_
              else nrow(unique(round(as.matrix(st$theta_u), 10)))

    acc <- suppressWarnings(sum(as.numeric(st$accept),       na.rm = TRUE))
    att <- suppressWarnings(sum(as.numeric(st$move_attempt), na.rm = TRUE))

    # median inner-PF ESS, deterministic subsample (no RNG consumed)
    inner_ess_pct <- NA_real_
    if (!is.null(st$inner) && length(st$inner)) {
      js <- unique(round(seq(1, length(st$inner), length.out = min(10L, length(st$inner)))))
      vals <- unlist(lapply(js, function(j) {
        inn <- st$inner[[j]]
        if (!length(inn)) return(NULL)
        ks <- unique(round(seq(1, length(inn), length.out = min(5L, length(inn)))))
        vapply(ks, function(k) {
          w <- inn[[k]]$w
          if (is.null(w) || !sum(w, na.rm = TRUE)) return(NA_real_)
          wn <- w / sum(w)
          100 / sum(wn^2) / length(wn)
        }, numeric(1))
      }))
      inner_ess_pct <- suppressWarnings(median(vals, na.rm = TRUE))
    }

  data.frame(
    dataset = run_info$dataset,
    score = run_info$score,
    model = run_info$mdl_name,
    model_family = run_info$model_family,
    Patient = NA_integer_,
    iter = as.integer(it),
    run_time = compute_time,
    compute_time = compute_time,
    update_time = update_result$update_time,
    pmmh_time = update_result$pmmh_time,
    non_pmmh_time = update_result$non_pmmh_time,
    did_resample = update_result$did_resample,
    forecast_time = forecast_time,
    score_time = score_time,
    n_test_patients = length(patient_outputs),
    n_test_rows = n_test_rows,

    seed            = run_info$seed,
    n_theta         = as.numeric(np),
    n_x             = as.numeric(st$Nx %||% NA),
    np_sim          = as.numeric(run_info$n_particles),
    ess_before      = as.numeric(update_result$ess_before %||% NA),
    ess_after       = ess_after,
    ess_pct         = 100 * ess_after / np,
    maxW            = if (is.null(W)) NA_real_ else max(W, na.rm = TRUE),
    inner_ess_pct   = inner_ess_pct,
    unique_theta    = as.numeric(n_uniq),
    unique_pct      = 100 * n_uniq / np,
    accept_cum      = acc,
    attempt_cum     = att,
    accept_pct      = 100 * acc / pmax(att, 1),
    pmmh_window     = as.numeric(st$pmmh_window_n %||% NA),
    ess_theta_thr   = as.numeric(st$ess_theta %||% NA),
    ess_x_thr       = as.numeric(st$ess_x %||% NA),
    proposal_sd     = if (is.null(st$proposal_sd)) NA_character_
                      else paste(signif(st$proposal_sd, 4), collapse = "; "),
    
    stringsAsFactors = FALSE
  )
}
