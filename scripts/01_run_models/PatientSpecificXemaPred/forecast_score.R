# -----------------------------------------------------------------------
# Forecasting and scoring wrappers
# -----------------------------------------------------------------------

compute_rps <- function(y_true, samples, M_max) {
  y_true <- as.integer(y_true)
  samples <- as.integer(samples)

  y_true <- max(0L, min(M_max, y_true))
  samples <- pmin(M_max, pmax(0L, samples))

  n <- length(samples)
  if (n == 0L) return(NA_real_)

  counts <- tabulate(samples + 1L, nbins = M_max + 1L)
  pmf <- counts / n
  cdf_pred <- cumsum(pmf)

  k <- 0:M_max
  cdf_true <- as.numeric(k >= y_true)

  sum((cdf_pred - cdf_true)^2)
}

prepare_test_for_scoring <- function(test, run_info, t0_patient) {
  if (identical(run_info$model_family, "intensity")) {
    time_test <- as.integer(test$Time - t0_patient)
    test_vals <- as.integer(round(test$Score / run_info$reso))
    test_vals <- pmin(run_info$M_max, pmax(0L, test_vals))
    return(list(time_test = time_test, test_vals = test_vals, yc_test = test_vals + 1L))
  }

  if (identical(run_info$model_family, "subjective")) {
    time_test <- as.integer(test$Time - t0_patient)
    test_vals <- pmin(100L, pmax(0L, as.integer(round(test$Score * 10))))
    return(list(time_test = time_test, test_vals = test_vals))
  }

  if (identical(run_info$model_family, "extent")) {
    day_test <- as.integer(test$Time - t0_patient)
    test_vals <- pmin(100L, pmax(0L, as.integer(round(test$Score))))
    return(list(day_test = day_test, test_vals = test_vals))
  }

  stop("[ERROR] Unknown model family: ", run_info$model_family)
}

forecast_patient_iteration <- function(state, test_info, run_info) {
  if (identical(run_info$model_family, "intensity")) {
    return(pf_forecast_intensity(
      x_last = state$x,
      w_last = state$w,
      theta_last_u = state$theta_u,
      time_last = state$time_last,
      time_test = test_info$time_test,
      np_sim = run_info$n_particles,
      M_max = run_info$M_max
    ))
  }

  if (identical(run_info$model_family, "subjective")) {
    return(pf_forecast_subjective_from_state(
      state,
      time_test = test_info$time_test,
      np_sim = run_info$n_particles
    ))
  }

  if (identical(run_info$model_family, "extent")) {
    return(pf_forecast_extent_from_state(
      state,
      day_test_vec = test_info$day_test,
      np_sim = run_info$n_particles
    ))
  }

  stop("[ERROR] Unknown model family: ", run_info$model_family)
}

score_patient_forecast <- function(fc, test, test_info, run_info) {
  nrow_test <- nrow(test)
  lpd_vec <- numeric(nrow_test)
  rps_vec <- numeric(nrow_test)
  samples_list <- vector("list", nrow_test)
  y_pred_vec <- numeric(nrow_test)

  if (identical(run_info$model_family, "intensity")) {
    for (i in seq_len(nrow_test)) {
      draws0 <- fc$y_forecast[, i] - 1L
      samples_list[[i]] <- draws0
      lpd_vec[i] <- compute_lpd_intensity(
        yc_true = test_info$yc_test[i],
        x_particles = fc$x_forecast[, i],
        theta_u_mat = fc$theta_forecast_u,
        M_max = run_info$M_max
      )
      rps_vec[i] <- compute_rps(test_info$test_vals[i], draws0, run_info$M_max)
      y_pred_vec[i] <- median(draws0)
    }

    samples_list <- lapply(samples_list, function(v) v * run_info$reso)
    y_pred_vec <- y_pred_vec * run_info$reso

    return(build_perf_df(test, lpd_vec, rps_vec, y_pred_vec, samples_list))
  }

  if (identical(run_info$model_family, "subjective")) {
    for (i in seq_len(nrow_test)) {
      draws <- fc$y_forecast[, i]
      samples_list[[i]] <- draws / 10
      lpd_vec[i] <- compute_lpd_subjective(test_info$test_vals[i], fc$x_forecast[, i], M = 100L)
      rps_vec[i] <- compute_rps(test_info$test_vals[i], draws, 100L)
      y_pred_vec[i] <- median(draws) / 10
    }

    return(build_perf_df(test, lpd_vec, rps_vec, y_pred_vec, samples_list))
  }

  if (identical(run_info$model_family, "extent")) {
    for (i in seq_len(nrow_test)) {
      draws <- fc$y_forecast[, i]
      samples_list[[i]] <- draws
      lpd_vec[i] <- compute_lpd_extent(test_info$test_vals[i], fc$x_forecast[, i])
      rps_vec[i] <- compute_rps(test_info$test_vals[i], draws, 100L)
      y_pred_vec[i] <- median(draws)
    }

    return(build_perf_df(test, lpd_vec, rps_vec, y_pred_vec, samples_list))
  }

  stop("[ERROR] Unknown model family: ", run_info$model_family)
}

build_perf_df <- function(test, lpd_vec, rps_vec, y_pred_vec, samples_list) {
  out <- test %>%
    dplyr::mutate(
      lpd = lpd_vec,
      RPS = rps_vec,
      y_pred = y_pred_vec,
      Samples = samples_list
    )

  cols_to_drop <- intersect(c("LastTime", "LastScore"), colnames(out))
  if (length(cols_to_drop)) out <- dplyr::select(out, -dplyr::all_of(cols_to_drop))

  out
}
