# -----------------------------------------------------------------------
# Patient-specific cohort validation loop
# -----------------------------------------------------------------------

run_one_patient <- function(pid, data_obj, run_info, paths) {
  pred_list <- list()
  diag_list <- list()

  df <- prepare_patient_dataframe(data_obj$df_all, pid, run_info$t_horizon)
  if (nrow(df) == 0L) {
    cat(glue::glue("[WARN] Patient {pid}: no rows. Skipping.\n"))
    return(list(pred = empty_prediction_df(), diag = empty_diagnostic_df()))
  }

  t0_patient <- min(df$Time)
  train_it <- get_fc_training_iteration(df[["Iteration"]])

  if (length(train_it) == 0L) {
    cat(glue::glue("[WARN] Patient {pid}: no forward-chaining iterations. Skipping.\n"))
    return(list(pred = empty_prediction_df(), diag = empty_diagnostic_df()))
  }

  iters_p <- patient_iters_dir(paths, pid)
  diag_p <- patient_diag_dir(paths, pid)

  state <- load_patient_state(paths, pid, run_info)
  missing_its <- get_missing_iterations(train_it, iters_p, run_info$save_mode, state)

  if (length(missing_its) == 0L) {
    cat(glue::glue("[INFO] Patient {pid}: nothing to do (already complete).\n"))
  } else {
    cat(glue::glue("[INFO] Patient {pid}: running {length(missing_its)} missing iterations.\n"))
  }

  for (it in sort(missing_its, decreasing = FALSE)) {
    set.seed(run_info$seed + 100000L * as.integer(pid) + as.integer(it))

    split <- split_fc_dataset(df, it)
    train <- split$Training
    test <- split$Testing

    if (nrow(train) == 0L || nrow(test) == 0L) next

    t_start_total <- Sys.time()

    t_fit0 <- Sys.time()
    if (!is.null(state)) { state$ess_min <- NULL; state$ess_sum <- 0; state$n_check <- 0L; state$n_rejuvenated <- NA_integer_ }
    state <- update_patient_state(state, train, run_info, t0_patient)
    t_fit1 <- Sys.time()

    test_info <- prepare_test_for_scoring(test, run_info, t0_patient)

    t_fc0 <- Sys.time()
    fc <- forecast_patient_iteration(state, test_info, run_info)
    t_fc1 <- Sys.time()

    t_sc0 <- Sys.time()
    perf <- score_patient_forecast(fc, test, test_info, run_info)
    t_sc1 <- Sys.time()

    pred_list[[length(pred_list) + 1L]] <- perf

    t_io0 <- Sys.time()
    if (identical(run_info$save_mode, "all")) {
      saveRDS(perf, file = file.path(iters_p, sprintf("iter-%04d.rds", it)))
    }
    t_io1 <- Sys.time()

    t_end_total <- Sys.time()

    compute_time <- elapsed_seconds(t_fit0, t_fit1)
    forecast_time <- elapsed_seconds(t_fc0, t_fc1)
    score_time <- elapsed_seconds(t_sc0, t_sc1)
    io_time <- elapsed_seconds(t_io0, t_io1)
    run_time_total <- elapsed_seconds(t_start_total, t_end_total)

    W   <- state$w
    np  <- length(W)

    ess_after <- if (is.null(W) || !sum(W, na.rm = TRUE)) NA_real_
                 else { wn <- W / sum(W); 1 / sum(wn^2) }

    n_uniq <- if (is.null(state$theta_u)) NA_integer_
              else nrow(unique(round(as.matrix(state$theta_u), 10)))

    acc <- suppressWarnings(sum(as.numeric(state$accept),       na.rm = TRUE))
    att <- suppressWarnings(sum(as.numeric(state$move_attempt), na.rm = TRUE))

    diag_vi <- data.frame(
      dataset = run_info$dataset,
      score = run_info$score,
      model = run_info$model,
      model_family = run_info$model_family,
      Patient = pid,
      iter = it,
      run_time = compute_time,
      compute_time = compute_time,
      forecast_time = forecast_time,
      score_time = score_time,
      io_time = io_time,
      run_time_total = run_time_total,

      seed            = run_info$seed,
      n_particles     = as.numeric(np),
      pmmh_nx         = as.numeric(state$pmmh_nx %||% NA),
      np_sim          = as.numeric(run_info$n_particles),

      ess_after       = ess_after,
      ess_pct         = 100 * ess_after / np,
      ess_min         = as.numeric(state$ess_min %||% NA),
      ess_min_pct     = 100 * as.numeric(state$ess_min %||% NA) / np,
      ess_mean        = if ((state$n_check %||% 0L) > 0L)
                          state$ess_sum / state$n_check else NA_real_,
      n_ess_checks    = as.numeric(state$n_check %||% NA),
      maxW            = if (is.null(W)) NA_real_ else max(W, na.rm = TRUE),

      unique_theta    = as.numeric(n_uniq),
      unique_pct      = 100 * n_uniq / np,

      accept_cum      = acc,
      attempt_cum     = att,
      accept_pct      = 100 * acc / pmax(att, 1),

      n_resample_cum  = as.numeric(state$n_resample %||% 0),
      n_moves_cum     = as.numeric(state$n_moves_cum %||% 0),
      n_rejuvenated   = as.numeric(state$n_rejuvenated %||% NA),

      pmmh_window     = as.numeric(state$pmmh_window_n %||% NA),
      ess_frac_thr    = as.numeric(state$ess_frac %||% NA),
      proposal_sd     = if (is.null(state$proposal_sd)) NA_character_
                        else paste(signif(state$proposal_sd, 4), collapse = "; "),

      prior_mode = run_info$prior_mode,
      prior_power = run_info$prior_power,
      prior_source_dataset = run_info$prior_source_dataset %||% NA_character_,
      prior_source_iteration = run_info$population_prior$source_iteration %||% NA_integer_,
      prior_source_training_day = run_info$population_prior$source_training_day %||% NA_real_,
      stringsAsFactors = FALSE
    )

    diag_list[[length(diag_list) + 1L]] <- diag_vi

    if (identical(run_info$save_mode, "all")) {
      saveRDS(diag_vi, file = file.path(diag_p, sprintf("diag-%04d.rds", it)))
    }

    save_patient_state(state, paths, pid)

    rm(list = intersect(c("test_info", "fc", "perf"), ls()))
    gc()
  }

  cat(glue::glue("[INFO] Patient {pid}: done.\n"))

  list(
    pred = dplyr::bind_rows(pred_list),
    diag = dplyr::bind_rows(diag_list)
  )
}

run_patients_parallel <- function(data_obj, run_info, paths) {
  patients <- data_obj$patients

  cl <- make_patient_cluster(run_info$n_cluster)
  on.exit(stop_patient_cluster(cl), add = TRUE)

  if (!is.null(cl)) {
    doRNG::registerDoRNG(run_info$seed)

    res_list <- foreach::foreach(
      pid = patients,
      .packages = c("dplyr", "tidyr", "purrr", "glue", "here", "yaml", "jsonlite"),
      .export = get_parallel_exports()
    ) %dopar% {
      tryCatch(
        run_one_patient(pid, data_obj, run_info, paths),
        error = function(e) {
          msg <- paste0("[ERROR] Patient ", pid, " failed: ", conditionMessage(e))
          message(msg)

          if (identical(run_info$save_mode, "all")) {
            final_p <- file.path(paths$final_dir, glue::glue("patient_{pid}"))
            dir.create(final_p, recursive = TRUE, showWarnings = FALSE)
            saveRDS(list(Patient = pid, error = conditionMessage(e)), file = file.path(final_p, "ERROR.rds"))
          }

          list(
            pred = error_prediction_df(pid),
            diag = error_diagnostic_df(pid, run_info)
          )
        }
      )
    }
  } else {
    res_list <- foreach::foreach(
      pid = patients,
      .packages = c("dplyr", "tidyr", "purrr", "glue", "here", "yaml", "jsonlite")
    ) %do% {
      run_one_patient(pid, data_obj, run_info, paths)
    }
  }

  cat("[INFO] All patients complete.\n")
  res_list
}
