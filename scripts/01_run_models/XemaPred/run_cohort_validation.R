# -----------------------------------------------------------------------
# Run population-level XemaPred validation
#
# Purpose:
# - Run the full forward-chaining loop for population-level XemaPred.
# - Resume from existing iteration-level files where possible.
# - Use the cached SMC2 state in fits/smc2/pf_state.rds.
#
# Output layout:
# - iters/iter-XXXX.rds contains all patient forecasts for one iteration.
# - diag/diag-XXXX.rds contains one diagnostic row for one iteration.
# -----------------------------------------------------------------------

xemapred_iter_file <- function(paths, it) {
  file.path(paths$iters_dir, sprintf("iter-%04d.rds", as.integer(it)))
}

xemapred_diag_file <- function(paths, it) {
  file.path(paths$diag_dir, sprintf("diag-%04d.rds", as.integer(it)))
}

xemapred_iteration_complete <- function(paths, it) {
  file.exists(xemapred_iter_file(paths, it)) &&
    file.exists(xemapred_diag_file(paths, it)) &&
    file.exists(xemapred_posterior_file(paths, it))
}

run_xemapred_cohort_validation <- function(data_obj, run_info, paths) {
  cat(
    glue::glue(
      "[INFO] Running {toupper(run_info$model_family)} cohort XemaPred validation.\n"
    )
  )

  state <- get_or_initialise_xemapred_state(
    run_info = run_info,
    data_obj = data_obj,
    paths = paths
  )

  for (it in data_obj$iters_all) {
    it <- as.integer(it)

    if (xemapred_iteration_complete(paths, it)) {
      cat(glue::glue("[INFO] Iteration {it} already complete. Skipping.\n"))
      next
    }

    set.seed(run_info$seed + 100000L + it)

    cat("\n")
    cat(glue::glue("[INFO] Starting XemaPred cohort iteration {it}\n"))

    obs_df <- build_xemapred_cohort_obs(
      it = it,
      state = state,
      data_obj = data_obj,
      run_info = run_info
    )

    update_result <- update_xemapred_state_one_iter(
      it = it,
      state = state,
      obs_df = obs_df,
      run_info = run_info,
      paths = paths
    )

    state <- update_result$state

    save_xemapred_posterior(
      state = state,
      it = it,
      data_obj = data_obj,
      run_info = run_info,
      paths = paths
    )

    patient_outputs <- vector("list", length(data_obj$patients))
    names(patient_outputs) <- as.character(data_obj$patients)

    for (pid in data_obj$patients) {
      patient_outputs[[as.character(pid)]] <- forecast_score_xemapred_patient_iteration(
        pid = pid,
        it = it,
        state = state,
        data_obj = data_obj,
        update_result = update_result,
        run_info = run_info,
        paths = paths
      )
    }

    patient_outputs <- purrr::compact(patient_outputs)

    iter_perf <- if (length(patient_outputs)) {
      purrr::map(patient_outputs, "perf") |>
        dplyr::bind_rows()
    } else {
      tibble::tibble()
    }

    iter_diag <- make_xemapred_iteration_diag(
      run_info = run_info,
      it = it,
      update_result = update_result,
      patient_outputs = patient_outputs
    )

    saveRDS(
      iter_perf,
      file = xemapred_iter_file(paths, it)
    )

    saveRDS(
      iter_diag,
      file = xemapred_diag_file(paths, it)
    )

    cat(
      glue::glue(
        "[INFO] Saved iteration {it}: {nrow(iter_perf)} forecast rows, {nrow(iter_diag)} diagnostic rows.\n"
      )
    )

    cat(glue::glue("[INFO] Finished XemaPred cohort iteration {it}\n"))

    gc()
  }

  cat("[INFO] XemaPred cohort validation loop complete.\n")

  invisible(TRUE)
}
