# -----------------------------------------------------------------------
# XemaPred cohort state update
#
# Purpose:
# - Assimilate one iteration's new observations into the cached SMC2 state.
# - Dispatch to the correct model-family-specific SMC2 update function.
# - Return the updated state and runtime diagnostics.
# -----------------------------------------------------------------------

update_xemapred_state_one_iter <- function(it, state, obs_df, run_info, paths) {
  update_time <- NA_real_
  pmmh_time <- NA_real_
  non_pmmh_time <- NA_real_
  did_resample <- FALSE

  if (is.null(obs_df) || nrow(obs_df) == 0) {
    t_fit0 <- Sys.time()
    t_fit1 <- t_fit0

    return(
      list(
        state = state,
        t_fit0 = t_fit0,
        t_fit1 = t_fit1,
        update_time = update_time,
        pmmh_time = pmmh_time,
        non_pmmh_time = non_pmmh_time,
        did_resample = did_resample
      )
    )
  }

  family <- run_info$model_family

  t_fit0 <- Sys.time()

  if (family == "intensity") {
    t0 <- Sys.time()

    rr <- smc2_intensity_update_many(
      st = state,
      obs_df = obs_df,
      n_workers = run_info$n_cluster,
      seed = run_info$seed + 50000L + as.integer(it)
    )

    t1 <- Sys.time()

    cat(
      glue::glue(
        "[INFO] Intensity update took {round(difftime(t1, t0, units = 'secs'), 2)} seconds\n"
      )
    )

    state <- rr$st
  } else if (family == "subjective") {
    t0 <- Sys.time()

    rr <- smc2_subjective_update_many(
      st = state,
      obs_df = obs_df,
      n_workers = run_info$n_cluster,
      seed = run_info$seed + 50000L + as.integer(it)
    )

    t1 <- Sys.time()

    cat(
      glue::glue(
        "[INFO] Subjective update took {round(difftime(t1, t0, units = 'secs'), 2)} seconds\n"
      )
    )

    state <- rr$st
  } else if (family == "extent") {
    t0 <- Sys.time()

    rr <- smc2_extent_update_many(
      st = state,
      obs_df = obs_df,
      n_workers = run_info$n_cluster,
      seed = run_info$seed + 50000L + as.integer(it)
    )

    t1 <- Sys.time()

    cat(
      glue::glue(
        "[INFO] Extent update took {round(difftime(t1, t0, units = 'secs'), 2)} seconds\n"
      )
    )

    state <- rr$st
  } else {
    stop("[ERROR] Unknown XemaPred model family: ", family)
  }

  update_time <- rr$update_time %||% NA_real_
  pmmh_time <- rr$pmmh_time %||% NA_real_
  did_resample <- rr$did_resample %||% FALSE
  ess_before <- rr$ess_before %||% NA_real_

  non_pmmh_time <- ifelse(
    is.finite(update_time) & is.finite(pmmh_time),
    update_time - pmmh_time,
    NA_real_
  )

  t_fit1 <- Sys.time()

  save_xemapred_state(state, paths)

  list(
    state = state,
    t_fit0 = t_fit0,
    t_fit1 = t_fit1,
    update_time = update_time,
    pmmh_time = pmmh_time,
    non_pmmh_time = non_pmmh_time,
    did_resample = did_resample,
    ess_before = ess_before
  )
}