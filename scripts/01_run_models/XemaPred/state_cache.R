# -----------------------------------------------------------------------
# XemaPred SMC2 state-cache utilities
#
# Purpose:
# - Select model-family-specific priors.
# - Load an existing cached SMC2 state if available.
# - Validate that the cached state matches the current score/model.
# - Initialise a new SMC2 state when no valid cache exists.
# - Refresh runtime settings from the YAML config after restart.
#
# Used by:
# - XemaPred/run_validation.R
# - XemaPred/run_cohort_validation.R
# -----------------------------------------------------------------------

get_xemapred_priors <- function(run_info) {
  family <- run_info$model_family

  if (family == "intensity") {
    return(get_priors_intensity(run_info$M_max))
  }

  if (family == "subjective") {
    return(get_priors_subjective_smc2())
  }

  if (family == "extent") {
    return(get_priors_extent_hyper())
  }

  stop("[ERROR] Unknown XemaPred model family: ", family)
}

load_xemapred_state <- function(paths) {
  if (!file.exists(paths$cohort_state_file)) {
    cat("[INFO] No cached SMC2 state found.\n")
    return(NULL)
  }

  cat("[INFO] Loading cached SMC2 state:\n")
  cat("  ", paths$cohort_state_file, "\n", sep = "")

  state <- tryCatch(
    readRDS(paths$cohort_state_file),
    error = function(e) {
      warning(
        "[WARN] Failed to read cached SMC2 state: ",
        conditionMessage(e)
      )
      NULL
    }
  )

  state
}

validate_xemapred_state <- function(state, run_info) {
  if (is.null(state)) {
    return(FALSE)
  }

  family <- run_info$model_family

  if (family == "intensity") {
    return(
      is.list(state) &&
        !is.null(state$J) &&
        !is.null(state$Nx) &&
        !is.null(state$theta_u) &&
        nrow(state$theta_u) == state$J &&
        !is.null(state$W) &&
        length(state$W) == state$J &&
        !is.null(state$inner) &&
        length(state$inner) == state$J &&
        !is.null(state$M_max) &&
        as.integer(state$M_max) == as.integer(run_info$M_max)
    )
  }

  if (family == "subjective") {
    M_binom <- as.integer(run_info$M_max)

    return(
      is.list(state) &&
        !is.null(state$J) &&
        !is.null(state$Nx) &&
        !is.null(state$theta_u) &&
        nrow(state$theta_u) == state$J &&
        !is.null(state$W) &&
        length(state$W) == state$J &&
        !is.null(state$inner) &&
        length(state$inner) == state$J &&
        !is.null(state$M) &&
        as.integer(state$M) == M_binom
    )
  }

  if (family == "extent") {
    return(
      is.list(state) &&
        !is.null(state$J) &&
        length(state$J) == 1L &&
        !is.null(state$Nx) &&
        length(state$Nx) == 1L &&
        !is.null(state$theta_u) &&
        nrow(state$theta_u) == state$J &&
        !is.null(state$W) &&
        length(state$W) == state$J &&
        !is.null(state$inner) &&
        length(state$inner) == state$J
    )
  }

  stop("[ERROR] Unknown XemaPred model family: ", family)
}

initialise_xemapred_state <- function(run_info, data_obj, priors) {
  family <- run_info$model_family

  cat("[INFO] Initialising new SMC2 state for family: ", family, "\n", sep = "")

  if (family == "intensity") {
    return(
      smc2_intensity_init(
        J = run_info$n_theta,
        Nx = run_info$n_x,
        priors = priors,
        M_max = run_info$M_max,
        patient_ids = data_obj$patients,
        ess_theta = run_info$ess_theta,
        a_theta = run_info$a_theta,
        ess_x = run_info$ess_x,
        a_x = run_info$a_x
      )
    )
  }

  if (family == "subjective") {
    M_binom <- as.integer(run_info$M_max)

    return(
      smc2_subjective_init(
        J = run_info$n_theta,
        Nx = run_info$n_x,
        priors = priors,
        M = M_binom,
        patient_ids = data_obj$patients,
        ess_theta = run_info$ess_theta,
        a_theta = run_info$a_theta,
        ess_x = run_info$ess_x,
        a_x = run_info$a_x
      )
    )
  }

  if (family == "extent") {
    return(
      smc2_extent_init(
        J = run_info$n_theta,
        Nx = run_info$n_x,
        priors_h = priors,
        patient_ids = data_obj$patients,
        ess_theta = run_info$ess_theta,
        a_theta = run_info$a_theta,
        ess_x = run_info$ess_x,
        a_x = run_info$a_x
      )
    )
  }

  stop("[ERROR] Unknown XemaPred model family: ", family)
}

refresh_xemapred_state_settings <- function(state, run_info) {
  state$ess_theta <- run_info$ess_theta
  state$a_theta <- run_info$a_theta

  state$ess_x <- run_info$ess_x
  state$a_x <- run_info$a_x

  state$pmmh_window_n <- run_info$pmmh_window_n

  state
}

get_or_initialise_xemapred_state <- function(run_info, data_obj, paths) {
  priors <- get_xemapred_priors(run_info)

  state <- load_xemapred_state(paths)

  if (!validate_xemapred_state(state, run_info)) {
    if (!is.null(state)) {
      cat("[WARN] Cached SMC2 state is incompatible with current config. Reinitialising.\n")
    }

    state <- initialise_xemapred_state(
      run_info = run_info,
      data_obj = data_obj,
      priors = priors
    )
  } else {
    cat("[INFO] Cached SMC2 state is valid. Resuming from cache.\n")
  }

  state <- refresh_xemapred_state_settings(
    state = state,
    run_info = run_info
  )

  state
}

save_xemapred_state <- function(state, paths) {
  dir.create(
    dirname(paths$cohort_state_file),
    recursive = TRUE,
    showWarnings = FALSE
  )

  saveRDS(state, paths$cohort_state_file)

  invisible(paths$cohort_state_file)
}