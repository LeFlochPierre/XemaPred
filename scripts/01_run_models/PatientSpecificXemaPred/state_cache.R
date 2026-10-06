# -----------------------------------------------------------------------
# Patient-specific PF state cache
# -----------------------------------------------------------------------

load_patient_state <- function(
    paths,
    pid,
    run_info
) {

  f <- patient_state_file(
    paths,
    pid
  )

  if (!file.exists(f)) {
    return(NULL)
  }

  state <- safe_read_rds(f)

  if (is.null(state)) {
    return(NULL)
  }

  if (!validate_patient_state(
    state,
    run_info
  )) {

    cat(
      glue::glue(
        "[WARN] Patient {pid}: ",
        "cache state incompatible; resetting.\n"
      )
    )

    return(NULL)
  }

  # ------------------------------------------------------------
  # Refresh runtime settings
  # ------------------------------------------------------------

  state$prior_mode <-
    run_info$prior_mode

  state$prior_combination <-
    run_info$prior_combination

  state$prior_power <-
    run_info$prior_power

  state$population_prior <-
    run_info$population_prior

  state$prior_id <-
    run_info$population_prior$prior_id %||%
    "base"

  state$pmmh_prior_mode <-
    run_info$pmmh_prior_mode

  state$ess_frac <-
    run_info$ess_frac

  state$pmmh_window_n <-
    run_info$pmmh_window_n

  state
}


save_patient_state <- function(
    state,
    paths,
    pid
) {

  saveRDS(
    state,
    file = patient_state_file(
      paths,
      pid
    )
  )
}


validate_patient_state <- function(
    state,
    run_info
) {

  if (!is.list(state)) {
    return(FALSE)
  }

  model_family <-
    run_info$model_family

  np <-
    run_info$n_particles

  # ------------------------------------------------------------
  # Validate prior configuration
  #
  # Old cache states without these fields are interpreted as
  # having been generated using the original broad prior.
  # ------------------------------------------------------------

  state_prior_mode <-
    state$prior_mode %||%
    "base"

  run_prior_mode <-
    run_info$prior_mode %||%
    "base"


  # ------------------------------------------------------------
  # Prior combination
  #
  # Old population caches were generated using the original
  # geometric/power formulation.
  # ------------------------------------------------------------

  state_prior_combination <-
    state$prior_combination %||%
    if (identical(
      state_prior_mode,
      "population"
    )) {
      "power"
    } else {
      "base"
    }

  run_prior_combination <-
    run_info$prior_combination %||%
    if (identical(
      run_prior_mode,
      "population"
    )) {
      "power"
    } else {
      "base"
    }


  # ------------------------------------------------------------
  # Prior strength
  # ------------------------------------------------------------

  state_prior_power <-
    state$prior_power %||%
    0

  run_prior_power <-
    run_info$prior_power %||%
    0


  # ------------------------------------------------------------
  # Population-prior identity
  # ------------------------------------------------------------

  state_prior_id <-
    state$prior_id %||%
    "base"

  run_prior_id <-
    run_info$population_prior$prior_id %||%
    "base"


  # ------------------------------------------------------------
  # PMMH prior target
  # ------------------------------------------------------------

  state_pmmh_prior_mode <-
    state$pmmh_prior_mode %||%
    "matched"

  run_pmmh_prior_mode <-
    run_info$pmmh_prior_mode %||%
    "matched"


  # ------------------------------------------------------------
  # Cache is reusable only if ALL prior settings match
  # ------------------------------------------------------------

  prior_ok <-
    identical(
      state_prior_mode,
      run_prior_mode
    ) &&

    identical(
      state_prior_combination,
      run_prior_combination
    ) &&

    isTRUE(
      all.equal(
        state_prior_power,
        run_prior_power
      )
    ) &&

    identical(
      state_prior_id,
      run_prior_id
    ) &&

    identical(
      state_pmmh_prior_mode,
      run_pmmh_prior_mode
    )

  if (!prior_ok) {
    return(FALSE)
  }

  # ------------------------------------------------------------
  # Intensity
  # ------------------------------------------------------------

  if (identical(
    model_family,
    "intensity"
  )) {

    return(
      !is.null(state$M_max) &&
        state$M_max == run_info$M_max &&

        !is.null(state$np) &&
        state$np == np &&

        !is.null(state$x) &&
        length(state$x) == np &&

        !is.null(state$w) &&
        length(state$w) == np &&

        !is.null(state$theta_u) &&
        nrow(state$theta_u) == np &&

        !is.null(state$time_last) &&
        !is.null(state$train_cache) &&
        !is.null(state$priors) &&
        !is.null(state$pmmh_nx) &&
        !is.null(state$proposal_sd) &&

        !is.null(state$accept) &&
        length(state$accept) == np &&

        !is.null(state$move_attempt) &&
        length(state$move_attempt) == np
    )
  }

  # ------------------------------------------------------------
  # Extent
  # ------------------------------------------------------------

  if (identical(
    model_family,
    "extent"
  )) {

    return(
      !is.null(state$np) &&
        state$np == np &&

        !is.null(state$y_lat) &&
        length(state$y_lat) == np &&

        !is.null(state$logit_tss1) &&
        length(state$logit_tss1) == np &&

        !is.null(state$w) &&
        length(state$w) == np &&

        !is.null(state$theta_u) &&
        nrow(state$theta_u) == np &&

        !is.null(state$time_last) &&
        !is.null(state$train_cache) &&
        !is.null(state$priors) &&
        !is.null(state$pmmh_nx) &&
        !is.null(state$proposal_sd) &&

        !is.null(state$accept) &&
        length(state$accept) == np &&

        !is.null(state$move_attempt) &&
        length(state$move_attempt) == np
    )
  }

  # ------------------------------------------------------------
  # Subjective
  # ------------------------------------------------------------

  if (identical(
    model_family,
    "subjective"
  )) {

    return(
      !is.null(state$np) &&
        state$np == np &&

        !is.null(state$x) &&
        length(state$x) == np &&

        !is.null(state$w) &&
        length(state$w) == np &&

        !is.null(state$theta_u) &&
        nrow(state$theta_u) == np &&

        !is.null(state$time_last) &&
        !is.null(state$M) &&
        !is.null(state$train_cache) &&
        !is.null(state$priors) &&
        !is.null(state$pmmh_nx) &&
        !is.null(state$proposal_sd) &&

        !is.null(state$accept) &&
        length(state$accept) == np &&

        !is.null(state$move_attempt) &&
        length(state$move_attempt) == np
    )
  }

  FALSE
}