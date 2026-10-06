# -----------------------------------------------------------------------
# Family-specific state initialization and update wrappers
# -----------------------------------------------------------------------

init_patient_state <- function(run_info) {

  if (identical(run_info$model_family, "intensity")) {

    return(
      pf_intensity_init(
        np = run_info$n_particles,
        priors = get_priors_intensity(run_info$M_max),
        M_max = run_info$M_max,
        ess_frac = run_info$ess_frac,
        pmmh_nx = run_info$pmmh_nx,
        pmmh_window_n = run_info$pmmh_window_n,
        proposal_sd = c(0.074, 0.074, 0.25, 0.074, 0.074),

        prior_mode = run_info$prior_mode,
        prior_combination = run_info$prior_combination,
        prior_power = run_info$prior_power,
        population_prior = run_info$population_prior,
        pmmh_prior_mode = run_info$pmmh_prior_mode
      )
    )
  }


  if (identical(run_info$model_family, "subjective")) {

    return(
      pf_subjective_init(
        np = run_info$n_particles,
        priors = get_priors_subjective(),
        ess_frac = run_info$ess_frac,
        M = 100L,
        pmmh_nx = run_info$pmmh_nx,
        pmmh_window_n = run_info$pmmh_window_n,
        proposal_sd = c(0.03, 0.10, 0.03),

        prior_mode = run_info$prior_mode,
        prior_combination = run_info$prior_combination,
        prior_power = run_info$prior_power,
        population_prior = run_info$population_prior,
        pmmh_prior_mode = run_info$pmmh_prior_mode
      )
    )
  }


  if (identical(run_info$model_family, "extent")) {

    return(
      pf_extent_init(
        np = run_info$n_particles,
        priors = get_priors_extent(),
        ess_frac = run_info$ess_frac,
        pmmh_nx = run_info$pmmh_nx,
        pmmh_window_n = run_info$pmmh_window_n,
        proposal_sd = c(0.05, 0.10),

        prior_mode = run_info$prior_mode,
        prior_combination = run_info$prior_combination,
        prior_power = run_info$prior_power,
        population_prior = run_info$population_prior,
        pmmh_prior_mode = run_info$pmmh_prior_mode
      )
    )
  }


  stop(
    "[ERROR] Unknown model family: ",
    run_info$model_family
  )
}

update_patient_state <- function(state, train, run_info, t0_patient) {
  if (is.null(state)) state <- init_patient_state(run_info)

  if (identical(run_info$model_family, "intensity")) {
    time_train <- as.integer(train$Time - t0_patient)
    train_vals <- as.integer(round(train$Score / run_info$reso))
    train_vals <- pmin(run_info$M_max, pmax(0L, train_vals))
    yc_train <- train_vals + 1L

    new_idx <- which(time_train > state$time_last)
    if (length(new_idx) > 0L) {
      state <- pf_intensity_update_many(
        state,
        yc_vec = yc_train[new_idx],
        time_vec = time_train[new_idx]
      )
    }

    return(state)
  }

  if (identical(run_info$model_family, "subjective")) {
    time_train <- as.integer(train$Time - t0_patient)
    train_vals <- pmin(100L, pmax(0L, as.integer(round(train$Score * 10))))

    new_idx <- which(time_train > state$time_last)
    if (length(new_idx) > 0L) {
      state <- pf_subjective_update_many(
        state,
        y_vec = train_vals[new_idx],
        time_vec = time_train[new_idx]
      )
    }

    return(state)
  }

  if (identical(run_info$model_family, "extent")) {
    day_train <- as.integer(train$Time - t0_patient)
    train_vals <- pmin(100L, pmax(0L, as.integer(round(train$Score))))

    new_idx <- which(day_train > state$time_last)
    if (length(new_idx) > 0L) {
      state <- pf_extent_update_many(
        state,
        y_vec = train_vals[new_idx],
        day_vec = day_train[new_idx]
      )
    }

    return(state)
  }

  stop("[ERROR] Unknown model family: ", run_info$model_family)
}
