# -----------------------------------------------------------------------
# Reference-model fitting functions
#
# Purpose:
# - Fit all reference/comparison models used in validation.
# - Keep reference-model logic separate from EczemaPred item-level models.
# - Provide one dispatcher, fit_reference_iteration(), called by the
#   reference-model runner.
#
# Models handled here:
# - uniform
# - historical
# - MC
# - RW
# - AR1
# - MixedAR1
# - Smoothing
#
# Notes:
# - RW, AR1, MixedAR1, and Smoothing use the util Stan iteration fitter:
#   utils/fit_stan_iteration.R
# - uniform and historical are non-Stan baselines.
# - MC is a special Stan model with one-indexed score transitions.
# -----------------------------------------------------------------------

# -----------------------------------------------------------------------
# AR(1) reference model
# -----------------------------------------------------------------------

fit_ar1_iteration <- function(it, data_obj, run_info, paths) {
  fit_stan_validation_iteration(it, data_obj, run_info, paths)
}

# -----------------------------------------------------------------------
# Mixed AR(1) reference model
# -----------------------------------------------------------------------

fit_mixed_ar1_iteration <- function(it, data_obj, run_info, paths) {
  fit_stan_validation_iteration(it, data_obj, run_info, paths)
}

# -----------------------------------------------------------------------
# Random-walk reference model
# -----------------------------------------------------------------------

fit_rw_iteration <- function(it, data_obj, run_info, paths) {
  fit_stan_validation_iteration(it, data_obj, run_info, paths)
}

# -----------------------------------------------------------------------
# Smoothing reference model
# -----------------------------------------------------------------------

fit_smoothing_iteration <- function(it, data_obj, run_info, paths) {
  fit_stan_validation_iteration(it, data_obj, run_info, paths)
}

# -----------------------------------------------------------------------
# Historical reference model
#
# Purpose:
# - Generate predictions from the empirical distribution of previous scores.
# - Used as a simple data-driven baseline.
# -----------------------------------------------------------------------

fit_historical_iteration <- function(it, data_obj, run_info, paths) {
  cat(glue::glue("[INFO] Starting iteration {it} for historical model...\n"))

  split <- split_iteration_data(data_obj, it)

  train_tmp <- rescale_for_model(split$Training, run_info)
  test_tmp <- rescale_for_model(split$Testing, run_info)

  start <- Sys.time()

  perf <- test_tmp %>%
    add_historical_pred(
      train = train_tmp,
      max_score = run_info$M,
      discrete = !run_info$is_continuous,
      add_uniform = TRUE,
      include_samples = run_info$is_continuous
    )

  end <- Sys.time()

  if (!run_info$is_continuous) {
    perf <- perf %>%
      dplyr::mutate(Score = .data$Score * run_info$reso)
  }

  perf <- perf %>%
    remove_last_observation_columns()

  diag_vi <- make_runtime_diag(run_info, it, start, end)

  saveRDS(
    perf,
    file = file.path(paths$iters_dir, sprintf("iter-%04d.rds", it))
  )

  saveRDS(
    diag_vi,
    file = file.path(paths$diag_dir, sprintf("diag-%04d.rds", it))
  )

  cat(glue::glue("[INFO] Iteration {it} complete.\n"))

  invisible(NULL)
}

# -----------------------------------------------------------------------
# Markov-chain reference model
#
# Purpose:
# - Fit the MC reference model for intensity signs.
# - Convert scores to the one-indexed format expected by the MC model.
# -----------------------------------------------------------------------

fit_mc_iteration <- function(it, data_obj, run_info, paths) {
  cat(glue::glue("[INFO] Starting iteration {it} for MC model...\n"))

  split <- split_iteration_data(data_obj, it)

  train <- split$Training
  test <- split$Testing

  train_MC <- train %>%
    dplyr::rename(y0 = Score) %>%
    dplyr::group_by(Patient) %>%
    dplyr::mutate(
      y0 = .data$y0 + 1,
      y1 = dplyr::lead(.data$y0),
      dt = dplyr::lead(.data$Time) - .data$Time
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(y0, y1, dt) %>%
    tidyr::drop_na()

  test_MC <- test %>%
    dplyr::rename(
      y0 = LastScore,
      y1 = Score,
      dt = Horizon
    ) %>%
    dplyr::mutate(
      y0 = .data$y0 + 1,
      y1 = .data$y1 + 1
    ) %>%
    dplyr::select(y0, y1, dt)

  model <- EczemaModel("MC", K = run_info$M + 1)

  start <- Sys.time()

  fit <- EczemaFit(
    model,
    train_MC,
    test_MC,
    iter = run_info$n_it,
    chains = run_info$n_chains,
    init = 0,
    pars = run_info$param,
    refresh = 0
  )

  end <- Sys.time()

  check_stan_fit_valid(fit)

  perf <- test %>%
    add_predictions(
      fit = fit,
      discrete = TRUE,
      include_samples = TRUE
    ) %>%
    dplyr::mutate(
      Samples = purrr::map(.data$Samples, ~ .x - 1)
    ) %>%
    remove_last_observation_columns()

  diag_vi <- make_runtime_diag(run_info, it, start, end)

  saveRDS(
    perf,
    file = file.path(paths$iters_dir, sprintf("iter-%04d.rds", it))
  )

  saveRDS(
    diag_vi,
    file = file.path(paths$diag_dir, sprintf("diag-%04d.rds", it))
  )

  cat(glue::glue("[INFO] Iteration {it} complete.\n"))

  rm(list = intersect(
    c("fit", "perf", "train", "test", "train_MC", "test_MC", "model"),
    ls()
  ))

  gc()

  invisible(NULL)
}

# -----------------------------------------------------------------------
# Uniform reference model
#
# Purpose:
# - Generate predictions from a uniform distribution over possible scores.
# - Used as a simple non-informative baseline.
# -----------------------------------------------------------------------

fit_uniform_iteration <- function(it, data_obj, run_info, paths) {
  cat(glue::glue("[INFO] Starting iteration {it} for uniform model...\n"))

  split <- split_iteration_data(data_obj, it)

  test_tmp <- rescale_for_model(split$Testing, run_info)

  start <- Sys.time()

  perf <- test_tmp %>%
    add_uniform_pred(
      max_score = run_info$M,
      discrete = !run_info$is_continuous,
      include_samples = run_info$is_continuous,
      n_samples = 2 * run_info$M
    )

  end <- Sys.time()

  if (!run_info$is_continuous) {
    perf <- perf %>%
      dplyr::mutate(Score = .data$Score * run_info$reso)
  }

  perf <- perf %>%
    remove_last_observation_columns()

  diag_vi <- make_runtime_diag(run_info, it, start, end)

  saveRDS(
    perf,
    file = file.path(paths$iters_dir, sprintf("iter-%04d.rds", it))
  )

  saveRDS(
    diag_vi,
    file = file.path(paths$diag_dir, sprintf("diag-%04d.rds", it))
  )

  cat(glue::glue("[INFO] Iteration {it} complete.\n"))

  invisible(NULL)
}

# -----------------------------------------------------------------------
# Reference-model dispatcher
#
# Purpose:
# - Route one forward-chaining iteration to the correct reference model
#   implementation based on run_info$mdl_name.
# -----------------------------------------------------------------------

fit_reference_iteration <- function(it, data_obj, run_info, paths) {
  cat(glue::glue(
    "[INFO] Starting reference iteration {it} for {run_info$mdl_name}...\n"
  ))

  switch(
    run_info$mdl_name,
    uniform = fit_uniform_iteration(it, data_obj, run_info, paths),
    historical = fit_historical_iteration(it, data_obj, run_info, paths),
    MC = fit_mc_iteration(it, data_obj, run_info, paths),
    RW = fit_rw_iteration(it, data_obj, run_info, paths),
    AR1 = fit_ar1_iteration(it, data_obj, run_info, paths),
    MixedAR1 = fit_mixed_ar1_iteration(it, data_obj, run_info, paths),
    Smoothing = fit_smoothing_iteration(it, data_obj, run_info, paths),
    stop("Unknown reference model: ", run_info$mdl_name)
  )
}