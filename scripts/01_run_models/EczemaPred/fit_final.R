# -----------------------------------------------------------------------
# Final fit for EczemaPred item-level models
#
# Purpose:
# - Fit the model at the final forward-chaining iteration.
# - Extract posterior summaries that can be reused as informative priors.
# - Save the final Stan fit and the corresponding predictions/diagnostics.
#
# Notes:
# - This is not the normal validation path.
# - Normal validation uses fit_eczemapred_item_iteration().
# - This function is only used when fit_only: true in the YAML config.
# -----------------------------------------------------------------------

run_final_fit_for_power_priors <- function(data_obj, run_info, paths) {
  cat("[INFO] Running final fit for posterior summaries / power priors...\n")

  it <- max(data_obj$train_it)

  split <- split_iteration_data(data_obj, it)
  split <- apply_debug_one_patient(split, run_info)

  train <- split$Training
  test <- split$Testing

  if (nrow(train) == 0 || nrow(test) == 0) {
    stop("[ERROR] No valid training or test data for final fit.")
  }

  train_tmp <- rescale_for_model(train, run_info)
  test_tmp <- rescale_for_model(test, run_info)

  model <- EczemaModel(
    run_info$mdl_name,
    max_score = run_info$M,
    discrete = !run_info$is_continuous
  )

  start <- Sys.time()

  # final_seed <- as.integer(run_info$seed + 900000L + as.integer(it))
  # set.seed(final_seed)

  # fit <- EczemaFit(
  #   model,
  #   train_tmp,
  #   test_tmp,
  #   iter = run_info$n_it,
  #   chains = run_info$n_chains,
  #   refresh = 0,
  #   seed = final_seed
  # )

  fit <- EczemaFit(
    model,
    train_tmp,
    test_tmp,
    iter = run_info$n_it,
    chains = run_info$n_chains,
    refresh = 0
  )

  end <- Sys.time()

  check_stan_fit_valid(fit)
  check_stan_diagnostics(fit)

  power_prior <- extract_power_prior(fit, run_info$mdl_name)

  cat("[DEBUG] Power prior contents before saving:\n")
  print(power_prior)

  fit_final_path <- file.path(paths$fits_dir, "fit-final.rds")

  saveRDS(
    list(
      fit = fit,
      power_prior = power_prior,
      model = run_info$mdl_name,
      score = run_info$score,
      dataset = run_info$dataset,
      iteration = it
    ),
    file = fit_final_path
  )

  cat("[INFO] Saved fit + power_prior to: ", fit_final_path, "\n", sep = "")

  perf <- test %>%
    add_predictions(
      fit = fit,
      discrete = !run_info$is_continuous,
      include_samples = TRUE
    ) %>%
    restore_score_scale(run_info) %>%
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

  invisible(list(
    fit = fit,
    power_prior = power_prior,
    perf = perf,
    diag = diag_vi
  ))
}

# -----------------------------------------------------------------------
# Backward-compatible alias
#
# Purpose:
# - Keep old runner code working if it still calls run_final_fit().
# - New code should call run_final_fit_for_power_priors().
# -----------------------------------------------------------------------

run_final_fit <- function(data_obj, run_info, paths) {
  run_final_fit_for_power_priors(data_obj, run_info, paths)
}