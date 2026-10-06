# -----------------------------------------------------------------------
# Shared Stan iteration fitter
#
# Purpose:
# - Fit one Stan-based EczemaPred model for one forward-chaining iteration.
# - Save iteration-level predictions and runtime diagnostics.
# - Reused by both EczemaPred item models and Stan-based reference models.
#
# Used by:
# - eczemapred_item_models/fit_iteration.R
# - reference_models/fit_reference_models.R
# -----------------------------------------------------------------------

get_key_posterior_parameters <- function(model) {
  switch(
    model,
    "BinMC" = c("mu_logit_p10", "sigma_logit_p10", "sigma"),
    "OrderedRW" = c("sigma_lat", "sigma_meas", "mu_y0", "sigma_y0","ct"),
    "BinRW" = c("sigma", "mu_logit_y0", "sigma_logit_y0"),
    stop("Unknown EczemaPred model: ", model)
  )
}

fit_stan_validation_iteration <- function(it, data_obj, run_info, paths) {
  cat(glue::glue(
    "[INFO] Starting iteration {it} for {run_info$mdl_name}...\n"
  ))

  split <- split_iteration_data(data_obj, it)

  train <- split$Training
  test <- split$Testing

  train_tmp <- rescale_for_model(train, run_info)
  test_tmp <- rescale_for_model(test, run_info)

  model <- EczemaModel(
    run_info$mdl_name,
    max_score = run_info$M,
    discrete = !run_info$is_continuous
  )

  # Posterior saving only for EczemaPred item-level models.
  save_posterior <- identical(
    run_info$model_family,
    "eczemapred_item_models"
  )

  posterior_pars <- if (save_posterior) {
    get_key_posterior_parameters(run_info$mdl_name)
  } else {
    character(0)
  }

  pars_to_save <- unique(c(run_info$param, posterior_pars))

  iter_seed <- as.integer(
    run_info$seed + 100000L + as.integer(it)
  )

  start <- Sys.time()

  fit <- EczemaFit(
    model,
    train_tmp,
    test_tmp,
    iter = run_info$n_it,
    chains = run_info$n_chains,
    pars = pars_to_save,
    refresh = 0,
    seed = iter_seed
  )

  end <- Sys.time()

  check_stan_fit_valid(fit)

  # Save compact HMC posterior for EczemaPred.
  if (save_posterior) {
    posterior_draws <- rstan::extract(
      fit,
      pars = posterior_pars,
      permuted = TRUE
    )

    posterior_summary <- as.data.frame(
      summary(fit, pars = posterior_pars)$summary
    )

    posterior_summary$parameter <- rownames(posterior_summary)
    rownames(posterior_summary) <- NULL

    posterior_result <- list(
      dataset = run_info$dataset,
      score = run_info$score,
      model = run_info$mdl_name,
      iteration = it,
      seed = iter_seed,
      max_training_day = max(train$Time, na.rm = TRUE),
      n_training_observations = nrow(train),
      n_training_patients = dplyr::n_distinct(train$Patient),
      parameters = posterior_pars,
      draws = posterior_draws,
      hmc_summary = posterior_summary
    )

    saveRDS(
      posterior_result,
      file.path(
        paths$posterior_dir,
        sprintf("posterior-%04d.rds", it)
      )
    )
  }

  perf <- test %>%
    add_predictions(
      fit = fit,
      discrete = !run_info$is_continuous,
      include_samples = TRUE
    ) %>%
    restore_score_scale(run_info) %>%
    remove_last_observation_columns()

  diag_vi <- make_runtime_diag(
    run_info,
    it,
    start,
    end
  )

  saveRDS(
    perf,
    file.path(paths$iters_dir, sprintf("iter-%04d.rds", it))
  )

  saveRDS(
    diag_vi,
    file.path(paths$diag_dir, sprintf("diag-%04d.rds", it))
  )

  cat(glue::glue("[INFO] Iteration {it} complete.\n"))

  rm(list = intersect(
    c(
      "fit", "perf", "train", "test", "train_tmp", "test_tmp", "model",
      "posterior_draws", "posterior_summary", "posterior_result"
    ),
    ls()
  ))

  gc()
  invisible(NULL)
}