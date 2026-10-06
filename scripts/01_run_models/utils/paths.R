# -----------------------------------------------------------------------
# Shared output-path utilities
#
# Purpose:
# - Build the canonical output directory for each dataset, score, model,
#   and forecasting horizon.
# - Save the exact YAML config and run metadata next to the model outputs.
#
# Output layout:
# results/<dataset>/<score>/<model>/H<horizon>/
#
# Example:
# results/Derexyl/dryness/OrderedRW/H4/
# results/PFDC/SCORAD/MixedAR1/H4/
# -----------------------------------------------------------------------

make_output_paths <- function(run_info) {
  base_dir <- here::here(
    "results",
    run_info$dataset,
    run_info$score,
    run_info$mdl_name,
    paste0("H", run_info$t_horizon)
  )

  paths <- list(
    base_dir = base_dir,
    iters_dir = file.path(base_dir, "iters"),
    diag_dir = file.path(base_dir, "diag"),
    fits_dir = file.path(base_dir, "fits"),
    posterior_dir = file.path(base_dir, "posterior"),
    final_dir = file.path(base_dir, "final"),
    config_out = file.path(base_dir, "cfg.yaml"),
    meta_out = file.path(base_dir, "meta.json"),
    predictions_out = file.path(base_dir, "final", "predictions.rds"),
    diagnostics_out = file.path(base_dir, "final", "diagnostics.rds")
  )

  invisible(lapply(
    c(
      paths$base_dir,
      paths$iters_dir,
      paths$diag_dir,
      paths$fits_dir,
      paths$posterior_dir,
      paths$final_dir
    ),
    dir.create,
    recursive = TRUE,
    showWarnings = FALSE
  ))

  cat(
    "[INFO] Output directories:\n",
    "  base:      ", paths$base_dir, "\n",
    "  iters:     ", paths$iters_dir, "\n",
    "  diag:      ", paths$diag_dir, "\n",
    "  fits:      ", paths$fits_dir, "\n",
    "  posterior: ", paths$posterior_dir, "\n",
    "  final:     ", paths$final_dir, "\n",
    sep = ""
  )

  paths
}

save_run_metadata <- function(cfg, run_info, paths) {
  yaml::write_yaml(cfg, paths$config_out)

  meta <- list(
    dataset = run_info$dataset,
    score = run_info$score,
    score_type = run_info$score_type,
    model = run_info$mdl_name,
    model_family = run_info$model_family,
    t_horizon = run_info$t_horizon,
    n_chains = run_info$n_chains,
    n_iter = run_info$n_it,
    n_cluster = run_info$n_cluster,
    run = run_info$run,
    fit_only = run_info$fit_only,
    debug_one_patient = run_info$debug_one_patient,
    seed = run_info$seed,
    config_path = run_info$config_path
  )

  jsonlite::write_json(
    meta,
    paths$meta_out,
    auto_unbox = TRUE,
    pretty = TRUE
  )

  invisible(meta)
}