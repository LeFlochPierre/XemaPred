# -----------------------------------------------------------------------
# XemaPred output-path utilities
#
# Purpose:
# - Build canonical output directories for XemaPred runs.
# - Save the exact YAML config and run metadata for reproducibility.
#
# Output layout:
# results_XemaPred/<dataset>/<score>/<model>-XemaPred[-run_suffix]/H<horizon>/
# -----------------------------------------------------------------------

format_xemapred_run_suffix <- function(run_suffix) {
  if (is.null(run_suffix) || run_suffix == "" || toupper(run_suffix) == "NONE") {
    return("")
  }

  paste0("-", run_suffix)
}

make_xemapred_paths <- function(run_info) {
  # XemaPred runner is always the population-level SMC2 model.
  # Output layout:
  # results_XemaPred/<dataset>/<score>/<model>-XemaPred[-run_suffix]/H<horizon>/
  run_suffix <- format_xemapred_run_suffix(run_info$run_suffix)

  model_dir <- paste0(
    run_info$mdl_name,
    "-XemaPred",
    run_suffix
  )

  base_dir <- here::here(
    "results",
    run_info$dataset,
    run_info$score,
    model_dir,
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

    # Cached population-level SMC2 state used for sequential resuming.
    cohort_state_file = file.path(base_dir, "fits", "smc2", "pf_state.rds"),

    # EczemaPred-style final outputs.
    predictions_out = file.path(base_dir, "final", "predictions.rds"),
    diagnostics_out = file.path(base_dir, "final", "diagnostics.rds"),

    # Backward-compatible aliases for older plotting scripts.
    predictions_all = file.path(base_dir, "final", "predictions_ALL_PATIENTS.rds"),
    diagnostics_all = file.path(base_dir, "final", "diagnostics_ALL_PATIENTS.rds")
  )

  invisible(lapply(
    c(
      paths$base_dir,
      paths$iters_dir,
      paths$diag_dir,
      paths$fits_dir,
      paths$posterior_dir,
      paths$final_dir,
      dirname(paths$cohort_state_file)
    ),
    dir.create,
    recursive = TRUE,
    showWarnings = FALSE
  ))

  cat(
    "[INFO] XemaPred output directories:\n",
    "  base:  ", paths$base_dir, "\n",
    "  iters: ", paths$iters_dir, "\n",
    "  diag:  ", paths$diag_dir, "\n",
    "  fits:  ", paths$fits_dir, "\n",
    "  final: ", paths$final_dir, "\n",
    sep = ""
  )

  paths
}

save_xemapred_metadata <- function(cfg, run_info, paths) {
  yaml::write_yaml(cfg, paths$config_out)

  meta <- list(
    dataset = run_info$dataset,
    score = run_info$score,
    model = run_info$mdl_name,
    model_family = run_info$model_family,
    t_horizon = run_info$t_horizon,
    n_particles = run_info$n_particles,
    ess_frac = run_info$ess_frac,
    lw_a = run_info$lw_a,
    n_cluster = run_info$n_cluster,
    run = run_info$run,
    n_theta = run_info$n_theta,
    n_x = run_info$n_x,
    ess_theta = run_info$ess_theta,
    ess_x = run_info$ess_x,
    a_theta = run_info$a_theta,
    a_x = run_info$a_x,
    pmmh_window_n = if (is.infinite(run_info$pmmh_window_n)) {
      "Inf"
    } else {
      run_info$pmmh_window_n
    },
    seed = run_info$seed,
    rng_kind = RNGkind()[1],
    M_max = run_info$M_max,
    K = run_info$K,
    reso = run_info$reso,
    max_score = run_info$max_score,
    run_suffix = run_info$run_suffix,
    config_path = run_info$config_path,
    output_layout = "iteration_level"
  )

  jsonlite::write_json(
    meta,
    paths$meta_out,
    auto_unbox = TRUE,
    pretty = TRUE
  )

  invisible(meta)
}
