# -----------------------------------------------------------------------
# Output paths for PatientSpecificXemaPred
# 
# Output layout:
# results/<dataset>/<score>/<model>-PatientSpecificXemaPred[-run_suffix]/H<horizon>/
# -----------------------------------------------------------------------

make_output_paths <- function(run_info) {
  if (!is.null(run_info$output_dir) && nzchar(run_info$output_dir)) {
    base_dir <- run_info$output_dir
  } else {
    model_output_name <- paste0(
      run_info$model,
      "-",
      run_info$output_label,
      run_info$run_suffix
    )

    base_dir <- here::here(
      run_info$result_root,
      run_info$dataset,
      run_info$score,
      model_output_name,
      paste0("H", run_info$t_horizon)
    )
  }

  final_dir <- file.path(base_dir, "final")
  fits_dir <- file.path(base_dir, "fits")
  iters_dir <- if (identical(run_info$save_mode, "all")) file.path(base_dir, "iters") else NULL
  diag_dir <- if (identical(run_info$save_mode, "all")) file.path(base_dir, "diag") else NULL

  dirs <- c(base_dir, final_dir, fits_dir, iters_dir, diag_dir)
  dirs <- dirs[!vapply(dirs, is.null, logical(1))]
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

  if (!dir.exists(base_dir)) stop("[ERROR] Could not create base_dir: ", base_dir)
  if (!dir.exists(final_dir)) stop("[ERROR] Could not create final_dir: ", final_dir)
  if (!dir.exists(fits_dir)) stop("[ERROR] Could not create fits_dir: ", fits_dir)

  list(
    base_dir = base_dir,
    final_dir = final_dir,
    fits_dir = fits_dir,
    iters_dir = iters_dir,
    diag_dir = diag_dir
  )
}

save_run_metadata <- function(cfg, run_info, paths) {
  yaml::write_yaml(cfg, file.path(paths$base_dir, "cfg.yaml"))

  if (!exists(".Random.seed", envir = .GlobalEnv)) runif(1)

  meta <- list(
    family = run_info$family,
    score = run_info$score,
    dataset = run_info$dataset,
    model = run_info$model,
    model_family = run_info$model_family,
    t_horizon = run_info$t_horizon,
    n_particles = run_info$n_particles,
    ess_frac = run_info$ess_frac,
    n_cluster = run_info$n_cluster,
    run = run_info$run,
    seed = run_info$seed,
    rng_kind = RNGkind()[1],
    M_max = run_info$M_max,
    K = run_info$K,
    reso = run_info$reso,
    max_score = run_info$max_score,
    save_mode = run_info$save_mode,
    pmmh_nx = run_info$pmmh_nx,
    pmmh_window_n = run_info$pmmh_window_n,
    prior_mode = run_info$prior_mode,
    prior_combination = run_info$prior_combination,
    prior_power = run_info$prior_power,
    prior_source_dataset = run_info$prior_source_dataset,
    population_prior_file = run_info$population_prior_file,
    population_prior_id = run_info$population_prior$prior_id %||% NULL,
    population_prior_iteration = run_info$population_prior$source_iteration %||% NULL,
    population_prior_training_day = run_info$population_prior$source_training_day %||% NULL,
    pmmh_prior_mode = run_info$pmmh_prior_mode,
    n_patient_folds = run_info$n_patient_folds,
    patient_fold = run_info$patient_fold,
    prior_source_fold = run_info$prior_source_fold
  )

  jsonlite::write_json(meta, file.path(paths$base_dir, "meta.json"), auto_unbox = TRUE, pretty = TRUE)
  cat("[INFO] Output base:\n", "  ", paths$base_dir, "\n", sep = "")
}

patient_iters_dir <- function(paths, pid) {
  if (is.null(paths$iters_dir)) return(NULL)
  out <- file.path(paths$iters_dir, glue::glue("patient_{pid}"))
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  out
}

patient_diag_dir <- function(paths, pid) {
  if (is.null(paths$diag_dir)) return(NULL)
  out <- file.path(paths$diag_dir, glue::glue("patient_{pid}"))
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  out
}

patient_state_file <- function(paths, pid) {
  out <- file.path(paths$fits_dir, glue::glue("patient_{pid}"), "pf_state.rds")
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  out
}
