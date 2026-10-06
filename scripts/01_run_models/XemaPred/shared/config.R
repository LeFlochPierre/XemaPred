# -----------------------------------------------------------------------
# XemaPred config utilities
#
# Purpose:
# - Read and validate YAML configs for population-level XemaPred.
# - Convert the raw YAML list into a normalised run_info object.
# - Keep config parsing out of the main runner.
# -----------------------------------------------------------------------

read_xemapred_config <- function(config_path = NULL) {
  if (is.null(config_path)) {
    args <- commandArgs(trailingOnly = TRUE)

    if (length(args) < 1) {
      stop("Usage: Rscript run_validation.R config.yaml")
    }

    config_path <- args[1]
  }

  cat(glue::glue("[INFO] Loading config: {config_path}\n"))

  cfg <- yaml::read_yaml(config_path)
  attr(cfg, "config_path") <- config_path

  cfg
}

parse_inf_integer <- function(x, default = Inf) {
  if (is.null(x)) {
    return(default)
  }

  if (is.character(x) && tolower(x) %in% c("inf", "infinity")) {
    return(Inf)
  }

  if (is.infinite(x)) {
    return(Inf)
  }

  as.integer(x)
}

validate_xemapred_config <- function(cfg) {
  config_path <- attr(cfg, "config_path") %||% NA_character_

  dataset <- cfg$dataset %||% stop("[ERROR] Missing config field: dataset")
  score <- cfg$score %||% stop("[ERROR] Missing config field: score")
  mdl_name <- cfg$model %||% get_expected_xemapred_item_model(score)

  dataset <- match.arg(dataset, c("Derexyl", "PFDC", "Fake"))

  score_info <- get_xemapred_score_info(score)

  expected_model <- score_info$expected_model

  if (!identical(mdl_name, expected_model)) {
    stop(
      "[ERROR] Model mismatch for score '",
      score,
      "'. Expected '",
      expected_model,
      "' but got '",
      mdl_name,
      "'."
    )
  }

  t_horizon <- as.integer(cfg$t_horizon %||% 1L)
  n_particles <- as.integer(cfg$n_particles %||% 2000L)
  ess_frac <- as.numeric(cfg$ess_frac %||% 0.5)
  lw_a <- as.numeric(cfg$lw_a %||% 0.98)
  n_cluster <- as.integer(cfg$n_cluster %||% 1L)

  run <- if (!is.null(cfg$run)) {
    isTRUE(cfg$run)
  } else {
    TRUE
  }

  n_theta <- as.integer(cfg$n_theta %||% 200L)
  n_x <- as.integer(cfg$n_x %||% 200L)

  ess_theta <- as.numeric(cfg$ess_theta %||% 0.5)
  ess_x <- as.numeric(cfg$ess_x %||% ess_frac)

  a_theta <- as.numeric(cfg$a_theta %||% lw_a)
  a_x <- as.numeric(cfg$a_x %||% lw_a)

  pmmh_window_n <- parse_inf_integer(cfg$pmmh_window_n, default = Inf)

  seed <- as.integer(cfg$seed %||% 1744834965L)

  run_suffix <- cfg$run_suffix %||% ""

  # -----------------------------------------------------------------------
  # Optional patient-fold selection
  #
  # Default:
  #   n_patient_folds = 1
  #   patient_fold = 1
  #
  # means use the complete cohort, preserving the original behaviour.
  # -----------------------------------------------------------------------

  n_patient_folds <- as.integer(
    cfg$n_patient_folds %||% 1L
  )

  patient_fold <- as.integer(
    cfg$patient_fold %||% 1L
  )

  if (n_patient_folds < 1L) {
    stop("[ERROR] n_patient_folds must be >= 1.")
  }

  if (patient_fold < 1L || patient_fold > n_patient_folds) {
    stop(
      "[ERROR] patient_fold must be between 1 and n_patient_folds. ",
      "Got patient_fold = ",
      patient_fold,
      " and n_patient_folds = ",
      n_patient_folds,
      "."
    )
  }

  stopifnot(
    is.numeric(t_horizon), t_horizon > 0,
    is.numeric(n_particles), n_particles > 0,
    is.numeric(ess_frac), ess_frac > 0, ess_frac <= 1,
    is.numeric(lw_a), lw_a > 0, lw_a < 1,
    is.numeric(n_cluster), n_cluster >= 1,
    is.numeric(n_theta), n_theta > 0,
    is.numeric(n_x), n_x > 0,
    is.numeric(ess_theta), ess_theta > 0, ess_theta <= 1,
    is.numeric(ess_x), ess_x > 0, ess_x <= 1,
    is.numeric(a_theta), a_theta > 0, a_theta < 1,
    is.numeric(a_x), a_x > 0, a_x < 1
  )

  run_info <- list(
    config_path = config_path,
    dataset = dataset,
    score = score,
    mdl_name = mdl_name,
    model_family = score_info$model_family,
    item_label = score_info$item_label,
    max_score = score_info$max_score,
    reso = score_info$reso,
    M_max = score_info$M_max,
    K = score_info$K,
    t_horizon = t_horizon,
    n_particles = n_particles,
    ess_frac = ess_frac,
    lw_a = lw_a,
    n_cluster = n_cluster,
    run = run,
    n_theta = n_theta,
    n_x = n_x,
    ess_theta = ess_theta,
    ess_x = ess_x,
    a_theta = a_theta,
    a_x = a_x,
    pmmh_window_n = pmmh_window_n,
    n_patient_folds = n_patient_folds,
    patient_fold = patient_fold,
    seed = seed,
    run_suffix = run_suffix
  )

  print_xemapred_run_info(run_info)

  run_info
}

print_xemapred_run_info <- function(run_info) {
  cat(
    glue::glue(
      "[INFO] XemaPred configuration:\n",
      "  Dataset: {run_info$dataset}\n",
      "  Score: {run_info$score}\n",
      "  Model: {run_info$mdl_name}\n",
      "  Family: {run_info$model_family}\n",
      "  Horizon: {run_info$t_horizon}\n",
      "  Particles: {run_info$n_particles}\n",
      "  ESS fraction: {run_info$ess_frac}\n",
      "  Liu-West a: {run_info$lw_a}\n",
      "  Workers: {run_info$n_cluster}\n",
      "  Run: {run_info$run}\n",
      "  n_theta: {run_info$n_theta}\n",
      "  n_x: {run_info$n_x}\n",
      "  ess_theta: {run_info$ess_theta}\n",
      "  ess_x: {run_info$ess_x}\n",
      "  a_theta: {run_info$a_theta}\n",
      "  a_x: {run_info$a_x}\n",
      "  Patient folds: {run_info$n_patient_folds}\n",
      "  Patient fold selected: {run_info$patient_fold}\n",
      "  PMMH window: {if (is.infinite(run_info$pmmh_window_n)) 'Inf' else run_info$pmmh_window_n}\n",
      "  Seed: {run_info$seed}\n",
      "  Run suffix: {run_info$run_suffix}\n"
    )
  )

  invisible(run_info)
}