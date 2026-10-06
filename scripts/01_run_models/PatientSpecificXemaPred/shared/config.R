# -----------------------------------------------------------------------
# YAML configuration parsing for PatientSpecificXemaPred
# -----------------------------------------------------------------------

read_run_config <- function(args = commandArgs(trailingOnly = TRUE)) {

  if (length(args) < 1L) {
    stop(
      paste0(
        "Usage: Rscript ",
        "scripts/01_run_models/PatientSpecificXemaPred/run_validation.R ",
        "config.yaml"
      )
    )
  }

  config_path <- args[1]

  cat(
    glue::glue(
      "[INFO] Loading config: {config_path}\n"
    )
  )

  cfg <- yaml::read_yaml(config_path)
  cfg$config_path <- config_path

  cfg
}


validate_run_config <- function(
    cfg,
    expected_family = "patient_specific_xemapred"
) {

  # ---------------------------------------------------------------------
  # Basic run information
  # ---------------------------------------------------------------------

  dataset <- match.arg(
    cfg$dataset,
    c("Derexyl", "PFDC", "Fake")
  )

  score <- as.character(cfg$score)

  mapping <- get_score_mapping(score)


  # ---------------------------------------------------------------------
  # Save mode
  # ---------------------------------------------------------------------

  save_mode <- tolower(
    cfg$save_mode %||% "all"
  )

  if (!save_mode %in% c("all", "cache_only")) {
    stop(
      "[ERROR] save_mode must be 'all' or 'cache_only'."
    )
  }


  # ---------------------------------------------------------------------
  # Run suffix
  # ---------------------------------------------------------------------

  run_suffix <- cfg$run_suffix %||% ""

  if (
    is.null(run_suffix) ||
    run_suffix == "" ||
    toupper(run_suffix) == "NONE"
  ) {

    run_suffix <- ""

  } else if (!startsWith(run_suffix, "-")) {

    run_suffix <- paste0(
      "-",
      run_suffix
    )
  }


  # ---------------------------------------------------------------------
  # Prior configuration
  # ---------------------------------------------------------------------

  prior_mode <- tolower(
    cfg$prior_mode %||% "base"
  )

  prior_combination <- tolower(
    cfg$prior_combination %||%
      if (identical(prior_mode, "population")) {
        "power"
      } else {
        "base"
      }
  )

  if (
    identical(prior_mode, "population") &&
    !prior_combination %in% c("power", "mixture")
  ) {
    stop(
      "[ERROR] prior_combination must be 'power' or 'mixture'."
    )
  }

  if (identical(prior_mode, "base")) {
    prior_combination <- "base"
  }

  if (!prior_mode %in% c("base", "population")) {
    stop(
      "[ERROR] prior_mode must be 'base' or 'population'."
    )
  }


  prior_power <- as.numeric(
    cfg$prior_power %||%
      if (
        identical(
          prior_mode,
          "population"
        )
      ) {
        1.0
      } else {
        0.0
      }
  )


  if (identical(prior_mode, "population")) {

    if (identical(prior_combination, "power")) {

      if (
        !is.finite(prior_power) ||
        prior_power <= 0 ||
        prior_power > 1
      ) {
        stop(
          paste0(
            "[ERROR] For prior_combination = 'power', ",
            "prior_power must satisfy 0 < prior_power <= 1."
          )
        )
      }

    } else if (identical(prior_combination, "mixture")) {

      if (
        !is.finite(prior_power) ||
        prior_power < 0 ||
        prior_power > 1
      ) {
        stop(
          paste0(
            "[ERROR] For prior_combination = 'mixture', ",
            "prior_power must satisfy 0 <= prior_power <= 1."
          )
        )
      }
    }
  }


  pmmh_prior_mode <- tolower(
    cfg$pmmh_prior_mode %||% "matched"
  )

  if (!pmmh_prior_mode %in% c("matched", "base")) {
    stop(
      "[ERROR] pmmh_prior_mode must be 'matched' or 'base'."
    )
  }


  allow_same_dataset_prior <- isTRUE(
    cfg$allow_same_dataset_prior %||% FALSE
  )


  # ---------------------------------------------------------------------
  # Optional patient-fold selection
  #
  # patient_fold:
  #   fold of patients that will actually be fitted / forecast.
  #
  # prior_source_fold:
  #   fold whose population posterior is used as the informative prior.
  #
  # Example:
  #
  #   patient_fold:      2
  #   prior_source_fold: 1
  #
  # means:
  #
  #   population Fold 1 -> patient-specific Fold 2
  # ---------------------------------------------------------------------

  n_patient_folds <- as.integer(
    cfg$n_patient_folds %||% 1L
  )

  patient_fold <- as.integer(
    cfg$patient_fold %||% 1L
  )

  prior_source_fold <- if (
    !is.null(
      cfg$prior_source_fold
    )
  ) {

    as.integer(
      cfg$prior_source_fold
    )

  } else {

    NA_integer_
  }


  if (
    is.na(n_patient_folds) ||
    n_patient_folds < 1L
  ) {

    stop(
      "[ERROR] n_patient_folds must be an integer >= 1."
    )
  }


  if (
    is.na(patient_fold) ||
    patient_fold < 1L ||
    patient_fold > n_patient_folds
  ) {

    stop(
      paste0(
        "[ERROR] patient_fold must be between 1 and ",
        "n_patient_folds."
      )
    )
  }


  # ---------------------------------------------------------------------
  # Population-prior source dataset
  # ---------------------------------------------------------------------

  prior_root <-
    cfg$prior_root %||%
    file.path(
      "priors",
      "population_xemapred"
    )


  prior_source_dataset <-
    cfg$prior_source_dataset %||%
    if (
      identical(
        dataset,
        "Derexyl"
      )
    ) {

      "PFDC"

    } else if (
      identical(
        dataset,
        "PFDC"
      )
    ) {

      "Derexyl"

    } else {

      NULL
    }


  # ---------------------------------------------------------------------
  # Validate same-dataset cross-fold experiment
  # ---------------------------------------------------------------------

  same_dataset_prior <-
    identical(
      prior_mode,
      "population"
    ) &&
    !is.null(
      prior_source_dataset
    ) &&
    identical(
      prior_source_dataset,
      dataset
    )


  if (same_dataset_prior) {

    # ---------------------------------------------------------------
    # Same-dataset prior must be explicitly enabled.
    # ---------------------------------------------------------------

    if (!allow_same_dataset_prior) {

      stop(
        paste0(
          "[ERROR] Same-dataset population prior requested, ",
          "but allow_same_dataset_prior is FALSE. ",
          "Set `allow_same_dataset_prior: true` for an explicit ",
          "held-out cross-fold experiment."
        )
      )
    }


    # ---------------------------------------------------------------
    # For the experiment considered here, same-dataset informative
    # priors must be generated from a held-out patient fold.
    # ---------------------------------------------------------------

    if (n_patient_folds < 2L) {

      stop(
        paste0(
          "[ERROR] Same-dataset population prior requires ",
          "n_patient_folds >= 2 to avoid patient leakage."
        )
      )
    }


    if (is.na(prior_source_fold)) {

      stop(
        paste0(
          "[ERROR] Same-dataset cross-fold population prior requires ",
          "`prior_source_fold` in the YAML."
        )
      )
    }


    if (
      prior_source_fold < 1L ||
      prior_source_fold > n_patient_folds
    ) {

      stop(
        paste0(
          "[ERROR] prior_source_fold must be between 1 and ",
          "n_patient_folds."
        )
      )
    }


    # ---------------------------------------------------------------
    # Critical leakage check.
    # ---------------------------------------------------------------

    if (
      prior_source_fold ==
      patient_fold
    ) {

      stop(
        paste0(
          "[ERROR] Leakage detected: prior_source_fold cannot equal ",
          "patient_fold in the held-out population-prior experiment."
        )
      )
    }


    # ---------------------------------------------------------------
    # Require an explicit fold-specific prior file.
    #
    # Do NOT silently use:
    #
    #   priors/population_xemapred/PFDC/dryness.rds
    #
    # because that may be a posterior learned from the full cohort.
    # ---------------------------------------------------------------

    if (
      is.null(
        cfg$population_prior_file
      ) ||
      !nzchar(
        as.character(
          cfg$population_prior_file
        )
      )
    ) {

      stop(
        paste0(
          "[ERROR] Same-dataset cross-fold prior requires an explicit ",
          "`population_prior_file` pointing to the correct ",
          "fold-specific population posterior."
        )
      )
    }
  }


  # ---------------------------------------------------------------------
  # Population-prior file
  # ---------------------------------------------------------------------

  population_prior_file <-
    cfg$population_prior_file %||%
    if (
      identical(
        prior_mode,
        "population"
      ) &&
      !is.null(
        prior_source_dataset
      )
    ) {

      file.path(
        prior_root,
        prior_source_dataset,
        paste0(
          score,
          ".rds"
        )
      )

    } else {

      NULL
    }


  # ---------------------------------------------------------------------
  # Validate population prior
  # ---------------------------------------------------------------------

  if (
    identical(
      prior_mode,
      "population"
    )
  ) {

    if (
      is.null(
        prior_source_dataset
      )
    ) {

      stop(
        "[ERROR] Could not determine prior_source_dataset."
      )
    }


    if (
      is.null(
        population_prior_file
      ) ||
      !file.exists(
        population_prior_file
      )
    ) {

      stop(
        paste0(
          "[ERROR] Population prior file not found: ",
          population_prior_file
        )
      )
    }
  }


  # ---------------------------------------------------------------------
  # Build validated run information
  # ---------------------------------------------------------------------

  run_info <- list(

    family =
      expected_family,

    config_path =
      cfg$config_path,

    cfg =
      cfg,

    dataset =
      dataset,

    score =
      score,

    model =
      cfg$model %||%
      "SMC",

    t_horizon =
      as.integer(
        cfg$t_horizon %||% 1L
      ),

    n_particles =
      as.integer(
        cfg$n_particles %||% 2000L
      ),

    ess_frac =
      as.numeric(
        cfg$ess_frac %||% 0.5
      ),

    n_cluster =
      as.integer(
        cfg$n_cluster %||% 1L
      ),

    run =
      if (
        !is.null(
          cfg$run
        )
      ) {
        isTRUE(
          cfg$run
        )
      } else {
        TRUE
      },

    pmmh_window_n =
      if (
        !is.null(
          cfg$pmmh_window_n
        )
      ) {
        as.integer(
          cfg$pmmh_window_n
        )
      } else {
        20L
      },

    pmmh_nx =
      as.integer(
        cfg$pmmh_nx %||% 100L
      ),

    save_mode =
      save_mode,

    seed =
      as.integer(
        cfg$seed %||% 1744834965L
      ),

    output_dir =
      cfg$output_dir %||%
      NULL,

    result_root =
      cfg$result_root %||%
      "results",

    output_label =
      cfg$output_label %||%
      "PatientSpecificXemaPred",

    run_suffix =
      run_suffix,

    item_label =
      mapping$item_label,

    max_score =
      mapping$max_score,

    reso =
      mapping$reso,

    M_max =
      mapping$M_max,

    K =
      mapping$K,

    model_family =
      mapping$model_family,

    prior_mode =
      prior_mode,

    prior_combination =
      prior_combination,

    prior_power =
      prior_power,

    prior_root =
      prior_root,

    prior_source_dataset =
      prior_source_dataset,

    population_prior_file =
      population_prior_file,

    allow_same_dataset_prior =
      allow_same_dataset_prior,

    pmmh_prior_mode =
      pmmh_prior_mode,

    n_patient_folds =
      n_patient_folds,

    patient_fold =
      patient_fold,

    prior_source_fold =
      prior_source_fold
  )


  # ---------------------------------------------------------------------
  # Final sanity checks
  # ---------------------------------------------------------------------

  if (
    run_info$t_horizon < 1L
  ) {

    stop(
      "[ERROR] t_horizon must be >= 1."
    )
  }


  if (
    run_info$n_particles < 1L
  ) {

    stop(
      "[ERROR] n_particles must be >= 1."
    )
  }


  if (
    run_info$n_cluster < 1L
  ) {

    stop(
      "[ERROR] n_cluster must be >= 1."
    )
  }


  run_info
}


print_run_config <- function(run_info) {

  cat(
    glue::glue(
      "[INFO] Configuration:\n",

      "  Dataset: {run_info$dataset}\n",
      "  Score: {run_info$score}\n",
      "  Model: {run_info$model}\n",
      "  Family: {run_info$model_family}\n",

      "  Horizon: {run_info$t_horizon}\n",
      "  Particles: {run_info$n_particles}\n",
      "  ESS fraction: {run_info$ess_frac}\n",
      "  Workers: {run_info$n_cluster}\n",
      "  Run: {run_info$run}\n",

      "  PMMH nx: {run_info$pmmh_nx}\n",
      "  PMMH window: ",
      "{if (is.infinite(run_info$pmmh_window_n)) ",
      "'Inf (all data)' else run_info$pmmh_window_n}\n",

      "  Prior mode: {run_info$prior_mode}\n",
      "  Prior combination: {run_info$prior_combination}\n",
      "  Prior power/weight: {run_info$prior_power}\n",
      "  Prior root: {run_info$prior_root}\n",

      "  Prior source dataset: ",
      "{run_info$prior_source_dataset}\n",

      "  Allow same dataset prior: ",
      "{run_info$allow_same_dataset_prior}\n",

      "  Population prior file: ",
      "{run_info$population_prior_file}\n",

      "  PMMH prior mode: ",
      "{run_info$pmmh_prior_mode}\n",

      "  Save mode: {run_info$save_mode}\n",
      "  Output dir: {run_info$output_dir}\n",
      "  Result root: {run_info$result_root}\n",
      "  Output label: {run_info$output_label}\n",
      "  Run suffix: {run_info$run_suffix}\n",

      "  Patient folds: {run_info$n_patient_folds}\n",
      "  Target patient fold: {run_info$patient_fold}\n",
      "  Prior source fold: {run_info$prior_source_fold}\n",

      .null = "NULL",
      .trim = FALSE
    )
  )
}