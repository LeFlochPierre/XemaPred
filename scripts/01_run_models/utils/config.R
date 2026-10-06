# -----------------------------------------------------------------------
# Shared configuration utilities
#
# Purpose:
# - Read the YAML file passed to a validation script.
# - Validate dataset, score, model, horizon, chains, and cluster settings.
# - Classify each model as either an EczemaPred item-level model or a
#   reference/comparison model.
#
# Used by:
# - eczemapred_item_models/run_validation.R
# - reference_models/run_validation.R
# -----------------------------------------------------------------------

read_run_config <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) < 1) {
    stop("Please provide a YAML config file, e.g. Rscript run_validation.R configs/Derexyl/dryness_OrderedRW.yaml")
  }

  config_path <- args[[1]]
  cat(glue::glue("[INFO] Loading config: {config_path}\n"))

  cfg <- yaml::read_yaml(config_path)
  cfg$config_path <- config_path

  cfg
}

get_score_type <- function(score) {
  intensity_scores <- c("dryness", "redness", "swelling", "oozing", "thickening", "scratching")
  subjective_scores <- c("itching", "sleep")

  if (score %in% c("SCORAD", "oSCORAD")) {
    return("poscorad")
  }

  if (score == "extent") {
    return("extent")
  }

  if (score %in% intensity_scores) {
    return("intensity")
  }

  if (score %in% subjective_scores) {
    return("subjective")
  }

  return("unknown")
}

get_model_family <- function(mdl_name) {
  eczemapred_item_models <- c("BinMC", "OrderedRW", "BinRW")
  reference_models <- c("uniform", "historical", "RW", "AR1", "MixedAR1", "Smoothing", "MC")

  if (mdl_name %in% eczemapred_item_models) {
    return("eczemapred_item_models")
  }

  if (mdl_name %in% reference_models) {
    return("reference_models")
  }

  stop("Unknown model name: ", mdl_name)
}

validate_item_model_combination <- function(score_type, mdl_name) {
  valid <- list(
    extent = c("BinMC"),
    intensity = c("OrderedRW"),
    subjective = c("BinRW")
  )

  if (!score_type %in% names(valid)) {
    stop("EczemaPred item models are only defined for extent, intensity, and subjective scores.")
  }

  if (!mdl_name %in% valid[[score_type]]) {
    stop(
      "Invalid EczemaPred item-model combination: score_type = ", score_type,
      ", model = ", mdl_name,
      ". Expected one of: ", paste(valid[[score_type]], collapse = ", ")
    )
  }

  invisible(TRUE)
}

validate_reference_model_combination <- function(score_type, mdl_name) {
  valid <- list(
    poscorad = c("uniform", "historical", "RW", "AR1", "MixedAR1", "Smoothing"),
    extent = c("uniform", "historical", "RW"),
    intensity = c("uniform", "historical", "MC"),
    subjective = c("uniform", "historical", "RW")
  )

  if (!score_type %in% names(valid)) {
    stop("Unknown score type for reference model validation: ", score_type)
  }

  if (!mdl_name %in% valid[[score_type]]) {
    stop(
      "Invalid reference-model combination: score_type = ", score_type,
      ", model = ", mdl_name,
      ". Expected one of: ", paste(valid[[score_type]], collapse = ", ")
    )
  }

  invisible(TRUE)
}

is_one_logical <- function(x) {
  is.logical(x) && length(x) == 1 && !is.na(x)
}

is_one_whole_number <- function(x) {
  is.numeric(x) && length(x) == 1 && !is.na(x) && x == as.integer(x)
}

validate_run_config <- function(cfg, expected_family = NULL) {
  required_fields <- c("dataset", "score", "model", "t_horizon", "n_chains", "n_iter", "n_cluster", "run")
  missing_fields <- setdiff(required_fields, names(cfg))

  if (length(missing_fields) > 0) {
    stop("Missing required config fields: ", paste(missing_fields, collapse = ", "))
  }

  item_dict_all <- detail_POSCORAD()

  dataset <- match.arg(cfg$dataset, c("Derexyl", "PFDC", "Fake"))
  score <- match.arg(cfg$score, item_dict_all[["Name"]])
  mdl_name <- cfg$model

  score_type <- get_score_type(score)
  model_family <- get_model_family(mdl_name)

  if (!is.null(expected_family) && model_family != expected_family) {
    stop(
      "This runner is for ", expected_family, ", but config model '", mdl_name,
      "' belongs to ", model_family, "."
    )
  }

  if (model_family == "eczemapred_item_models") {
    validate_item_model_combination(score_type, mdl_name)
  }

  if (model_family == "reference_models") {
    validate_reference_model_combination(score_type, mdl_name)
  }

  t_horizon <- cfg$t_horizon
  n_chains <- cfg$n_chains
  n_it <- cfg$n_iter
  n_cluster <- cfg$n_cluster
  run <- cfg$run

  fit_only <- if (!is.null(cfg$fit_only)) cfg$fit_only else FALSE
  debug_one_patient <- if (!is.null(cfg$debug_one_patient)) cfg$debug_one_patient else FALSE

  seed <- as.integer(cfg$seed %||% 1744834965L)

  if (!is_one_logical(run)) stop("`run` must be TRUE or FALSE.")
  if (!is_one_logical(fit_only)) stop("`fit_only` must be TRUE or FALSE.")
  if (!is_one_logical(debug_one_patient)) stop("`debug_one_patient` must be TRUE or FALSE.")

  if (!is_one_whole_number(n_chains) || n_chains <= 0) stop("`n_chains` must be a positive whole number.")
  if (!is_one_whole_number(n_it) || n_it <= 0) stop("`n_iter` must be a positive whole number.")
  if (!is_one_whole_number(t_horizon) || t_horizon <= 0) stop("`t_horizon` must be a positive whole number.")
  if (!is_one_whole_number(n_cluster) || n_cluster <= 0) stop("`n_cluster` must be a positive whole number.")

  max_cluster <- floor((parallel::detectCores() - 2) / n_chains)
  if (n_cluster > max_cluster) {
    stop(
      "`n_cluster` is too high for the requested number of chains. ",
      "Requested n_cluster = ", n_cluster,
      ", maximum recommended = ", max_cluster, "."
    )
  }

  item_dict <- item_dict_all %>% dplyr::filter(.data$Name == score)

  item_lbl <- as.character(item_dict[["Label"]])
  max_score <- item_dict[["Maximum"]]

  is_continuous <- score %in% c("SCORAD", "oSCORAD")
  reso <- if (is_continuous) 1 else item_dict[["Resolution"]]
  M <- round(max_score / reso)

  param <- if (is_continuous) {
    c("lpd", "y_pred")
  } else {
    c("lpd", "cum_err", "y_pred")
  }

  run_info <- list(
    config_path = cfg$config_path,
    dataset = dataset,
    score = score,
    score_type = score_type,
    mdl_name = mdl_name,
    model_family = model_family,
    t_horizon = t_horizon,
    n_chains = n_chains,
    n_it = n_it,
    n_cluster = n_cluster,
    run = run,
    fit_only = fit_only,
    debug_one_patient = debug_one_patient,
    seed = seed,
    item_lbl = item_lbl,
    max_score = max_score,
    reso = reso,
    M = M,
    is_continuous = is_continuous,
    param = param
  )

  cat(glue::glue(
    "[INFO] Configuration:\n",
    "  Dataset: {dataset}\n",
    "  Score: {score}\n",
    "  Score type: {score_type}\n",
    "  Model: {mdl_name}\n",
    "  Model family: {model_family}\n",
    "  Horizon: {t_horizon}\n",
    "  Chains: {n_chains}\n",
    "  Iterations: {n_it}\n",
    "  Clusters: {n_cluster}\n",
    "  Fit only: {fit_only}\n",
    "  Debug one patient: {debug_one_patient}\n",
    "  Seed: {seed}\n\n"
  ))

  run_info
}
