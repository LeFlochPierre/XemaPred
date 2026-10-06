#!/usr/bin/env Rscript

# ======================================================================
# process_patient_specific_xemapred.R
#
# Purpose:
# - Process patient-specific XemaPred outputs.
# - No population-level XemaPred.
# - No EczemaPred.
# - No reference models.
# - Read predictions_ALL_PATIENTS.rds and diagnostics_ALL_PATIENTS.rds.
# - Recompute item-level and PO-SCORAD performance.
#
# Expected raw inputs:
#   results/<Dataset>/<Item>/<Model>-PatientSpecificXemaPred[-run_suffix]/H<horizon>/final/predictions_ALL_PATIENTS.rds
#   results/<Dataset>/<Item>/<Model>-PatientSpecificXemaPred[-run_suffix]/H1/final/diagnostics_ALL_PATIENTS.rds
#
# Main patient-specific XemaPred item models:
#   extent                     -> BinMC-PatientSpecificXemaPred
#   itching, sleep             -> BinRW-PatientSpecificXemaPred
#   intensity signs            -> OrderedRW-PatientSpecificXemaPred
#
# Outputs:
#   results/<Dataset>/ALL_MODELS/ITEMS/<Metric>/perf_PatientSpecificXemaPred_<Dataset>_<Metric>.RData
#   results/<Dataset>/ALL_MODELS/POSCORAD/<Metric>/perf_PatientSpecificXemaPred_<Dataset>_<Metric>.RData
#   results/<Dataset>/ALL_MODELS/comp_time/comp_time_PatientSpecificXemaPred_everyday_<Dataset>.RData
#
# Usage:
#   Rscript scripts/02_process_results/process_patient_specific_xemapred.R \
#     --datasets PFDC,Derexyl \
#     --metrics lpd,CRPS,Accuracy_median,Accuracy_map,Accuracy_prob \
#     --perf_horizon 4 \
#     --runtime_horizon 1 \
#     --overwrite 1 \
#     --debug_paths 1
# ======================================================================

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(glue)
  library(tibble)
  library(stringr)
  library(scoringRules)
})

source(here::here("scripts", "00_setup", "00_init.R"))

# Change this line if your shared helper file has another name.
source(here::here("scripts", "02_process_results", "utils", "process_common.R"))

# ----------------------------------------------------------------------
# Patient-specific XemaPred constants
# ----------------------------------------------------------------------

# Base name used in the raw result directories.
INPUT_MODEL_TYPE <- "PatientSpecificXemaPred"

INFERENCE <- "SMC"

# ----------------------------------------------------------------------
# Patient-specific output naming
# ----------------------------------------------------------------------

format_patient_specific_run_suffix <- function(run_suffix) {
  if (is.null(run_suffix) || run_suffix == "" || toupper(run_suffix) == "NONE") {
    return("")
  }

  run_suffix <- as.character(run_suffix)

  if (startsWith(run_suffix, "-")) {
    return(run_suffix)
  }

  paste0("-", run_suffix)
}

patient_specific_output_model_type <- function(run_suffix = "") {
  suffix <- format_patient_specific_run_suffix(run_suffix)

  if (!nzchar(suffix)) {
    return(INPUT_MODEL_TYPE)
  }

  version <- sub("^-", "", suffix)

  paste0(
    INPUT_MODEL_TYPE,
    version
  )
}

patient_specific_model_dir_for_item <- function(item, run_suffix = "") {
  base_model <- eczemapred_model_for_item(item)
  suffix <- format_patient_specific_run_suffix(run_suffix)

  paste0(
    base_model,
    "-",
    INPUT_MODEL_TYPE,
    suffix
  )
}

patient_specific_prediction_file <- function(ds, item, horizon, run_suffix = "") {
  model_dir <- patient_specific_model_dir_for_item(item, run_suffix = run_suffix)

  file.path(
    "results",
    ds,
    item,
    model_dir,
    paste0("H", horizon),
    "final",
    "predictions_ALL_PATIENTS.rds"
  )
}

patient_specific_diag_file <- function(ds, item, horizon, run_suffix = "") {
  model_dir <- patient_specific_model_dir_for_item(item, run_suffix = run_suffix)

  file.path(
    "results",
    ds,
    item,
    model_dir,
    paste0("H", horizon),
    "final",
    "diagnostics_ALL_PATIENTS.rds"
  )
}

# ----------------------------------------------------------------------
# File grids
# ----------------------------------------------------------------------

build_patient_specific_item_grid <- function(
    ds,
    horizon,
    run_suffix = ""
) {
  tibble(
    Dataset = ds,
    Item = eczemapred_items()
  ) %>%
    mutate(
      Model = vapply(
        .data$Item,
        eczemapred_model_for_item,
        character(1)
      ),

      Model_dir = vapply(
        .data$Item,
        function(item) {
          patient_specific_model_dir_for_item(
            item,
            run_suffix = run_suffix
          )
        },
        character(1)
      ),

      Model_type = patient_specific_output_model_type(
        run_suffix
      ),

      Component = case_when(
        .data$Item == "extent" ~ "Extent",
        .data$Item %in% c("itching", "sleep") ~
          "Subjective symptoms",
        TRUE ~ "Intensity signs"
      ),

      Inference = INFERENCE,
      Label = "Main model",
      Section = "ITEMS",

      File = vapply(
        .data$Item,
        function(item) {
          patient_specific_prediction_file(
            ds = ds,
            item = item,
            horizon = horizon,
            run_suffix = run_suffix
          )
        },
        character(1)
      )
    )
}

build_patient_specific_diag_grid <- function(ds, horizon, run_suffix = "") {
  build_patient_specific_item_grid(ds, horizon, run_suffix = run_suffix) %>%
    mutate(
      Section = "DIAGNOSTICS",
      File = vapply(
        .data$Item,
        function(item) patient_specific_diag_file(
          ds = ds,
          item = item,
          horizon = horizon,
          run_suffix = run_suffix
        ),
        character(1)
      )
    )
}

print_grid_availability <- function(grid, title) {
  log_section(title)

  out <- grid %>%
    mutate(exists = file.exists(.data$File)) %>%
    select(
      Dataset,
      Section,
      Item,
      Model,
      Model_dir,
      File,
      exists
    )

  print(out, n = Inf)

  invisible(out)
}

# ----------------------------------------------------------------------
# Item-level patient-specific XemaPred performance
# ----------------------------------------------------------------------

compute_patient_specific_items_perf <- function(
    ds,
    metric,
    perf_horizon,
    dict_datasets,
    run_suffix = "",
    debug_paths = FALSE,
    strict_files = FALSE
) {
  log_section(glue("[PatientSpecificXemaPred ITEMS] dataset={ds} metric={metric} horizon=H{perf_horizon}"))

  grid <- build_patient_specific_item_grid(
    ds = ds,
    horizon = perf_horizon,
    run_suffix = run_suffix
  )

  if (isTRUE(debug_paths)) {
    print_grid_availability(grid, "PatientSpecificXemaPred item prediction files")
  }

  missing <- grid %>% filter(!file.exists(.data$File))

  if (nrow(missing)) {
    log_warn(nrow(missing), " PatientSpecificXemaPred item prediction file(s) missing.")

    print(
      missing %>% select(Dataset, Item, Model, Model_dir, File),
      n = Inf
    )

    if (isTRUE(strict_files)) {
      stop("Stopping because strict_files=TRUE and some PatientSpecificXemaPred item files are missing.")
    }
  }

  grid <- grid %>% filter(file.exists(.data$File))

  if (!nrow(grid)) {
    log_warn("No PatientSpecificXemaPred item prediction files available for dataset=", ds)
    return(tibble())
  }

  compute_perf_from_grid(
    grid = grid,
    ds = ds,
    metric = metric,
    perf_horizon = perf_horizon,
    dict_datasets = dict_datasets,
    allow_no_samples = FALSE,
    filter_rw_iteration_zero = FALSE,
    adjust_horizon_fn = function(model) TRUE
  )
}

# ----------------------------------------------------------------------
# PO-SCORAD reconstruction from patient-specific item predictions
# ----------------------------------------------------------------------

attach_patient_specific_item_predictions <- function(
    ds,
    pred_sc,
    perf_horizon,
    run_suffix = "",
    debug_paths = FALSE
) {
  for (severity_item in eczemapred_items()) {
    model <- eczemapred_model_for_item(severity_item)
    model_dir <- patient_specific_model_dir_for_item(severity_item, run_suffix = run_suffix)

    f <- patient_specific_prediction_file(
      ds = ds,
      item = severity_item,
      horizon = perf_horizon,
      run_suffix = run_suffix
    )

    if (isTRUE(debug_paths)) {
      cat(
        "[DEBUG][POSCORAD-PatientSpecificXemaPred] Item prediction file:\n",
        "  dataset:   ", ds, "\n",
        "  item:      ", severity_item, "\n",
        "  model:     ", model, "\n",
        "  model_dir: ", model_dir, "\n",
        "  expected:  ", f, "\n",
        "  exists:    ", file.exists(f), "\n\n",
        sep = ""
      )
    }

    if (!file.exists(f)) {
      stop(
        "[POSCORAD-PatientSpecificXemaPred] Missing item prediction file\n",
        "  dataset:   ", ds, "\n",
        "  item:      ", severity_item, "\n",
        "  model:     ", model, "\n",
        "  model_dir: ", model_dir, "\n",
        "  expected:  ", f
      )
    }

    res_item <- readRDS(f) %>%
      select(Patient, Time, Horizon, Iteration, Samples) %>%
      rename(!!paste0(severity_item, "_pred") := Samples)

    pred_sc <- left_join(
      pred_sc,
      res_item,
      by = c("Patient", "Time", "Horizon", "Iteration")
    )
  }

  pred_sc %>% drop_na()
}

compute_patient_specific_poscorad_perf <- function(
    ds,
    metric,
    perf_horizon,
    dict_datasets,
    run_suffix = "",
    debug_paths = FALSE
) {
  log_section(glue("[PatientSpecificXemaPred POSCORAD] dataset={ds} metric={metric} horizon=H{perf_horizon}"))

  output_model_type <-
    patient_specific_output_model_type(
      run_suffix
    )

  tmp <- build_poscorad_testing_frame(
    ds = ds,
    perf_horizon = perf_horizon,
    dict_datasets = dict_datasets
  )

  pred_sc <- tmp$pred
  POSCORAD <- tmp$POSCORAD

  pred_sc <- attach_patient_specific_item_predictions(
    ds = ds,
    pred_sc = pred_sc,
    perf_horizon = perf_horizon,
    run_suffix = run_suffix,
    debug_paths = debug_paths
  )

  res_main <- build_scorad_samples(pred_sc)

  res_main <- ensure_metric_column(
    res = res_main,
    item = "SCORAD",
    metric = metric,
    allow_no_samples = FALSE
  )

  fc_it_sc <- detail_fc_training(
    POSCORAD %>% rename(Time = Day),
    perf_horizon
  )

  estimate_performance(
    metric,
    res_main,
    fc_it_sc,
    adjust_horizon = TRUE
  ) %>%
    filter(.data$Variable == "Fit") %>%
    mutate(
      Model_type = output_model_type,
      Label = "Main model",
      Metric = metric,
      Dataset = ds,
      Item = "SCORAD",

      Model = paste0(
        output_model_type,
        "_rebuilt_from_items"
      ),

      Model_dir = paste0(
        "SCORAD-",
        INPUT_MODEL_TYPE,
        format_patient_specific_run_suffix(
          run_suffix
        )
      ),

      Component = "SCORAD",
      Inference = INFERENCE
    )
}

# ----------------------------------------------------------------------
# Patient-specific timing / computation time
# ----------------------------------------------------------------------

normalise_patient_specific_diag <- function(
    diag_df,
    ds,
    item,
    model,
    model_dir,
    timing_horizon,
    run_suffix = ""
) {
  step_days <- if (as.integer(timing_horizon) == 1L) {
    1L
  } else {
    as.integer(timing_horizon)
  }

  diag_df <- diag_df %>%
    mutate(
      Dataset = ds,
      Item = item,
      Model = model,
      Model_dir = model_dir,
      Inference = INFERENCE,
      Model_type = patient_specific_output_model_type(run_suffix),
      timing_horizon = as.integer(timing_horizon),
      timing_step_days = as.integer(step_days)
    )

  if ("compute_time" %in% names(diag_df) && !"run_time" %in% names(diag_df)) {
    diag_df <- diag_df %>% mutate(run_time = .data$compute_time)
  }

  if ("run_time" %in% names(diag_df) && !"compute_time" %in% names(diag_df)) {
    diag_df <- diag_df %>% mutate(compute_time = .data$run_time)
  }

  if ("iter" %in% names(diag_df)) {
    diag_df <- diag_df %>%
      mutate(iter_day = as.integer(.data$iter) * as.integer(.data$timing_step_days))
  }

  diag_df
}

compute_patient_specific_comp_time <- function(
    ds,
    runtime_horizon,
    run_suffix = "",
    debug_paths = FALSE,
    strict_files = FALSE
) {
  log_section(glue("[PatientSpecificXemaPred TIMING] dataset={ds} horizon=H{runtime_horizon}"))

  grid <- build_patient_specific_diag_grid(
    ds = ds,
    horizon = runtime_horizon,
    run_suffix = run_suffix
  )

  if (isTRUE(debug_paths)) {
    print_grid_availability(grid, "PatientSpecificXemaPred diagnostic files")
  }

  missing <- grid %>% filter(!file.exists(.data$File))

  if (nrow(missing)) {
    log_warn(nrow(missing), " PatientSpecificXemaPred diagnostic file(s) missing.")

    print(
      missing %>% select(Dataset, Item, Model, Model_dir, File),
      n = Inf
    )

    if (isTRUE(strict_files)) {
      stop("Stopping because strict_files=TRUE and some PatientSpecificXemaPred diagnostic files are missing.")
    }
  }

  grid <- grid %>% filter(file.exists(.data$File))

  if (!nrow(grid)) {
    log_warn("No PatientSpecificXemaPred diagnostic files available for dataset=", ds)
    return(tibble())
  }

  map_dfr(seq_len(nrow(grid)), function(i) {
    f <- grid$File[i]

    readRDS(f) %>%
      normalise_patient_specific_diag(
        ds = ds,
        item = grid$Item[i],
        model = grid$Model[i],
        model_dir = grid$Model_dir[i],
        timing_horizon = runtime_horizon,
        run_suffix = run_suffix
      )
  })
}

save_patient_specific_comp_time_once_per_dataset <- function(
    ds,
    runtime_horizon,
    run_suffix = "",
    overwrite = TRUE,
    debug_paths = FALSE,
    strict_files = FALSE
) {
  out_dir <- comp_time_dir(ds)
  dir_create(out_dir)

  tag <- if (as.integer(runtime_horizon) == 1L) {
    "everyday"
  } else {
    paste0("every", as.integer(runtime_horizon), "days")
  }

  output_model_type <-
    patient_specific_output_model_type(
      run_suffix
    )

  out_path <- file.path(
    out_dir,
    paste0(
      "comp_time_",
      output_model_type,
      "_",
      tag,
      "_",
      ds,
      ".RData"
    )
  )

  if (file.exists(out_path) && !isTRUE(overwrite)) {
    log_info("Existing PatientSpecificXemaPred comp_time kept: ", out_path)
    return(invisible(TRUE))
  }

  comp_time <- compute_patient_specific_comp_time(
    ds = ds,
    runtime_horizon = runtime_horizon,
    run_suffix = run_suffix,
    debug_paths = debug_paths,
    strict_files = strict_files
  )

  if (!nrow(comp_time)) {
    log_warn("No PatientSpecificXemaPred comp_time rows to save for dataset=", ds)
    return(invisible(FALSE))
  }

  save(comp_time, file = out_path)
  log_save(out_path)

  invisible(TRUE)
}

# ----------------------------------------------------------------------
# Driver
# ----------------------------------------------------------------------

run_one_dataset_metric_patient_specific <- function(
    ds,
    metric,
    perf_horizon,
    overwrite,
    dict_datasets,
    run_suffix = "",
    process_items = TRUE,
    process_poscorad = TRUE,
    debug_paths = FALSE,
    strict_files = FALSE
) {
  ok <- TRUE

  if (isTRUE(process_items)) {
    perf_items <- compute_patient_specific_items_perf(
      ds = ds,
      metric = metric,
      perf_horizon = perf_horizon,
      dict_datasets = dict_datasets,
      run_suffix = run_suffix,
      debug_paths = debug_paths,
      strict_files = strict_files
    )

    if (nrow(perf_items)) {
      ok <- isTRUE(write_split_perf_by_model(
        out_root = items_metric_dir(ds, metric),
        perf_df = perf_items,
        ds = ds,
        metric = metric,
        overwrite = overwrite
      )) && ok
    } else {
      log_warn("No PatientSpecificXemaPred item performance rows for dataset=", ds, ", metric=", metric)
      ok <- FALSE
    }
  }

  if (isTRUE(process_poscorad)) {
    perf_poscorad <- compute_patient_specific_poscorad_perf(
      ds = ds,
      metric = metric,
      perf_horizon = perf_horizon,
      dict_datasets = dict_datasets,
      run_suffix = run_suffix,
      debug_paths = debug_paths
    )

    if (nrow(perf_poscorad)) {
      ok <- isTRUE(write_split_perf_by_model(
        out_root = poscorad_metric_dir(ds, metric),
        perf_df = perf_poscorad,
        ds = ds,
        metric = metric,
        overwrite = overwrite
      )) && ok
    } else {
      log_warn("No PatientSpecificXemaPred POSCORAD performance rows for dataset=", ds, ", metric=", metric)
      ok <- FALSE
    }
  }

  invisible(ok)
}

# ----------------------------------------------------------------------
# Cross-fitted population-prior run identification
#
# PopPrior-fold-1:
#   Target Fold 1 <- Population Source Fold 2
#
# PopPrior-fold-2:
#   Target Fold 2 <- Population Source Fold 1
#
# These runs are the renamed W100 experiment:
#   100% population prior
# ----------------------------------------------------------------------

TARGET_FOLDS <- c(1L, 2L)

FOLD_RUN_SUFFIX <- c(
  `1` = "PopPrior-fold-1",
  `2` = "PopPrior-fold-2"
)

# Combined full held-out cohort.
COMBINED_RUN_SUFFIX <- "PopPrior"

# Metadata only: renamed experiment is mathematically W100.
MIXTURE_WEIGHT_TAG <- "W100"
MIXTURE_POPULATION_WEIGHT <- 1.0


source_fold_for_target <- function(target_fold) {
  target_fold <- as.integer(target_fold)

  if (target_fold == 1L) return(2L)
  if (target_fold == 2L) return(1L)

  stop("[ERROR] Only target folds 1 and 2 are supported.")
}


fold_suffix_for_target <- function(target_fold) {
  key <- as.character(as.integer(target_fold))

  suffix <- FOLD_RUN_SUFFIX[[key]]

  if (is.null(suffix)) {
    stop(
      "[ERROR] No run suffix defined for target fold ",
      target_fold,
      "."
    )
  }

  suffix
}


read_required_rds <- function(file, strict_files = TRUE) {
  if (!file.exists(file)) {
    msg <- paste0("[ERROR] Missing required file: ", file)

    if (isTRUE(strict_files)) {
      stop(msg)
    }

    log_warn(msg)
    return(NULL)
  }

  readRDS(file)
}


assert_required_columns_2fold <- function(df, required, file) {
  missing <- setdiff(required, names(df))

  if (length(missing)) {
    stop(
      "[ERROR] Missing required column(s) in ",
      file,
      ": ",
      paste(missing, collapse = ", ")
    )
  }

  invisible(TRUE)
}


assert_no_patient_overlap_2fold <- function(
    fold_objects,
    ds,
    item,
    kind
) {
  if (
    is.null(fold_objects[["1"]]) ||
    is.null(fold_objects[["2"]])
  ) {
    return(invisible(TRUE))
  }

  p1 <- sort(unique(fold_objects[["1"]]$Patient))
  p2 <- sort(unique(fold_objects[["2"]]$Patient))

  overlap <- intersect(p1, p2)

  if (length(overlap)) {
    stop(
      "[ERROR] Patient overlap between target folds while combining ",
      kind,
      " for dataset=",
      ds,
      ", item=",
      item,
      ": ",
      paste(overlap, collapse = ", ")
    )
  }

  invisible(TRUE)
}


assert_no_duplicate_prediction_keys_2fold <- function(df, ds, item) {
  key_cols <- c(
    "Patient",
    "Time",
    "Horizon",
    "Iteration"
  )

  assert_required_columns_2fold(
    df,
    key_cols,
    paste0(ds, "/", item, " combined predictions")
  )

  dup <- df %>%
    count(
      across(all_of(key_cols)),
      name = "n"
    ) %>%
    filter(.data$n > 1L)

  if (nrow(dup)) {
    print(
      dup,
      n = min(50L, nrow(dup))
    )

    stop(
      "[ERROR] Duplicate prediction keys after combining folds for ",
      "dataset=",
      ds,
      ", item=",
      item,
      "."
    )
  }

  invisible(TRUE)
}


assert_no_duplicate_diagnostic_keys_2fold <- function(df, ds, item) {
  key_cols <- c(
    "Patient",
    "iter"
  )

  if (!all(key_cols %in% names(df))) {
    return(invisible(TRUE))
  }

  dup <- df %>%
    count(
      across(all_of(key_cols)),
      name = "n"
    ) %>%
    filter(.data$n > 1L)

  if (nrow(dup)) {
    print(
      dup,
      n = min(50L, nrow(dup))
    )

    stop(
      "[ERROR] Duplicate diagnostic keys after combining folds for ",
      "dataset=",
      ds,
      ", item=",
      item,
      "."
    )
  }

  invisible(TRUE)
}


# ----------------------------------------------------------------------
# Combine prediction files across the two held-out folds
#
# Used for both H1 and H4 combined prediction outputs.
# Standard manuscript forecasting performance remains based on H4.
# ----------------------------------------------------------------------

combine_prediction_folds_2fold <- function(
    ds,
    item,
    perf_horizon = 4L,
    overwrite = TRUE,
    debug_paths = FALSE,
    strict_files = TRUE
) {
  fold_objects <- list()

  for (target_fold in TARGET_FOLDS) {

    source_fold <- source_fold_for_target(
      target_fold
    )

    run_suffix <- fold_suffix_for_target(
      target_fold
    )

    f <- patient_specific_prediction_file(
      ds = ds,
      item = item,
      horizon = perf_horizon,
      run_suffix = run_suffix
    )

    if (isTRUE(debug_paths)) {
      cat(
        "[DEBUG][2FOLD-PRED]\n",
        "  dataset:      ", ds, "\n",
        "  item:         ", item, "\n",
        "  horizon:      H", perf_horizon, "\n",
        "  target fold:  ", target_fold, "\n",
        "  source fold:  ", source_fold, "\n",
        "  run suffix:   ", run_suffix, "\n",
        "  file:         ", f, "\n",
        "  exists:       ", file.exists(f), "\n\n",
        sep = ""
      )
    }

    res <- read_required_rds(
      f,
      strict_files = strict_files
    )

    if (is.null(res)) {
      next
    }

    assert_required_columns_2fold(
      res,
      c(
        "Patient",
        "Time",
        "Horizon",
        "Iteration"
      ),
      f
    )

    res <- res %>%
      mutate(
        CrossFit_TargetFold = as.integer(target_fold),
        CrossFit_SourceFold = as.integer(source_fold),
        PriorCombination = "mixture",
        MixtureWeightTag = MIXTURE_WEIGHT_TAG,
        MixturePopulationWeight = MIXTURE_POPULATION_WEIGHT
      )

    fold_objects[[as.character(target_fold)]] <- res
  }


  if (!length(fold_objects)) {
    log_warn(
      "No fold prediction files available for dataset=",
      ds,
      ", item=",
      item
    )

    return(FALSE)
  }


  assert_no_patient_overlap_2fold(
    fold_objects = fold_objects,
    ds = ds,
    item = item,
    kind = "predictions"
  )


  combined <- bind_rows(
    fold_objects
  ) %>%
    arrange(
      .data$Patient,
      .data$Iteration,
      .data$Time,
      .data$Horizon
    )


  assert_no_duplicate_prediction_keys_2fold(
    combined,
    ds = ds,
    item = item
  )


  out_file <- patient_specific_prediction_file(
    ds = ds,
    item = item,
    horizon = perf_horizon,
    run_suffix = COMBINED_RUN_SUFFIX
  )


  dir_create(
    dirname(out_file)
  )


  if (
    file.exists(out_file) &&
    !isTRUE(overwrite)
  ) {
    log_info(
      "Existing combined 2-fold prediction file kept: ",
      out_file
    )

    return(TRUE)
  }


  saveRDS(
    combined,
    out_file
  )


  log_save(
    out_file
  )

  log_info(
    "Combined 2-fold H",
    perf_horizon,
    " predictions: dataset=",
    ds,
    " item=",
    item,
    " | rows=",
    nrow(combined),
    " | patients=",
    n_distinct(combined$Patient)
  )


  invisible(TRUE)
}


# ----------------------------------------------------------------------
# Combine H1 diagnostics
#
# These are used ONLY for timing. No H1 predictive performance is
# estimated anywhere in this script.
# ----------------------------------------------------------------------

combine_diagnostic_folds_2fold <- function(
    ds,
    item,
    runtime_horizon = 1L,
    overwrite = TRUE,
    debug_paths = FALSE,
    strict_files = TRUE
) {
  fold_objects <- list()

  for (target_fold in TARGET_FOLDS) {

    source_fold <- source_fold_for_target(
      target_fold
    )

    run_suffix <- fold_suffix_for_target(
      target_fold
    )

    f <- patient_specific_diag_file(
      ds = ds,
      item = item,
      horizon = runtime_horizon,
      run_suffix = run_suffix
    )

    if (isTRUE(debug_paths)) {
      cat(
        "[DEBUG][2FOLD-DIAG]\n",
        "  dataset:      ", ds, "\n",
        "  item:         ", item, "\n",
        "  horizon:      H", runtime_horizon, "\n",
        "  target fold:  ", target_fold, "\n",
        "  source fold:  ", source_fold, "\n",
        "  run suffix:   ", run_suffix, "\n",
        "  file:         ", f, "\n",
        "  exists:       ", file.exists(f), "\n\n",
        sep = ""
      )
    }

    res <- read_required_rds(
      f,
      strict_files = strict_files
    )

    if (is.null(res)) {
      next
    }

    assert_required_columns_2fold(
      res,
      "Patient",
      f
    )

    res <- res %>%
      mutate(
        CrossFit_TargetFold = as.integer(target_fold),
        CrossFit_SourceFold = as.integer(source_fold),
        PriorCombination = "mixture",
        MixtureWeightTag = MIXTURE_WEIGHT_TAG,
        MixturePopulationWeight = MIXTURE_POPULATION_WEIGHT
      )

    fold_objects[[as.character(target_fold)]] <- res
  }


  if (!length(fold_objects)) {
    log_warn(
      "No fold diagnostic files available for dataset=",
      ds,
      ", item=",
      item
    )

    return(FALSE)
  }


  assert_no_patient_overlap_2fold(
    fold_objects = fold_objects,
    ds = ds,
    item = item,
    kind = "diagnostics"
  )


  combined <- bind_rows(
    fold_objects
  )


  if ("iter" %in% names(combined)) {
    combined <- combined %>%
      arrange(
        .data$Patient,
        .data$iter
      )
  } else {
    combined <- combined %>%
      arrange(
        .data$Patient
      )
  }


  assert_no_duplicate_diagnostic_keys_2fold(
    combined,
    ds = ds,
    item = item
  )


  out_file <- patient_specific_diag_file(
    ds = ds,
    item = item,
    horizon = runtime_horizon,
    run_suffix = COMBINED_RUN_SUFFIX
  )


  dir_create(
    dirname(out_file)
  )


  if (
    file.exists(out_file) &&
    !isTRUE(overwrite)
  ) {
    log_info(
      "Existing combined 2-fold diagnostic file kept: ",
      out_file
    )

    return(TRUE)
  }


  saveRDS(
    combined,
    out_file
  )


  log_save(
    out_file
  )

  log_info(
    "Combined 2-fold H",
    runtime_horizon,
    " diagnostics for TIMING ONLY: dataset=",
    ds,
    " item=",
    item,
    " | rows=",
    nrow(combined),
    " | patients=",
    n_distinct(combined$Patient)
  )


  invisible(TRUE)
}


# ----------------------------------------------------------------------
# Combine all 9 items for one dataset
# ----------------------------------------------------------------------

combine_one_dataset_2fold <- function(
    ds,
    perf_horizon = 4L,
    runtime_horizon = 1L,
    overwrite = TRUE,
    debug_paths = FALSE,
    strict_files = TRUE,
    record_timing = TRUE
) {
  log_header(
    glue(
      "[2-FOLD COMBINE] {ds}"
    )
  )


  ok <- TRUE


  for (item in eczemapred_items()) {

    # --------------------------------------------------------------
    # H4 predictions
    #
    # Used for the standard forecasting-performance analysis.
    # --------------------------------------------------------------

    ok_pred_h4 <- combine_prediction_folds_2fold(
      ds = ds,
      item = item,
      perf_horizon = perf_horizon,
      overwrite = overwrite,
      debug_paths = debug_paths,
      strict_files = strict_files
    )

    ok <- isTRUE(ok_pred_h4) && ok


    # --------------------------------------------------------------
    # H1 predictions
    #
    # Save the combined held-out predictions as well.
    # These can subsequently be used for H1-vs-H4 predictive
    # comparisons without changing the standard H4 analysis.
    # --------------------------------------------------------------

    ok_pred_h1 <- combine_prediction_folds_2fold(
      ds = ds,
      item = item,
      perf_horizon = runtime_horizon,
      overwrite = overwrite,
      debug_paths = debug_paths,
      strict_files = strict_files
    )

    ok <- isTRUE(ok_pred_h1) && ok


    # --------------------------------------------------------------
    # H1 diagnostics
    #
    # Used for computation-time analysis.
    # --------------------------------------------------------------

    if (isTRUE(record_timing)) {

      ok_diag <- combine_diagnostic_folds_2fold(
        ds = ds,
        item = item,
        runtime_horizon = runtime_horizon,
        overwrite = overwrite,
        debug_paths = debug_paths,
        strict_files = strict_files
      )

      ok <- isTRUE(ok_diag) && ok
    }
  }


  invisible(ok)
}


# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------

args <- commandArgs(
  trailingOnly = TRUE
)


datasets <- c(
  "Derexyl"
)

metrics <- c(
  "lpd",
  "Accuracy_map",
  "Accuracy_prob",
  "CRPS",
  "Accuracy_median"
)

# Predictive performance is ALWAYS based on H4 by default.
perf_horizon <- 4L
# computation time.
runtime_horizon <- 1L

overwrite <- TRUE

process_items <- TRUE
process_poscorad <- TRUE
record_timing <- TRUE

debug_paths <- FALSE
strict_files <- TRUE


dict_datasets <- tibble(
  Dataset = c(
    "PFDC",
    "Derexyl"
  ),
  Max_train_day = c(
    80,
    115
  )
)


i <- 1L

while (i <= length(args)) {

  key <- args[[i]]


  if (key %in% c("-h", "--help")) {

    cat(
"process_patient_specific_xemapred_2fold.R
==========================================

Purpose:
  Combine held-out 2-fold PatientSpecificXemaPred outputs and calculate
  performance directly, without calling another processing script.

Cross-fitting:
  Target Fold 1 <- Population Source Fold 2
  Target Fold 2 <- Population Source Fold 1

Population-prior inputs:
  PopPrior-fold-1 = Target Fold 1 <- Population Source Fold 2
  PopPrior-fold-2 = Target Fold 2 <- Population Source Fold 1

Combined output:
  PopPrior

Usage:
  Rscript scripts/02_process_results/process_patient_specific_xemapred_2fold.R
    [--datasets PFDC,Derexyl]
    [--metrics lpd,CRPS,Accuracy_median,Accuracy_map,Accuracy_prob]
    [--overwrite 0|1]
    [--items 0|1]
    [--poscorad 0|1]
    [--record_timing 0|1]
    [--debug_paths 0|1]
    [--strict_files 0|1]

Fixed horizons:
  H4 = predictions / item performance / PO-SCORAD performance
  H1 = diagnostics / computation time ONLY
"
    )

    quit(
      status = 0
    )


  } else if (key == "--datasets") {

    datasets <- parse_csv_arg(
      args[[i + 1L]]
    )

    i <- i + 2L


  } else if (key == "--metrics") {

    metrics <- parse_csv_arg(
      args[[i + 1L]]
    )

    i <- i + 2L
  } else if (key == "--overwrite") {

    overwrite <- bool_arg(
      args[[i + 1L]]
    )

    i <- i + 2L


  } else if (key == "--items") {

    process_items <- bool_arg(
      args[[i + 1L]]
    )

    i <- i + 2L


  } else if (key == "--poscorad") {

    process_poscorad <- bool_arg(
      args[[i + 1L]]
    )

    i <- i + 2L


  } else if (key == "--record_timing") {

    record_timing <- bool_arg(
      args[[i + 1L]]
    )

    i <- i + 2L


  } else if (key == "--debug_paths") {

    debug_paths <- bool_arg(
      args[[i + 1L]]
    )

    i <- i + 2L


  } else if (key == "--strict_files") {

    strict_files <- bool_arg(
      args[[i + 1L]]
    )

    i <- i + 2L


  } else {

    stop(
      "Unknown argument: ",
      key
    )
  }
}


bad_datasets <- setdiff(
  datasets,
  c(
    "PFDC",
    "Derexyl"
  )
)


if (length(bad_datasets)) {
  stop(
    "[ERROR] Unsupported dataset(s): ",
    paste(
      bad_datasets,
      collapse = ", "
    )
  )
}


# ----------------------------------------------------------------------
# Initial logging
# ----------------------------------------------------------------------

log_header(
  "[INIT] Process held-out 2-fold population-prior PatientSpecificXemaPred"
)

log_info(
  "datasets               : ",
  paste(
    datasets,
    collapse = ","
  )
)

log_info(
  "metrics                : ",
  paste(
    metrics,
    collapse = ","
  )
)

log_info(
  "mixture weight tag     : ",
  MIXTURE_WEIGHT_TAG
)

log_info(
  "population weight      : ",
  MIXTURE_POPULATION_WEIGHT
)


log_info(
  "PERFORMANCE horizon    : H",
  perf_horizon,
  " (predictions only)"
)

log_info(
  "TIMING horizon         : H",
  runtime_horizon,
  " (diagnostics only)"
)

log_info(
  "target F1 source suffix: ",
  FOLD_RUN_SUFFIX[["1"]]
)

log_info(
  "target F2 source suffix: ",
  FOLD_RUN_SUFFIX[["2"]]
)

log_info(
  "combined run suffix    : ",
  COMBINED_RUN_SUFFIX
)

log_info(
  "overwrite              : ",
  as.integer(overwrite)
)

log_info(
  "process_items          : ",
  as.integer(process_items)
)

log_info(
  "process_poscorad       : ",
  as.integer(process_poscorad)
)

log_info(
  "record_timing          : ",
  as.integer(record_timing)
)

log_info(
  "debug_paths            : ",
  as.integer(debug_paths)
)

log_info(
  "strict_files           : ",
  as.integer(strict_files)
)


# ----------------------------------------------------------------------
# 1. Combine raw held-out folds
# ----------------------------------------------------------------------

ok <- TRUE


for (ds in datasets) {

  ok_ds <- combine_one_dataset_2fold(
    ds = ds,
    perf_horizon = perf_horizon,
    runtime_horizon = runtime_horizon,
    overwrite = overwrite,
    debug_paths = debug_paths,
    strict_files = strict_files,
    record_timing = record_timing
  )

  ok <- isTRUE(ok_ds) && ok
}


if (
  !ok &&
  isTRUE(strict_files)
) {
  stop(
    "[ERROR] One or more 2-fold combination steps failed."
  )
}


# ----------------------------------------------------------------------
# 2. Directly process the combined held-out cohort
#
# IMPORTANT:
#
# Performance functions receive perf_horizon = 4 ONLY.
# Timing functions receive runtime_horizon = 1 ONLY.
# ----------------------------------------------------------------------

for (ds in datasets) {

  log_header(
    glue(
      "[PROCESS COMBINED 2-FOLD DATASET] {ds}"
    )
  )


  # --------------------------------------------------------------------
  # H1 DIAGNOSTICS -> TIMING ONLY
  # --------------------------------------------------------------------

  if (isTRUE(record_timing)) {

    ok <- isTRUE(
      save_patient_specific_comp_time_once_per_dataset(
        ds = ds,
        runtime_horizon = runtime_horizon,
        run_suffix = COMBINED_RUN_SUFFIX,
        overwrite = overwrite,
        debug_paths = debug_paths,
        strict_files = strict_files
      )
    ) && ok

  } else {

    log_info(
      "record_timing=0, skipping H1 timing for dataset=",
      ds
    )
  }


  # --------------------------------------------------------------------
  # H4 PREDICTIONS -> ITEM + PO-SCORAD PERFORMANCE ONLY
  # --------------------------------------------------------------------

  for (met in metrics) {

    ok <- isTRUE(
      run_one_dataset_metric_patient_specific(
        ds = ds,
        metric = met,
        perf_horizon = perf_horizon,
        overwrite = overwrite,
        dict_datasets = dict_datasets,
        run_suffix = COMBINED_RUN_SUFFIX,
        process_items = process_items,
        process_poscorad = process_poscorad,
        debug_paths = debug_paths,
        strict_files = strict_files
      )
    ) && ok
  }
}


log_header(
  if (ok) {
    "[OK] Held-out 2-fold population-prior PatientSpecificXemaPred processing complete."
  } else {
    "[WARN] Held-out 2-fold population-prior PatientSpecificXemaPred processing completed with issues."
  }
)


quit(
  status = if (ok) 0 else 1
)
