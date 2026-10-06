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
# CLI
# ----------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

datasets <- c("PFDC", "Derexyl")
metrics <- c("lpd", "CRPS", "Accuracy_median", "Accuracy_map", "Accuracy_prob")

perf_horizon <- 4L
runtime_horizon <- 1L
overwrite <- TRUE

run_suffix <- ""
process_items <- TRUE
process_poscorad <- TRUE
record_timing <- TRUE
debug_paths <- FALSE
strict_files <- FALSE

dict_datasets <- tibble(
  Dataset = c("PFDC", "Derexyl"),
  Max_train_day = c(80, 115)
)

i <- 1
while (i <= length(args)) {
  key <- args[[i]]

  if (key %in% c("-h", "--help")) {
    cat(
"process_patient_specific_xemapred.R
========================================

Usage:
  Rscript scripts/02_process_results/process_patient_specific_xemapred.R
    [--datasets PFDC,Derexyl]
    [--metrics lpd,CRPS,Accuracy_median,Accuracy_map,Accuracy_prob]
    [--perf_horizon 4]
    [--runtime_horizon 1]
    [--run_suffix NONE]
    [--overwrite 0|1]
    [--items 0|1]
    [--poscorad 0|1]
    [--record_timing 0|1]
    [--debug_paths 0|1]
    [--strict_files 0|1]

Notes:
  - PatientSpecificXemaPred is processed.
  - No population-level XemaPred.
  - No EczemaPred.
  - No reference models.
  - Performance uses perf_horizon, usually H4.
  - Timing uses runtime_horizon, usually H1.
  - Expected model folders are:
      BinMC-PatientSpecificXemaPred[-run_suffix]
      BinRW-PatientSpecificXemaPred[-run_suffix]
      OrderedRW-PatientSpecificXemaPred[-run_suffix]
"
    )
    quit(status = 0)

  } else if (key == "--datasets") {
    datasets <- parse_csv_arg(args[[i + 1]])
    i <- i + 2

  } else if (key == "--metrics") {
    metrics <- parse_csv_arg(args[[i + 1]])
    i <- i + 2

  } else if (key == "--perf_horizon") {
    perf_horizon <- as.integer(args[[i + 1]])
    i <- i + 2

  } else if (key == "--t_horizon") {
    perf_horizon <- as.integer(args[[i + 1]])
    i <- i + 2

  } else if (key == "--runtime_horizon") {
    runtime_horizon <- as.integer(args[[i + 1]])
    i <- i + 2

  } else if (key == "--run_suffix") {
    run_suffix <- args[[i + 1]]
    if (toupper(run_suffix) == "NONE") {
      run_suffix <- ""
    }
    i <- i + 2

  } else if (key == "--overwrite") {
    overwrite <- bool_arg(args[[i + 1]])
    i <- i + 2

  } else if (key == "--items") {
    process_items <- bool_arg(args[[i + 1]])
    i <- i + 2

  } else if (key == "--poscorad") {
    process_poscorad <- bool_arg(args[[i + 1]])
    i <- i + 2

  } else if (key == "--record_timing") {
    record_timing <- bool_arg(args[[i + 1]])
    i <- i + 2

  } else if (key == "--debug_paths") {
    debug_paths <- bool_arg(args[[i + 1]])
    i <- i + 2

  } else if (key == "--strict_files") {
    strict_files <- bool_arg(args[[i + 1]])
    i <- i + 2

  } else {
    stop("Unknown argument: ", key)
  }
}

log_header("[INIT] Recompute PatientSpecificXemaPred performance")

log_info("datasets        : ", paste(datasets, collapse = ","))
log_info("metrics         : ", paste(metrics, collapse = ","))
log_info("perf_horizon    : ", perf_horizon)
log_info("runtime_horizon : ", runtime_horizon)
log_info("run_suffix      : ", ifelse(nzchar(run_suffix), run_suffix, "NONE"))
log_info("overwrite       : ", as.integer(overwrite))
log_info("process_items   : ", as.integer(process_items))
log_info("process_poscorad: ", as.integer(process_poscorad))
log_info("record_timing   : ", as.integer(record_timing))
log_info("debug_paths     : ", as.integer(debug_paths))
log_info("strict_files    : ", as.integer(strict_files))

ok <- TRUE

for (ds in datasets) {
  log_header(glue("[DATASET] {ds}"))

  if (isTRUE(record_timing)) {
    ok <- isTRUE(save_patient_specific_comp_time_once_per_dataset(
      ds = ds,
      runtime_horizon = runtime_horizon,
      run_suffix = run_suffix,
      overwrite = overwrite,
      debug_paths = debug_paths,
      strict_files = strict_files
    )) && ok
  } else {
    log_info("record_timing=0, skipping PatientSpecificXemaPred comp_time for dataset=", ds)
  }

  for (met in metrics) {
    ok <- isTRUE(run_one_dataset_metric_patient_specific(
      ds = ds,
      metric = met,
      perf_horizon = perf_horizon,
      overwrite = overwrite,
      dict_datasets = dict_datasets,
      run_suffix = run_suffix,
      process_items = process_items,
      process_poscorad = process_poscorad,
      debug_paths = debug_paths,
      strict_files = strict_files
    )) && ok
  }
}

log_header(if (ok) "[OK] Done." else "[WARN] Done with issues.")

quit(status = if (ok) 0 else 1)