#!/usr/bin/env Rscript

# ======================================================================
# process_xemapred.R
#
# Purpose:
# - Process population-level XemaPred outputs.
# - No EczemaPred.
# - No reference models.
# - No PatientSpecificXemaPred.
# - Read predictions_ALL_PATIENTS.rds and diagnostics_ALL_PATIENTS.rds.
# - Recompute item-level and PO-SCORAD performance.
#
# Expected raw inputs:
#   results/<Dataset>/<Item>/<Model>-XemaPred[-run_suffix]/H<horizon>/final/predictions_ALL_PATIENTS.rds
#   results/<Dataset>/<Item>/<Model>-XemaPred[-run_suffix]/H1/final/diagnostics_ALL_PATIENTS.rds
#
# Main XemaPred item models:
#   extent                     -> BinMC-XemaPred
#   itching, sleep             -> BinRW-XemaPred
#   intensity signs            -> OrderedRW-XemaPred
#
# Outputs:
#   results/<Dataset>/ALL_MODELS/ITEMS/<Metric>/perf_XemaPred_<Dataset>_<Metric>.RData
#   results/<Dataset>/ALL_MODELS/POSCORAD/<Metric>/perf_XemaPred_<Dataset>_<Metric>.RData
#   results/<Dataset>/ALL_MODELS/comp_time/comp_time_XemaPred_everyday_<Dataset>.RData
#
# Usage:
#   Rscript scripts/02_process_results/process_xemapred.R \
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
source(here::here("scripts", "02_process_results", "utils", "process_common.R"))

# ----------------------------------------------------------------------
# XemaPred-specific constants
# ----------------------------------------------------------------------

MODEL_TYPE <- "XemaPred"
INFERENCE <- "SMC"

# ----------------------------------------------------------------------
# XemaPred output naming
# ----------------------------------------------------------------------

format_xemapred_run_suffix <- function(run_suffix) {
  if (is.null(run_suffix) || run_suffix == "" || toupper(run_suffix) == "NONE") {
    return("")
  }

  run_suffix <- as.character(run_suffix)

  if (startsWith(run_suffix, "-")) {
    return(run_suffix)
  }

  paste0("-", run_suffix)
}

xemapred_model_dir_for_item <- function(item, run_suffix = "") {
  base_model <- eczemapred_model_for_item(item)
  suffix <- format_xemapred_run_suffix(run_suffix)

  paste0(base_model, "-XemaPred", suffix)
}

xemapred_prediction_file <- function(ds, item, horizon, run_suffix = "") {
  model_dir <- xemapred_model_dir_for_item(item, run_suffix = run_suffix)

  file.path(
    "results",
    ds,
    item,
    model_dir,
    paste0("H", horizon),
    "final",
    "predictions.rds"
  )
}

xemapred_diag_file <- function(ds, item, horizon, run_suffix = "") {
  model_dir <- xemapred_model_dir_for_item(item, run_suffix = run_suffix)

  file.path(
    "results",
    ds,
    item,
    model_dir,
    paste0("H", horizon),
    "final",
    "diagnostics.rds"
  )
}

# ----------------------------------------------------------------------
# XemaPred file grids
# ----------------------------------------------------------------------

build_xemapred_item_grid <- function(ds, horizon, run_suffix = "") {
  tibble(
    Dataset = ds,
    Item = eczemapred_items()
  ) %>%
    mutate(
      Model = vapply(.data$Item, eczemapred_model_for_item, character(1)),
      Model_dir = vapply(
        .data$Item,
        function(item) xemapred_model_dir_for_item(item, run_suffix = run_suffix),
        character(1)
      ),
      Model_type = MODEL_TYPE,
      Component = case_when(
        .data$Item == "extent" ~ "Extent",
        .data$Item %in% c("itching", "sleep") ~ "Subjective symptoms",
        TRUE ~ "Intensity signs"
      ),
      Inference = INFERENCE,
      Label = "Main model",
      Section = "ITEMS",
      File = vapply(
        .data$Item,
        function(item) xemapred_prediction_file(
          ds = ds,
          item = item,
          horizon = horizon,
          run_suffix = run_suffix
        ),
        character(1)
      )
    )
}

build_xemapred_diag_grid <- function(ds, horizon, run_suffix = "") {
  build_xemapred_item_grid(ds, horizon, run_suffix = run_suffix) %>%
    mutate(
      Section = "DIAGNOSTICS",
      File = vapply(
        .data$Item,
        function(item) xemapred_diag_file(
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
# Item-level XemaPred performance
# ----------------------------------------------------------------------

compute_xemapred_items_perf <- function(
    ds,
    metric,
    perf_horizon,
    dict_datasets,
    run_suffix = "",
    debug_paths = FALSE,
    strict_files = FALSE
) {
  log_section(glue("[XemaPred ITEMS] dataset={ds} metric={metric} horizon=H{perf_horizon}"))

  grid <- build_xemapred_item_grid(
    ds = ds,
    horizon = perf_horizon,
    run_suffix = run_suffix
  )

  if (isTRUE(debug_paths)) {
    print_grid_availability(grid, "XemaPred item prediction files")
  }

  missing <- grid %>% filter(!file.exists(.data$File))

  if (nrow(missing)) {
    log_warn(nrow(missing), " XemaPred item prediction file(s) missing.")

    print(
      missing %>% select(Dataset, Item, Model, Model_dir, File),
      n = Inf
    )

    if (isTRUE(strict_files)) {
      stop("Stopping because strict_files=TRUE and some XemaPred item files are missing.")
    }
  }

  grid <- grid %>% filter(file.exists(.data$File))

  if (!nrow(grid)) {
    log_warn("No XemaPred item prediction files available for dataset=", ds)
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
# PO-SCORAD reconstruction from XemaPred item predictions
# ----------------------------------------------------------------------

attach_xemapred_item_predictions <- function(
    ds,
    pred_sc,
    perf_horizon,
    run_suffix = "",
    debug_paths = FALSE
) {
  for (severity_item in eczemapred_items()) {
    model <- eczemapred_model_for_item(severity_item)
    model_dir <- xemapred_model_dir_for_item(severity_item, run_suffix = run_suffix)

    f <- xemapred_prediction_file(
      ds = ds,
      item = severity_item,
      horizon = perf_horizon,
      run_suffix = run_suffix
    )

    if (isTRUE(debug_paths)) {
      cat(
        "[DEBUG][POSCORAD-XemaPred] Item prediction file:\n",
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
        "[POSCORAD-XemaPred] Missing XemaPred item prediction file\n",
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

compute_xemapred_poscorad_perf <- function(
    ds,
    metric,
    perf_horizon,
    dict_datasets,
    run_suffix = "",
    debug_paths = FALSE
) {
  log_section(glue("[XemaPred POSCORAD] dataset={ds} metric={metric} horizon=H{perf_horizon}"))

  tmp <- build_poscorad_testing_frame(
    ds = ds,
    perf_horizon = perf_horizon,
    dict_datasets = dict_datasets
  )

  pred_sc <- tmp$pred
  POSCORAD <- tmp$POSCORAD

  pred_sc <- attach_xemapred_item_predictions(
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
      Model_type = MODEL_TYPE,
      Label = "Main model",
      Metric = metric,
      Dataset = ds,
      Item = "SCORAD",
      Model = "XemaPred_rebuilt_from_items",
      Model_dir = paste0("SCORAD-", MODEL_TYPE, format_xemapred_run_suffix(run_suffix)),
      Component = "SCORAD",
      Inference = INFERENCE
    )
}

# ----------------------------------------------------------------------
# XemaPred timing / computation time
# ----------------------------------------------------------------------

normalise_xemapred_diag <- function(diag_df, ds, item, model, model_dir, timing_horizon) {
  diag_df <- diag_df %>%
    mutate(
      Dataset = ds,
      Item = item,
      Model = model,
      Model_dir = model_dir,
      Inference = INFERENCE,
      Model_type = MODEL_TYPE,
      timing_horizon = timing_horizon,
      timing_step_days = ifelse(as.integer(timing_horizon) == 1L, 1L, as.integer(timing_horizon))
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

dedupe_xemapred_timing <- function(comp_time) {
  if (!nrow(comp_time)) {
    return(comp_time)
  }

  if (!all(c("Item", "Model_dir", "iter") %in% names(comp_time))) {
    return(comp_time)
  }

  numeric_cols <- names(comp_time)[vapply(comp_time, is.numeric, logical(1))]
  numeric_cols <- setdiff(numeric_cols, c("iter"))

  other_cols <- setdiff(names(comp_time), c(numeric_cols, "iter"))

  comp_time %>%
    group_by(across(all_of(c("Dataset", "Item", "Model", "Model_dir", "iter")))) %>%
    summarise(
      across(all_of(numeric_cols), ~ dplyr::first(.x)),
      across(
        all_of(setdiff(other_cols, c("Dataset", "Item", "Model", "Model_dir"))),
        ~ dplyr::first(.x)
      ),
      .groups = "drop"
    )
}

compute_xemapred_comp_time <- function(
    ds,
    runtime_horizon,
    run_suffix = "",
    debug_paths = FALSE,
    strict_files = FALSE
) {
  log_section(glue("[XemaPred TIMING] dataset={ds} horizon=H{runtime_horizon}"))

  grid <- build_xemapred_diag_grid(
    ds = ds,
    horizon = runtime_horizon,
    run_suffix = run_suffix
  )

  if (isTRUE(debug_paths)) {
    print_grid_availability(grid, "XemaPred diagnostic files")
  }

  missing <- grid %>% filter(!file.exists(.data$File))

  if (nrow(missing)) {
    log_warn(nrow(missing), " XemaPred diagnostic file(s) missing.")

    print(
      missing %>% select(Dataset, Item, Model, Model_dir, File),
      n = Inf
    )

    if (isTRUE(strict_files)) {
      stop("Stopping because strict_files=TRUE and some XemaPred diagnostic files are missing.")
    }
  }

  grid <- grid %>% filter(file.exists(.data$File))

  if (!nrow(grid)) {
    log_warn("No XemaPred diagnostic files available for dataset=", ds)
    return(tibble())
  }

  comp_time <- map_dfr(seq_len(nrow(grid)), function(i) {
    f <- grid$File[i]

    readRDS(f) %>%
      normalise_xemapred_diag(
        ds = ds,
        item = grid$Item[i],
        model = grid$Model[i],
        model_dir = grid$Model_dir[i],
        timing_horizon = runtime_horizon
      )
  })

  dedupe_xemapred_timing(comp_time)
}

save_xemapred_comp_time_once_per_dataset <- function(
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

  out_path <- file.path(
      out_dir,
      paste0("comp_time_", MODEL_TYPE, "_", tag, "_", ds, ".RData")
    )

  if (file.exists(out_path) && !isTRUE(overwrite)) {
    log_info("Existing XemaPred comp_time kept: ", out_path)
    return(invisible(TRUE))
  }

  comp_time <- compute_xemapred_comp_time(
    ds = ds,
    runtime_horizon = runtime_horizon,
    run_suffix = run_suffix,
    debug_paths = debug_paths,
    strict_files = strict_files
  )

  if (!nrow(comp_time)) {
    log_warn("No XemaPred comp_time rows to save for dataset=", ds)
    return(invisible(FALSE))
  }

  save(comp_time, file = out_path)
  log_save(out_path)

  invisible(TRUE)
}

# ----------------------------------------------------------------------
# Driver
# ----------------------------------------------------------------------

run_one_dataset_metric_xemapred <- function(
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
    perf_items <- compute_xemapred_items_perf(
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
      log_warn("No XemaPred item performance rows for dataset=", ds, ", metric=", metric)
      ok <- FALSE
    }
  }

  if (isTRUE(process_poscorad)) {
    perf_poscorad <- compute_xemapred_poscorad_perf(
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
      log_warn("No XemaPred POSCORAD performance rows for dataset=", ds, ", metric=", metric)
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
"process_xemapred.R
=======================

Usage:
  Rscript scripts/02_process_results/process_xemapred.R
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
  - Population-level XemaPred is processed.
  - No EczemaPred.
  - No reference models.
  - No PatientSpecificXemaPred.
  - Performance uses perf_horizon, usually H4.
  - Timing uses runtime_horizon, usually H1.
  - Expected model folders are:
      BinMC-XemaPred[-run_suffix]
      BinRW-XemaPred[-run_suffix]
      OrderedRW-XemaPred[-run_suffix]
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

# Propagate run_suffix into the OUTPUT identity, not just the input paths.
# build_xemapred_item_grid(), normalise_xemapred_diag() and
# compute_xemapred_poscorad_perf() read MODEL_TYPE at call time from the
# global environment, so reassigning it here is sufficient.
if (nzchar(run_suffix)) {
  MODEL_TYPE <- paste0("XemaPred", gsub("^-", "", format_xemapred_run_suffix(run_suffix)))
}

log_header("[INIT] Recompute XemaPred performance")

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
    ok <- isTRUE(save_xemapred_comp_time_once_per_dataset(
      ds = ds,
      runtime_horizon = runtime_horizon,
      run_suffix = run_suffix,
      overwrite = overwrite,
      debug_paths = debug_paths,
      strict_files = strict_files
    )) && ok
  } else {
    log_info("record_timing=0, skipping XemaPred comp_time for dataset=", ds)
  }

  for (met in metrics) {
    ok <- isTRUE(run_one_dataset_metric_xemapred(
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