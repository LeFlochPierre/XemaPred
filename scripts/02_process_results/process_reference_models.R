#!/usr/bin/env Rscript

# ======================================================================
# process_reference_models.R
#
# Purpose:
# - Process reference-model outputs.
# - Recompute item-level reference performance.
# - Recompute PO-SCORAD reference performance from models run directly on SCORAD.
# - Save performance objects in results/<Dataset>/ALL_MODELS/.
# ======================================================================

suppressWarnings(suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(glue)
  library(tibble)
  library(stringr)
  library(scoringRules)
}))

quiet_source <- function(path) {
  suppressWarnings(suppressPackageStartupMessages(source(path)))
}

quiet_source(here::here("scripts", "00_setup", "00_init.R"))
quiet_source(here::here("scripts", "02_process_results", "utils", "process_common.R"))

# ----------------------------------------------------------------------
# Driver
# ----------------------------------------------------------------------

run_one_dataset_metric_reference <- function(
    ds,
    metric,
    perf_horizon,
    overwrite,
    dict_datasets,
    process_items,
    process_poscorad
) {
  ok <- TRUE

  if (isTRUE(process_items)) {
    item_grid <- reference_item_grid(ds, perf_horizon)

    perf_items <- compute_perf_from_grid(
      grid = item_grid,
      ds = ds,
      metric = metric,
      perf_horizon = perf_horizon,
      dict_datasets = dict_datasets,
      allow_no_samples = TRUE,
      filter_rw_iteration_zero = TRUE,
      adjust_horizon_fn = function(model) !(model %in% c("historical", "uniform"))
    )

    ok <- isTRUE(write_split_perf_by_model(
      out_root = items_metric_dir(ds, metric),
      perf_df = perf_items,
      ds = ds,
      metric = metric,
      overwrite = overwrite
    )) && ok
  }

  if (isTRUE(process_poscorad)) {
    pos_grid <- reference_poscorad_grid(ds, perf_horizon)

    perf_poscorad <- compute_perf_from_grid(
      grid = pos_grid,
      ds = ds,
      metric = metric,
      perf_horizon = perf_horizon,
      dict_datasets = dict_datasets,
      allow_no_samples = FALSE,
      filter_rw_iteration_zero = TRUE,
      adjust_horizon_fn = function(model) !(model %in% c("historical", "uniform"))
    )

    ok <- isTRUE(write_split_perf_by_model(
      out_root = poscorad_metric_dir(ds, metric),
      perf_df = perf_poscorad,
      ds = ds,
      metric = metric,
      overwrite = overwrite
    )) && ok
  }

  invisible(ok)
}

build_preflight_grids_reference <- function(
    datasets,
    perf_horizon,
    process_items,
    process_poscorad
) {
  grids <- list()

  for (ds in datasets) {
    if (isTRUE(process_items)) {
      grids[[length(grids) + 1]] <- reference_item_grid(ds, perf_horizon) %>%
        mutate(Section = paste0("Reference item predictions H", perf_horizon))
    }

    if (isTRUE(process_poscorad)) {
      grids[[length(grids) + 1]] <- reference_poscorad_grid(ds, perf_horizon) %>%
        mutate(Section = paste0("Reference PO-SCORAD predictions H", perf_horizon))
    }
  }

  grids
}

# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

datasets <- c("PFDC", "Derexyl")
metrics <- c("lpd", "CRPS", "Accuracy_median", "Accuracy_map", "Accuracy_prob")

perf_horizon <- 4L
overwrite <- TRUE
strict_files <- TRUE

process_items <- TRUE
process_poscorad <- TRUE

dict_datasets <- tibble(
  Dataset = c("PFDC", "Derexyl"),
  Max_train_day = c(80, 115)
)

i <- 1
while (i <= length(args)) {
  key <- args[[i]]

  if (key %in% c("-h", "--help")) {
    cat(
"process_reference_models.R
=====================================

Usage:
  Rscript scripts/02_process_results/process_reference_models.R
    [--datasets PFDC,Derexyl]
    [--metrics lpd,CRPS,Accuracy_median,Accuracy_map,Accuracy_prob]
    [--perf_horizon 4]
    [--overwrite 0|1]
    [--items 0|1]
    [--poscorad 0|1]
    [--strict_files 0|1]

Notes:
  - Reference models are processed.
  - PO-SCORAD reference models are read directly from score = SCORAD predictions.
  - The script always prints a compact input-file availability check.
"
    )
    quit(status = 0)

  } else if (key == "--datasets") {
    datasets <- parse_csv_arg(args[[i + 1]])
    i <- i + 2

  } else if (key == "--metrics") {
    metrics <- parse_csv_arg(args[[i + 1]])
    i <- i + 2

  } else if (key %in% c("--perf_horizon", "--t_horizon")) {
    perf_horizon <- as.integer(args[[i + 1]])
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

  } else if (key == "--strict_files") {
    strict_files <- bool_arg(args[[i + 1]])
    i <- i + 2

  } else {
    stop("Unknown argument: ", key)
  }
}

log_header("Recompute reference-model performance")
log_info("datasets        : ", paste(datasets, collapse = ","))
log_info("metrics         : ", paste(metrics, collapse = ","))
log_info("perf_horizon    : H", perf_horizon)
log_info("overwrite       : ", as.integer(overwrite))
log_info("process_items   : ", as.integer(process_items))
log_info("process_poscorad: ", as.integer(process_poscorad))
log_info("strict_files    : ", as.integer(strict_files))

preflight_files(
  build_preflight_grids_reference(
    datasets = datasets,
    perf_horizon = perf_horizon,
    process_items = process_items,
    process_poscorad = process_poscorad
  ),
  strict = strict_files
)

ok <- TRUE

for (ds in datasets) {
  log_section(glue("Dataset: {ds}"))

  for (met in metrics) {
    log_info("Processing ", ds, " / ", met)

    ok <- isTRUE(run_one_dataset_metric_reference(
      ds = ds,
      metric = met,
      perf_horizon = perf_horizon,
      overwrite = overwrite,
      dict_datasets = dict_datasets,
      process_items = process_items,
      process_poscorad = process_poscorad
    )) && ok
  }
}

log_header(if (ok) "[OK] Done." else "[WARN] Done with issues.")
quit(status = if (ok) 0 else 1)
