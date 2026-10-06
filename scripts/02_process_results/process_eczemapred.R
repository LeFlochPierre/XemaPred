#!/usr/bin/env Rscript

# ======================================================================
# process_eczemapred.R
#
# Purpose:
# - Process EczemaPred main-model outputs.
# - Recompute item-level performance.
# - Rebuild PO-SCORAD from the 9 item-level main-model predictions.
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
# EczemaPred computation time
# ----------------------------------------------------------------------

compute_comp_time_eczemapred <- function(ds, runtime_horizon) {
  diag_grid <- build_eczemapred_diag_grid(ds, runtime_horizon)

  map_dfr(seq_len(nrow(diag_grid)), function(i) {
    f <- diag_grid$File[i]

    if (!file.exists(f)) {
      stop("Missing EczemaPred diagnostics file after preflight: ", f)
    }

    readRDS(f) %>%
      mutate(
        Item = diag_grid$Item[i],
        Model = diag_grid$Model[i],
        Dataset = ds,
        Inference = "HMC",
        Model_type = "EczemaPred"
      )
  })
}

save_comp_time_eczemapred <- function(ds, runtime_horizon, overwrite = TRUE) {
  out_dir <- comp_time_dir(ds)
  dir_create(out_dir)

  out_path <- file.path(
    out_dir,
    paste0("comp_time_EczemaPred_everyday_", ds, ".RData")
  )

  if (file.exists(out_path) && !isTRUE(overwrite)) {
    log_info("Existing comp_time kept: ", out_path)
    return(invisible(TRUE))
  }

  comp_time <- compute_comp_time_eczemapred(ds, runtime_horizon)
  save(comp_time, file = out_path)
  log_save(out_path)

  invisible(TRUE)
}

# ----------------------------------------------------------------------
# EczemaPred PO-SCORAD rebuilt from item predictions
# ----------------------------------------------------------------------

compute_eczemapred_poscorad_perf <- function(ds, metric, perf_horizon, dict_datasets) {
  tmp <- build_poscorad_testing_frame(ds, perf_horizon, dict_datasets)

  pred_sc <- attach_eczemapred_item_predictions(
    ds = ds,
    pred_sc = tmp$pred,
    perf_horizon = perf_horizon
  )

  res_main <- build_scorad_samples(pred_sc) %>%
    ensure_metric_column(
      item = "SCORAD",
      metric = metric,
      allow_no_samples = FALSE
    )

  fc_it_sc <- detail_fc_training(
    tmp$POSCORAD %>% rename(Time = Day),
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
      Model_type = "EczemaPred",
      Label = "Main model",
      Metric = metric,
      Dataset = ds,
      Item = "SCORAD",
      Model = "EczemaPred_rebuilt_from_items",
      Component = "SCORAD",
      Inference = "HMC"
    )
}

# ----------------------------------------------------------------------
# Driver
# ----------------------------------------------------------------------

run_one_dataset_metric_eczemapred <- function(
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
    grid <- build_eczemapred_item_grid(ds, perf_horizon)

    perf_items <- compute_perf_from_grid(
      grid = grid,
      ds = ds,
      metric = metric,
      perf_horizon = perf_horizon,
      dict_datasets = dict_datasets,
      allow_no_samples = FALSE,
      filter_rw_iteration_zero = FALSE,
      adjust_horizon_fn = function(model) TRUE
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
    perf_poscorad <- compute_eczemapred_poscorad_perf(
      ds = ds,
      metric = metric,
      perf_horizon = perf_horizon,
      dict_datasets = dict_datasets
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

build_preflight_grids_eczemapred <- function(
    datasets,
    perf_horizon,
    runtime_horizon,
    process_items,
    process_poscorad
) {
  grids <- list()

  for (ds in datasets) {
    if (isTRUE(process_items) || isTRUE(process_poscorad)) {
      grids[[length(grids) + 1]] <- build_eczemapred_item_grid(ds, perf_horizon) %>%
        mutate(Section = paste0("EczemaPred predictions H", perf_horizon))
    }

    grids[[length(grids) + 1]] <- build_eczemapred_diag_grid(ds, runtime_horizon) %>%
      mutate(Section = paste0("EczemaPred diagnostics H", runtime_horizon))
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
runtime_horizon <- 1L
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
"process_eczemapred.R
================================

Usage:
  Rscript scripts/02_process_results/process_eczemapred.R
    [--datasets PFDC,Derexyl]
    [--metrics lpd,CRPS,Accuracy_median,Accuracy_map,Accuracy_prob]
    [--perf_horizon 4]
    [--runtime_horizon 1]
    [--overwrite 0|1]
    [--items 0|1]
    [--poscorad 0|1]
    [--strict_files 0|1]

Notes:
  - EczemaPred main models are processed.
  - PO-SCORAD is rebuilt from item-level EczemaPred predictions.
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

  } else if (key == "--runtime_horizon") {
    runtime_horizon <- as.integer(args[[i + 1]])
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

log_header("Recompute EczemaPred performance")
log_info("datasets        : ", paste(datasets, collapse = ","))
log_info("metrics         : ", paste(metrics, collapse = ","))
log_info("perf_horizon    : H", perf_horizon)
log_info("runtime_horizon : H", runtime_horizon)
log_info("overwrite       : ", as.integer(overwrite))
log_info("process_items   : ", as.integer(process_items))
log_info("process_poscorad: ", as.integer(process_poscorad))
log_info("strict_files    : ", as.integer(strict_files))

preflight_files(
  build_preflight_grids_eczemapred(
    datasets = datasets,
    perf_horizon = perf_horizon,
    runtime_horizon = runtime_horizon,
    process_items = process_items,
    process_poscorad = process_poscorad
  ),
  strict = strict_files
)

ok <- TRUE

for (ds in datasets) {
  log_section(glue("Dataset: {ds}"))

  ok <- isTRUE(save_comp_time_eczemapred(
    ds = ds,
    runtime_horizon = runtime_horizon,
    overwrite = overwrite
  )) && ok

  for (met in metrics) {
    log_info("Processing ", ds, " / ", met)

    ok <- isTRUE(run_one_dataset_metric_eczemapred(
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
