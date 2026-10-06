# -----------------------------------------------------------------------
# Main runner for EczemaPred item-level models
#
# Purpose:
# - Run validation for the three EczemaPred item-level model families:
#   BinMC, OrderedRW, and BinRW.
# - Resume from existing iteration files when jobs are interrupted.
# - Combine predictions and diagnostics at the end.
#
# Run example:
# Rscript scripts/01_run_models/EczemaPred/run_validation.R \
#   configs/Derexyl/dryness_OrderedRW.yaml
# -----------------------------------------------------------------------

rm(list = ls())

cat("[INFO] Starting EczemaPred item-model validation script...\n")

source(here::here("scripts", "00_setup", "00_init.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(glue)
  library(yaml)
  library(jsonlite)
  library(foreach)
  library(doParallel)
})

source(here::here("scripts", "01_run_models", "utils", "config.R"))
source(here::here("scripts", "01_run_models", "utils", "paths.R"))
source(here::here("scripts", "01_run_models", "utils", "data_prep.R"))
source(here::here("scripts", "01_run_models", "utils", "diagnostics.R"))
source(here::here("scripts", "01_run_models", "utils", "combine_results.R"))
source(here::here("scripts", "01_run_models", "utils", "parallel.R"))
source(here::here("scripts", "01_run_models", "utils", "fit_stan_iteration.R"))

source(here::here("scripts", "01_run_models", "EczemaPred", "fit_iteration.R"))
source(here::here("scripts", "01_run_models", "EczemaPred", "fit_final.R"))
source(here::here("scripts", "01_run_models", "EczemaPred", "power_priors.R"))

cfg <- read_run_config()
run_info <- validate_run_config(cfg, expected_family = "eczemapred_item_models")

# RNGkind("L'Ecuyer-CMRG")
set.seed(run_info$seed)
cat(glue("[INFO] RNG seed set to {run_info$seed}\n"))

paths <- make_output_paths(run_info)
save_run_metadata(cfg, run_info, paths)

data_obj <- prepare_forward_chaining_data(run_info)

if (!isTRUE(run_info$run)) {
  cat("[INFO] run = FALSE. Configuration checked, but no model was run.\n")
  q("no")
}

if (isTRUE(run_info$fit_only)) {
  run_final_fit(data_obj, run_info, paths)
  combine_iteration_results(paths, expected_iterations = data_obj$train_it)
  combine_diagnostics(paths)
  q("no")
}

worker_sources <- c(
  file.path("scripts", "01_run_models", "utils", "data_prep.R"),
  file.path("scripts", "01_run_models", "utils", "diagnostics.R"),
  file.path("scripts", "01_run_models", "utils", "fit_stan_iteration.R"),
  file.path("scripts", "01_run_models", "EczemaPred", "fit_iteration.R")
)

run_missing_iterations_parallel(
  data_obj = data_obj,
  run_info = run_info,
  paths = paths,
  fit_function_name = "fit_eczemapred_item_iteration",
  worker_sources = worker_sources
)

combine_iteration_results(paths, expected_iterations = data_obj$train_it)
combine_diagnostics(paths)

cat("[INFO] EczemaPred item-model validation complete.\n")
