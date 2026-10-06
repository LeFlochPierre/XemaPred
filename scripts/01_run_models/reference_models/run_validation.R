# -----------------------------------------------------------------------
# Main runner for reference/comparison models
#
# Purpose:
# - Run validation for uniform, historical, RW, AR1, MixedAR1, Smoothing,
#   and MC reference models.
# - Keep reference-model logic separate from the EczemaPred item models.
# - Resume missing iterations and combine final outputs.
#
# Run example:
# Rscript scripts/01_run_models/reference_models/run_validation.R \
#   configs/Derexyl/SCORAD_MixedAR1.yaml
# -----------------------------------------------------------------------

rm(list = ls())
set.seed(1744834965)

cat("[INFO] Starting reference-model validation script...\n")

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

# -----------------------------------------------------------------------
# Source shared utilities
# -----------------------------------------------------------------------

source(here::here("scripts", "01_run_models", "utils", "config.R"))
source(here::here("scripts", "01_run_models", "utils", "paths.R"))
source(here::here("scripts", "01_run_models", "utils", "data_prep.R"))
source(here::here("scripts", "01_run_models", "utils", "diagnostics.R"))
source(here::here("scripts", "01_run_models", "utils", "combine_results.R"))
source(here::here("scripts", "01_run_models", "utils", "parallel.R"))
source(here::here("scripts", "01_run_models", "utils", "fit_stan_iteration.R"))

# -----------------------------------------------------------------------
# Source reference-model fitting functions
# -----------------------------------------------------------------------

source(here::here(
  "scripts",
  "01_run_models",
  "reference_models",
  "fit_reference_models.R"
))

# -----------------------------------------------------------------------
# Read configuration and prepare run
# -----------------------------------------------------------------------

cfg <- read_run_config()

run_info <- validate_run_config(
  cfg,
  expected_family = "reference_models"
)

paths <- make_output_paths(run_info)

save_run_metadata(cfg, run_info, paths)

data_obj <- prepare_forward_chaining_data(run_info)

# -----------------------------------------------------------------------
# Stop after validation if run = FALSE
# -----------------------------------------------------------------------

if (!isTRUE(run_info$run)) {
  cat("[INFO] run = FALSE. Configuration checked, but no model was run.\n")
  q("no")
}

# -----------------------------------------------------------------------
# Run missing forward-chaining iterations
# -----------------------------------------------------------------------

worker_sources <- c(
  file.path("scripts", "01_run_models", "utils", "data_prep.R"),
  file.path("scripts", "01_run_models", "utils", "diagnostics.R"),
  file.path("scripts", "01_run_models", "utils", "fit_stan_iteration.R"),
  file.path("scripts", "01_run_models", "reference_models", "fit_reference_models.R")
)

run_missing_iterations_parallel(
  data_obj = data_obj,
  run_info = run_info,
  paths = paths,
  fit_function_name = "fit_reference_iteration",
  worker_sources = worker_sources
)

# -----------------------------------------------------------------------
# Combine outputs
# -----------------------------------------------------------------------

combine_iteration_results(
  paths,
  expected_iterations = data_obj$train_it
)

combine_diagnostics(paths)

cat("[INFO] Reference-model validation complete.\n")