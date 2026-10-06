# -----------------------------------------------------------------------
# Main runner for population-level XemaPred
#
# Purpose:
# - Run population-level/cohort XemaPred using SMC2.
# - Read one YAML config.
# - Prepare data, initialise/resume SMC2 state, run forward chaining,
#   and combine outputs.
#
# Run example:
# Rscript scripts/01_run_models/XemaPred/run_validation.R \
#   configs/PFDC/XemaPred/t_horizon_1/pfdc_dryness_orderedrw.yaml
# -----------------------------------------------------------------------

rm(list = ls())

cat("[INFO] Starting XemaPred validation script...\n")

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
  library(doRNG)
})

source(here::here("scripts", "01_run_models", "XemaPred", "shared", "utils.R"))
source(here::here("scripts", "01_run_models", "XemaPred", "shared", "particle_utils.R"))
source(here::here("scripts", "01_run_models", "XemaPred", "shared", "score_info.R"))
source(here::here("scripts", "01_run_models", "XemaPred", "shared", "config.R"))
source(here::here("scripts", "01_run_models", "XemaPred", "shared", "paths.R"))
source(here::here("scripts", "01_run_models", "XemaPred", "shared", "data_prep.R"))
source(here::here("scripts", "01_run_models", "XemaPred", "shared", "parallel.R"))

source(here::here("scripts", "01_run_models", "XemaPred", "source_model_files.R"))
source(here::here("scripts", "01_run_models", "XemaPred", "state_cache.R"))
source(here::here("scripts", "01_run_models", "XemaPred", "build_observations.R"))
source(here::here("scripts", "01_run_models", "XemaPred", "update_state.R"))
source(here::here("scripts", "01_run_models", "XemaPred", "forecast_score.R"))
source(here::here("scripts", "01_run_models", "XemaPred", "posterior_snapshots.R"))
source(here::here("scripts", "01_run_models", "XemaPred", "combine_results.R"))
source(here::here("scripts", "01_run_models", "XemaPred", "run_cohort_validation.R"))

source_xemapred_model_files()

cfg <- read_xemapred_config()
run_info <- validate_xemapred_config(cfg)

RNGkind("L'Ecuyer-CMRG")
set.seed(run_info$seed)

paths <- make_xemapred_paths(run_info)
save_xemapred_metadata(cfg, run_info, paths)

data_obj <- prepare_xemapred_data(run_info)

if (!isTRUE(run_info$run)) {
  cat("[INFO] run = FALSE. Configuration checked, but no model was run.\n")
  q("no")
}

cl <- setup_xemapred_parallel(run_info)
on.exit(stop_xemapred_parallel(cl), add = TRUE)

run_xemapred_cohort_validation(
  data_obj = data_obj,
  run_info = run_info,
  paths = paths
)

combine_xemapred_results(
  data_obj = data_obj,
  run_info = run_info,
  paths = paths
)

cat("[INFO] XemaPred validation complete.\n")