#!/usr/bin/env Rscript
# -----------------------------------------------------------------------
# Main runner for patient-specific XemaPred item models
#
# Purpose:
# - Run patient-specific XemaPred validation for extent, intensity signs,
#   and subjective symptoms.
# - Resume from per-patient cached PF states and existing iteration files.
# - Combine predictions and diagnostics at the end.
#
# Run example:
# Rscript scripts/01_run_models/PatientSpecificXemaPred/run_validation.R \
#   configs/PFDC/PatientSpecificXemaPred/t_horizon_1/pfdc_dryness_orderedrw.yaml
# -----------------------------------------------------------------------

rm(list = ls())

cat("[INFO] Starting patient-specific XemaPred validation script...\n")

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

root <- here::here("scripts", "01_run_models", "PatientSpecificXemaPred")

source(file.path(root, "shared", "utils.R"))
source(file.path(root, "shared", "particle_utils.R"))
source(file.path(root, "shared", "logging.R"))
source(file.path(root, "shared", "score_info.R"))
source(file.path(root, "shared", "config.R"))
source(file.path(root, "shared", "paths.R"))
source(file.path(root, "shared", "data_prep.R"))
source(file.path(root, "shared", "parallel.R"))
source(file.path(root, "shared", "population_prior.R"))

source(file.path(root, "source_model_files.R"))
source(file.path(root, "build_observations.R"))
source(file.path(root, "state_cache.R"))
source(file.path(root, "update_state.R"))
source(file.path(root, "forecast_score.R"))
source(file.path(root, "combine_results.R"))
source(file.path(root, "run_cohort_validation.R"))

cfg <- read_run_config()
run_info <- validate_run_config(cfg, expected_family = "patient_specific_xemapred")

RNGkind("L'Ecuyer-CMRG")
set.seed(run_info$seed)
cat(glue("[INFO] RNG seed set to {run_info$seed}\n"))

print_run_config(run_info)

paths <- make_output_paths(run_info)
source_patient_specific_model_files(run_info$model_family)

run_info$population_prior <-
  load_population_prior_for_run(
    run_info
  )

save_run_metadata(cfg, run_info, paths)

data_obj <- prepare_patient_specific_data(run_info)

if (!isTRUE(run_info$run)) {
  cat("[INFO] run = FALSE. Configuration checked, but no model was run.\n")
  q("no")
}

res_list <- run_patients_parallel(data_obj, run_info, paths)
combine_patient_specific_results(res_list, paths)

cat("[INFO] Patient-specific XemaPred validation complete.\n")
