#!/usr/bin/env Rscript

# -----------------------------------------------------------------------
# Scalability plots
# -----------------------------------------------------------------------
# Choose one comparison:
#   xemapred_vs_eczemapred
#   xemapred_vs_patient_specific
#   eczemapred_vs_patient_specific
# -----------------------------------------------------------------------

source(here::here("scripts", "03_plotting", "shared", "00_plot_setup.R"))
source(here::here("scripts", "03_plotting", "shared", "plot_config.R"))
source(here::here("scripts", "03_plotting", "shared", "plot_paths.R"))
source(here::here("scripts", "03_plotting", "shared", "plot_scales.R"))
source(here::here("scripts", "03_plotting", "shared", "plot_themes.R"))
source(here::here("scripts", "03_plotting", "shared", "save_plots.R"))
source(here::here("scripts", "03_plotting", "shared", "scalability_helpers.R"))
source(here::here("scripts", "03_plotting", "shared", "score_metadata.R"))
source(here::here("scripts", "03_plotting", "shared", "observed_data.R"))
source(here::here("scripts", "03_plotting", "shared", "runtime_speedup_helper.R"))
source(here::here("scripts", "03_plotting", "shared", "learning_curve_plot_helpers.R"))
# -----------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------

# COMPARISON <- "all_three_models" 
COMPARISON <- "xemapred_vs_eczemapred"
# COMPARISON <- "xemapred_vs_patient_specific"
# COMPARISON <- "eczemapred_vs_patient_specific"
# COMPARISON <- "all_three_models"

DATASETS <- DATASETS_DEFAULT
RUNTIME_HORIZON <- 1

ITEMS <- c(
  "dryness",
  "extent",
  "itching",
  "redness",
  "oozing",
  "sleep",
  "swelling",
  "thickening",
  "scratching"
)

RESULT_ROOT <- "results"
PLOT_ROOT <- "plots"
PLOT_FORMATS <- c("png", "pdf")

USE_LOG_SCALE <- TRUE
Y_AXIS_MODES <- c("fixed", "free")
SAVE_ITEM_LEVEL <- TRUE

# -----------------------------------------------------------------------
# Run
# -----------------------------------------------------------------------

run_scalability_analysis(
  comparison = COMPARISON,
  datasets = DATASETS,
  items = ITEMS,
  result_root = RESULT_ROOT,
  plot_root = PLOT_ROOT,
  horizon = RUNTIME_HORIZON,
  logy = USE_LOG_SCALE,
  y_axis_modes = Y_AXIS_MODES,
  save_item_level = SAVE_ITEM_LEVEL,
  formats = PLOT_FORMATS
)