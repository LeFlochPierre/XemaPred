#!/usr/bin/env Rscript
# -----------------------------------------------------------------------
# Runtime speed-up analysis runner
# -----------------------------------------------------------------------
# Runs:
#   1) EczemaPred vs XemaPred
#   2) population-level XemaPred vs patient-specific XemaPred
#   3) EczemaPred vs patient-specific XemaPred
#
# Output layout for each comparison:
#   plots/runtime_speedup/<comparison>/
#     combined/   -> combined Dataset 1 + Dataset 2 plots
#     Derexyl/    -> Dataset 1 plots
#     PFDC/       -> Dataset 2 plots
#     csv/        -> analysis tables
#
# For trajectory plots, the x-axis is the actual accumulated number of
# training observations. Per-dataset plots additionally show training day
# on the upper x-axis. Combined plots omit the upper day axis because the
# observation-count-to-day mapping differs between datasets.
# -----------------------------------------------------------------------

suppressPackageStartupMessages({
  library(here)
  library(purrr)
})

# -----------------------------------------------------------------------
# Shared setup
# -----------------------------------------------------------------------

source(here::here("scripts", "03_plotting", "shared", "00_plot_setup.R"))
source(here::here("scripts", "03_plotting", "shared", "plot_config.R"))
source(here::here("scripts", "03_plotting", "shared", "plot_paths.R"))
source(here::here("scripts", "03_plotting", "shared", "plot_themes.R"))
source(here::here("scripts", "03_plotting", "shared", "plot_scales.R"))
source(here::here("scripts", "03_plotting", "shared", "score_metadata.R"))
source(here::here("scripts", "03_plotting", "shared", "observed_data.R"))
source(here::here("scripts", "03_plotting", "shared", "save_plots.R"))

# Runtime/scalability loaders.
source(here::here("scripts", "03_plotting", "shared", "scalability_helpers.R"))

# Runtime speed-up helper.
source(here::here("scripts", "03_plotting", "shared", "runtime_speedup_helper.R"))

# for the corretc x axis scales in the learning curve plots
source(here::here("scripts", "03_plotting", "shared", "learning_curve_plot_helpers.R"))

# -----------------------------------------------------------------------
# Arguments
# -----------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

comparisons <- c(
  "xemapred_vs_eczemapred",
  "xemapred_vs_patient_specific",
  "eczemapred_vs_patient_specific"
)

if (length(args) >= 1 && nzchar(args[[1]])) {
  comparisons <- strsplit(gsub("\\s+", "", args[[1]]), ",")[[1]]
}

items <- if (exists("ITEM_ORDER_MAIN", inherits = TRUE)) {
  ITEM_ORDER_MAIN
} else {
  rs_default_items()
}

datasets <- if (exists("DATASETS_DEFAULT", inherits = TRUE)) {
  DATASETS_DEFAULT
} else {
  c("Derexyl", "PFDC")
}

# -----------------------------------------------------------------------
# Run analyses
# -----------------------------------------------------------------------

results <- purrr::map(
  comparisons,
  ~ run_runtime_speedup_analysis(
    comparison = .x,
    datasets = datasets,
    items = items,
    result_root = "results",
    plot_root = "plots",
    horizon_for_counts = 1L,
    smooth_k = 6L,
    formats = c("png", "pdf")
  )
)

names(results) <- comparisons

cat("\n============================================================\n")
cat("[ALL RUNTIME SPEED-UP ANALYSES DONE]\n")
cat("============================================================\n")

print(vapply(results, `[[`, character(1), "out_dir"))
