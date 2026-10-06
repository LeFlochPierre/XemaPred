# -----------------------------------------------------------------------
# XemaPred project setup
#
# Purpose:
# - Load the core packages used across model-running scripts.
# - Define small global helpers that are needed by several workflows.
# - Configure project-root detection and Stan runtime options.
# -----------------------------------------------------------------------

suppressPackageStartupMessages({
  library(EczemaPredPOSCORAD)
  library(EczemaPred)
  library(HuraultMisc)
  library(TanakaData)

  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(glue)
  library(here)

  library(rstan)
})

# -----------------------------------------------------------------------
# Project root
# -----------------------------------------------------------------------

# Make here::here() use the XemaPred project root.
# This assumes this file is located at scripts/00_setup/00_init.R.
try(
  here::i_am("scripts/00_setup/00_init.R"),
  silent = TRUE
)

# -----------------------------------------------------------------------
# Stan options
# -----------------------------------------------------------------------

# Avoid accidental oversubscription on HPC.
# The model runners already parallelise over iterations using n_cluster, and
# each Stan fit can also use multiple chains. Setting mc.cores to all detected
# cores inside every worker can therefore overload the node.
stan_cores <- as.integer(Sys.getenv("R_STAN_CORES", unset = "1"))
options(mc.cores = stan_cores)

# Save compiled Stan models where possible.
# This is useful when repeatedly fitting the same EczemaPred models.
rstan::rstan_options(auto_write = TRUE)

# -----------------------------------------------------------------------
# Data loading
# -----------------------------------------------------------------------

#' Load Derexyl, PFDC, or fake POSCORAD data.
#'
#' Important:
#' - Removes patients with fewer than five observations.
#' - Regenerates patient IDs after filtering.
#' - Keep this behaviour unchanged if you want consistency with old runs.
#'
#' @param dataset One of "Derexyl", "PFDC", or "Fake".
#'
#' @return POSCORAD data frame.
load_dataset <- function(dataset = c("Derexyl", "PFDC", "Fake")) {
  dataset <- match.arg(dataset)

  if (dataset == "Derexyl") {
    out <- TanakaData::POSCORAD_Derexyl
  }

  if (dataset == "PFDC") {
    out <- TanakaData::POSCORAD_PFDC
  }

  if (dataset == "Fake") {
    out <- EczemaPredPOSCORAD::FakeData$Data
  }

  out <- out %>%
    dplyr::group_by(Patient) %>%
    dplyr::filter(dplyr::n() >= 5) %>%
    dplyr::mutate(Patient = dplyr::cur_group_id()) %>%
    dplyr::ungroup()

  out
}

# -----------------------------------------------------------------------
# Dataset metadata
# -----------------------------------------------------------------------

# Canonical dataset cutoffs.
# These values should match the EczemaPred/XemaPred processing scripts.
dict_datasets <- tibble::tibble(
  Dataset = c("PFDC", "Derexyl"),
  Max_train_day = c(80L, 115L)
)

# -----------------------------------------------------------------------
# Logging helper
# -----------------------------------------------------------------------

# Print one line and flush the console.
# Useful for PBS jobs, where logs can otherwise appear with a delay.
log_line <- function(..., sep = "", .flush = TRUE) {
  cat(..., "\n", sep = sep)

  if (.flush) {
    flush.console()
  }

  invisible(NULL)
}
