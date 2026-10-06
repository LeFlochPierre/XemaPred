# -----------------------------------------------------------------------
# Source the model files required by population-level XemaPred
#
# Notes:
# - Base files contain priors, likelihoods, parameter packing/unpacking,
#   and model-specific scoring helpers.
# - SMC2 files contain inner PFs, outer particles, PMMH moves, updates,
#   and forecast functions.
# -----------------------------------------------------------------------

source_xemapred_model_files <- function(model_family = NULL) {
  root <- here::here("scripts", "01_run_models", "XemaPred")

  families <- if (is.null(model_family)) {
    c("extent", "intensity", "subjective")
  } else {
    match.arg(model_family, c("extent", "intensity", "subjective"))
  }

  model_files <- unlist(lapply(families, function(fam) {
    file.path(root, "models", fam, paste0(fam, c("_base.R", "_smc2.R")))
  }))

  cat("[INFO] Sourcing XemaPred model files (", paste(families, collapse = ", "), ")...\n", sep = "")

  for (f in model_files) {
    if (!file.exists(f)) stop("[ERROR] Missing model source file: ", f)
    source(f)
  }

  cat("[INFO] XemaPred model files sourced.\n")

  invisible(model_files)
}