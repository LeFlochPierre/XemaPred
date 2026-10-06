# -----------------------------------------------------------------------
# Source the model files required by PatientSpecificXemaPred
# -----------------------------------------------------------------------

source_patient_specific_model_files <- function(model_family = c("extent", "intensity", "subjective")) {
  model_family <- match.arg(model_family)

  root <- here::here("scripts", "01_run_models", "PatientSpecificXemaPred")

  model_files <- switch(
    model_family,
    extent = c(
      file.path(root, "models", "extent", "extent_base.R"),
      file.path(root, "models", "extent", "extent_cache.R")
    ),
    intensity = c(
      file.path(root, "models", "intensity", "intensity_base.R"),
      file.path(root, "models", "intensity", "intensity_cache.R")
    ),
    subjective = c(
      file.path(root, "models", "subjective", "subjective_base.R"),
      file.path(root, "models", "subjective", "subjective_cache.R")
    )
  )

  for (f in model_files) {
    if (!file.exists(f)) stop("[ERROR] Missing model source file: ", f)
    source(f)
  }

  invisible(model_files)
}
