# -----------------------------------------------------------------------
# Shared result-combination utilities
#
# Purpose:
# - Detect which iteration files already exist.
# - Combine per-iteration prediction files into final predictions.rds.
# - Combine per-iteration diagnostic files into final diagnostics.rds.
# -----------------------------------------------------------------------

get_existing_iterations <- function(paths, require_posterior = FALSE) {
  pred_files <- list.files(
    paths$iters_dir,
    pattern = "^iter-[0-9]{4}\\.rds$",
    full.names = FALSE
  )

  diag_files <- list.files(
    paths$diag_dir,
    pattern = "^diag-[0-9]{4}\\.rds$",
    full.names = FALSE
  )

  pred_its <- as.integer(gsub("iter-([0-9]{4})\\.rds", "\\1", pred_files))
  diag_its <- as.integer(gsub("diag-([0-9]{4})\\.rds", "\\1", diag_files))

  existing_its <- intersect(pred_its, diag_its)

  if (isTRUE(require_posterior)) {
    post_files <- list.files(
      paths$posterior_dir,
      pattern = "^posterior-[0-9]{4}\\.rds$",
      full.names = FALSE
    )

    post_its <- as.integer(
      gsub("posterior-([0-9]{4})\\.rds", "\\1", post_files)
    )

    existing_its <- intersect(existing_its, post_its)
  }

  sort(existing_its)
}

get_missing_iterations <- function(
    data_obj,
    paths,
    require_posterior = FALSE
) {
  existing_its <- get_existing_iterations(
    paths,
    require_posterior = require_posterior
  )

  setdiff(data_obj$train_it, existing_its)
}

combine_iteration_results <- function(paths, expected_iterations = NULL) {
  val_files <- list.files(
    paths$iters_dir,
    pattern = "^iter-[0-9]{4}\\.rds$",
    full.names = TRUE
  )

  if (length(val_files) == 0) {
    warning("[WARNING] No iteration files found in: ", paths$iters_dir)
    return(invisible(NULL))
  }

  if (!is.null(expected_iterations) &&
      length(val_files) < length(expected_iterations)) {
    warning(glue::glue(
      "[WARNING] Missing results: expected {length(expected_iterations)} iterations, ",
      "found {length(val_files)}."
    ))
  }

  res <- lapply(sort(val_files), readRDS) %>%
    dplyr::bind_rows()

  saveRDS(res, paths$predictions_out)

  cat(
    "[INFO] Combined predictions saved to: ",
    paths$predictions_out,
    "\n",
    sep = ""
  )

  invisible(res)
}

combine_diagnostics <- function(paths) {
  diag_files <- list.files(
    paths$diag_dir,
    pattern = "^diag-[0-9]{4}\\.rds$",
    full.names = TRUE
  )

  if (length(diag_files) == 0) {
    warning("[WARNING] No diagnostic files found in: ", paths$diag_dir)
    return(invisible(NULL))
  }

  diag <- lapply(sort(diag_files), readRDS) %>%
    dplyr::bind_rows()

  saveRDS(diag, paths$diagnostics_out)

  cat(
    "[INFO] Combined diagnostics saved to: ",
    paths$diagnostics_out,
    "\n",
    sep = ""
  )

  invisible(diag)
}