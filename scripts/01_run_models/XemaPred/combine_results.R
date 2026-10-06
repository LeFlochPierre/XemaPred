# -----------------------------------------------------------------------
# Combine population-level XemaPred outputs
#
# Purpose:
# - Combine EczemaPred-style iteration-level files into final outputs.
# - Save canonical final/predictions.rds and final/diagnostics.rds.
# - Also save *_ALL_PATIENTS aliases for backward compatibility.
# -----------------------------------------------------------------------

get_existing_xemapred_iterations <- function(paths) {
  existing_files <- list.files(
    paths$iters_dir,
    pattern = "^iter-[0-9]{4}\\.rds$",
    full.names = FALSE
  )

  as.integer(gsub("iter-([0-9]{4})\\.rds", "\\1", existing_files))
}

combine_xemapred_iteration_results <- function(paths, expected_iterations = NULL) {
  val_files <- list.files(
    paths$iters_dir,
    pattern = "^iter-[0-9]{4}\\.rds$",
    full.names = TRUE
  )

  if (length(val_files) == 0) {
    warning("[WARNING] No XemaPred iteration files found in: ", paths$iters_dir)
    return(invisible(NULL))
  }

  if (!is.null(expected_iterations) && length(val_files) < length(expected_iterations)) {
    warning(glue::glue(
      "[WARNING] Missing XemaPred results: expected {length(expected_iterations)} iterations, found {length(val_files)}."
    ))
  }

  res <- purrr::map(val_files, safe_read_rds) |>
    purrr::compact() |>
    dplyr::bind_rows()

  saveRDS(res, file = paths$predictions_out)

  # Backward-compatible alias for older plotting scripts.
  saveRDS(res, file = paths$predictions_all)

  cat("[INFO] Combined XemaPred predictions saved to: ", paths$predictions_out, "\n", sep = "")
  cat("[INFO] Backward-compatible predictions alias saved to: ", paths$predictions_all, "\n", sep = "")

  invisible(res)
}

combine_xemapred_diagnostics <- function(paths, expected_iterations = NULL) {
  diag_files <- list.files(
    paths$diag_dir,
    pattern = "^diag-[0-9]{4}\\.rds$",
    full.names = TRUE
  )

  if (length(diag_files) == 0) {
    warning("[WARNING] No XemaPred diagnostic files found in: ", paths$diag_dir)
    return(invisible(NULL))
  }

  if (!is.null(expected_iterations) && length(diag_files) < length(expected_iterations)) {
    warning(glue::glue(
      "[WARNING] Missing XemaPred diagnostics: expected {length(expected_iterations)} iterations, found {length(diag_files)}."
    ))
  }

  diag <- purrr::map(diag_files, safe_read_rds) |>
    purrr::compact() |>
    dplyr::bind_rows()

  saveRDS(diag, file = paths$diagnostics_out)

  # Backward-compatible alias for older plotting scripts.
  saveRDS(diag, file = paths$diagnostics_all)

  cat("[INFO] Combined XemaPred diagnostics saved to: ", paths$diagnostics_out, "\n", sep = "")
  cat("[INFO] Backward-compatible diagnostics alias saved to: ", paths$diagnostics_all, "\n", sep = "")

  invisible(diag)
}

combine_xemapred_results <- function(data_obj, run_info, paths) {
  cat("[INFO] Combining XemaPred iteration-level outputs...\n")

  combine_xemapred_iteration_results(
    paths = paths,
    expected_iterations = data_obj$iters_all
  )

  combine_xemapred_diagnostics(
    paths = paths,
    expected_iterations = data_obj$iters_all
  )

  cat("[INFO] XemaPred result combination complete.\n")

  invisible(TRUE)
}
