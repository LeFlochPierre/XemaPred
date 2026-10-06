# -----------------------------------------------------------------------
# Combine patient-specific outputs
# -----------------------------------------------------------------------

combine_patient_specific_results <- function(res_list, paths) {
  res_list <- res_list[!vapply(res_list, is.null, logical(1))]

  read_many_rds <- function(files) {
    objs <- purrr::map(files, safe_read_rds)
    objs <- objs[!vapply(objs, is.null, logical(1))]

    if (length(objs) == 0L) {
      return(tibble::tibble())
    }

    dplyr::bind_rows(objs)
  }

  pred_files <- character(0)
  diag_files <- character(0)

  if (!is.null(paths$iters_dir) && dir.exists(paths$iters_dir)) {
    pred_files <- list.files(
      paths$iters_dir,
      pattern = "^iter-[0-9]{4}\\.rds$",
      recursive = TRUE,
      full.names = TRUE
    )
  }

  if (!is.null(paths$diag_dir) && dir.exists(paths$diag_dir)) {
    diag_files <- list.files(
      paths$diag_dir,
      pattern = "^diag-[0-9]{4}\\.rds$",
      recursive = TRUE,
      full.names = TRUE
    )
  }

  if (length(pred_files) > 0L) {
    pred_all <- read_many_rds(pred_files)
  } else {
    pred_all <- purrr::map(res_list, "pred") %>% dplyr::bind_rows()
  }

  if (length(diag_files) > 0L) {
    diag_all <- read_many_rds(diag_files)
  } else {
    diag_all <- purrr::map(res_list, "diag") %>% dplyr::bind_rows()
  }

  saveRDS(pred_all, file = file.path(paths$final_dir, "predictions_ALL_PATIENTS.rds"))
  saveRDS(diag_all, file = file.path(paths$final_dir, "diagnostics_ALL_PATIENTS.rds"))

  cat("[INFO] Wrote ALL_PATIENTS finals.\n")
  cat("[INFO] Predictions rows: ", nrow(pred_all), "\n", sep = "")
  cat("[INFO] Diagnostics rows: ", nrow(diag_all), "\n", sep = "")

  invisible(list(pred = pred_all, diag = diag_all))
}

collect_all_patients_from_disk <- function(paths, filename, out_name, tries = 30L, sleep_sec = 2) {
  for (k in seq_len(tries)) {
    files <- list.files(paths$final_dir, pattern = paste0("^", filename, "$"), recursive = TRUE, full.names = TRUE)

    if (length(files) > 0L) {
      objs <- purrr::map(files, safe_read_rds)
      objs <- objs[!vapply(objs, is.null, logical(1))]

      if (length(objs) > 0L) {
        df <- dplyr::bind_rows(objs)
        out <- file.path(paths$final_dir, out_name)
        saveRDS(df, out)
        message("[INFO] Saved combined ", filename, " -> ", out, " (n_files=", length(objs), ")")
        return(invisible(out))
      }
    }

    message("[INFO] Collector waiting for ", filename, " (attempt ", k, "/", tries, ") ...")
    Sys.sleep(sleep_sec)
  }

  message("[WARN] Could not create ", out_name, " (no readable inputs).")
  invisible(NULL)
}
