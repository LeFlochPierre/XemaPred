# -----------------------------------------------------------------------
# XemaPred data-preparation utilities
#
# Purpose:
# - Load the requested PO-SCORAD dataset.
# - Extract one score into a standard Patient/Time/Score format.
# - Add absolute and relative time columns.
# - Build the forward-chaining iteration schedule.
#
# Notes:
# - TimeAbs is used for forward-chaining splits.
# - TimeRel is used for particle-filter dynamics.
# -----------------------------------------------------------------------

prepare_xemapred_data <- function(run_info) {
  POSCORAD <- load_dataset(run_info$dataset)

  cat("[INFO] Dataset loaded. Rows: ", nrow(POSCORAD), "\n", sep = "")

  df_all <- POSCORAD |>
    dplyr::rename(
      Time = dplyr::all_of("Day"),
      Score = dplyr::all_of(run_info$item_label)
    ) |>
    dplyr::select(
      dplyr::all_of(c("Patient", "Time", "Score"))
    ) |>
    tidyr::drop_na()

  if (nrow(df_all) == 0) {
    stop("[ERROR] No data after drop_na().")
  }

  # -----------------------------------------------------------------------
  # Patient-fold selection
  # -----------------------------------------------------------------------

  all_patients <- sort(
    unique(df_all$Patient)
  )

  n_patients_total <- length(all_patients)

  cat(
    "[INFO] Patients found before fold selection: ",
    n_patients_total,
    "\n",
    sep = ""
  )

  if (run_info$n_patient_folds > n_patients_total) {
    stop(
      "[ERROR] n_patient_folds (",
      run_info$n_patient_folds,
      ") exceeds the number of available patients (",
      n_patients_total,
      ")."
    )
  }

  if (run_info$n_patient_folds > 1L) {

    # Deterministically divide sorted patient IDs into approximately
    # equal, non-overlapping folds.
    #
    # For the datasets here:
    #   PFDC    : 16 patients -> 8 / 8
    #   Derexyl : 336 patients -> 168 / 168

    fold_assignment <- cut(
      seq_along(all_patients),
      breaks = run_info$n_patient_folds,
      labels = FALSE
    )

    keep_patients <- all_patients[
      fold_assignment == run_info$patient_fold
    ]

    cat(
      "[INFO] Patient-fold selection:\n",
      "  Number of folds: ",
      run_info$n_patient_folds,
      "\n",
      "  Selected fold: ",
      run_info$patient_fold,
      "\n",
      "  Patients retained: ",
      length(keep_patients),
      " / ",
      n_patients_total,
      "\n",
      "  Patient IDs: ",
      paste(keep_patients, collapse = ", "),
      "\n",
      sep = ""
    )

    df_all <- df_all |>
      dplyr::filter(
        .data$Patient %in% keep_patients
      )
  }

  patients <- sort(
    unique(df_all$Patient)
  )

  cat(
    "[INFO] Patients used by XemaPred: ",
    length(patients),
    "\n",
    sep = ""
  )

  t0_global <- min(df_all$Time)

  df_all2 <- df_all |>
    dplyr::mutate(
      TimeAbs = as.integer(.data$Time),
      TimeRel = as.integer(.data$TimeAbs - t0_global)
    ) |>
    dplyr::group_by(.data$Patient) |>
    dplyr::arrange(.data$TimeAbs, .by_group = TRUE) |>
    dplyr::mutate(
      Time = .data$TimeAbs,
      Iteration = get_fc_iteration(.data$Time, run_info$t_horizon)
    ) |>
    dplyr::ungroup()

  df_by_pid <- split(df_all2, df_all2$Patient)

  iters_all <- df_all2 |>
    dplyr::group_by(.data$Patient) |>
    dplyr::summarise(
      train_it = list(get_fc_training_iteration(.data$Iteration)),
      .groups = "drop"
    ) |>
    dplyr::pull(.data$train_it) |>
    unlist() |>
    unique() |>
    sort()

  if (!length(iters_all)) {
    stop("[ERROR] No forward-chaining iterations found.")
  }

  cat(
    "[INFO] Time range: ",
    min(df_all2$TimeAbs),
    " to ",
    max(df_all2$TimeAbs),
    "\n",
    sep = ""
  )

  cat(
    "[INFO] Forward-chaining iterations: ",
    length(iters_all),
    "\n",
    sep = ""
  )

  list(
    POSCORAD = POSCORAD,
    df_all = df_all,
    df_all2 = df_all2,
    df_by_pid = df_by_pid,
    patients = patients,
    t0_global = t0_global,
    iters_all = iters_all
  )
}

split_xemapred_patient_iteration <- function(data_obj, pid, it) {
  dfp <- data_obj$df_by_pid[[as.character(pid)]]

  if (is.null(dfp)) {
    stop("[ERROR] Unknown patient ID: ", pid)
  }

  split <- split_fc_dataset(dfp, it)

  list(
    train = split$Training |>
      dplyr::mutate(TimeRel = as.integer(.data$Time - data_obj$t0_global)),
    test = split$Testing |>
      dplyr::mutate(TimeRel = as.integer(.data$Time - data_obj$t0_global))
  )
}