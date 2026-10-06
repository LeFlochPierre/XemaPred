# -----------------------------------------------------------------------
# Data preparation for patient-specific forward chaining
# -----------------------------------------------------------------------

prepare_patient_specific_data <- function(run_info) {
  POSCORAD <- load_dataset(run_info$dataset)
  cat("[INFO] Dataset loaded. Rows:", nrow(POSCORAD), "\n")

  df_all <- POSCORAD %>%
    dplyr::rename(Time = Day, Score = dplyr::all_of(run_info$item_label)) %>%
    dplyr::select(.data$Patient, .data$Time, .data$Score) %>%
    tidyr::drop_na()

  if (nrow(df_all) == 0) stop("[ERROR] No data after drop_na().")

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
      "[ERROR] n_patient_folds exceeds the number of available patients."
    )
  }

  if (run_info$n_patient_folds > 1L) {

    # IMPORTANT:
    # Same deterministic fold definition as population-level XemaPred.
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
      "  Target fold: ",
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

    df_all <- df_all %>%
      dplyr::filter(
        .data$Patient %in% keep_patients
      )
  }

  patients <- sort(
    unique(df_all$Patient)
  )

  cat(
    "[INFO] Patients used by PatientSpecificXemaPred: ",
    length(patients),
    "\n",
    sep = ""
  )

  list(
    df_all = df_all,
    patients = patients,
    model_family = run_info$model_family
  )
}

prepare_patient_dataframe <- function(df_all, pid, t_horizon) {
  df <- df_all %>%
    dplyr::filter(.data$Patient == pid) %>%
    dplyr::arrange(.data$Time) %>%
    dplyr::mutate(Iteration = get_fc_iteration(.data$Time, t_horizon))

  if (nrow(df) == 0L) return(df)
  df
}
