# -----------------------------------------------------------------------
# Shared data-preparation utilities
#
# Purpose:
# - Load the requested POSCORAD dataset.
# - Select the requested score and convert it to the common Time/Score format.
# - Define forward-chaining iterations.
# - Provide helpers for score scaling and debug-one-patient mode.
# -----------------------------------------------------------------------

prepare_forward_chaining_data <- function(run_info) {
  POSCORAD <- load_dataset(run_info$dataset)
  cat("[INFO] Dataset loaded. Rows:", nrow(POSCORAD), "\n")

  df <- POSCORAD %>%
    dplyr::rename(
      Time = Day,
      Score = dplyr::all_of(run_info$item_lbl)
    ) %>%
    dplyr::select(Patient, Time, Score) %>%
    tidyr::drop_na()

  df <- df %>%
    dplyr::mutate(Iteration = get_fc_iteration(.data$Time, run_info$t_horizon))

  train_it <- get_fc_training_iteration(df[["Iteration"]])

  cadence <- if (run_info$t_horizon == 1) {
    "every day"
  } else {
    paste("every", run_info$t_horizon, "days")
  }

  cat(glue::glue("[INFO] Cadence: {cadence}\n"))
  cat("[INFO] Number of forward-chaining iterations:", length(train_it), "\n")

  list(
    POSCORAD = POSCORAD,
    df = df,
    train_it = train_it
  )
}

split_iteration_data <- function(data_obj, it) {
  split_fc_dataset(data_obj$df, it)
}

rescale_for_model <- function(dat, run_info) {
  if (run_info$is_continuous) {
    return(dat)
  }

  dat %>%
    dplyr::mutate(Score = round(.data$Score / run_info$reso))
}

restore_score_scale <- function(perf, run_info) {
  if (run_info$is_continuous) {
    return(perf)
  }

  if ("Samples" %in% names(perf)) {
    perf <- perf %>%
      dplyr::mutate(
        Samples = purrr::map(.data$Samples, ~ .x * run_info$reso)
      )
  }

  perf
}

remove_last_observation_columns <- function(perf) {
  cols_to_remove <- intersect(c("LastTime", "LastScore"), names(perf))
  perf %>% dplyr::select(-dplyr::all_of(cols_to_remove))
}

apply_debug_one_patient <- function(split, run_info) {
  if (!isTRUE(run_info$debug_one_patient)) {
    return(split)
  }

  train <- split$Training
  test <- split$Testing

  intersecting_patients <- intersect(unique(train$Patient), unique(test$Patient))

  if (length(intersecting_patients) == 0) {
    stop("[ERROR] No patient has both valid training and test data for debug mode.")
  }

  debug_patient <- intersecting_patients[[1]]

  cat("[DEBUG] Running on patient:", debug_patient, "\n")

  split$Training <- train %>%
    dplyr::filter(.data$Patient == debug_patient) %>%
    tidyr::drop_na(Score)

  split$Testing <- test %>%
    dplyr::filter(.data$Patient == debug_patient) %>%
    tidyr::drop_na(Score)

  if (nrow(split$Training) == 0 || nrow(split$Testing) == 0) {
    stop("[ERROR] Selected debug patient has no valid training or test data.")
  }

  split
}
