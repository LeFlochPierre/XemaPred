# -----------------------------------------------------------------------
# Shared prediction data builders
# -----------------------------------------------------------------------

summarise_samples <- function(samples_list) {
  tibble(
    pred_median = vapply(samples_list, median, numeric(1), na.rm = TRUE),
    q05 = vapply(samples_list, quantile, numeric(1), probs = 0.05, na.rm = TRUE),
    q10 = vapply(samples_list, quantile, numeric(1), probs = 0.10, na.rm = TRUE),
    q25 = vapply(samples_list, quantile, numeric(1), probs = 0.25, na.rm = TRUE),
    q75 = vapply(samples_list, quantile, numeric(1), probs = 0.75, na.rm = TRUE),
    q90 = vapply(samples_list, quantile, numeric(1), probs = 0.90, na.rm = TRUE),
    q95 = vapply(samples_list, quantile, numeric(1), probs = 0.95, na.rm = TRUE)
  )
}

has_samples <- function(x) {
  !is.null(x) && length(x) > 0 && !all(is.na(x))
}

load_patient_item_predictions <- function(
    dataset,
    patient_id,
    item,
    horizon,
    model_type,
    result_root = "results"
) {
  f <- get_prediction_path(
    dataset = dataset,
    item = item,
    horizon = horizon,
    model_type = model_type,
    result_root = result_root
  )

  pred <- readRDS(f)

  if (!"Patient" %in% names(pred)) {
    pred$Patient <- patient_id
  }

  if (!"Iteration" %in% names(pred)) {
    pred$Iteration <- NA_integer_
  }

  pred <- pred %>%
    mutate(
      Patient = as.integer(as.character(.data$Patient))
    ) %>%
    filter(
      .data$Patient == patient_id,
      .data$Horizon <= horizon
    ) %>%
    arrange(.data$Time) %>%
    group_by(.data$Patient, .data$Time, .data$Horizon) %>%
    slice_tail(n = 1) %>%
    ungroup() %>%
    select(
      .data$Patient,
      .data$Time,
      .data$Horizon,
      .data$Iteration,
      .data$Samples
    )

  observed <- load_patient_observed_long(
    dataset = dataset,
    patient_id = patient_id,
    items = item
  ) %>%
    transmute(
      Patient = as.integer(as.character(.data$Patient)),
      Time = .data$Time,
      Score = .data$Score
    )

  pred %>%
    left_join(observed, by = c("Patient", "Time")) %>%
    filter(
      !is.na(.data$Score),
      purrr::map_lgl(.data$Samples, has_samples)
    ) %>%
    bind_cols(summarise_samples(.$Samples)) %>%
    mutate(
      Item = item,
      Item_label = ITEM_LABELS[[item]]
    )
}

load_patient_all_item_predictions <- function(
    dataset,
    patient_id,
    horizon = 4,
    model_type = "XemaPred",
    result_root = "results"
) {
  item_dfs <- purrr::map_dfr(ITEM_ORDER_MAIN, function(item) {
    load_patient_item_predictions(
      dataset = dataset,
      patient_id = patient_id,
      item = item,
      horizon = horizon,
      model_type = model_type,
      result_root = result_root
    )
  })

  required_items <- ITEM_ORDER_MAIN

  wide_samples <- item_dfs %>%
    select(.data$Patient, .data$Time, .data$Horizon, .data$Iteration, .data$Item, .data$Samples) %>%
    tidyr::pivot_wider(names_from = .data$Item, values_from = .data$Samples)

  true_sc <- load_dataset_checked(dataset) %>%
    filter(.data$Patient == patient_id) %>%
    transmute(
      Patient = .data$Patient,
      Time = .data$Day,
      Score = as.numeric(.data$SCORAD)
    )

  sc_df <- wide_samples %>%
    left_join(true_sc, by = c("Patient", "Time"))

  for (item in required_items) {
    if (!item %in% names(sc_df)) {
      stop("Cannot reconstruct PO-SCORAD because item is missing from predictions: ", item)
    }
  }

  sc_df <- sc_df %>%
    filter(
      purrr::pmap_lgl(
        select(., all_of(required_items)),
        function(...) all(purrr::map_lgl(list(...), has_samples))
      )
    ) %>%
    mutate(
      Samples = purrr::pmap(
        list(
          .data$extent,
          .data$redness,
          .data$dryness,
          .data$swelling,
          .data$oozing,
          .data$scratching,
          .data$thickening,
          .data$itching,
          .data$sleep
        ),
        function(extent, redness, dryness, swelling, oozing, scratching, thickening, itching, sleep) {
          B <- redness + dryness + swelling + oozing + scratching + thickening
          C <- itching + sleep
          0.2 * extent + 3.5 * B + C
        }
      )
    ) %>%
    select(.data$Patient, .data$Time, .data$Horizon, .data$Iteration, .data$Score, .data$Samples) %>%
    bind_cols(summarise_samples(.$Samples)) %>%
    mutate(
      Item = "SCORAD",
      Item_label = ITEM_LABELS[["SCORAD"]]
    )

  bind_rows(item_dfs, sc_df) %>%
    mutate(
      Item_label = factor(
        .data$Item_label,
        levels = unname(ITEM_LABELS[ITEM_ORDER_WITH_POSCORAD])
      )
    )
}
