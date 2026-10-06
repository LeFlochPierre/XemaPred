# -----------------------------------------------------------------------
# Shared observed data builders
# -----------------------------------------------------------------------

load_dataset_checked <- function(dataset) {
  if (!exists("load_dataset", mode = "function")) {
    stop("load_dataset() was not found. Check 00_plot_setup.R sourcing.")
  }

  load_dataset(dataset)
}

get_available_patients <- function(dataset) {
  load_dataset_checked(dataset) %>%
    pull(.data$Patient) %>%
    unique() %>%
    sort()
}

load_patient_observed_long <- function(dataset, patient_id, items = ITEM_ORDER_WITH_POSCORAD) {
  dat <- load_dataset_checked(dataset) %>%
    filter(.data$Patient == patient_id) %>%
    arrange(.data$Day)

  if (!nrow(dat)) {
    stop("No data found for patient ", patient_id, " in dataset ", dataset)
  }

  purrr::map_dfr(items, function(item) {
    score_col <- get_item_score_column(item)

    if (!score_col %in% names(dat)) {
      warning(
        "Skipping item ", item,
        " because column '", score_col,
        "' was not found in dataset ", dataset
      )
      return(NULL)
    }

    dat %>%
      transmute(
        Patient = .data$Patient,
        Time = .data$Day,
        Item = item,
        Item_label = ITEM_LABELS[[item]],
        Score = as.numeric(.data[[score_col]])
      ) %>%
      filter(!is.na(.data$Score))
  }) %>%
    mutate(
      Item_label = factor(
        .data$Item_label,
        levels = unname(ITEM_LABELS[items])
      )
    ) %>%
    arrange(.data$Item_label, .data$Time)
}

load_observed_items_long <- function(datasets = DATASETS_DEFAULT, items = ITEM_ORDER_DISTRIBUTION) {
  purrr::map_dfr(datasets, function(dataset) {
    dat <- load_dataset_checked(dataset)

    purrr::map_dfr(items, function(item) {
      score_col <- get_item_score_column(item)

      if (!score_col %in% names(dat)) {
        warning(
          "Skipping item ", item,
          " because column '", score_col,
          "' was not found in dataset ", dataset
        )
        return(NULL)
      }

      tibble::tibble(
        Dataset = dataset,
        Dataset_display = dataset_label(dataset),
        Patient = dat$Patient,
        Day = dat$Day,
        Item = item,
        Item_display = ITEM_LABELS[[item]],
        Score = as.numeric(dat[[score_col]])
      ) %>%
        filter(!is.na(.data$Score))
    })
  }) %>%
    mutate(
      Dataset_display = factor(.data$Dataset_display, levels = unname(DATASET_LABELS[datasets])),
      Item_display = factor(.data$Item_display, levels = unname(ITEM_LABELS[items]))
    )
}

load_observed_poscorad_long <- function(datasets = DATASETS_DEFAULT) {
  purrr::map_dfr(datasets, function(dataset) {
    dat <- load_dataset_checked(dataset)

    if (!"SCORAD" %in% names(dat)) {
      stop("Column 'SCORAD' not found for dataset: ", dataset)
    }

    dat %>%
      transmute(
        Dataset = dataset,
        Dataset_display = dataset_label(dataset),
        Patient = .data$Patient,
        Day = .data$Day,
        Score = as.numeric(.data$SCORAD)
      ) %>%
      filter(!is.na(.data$Score))
  }) %>%
    mutate(
      Dataset_display = factor(.data$Dataset_display, levels = unname(DATASET_LABELS[datasets]))
    )
}
