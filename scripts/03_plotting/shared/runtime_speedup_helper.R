# -----------------------------------------------------------------------
# Shared runtime speed-up analysis helpers
# -----------------------------------------------------------------------
# Purpose:
#   Build and plot item-level runtime reduction factors for:
#     1) EczemaPred vs XemaPred
#     2) population-level XemaPred vs patient-specific XemaPred
#
# Main entry point:
#   run_runtime_speedup_analysis(comparison = "xemapred_vs_eczemapred")
#   run_runtime_speedup_analysis(comparison = "xemapred_vs_patient_specific")
# -----------------------------------------------------------------------

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(tibble)
  library(scales)
  library(patchwork)
  library(readr)
})

# -----------------------------------------------------------------------
# Defaults from shared plotting config, with safe fallbacks
# -----------------------------------------------------------------------

rs_default_items <- function() {
  if (exists("ITEM_ORDER_MAIN", inherits = TRUE)) {
    return(get("ITEM_ORDER_MAIN", inherits = TRUE))
  }

  if (exists("sc_default_items", mode = "function", inherits = TRUE)) {
    return(sc_default_items())
  }

  c(
    "extent", "redness", "dryness", "swelling", "oozing",
    "scratching", "thickening", "itching", "sleep"
  )
}

rs_dataset_levels <- function() {
  if (exists("DATASETS_DEFAULT", inherits = TRUE)) {
    return(get("DATASETS_DEFAULT", inherits = TRUE))
  }

  c("Derexyl", "PFDC")
}

rs_dataset_labels <- function() {
  if (exists("DATASET_LABELS", inherits = TRUE)) {
    return(get("DATASET_LABELS", inherits = TRUE))
  }

  if (exists("sc_dataset_labels", mode = "function", inherits = TRUE)) {
    return(sc_dataset_labels())
  }

  c(Derexyl = "Dataset 1", PFDC = "Dataset 2")
}

rs_dataset_cols <- function() {
  if (exists("DATASET_COLS", inherits = TRUE)) {
    return(get("DATASET_COLS", inherits = TRUE))
  }
  c("Dataset 1" = "#a594ff", "Dataset 2" = "#61ce93")
}

rs_item_labels <- function() {
  if (exists("ITEM_LABELS", inherits = TRUE)) {
    return(get("ITEM_LABELS", inherits = TRUE))
  }

  if (exists("sc_item_labels", mode = "function", inherits = TRUE)) {
    return(sc_item_labels())
  }

  c(
    extent = "Extent",
    redness = "Redness",
    dryness = "Dryness",
    swelling = "Swelling",
    oozing = "Oozing",
    scratching = "Scratching",
    thickening = "Thickening",
    itching = "Itch",
    sleep = "Sleep"
  )
}

rs_model_cols <- function() {
  c(
    "EczemaPred" = "#0072B2",
    "XemaPred" = "#d66309",
    "Patient-specific XemaPred" = "#556B2F"
  )
}

rs_dir_create <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

rs_plot_dir <- function(..., plot_root = "plots") {
  if (exists("get_plot_dir", mode = "function", inherits = TRUE)) {
    return(get_plot_dir(..., plot_root = plot_root))
  }

  out <- here::here(plot_root, ...)
  rs_dir_create(out)
  out
}

rs_max_plot_iter <- function(dataset) {
  if (dataset == "PFDC") return(82L)
  if (dataset == "Derexyl") return(115L)

  stop("Unknown dataset: ", dataset)
}

rs_cap_runtime_days <- function(df) {
  if (!nrow(df)) return(df)

  df %>%
    rowwise() %>%
    mutate(iter_cap = rs_max_plot_iter(as.character(.data$Dataset))) %>%
    ungroup() %>%
    filter(.data$iter_plot <= .data$iter_cap) %>%
    select(-iter_cap)
}

rs_rolling_median <- function(x, k = 6L) {
  n <- length(x)
  out <- rep(NA_real_, n)

  for (i in seq_len(n)) {
    lo <- max(1L, i - floor(k / 2))
    hi <- min(n, i + floor(k / 2))
    out[[i]] <- stats::median(x[lo:hi], na.rm = TRUE)
  }

  out
}

rs_rolling_mean <- function(x, k = 6L) {
  n <- length(x)
  out <- rep(NA_real_, n)

  for (i in seq_len(n)) {
    lo <- max(1L, i - floor(k / 2))
    hi <- min(n, i + floor(k / 2))
    out[[i]] <- mean(x[lo:hi], na.rm = TRUE)
  }

  out
}


rs_factor_dataset_label <- function(x) {
  labels <- rs_dataset_labels()
  factor(dplyr::recode(as.character(x), !!!labels), levels = unname(labels))
}

rs_score_label <- function(x, items = rs_default_items()) {
  labels <- rs_item_labels()
  factor(
    dplyr::recode(as.character(x), !!!labels),
    levels = unname(labels[items])
  )
}

rs_speedup_distribution_axis_breaks <- function(settings) {
  if (settings$comparison == "xemapred_vs_eczemapred") {
    return(c(0.1, 1, 10, 100, 1000, 10000))
  }

  if (settings$comparison == "xemapred_vs_patient_specific") {
    return(c(0.1, 1, 10, 100, 1000))
  }

  if (settings$comparison == "eczemapred_vs_patient_specific") {
    return(c(0.1, 1, 10, 100, 1000, 10000, 100000, 1000000))
  }

  NULL
}

rs_speedup_overlay_axis_breaks <- function(settings) {
  if (settings$comparison == "xemapred_vs_eczemapred") {
    return(c(1, 10, 100, 1000))
  }

  if (settings$comparison == "xemapred_vs_patient_specific") {
    return(c(0.1, 1, 10, 100, 1000))
  }

  if (settings$comparison == "eczemapred_vs_patient_specific") {
    return(c(0.1, 1, 10, 100, 1000, 10000, 100000, 1000000))
  }

  NULL
}

rs_runtime_axis_breaks <- function(values = NULL) {
  base_breaks <- c(0.01, 0.1, 1, 10, 100, 1000)

  if (is.null(values)) {
    return(base_breaks)
  }

  values <- values[is.finite(values) & values > 0]

  if (!length(values)) {
    return(c(0.1, 1, 10, 100, 1000))
  }

  out <- base_breaks[
    base_breaks >= min(values, na.rm = TRUE) / 1.5 &
      base_breaks <= max(values, na.rm = TRUE) * 1.5
  ]

  if (!length(out)) {
    out <- pretty(values, n = 5)
    out <- out[is.finite(out) & out > 0]
  }

  out
}

rs_runtime_axis_labels <- function(x) {
  vapply(x, function(z) {
    if (is.na(z) || !is.finite(z)) {
      return("")
    }

    if (z >= 1000) {
      return(scales::comma(z, accuracy = 1))
    }

    if (z >= 10) {
      return(format(round(z, 0), trim = TRUE, scientific = FALSE))
    }

    if (z >= 1) {
      return(format(round(z, 1), trim = TRUE, scientific = FALSE))
    }

    format(signif(z, 2), trim = TRUE, scientific = FALSE)
  }, character(1))
}
# -----------------------------------------------------------------------
# Comparison settings
# -----------------------------------------------------------------------

rs_speedup_settings <- function(comparison) {
  comparison <- match.arg(
    comparison,
    c(
      "xemapred_vs_eczemapred",
      "xemapred_vs_patient_specific",
      "eczemapred_vs_patient_specific"
    )
  )

  # ----------------------------------------------------------
  # Cohort runtime vs cohort runtime
  # ----------------------------------------------------------

  if (comparison == "xemapred_vs_eczemapred") {
    return(list(
      comparison = comparison,
      out_stub = "xemapred_vs_eczemapred",
      title = "XemaPred vs EczemaPred",

      model_a = "EczemaPred",
      model_b = "XemaPred",

      model_b_runtime_name = "XemaPred_runtime_sec",

      ratio_label =
        "Cohort-level runtime ratio",

      heatmap_label =
        "EczemaPred cohort runtime / XemaPred cohort runtime",

      high_colour =
        rs_model_cols()[["XemaPred"]],

      comparison_scale = "cohort",
      direct = TRUE
    ))
  }

  # ----------------------------------------------------------
  # Per-patient runtime vs measured patient-specific runtime
  # ----------------------------------------------------------

  if (comparison == "xemapred_vs_patient_specific") {
    return(list(
      comparison = comparison,

      out_stub =
        "xemapred_vs_patient_specific_patientwise",

      title =
        "Patient-specific XemaPred vs population-level XemaPred",

      model_a =
        "XemaPred per patient",

      model_b =
        "Patient-specific XemaPred",

      model_b_runtime_name =
        "PatientSpecific_mean_runtime_sec",

      ratio_label =
        "Per-patient runtime ratio",

      heatmap_label =
        paste(
          "XemaPred cohort runtime per patient /",
          "patient-specific XemaPred runtime"
        ),

      high_colour =
        rs_model_cols()[["Patient-specific XemaPred"]],

      comparison_scale = "per_patient",
      direct = FALSE
    ))
  }

  # ----------------------------------------------------------
  # EczemaPred per patient vs measured patient-specific runtime
  # ----------------------------------------------------------

  list(
    comparison = comparison,

    out_stub =
      "eczemapred_vs_patient_specific_patientwise",

    title =
      "Patient-specific XemaPred vs EczemaPred",

    model_a =
      "EczemaPred per patient",

    model_b =
      "Patient-specific XemaPred",

    model_b_runtime_name =
      "PatientSpecific_mean_runtime_sec",

    ratio_label =
      "Per-patient runtime ratio",

    heatmap_label =
      paste(
        "EczemaPred cohort runtime per patient /",
        "patient-specific XemaPred runtime"
      ),

    high_colour =
      rs_model_cols()[["Patient-specific XemaPred"]],

    comparison_scale = "per_patient",
    direct = FALSE
  )
}

# -----------------------------------------------------------------------
# Runtime loaders
# -----------------------------------------------------------------------

rs_load_cohort_runtime <- function(
    dataset,
    model_label,
    model_prefixes,
    result_root = "results",
    items = rs_default_items(),
    time_priority = c("compute_time", "run_time", "runtime", "elapsed")
) {
  if (exists("sc_load_cohort_runtime", mode = "function", inherits = TRUE)) {
    return(sc_load_cohort_runtime(
      ds = dataset,
      model_label = model_label,
      model_prefixes = model_prefixes,
      result_root = result_root,
      items = items,
      time_priority = time_priority
    ))
  }

  stop(
    "sc_load_cohort_runtime() was not found. Source scalability_helpers.R ",
    "before runtime_speedup_helpers.R."
  )
}

rs_load_xemapred_runtime <- function(dataset, result_root = "results", items = rs_default_items()) {
  rs_load_cohort_runtime(
    dataset = dataset,
    model_label = "XemaPred",
    model_prefixes = c("XemaPred", "Fast"),
    result_root = result_root,
    items = items,
    time_priority = c("compute_time", "run_time", "runtime", "elapsed")
  )
}

rs_load_eczemapred_runtime <- function(dataset, result_root = "results", items = rs_default_items()) {
  rs_load_cohort_runtime(
    dataset = dataset,
    model_label = "EczemaPred",
    model_prefixes = c("EczemaPred"),
    result_root = result_root,
    items = items,
    time_priority = c("run_time", "compute_time", "runtime", "elapsed")
  )
}

rs_load_patient_specific_runtime <- function(
    dataset,
    result_root = "results",
    items = rs_default_items()
) {
  if (exists("sc_load_patient_specific_runtime", mode = "function", inherits = TRUE)) {
    return(sc_load_patient_specific_runtime(
      ds = dataset,
      result_root = result_root,
      items = items,
      model_prefixes = c("SoloFast", "PatientSpecificXemaPred"),
      model_label = "Patient-specific XemaPred"
    ))
  }

  stop(
    "sc_load_patient_specific_runtime() was not found. Source scalability_helpers.R ",
    "before runtime_speedup_helpers.R."
  )
}

# -----------------------------------------------------------------------
# Active-patient counts for the patient-wise comparison
# -----------------------------------------------------------------------

rs_check_forward_chaining_helpers <- function() {
  required <- c("get_fc_iteration", "get_fc_training_iteration", "split_fc_dataset")
  missing <- required[!vapply(required, exists, logical(1), mode = "function", inherits = TRUE)]

  if (length(missing)) {
    stop(
      "Missing forward-chaining helper(s): ", paste(missing, collapse = ", "),
      ". Source your project 00_init.R / functions.R before running the runtime speed-up analysis."
    )
  }
}

rs_load_dataset_checked <- function(dataset) {
  if (exists("load_dataset_checked", mode = "function", inherits = TRUE)) {
    return(load_dataset_checked(dataset))
  }

  if (exists("load_dataset", mode = "function", inherits = TRUE)) {
    return(load_dataset(dataset))
  }

  stop("load_dataset() was not found. Source your project setup first.")
}

rs_get_item_score_column <- function(item) {
  if (exists("get_item_score_column", mode = "function", inherits = TRUE)) {
    return(get_item_score_column(item))
  }

  if (exists("detail_POSCORAD", mode = "function", inherits = TRUE)) {
    item_dict <- detail_POSCORAD()
    item_row <- item_dict %>% filter(.data$Name == item)
    if (!nrow(item_row)) stop("Unknown score: ", item)
    return(as.character(item_row$Label[[1]]))
  }

  stop("Cannot resolve score column for item '", item, "'. Source score_metadata.R or project setup.")
}

rs_get_active_patients_per_iter_score <- function(
    dataset,
    score_name,
    horizon_for_counts = 1L
) {
  rs_check_forward_chaining_helpers()

  score_col <- rs_get_item_score_column(score_name)
  poscorad <- rs_load_dataset_checked(dataset)

  if (!score_col %in% names(poscorad)) {
    warning(
      "Column '", score_col, "' not found in ", dataset,
      "; skipping active-patient counts."
    )
    return(tibble())
  }

  df <- poscorad %>%
    dplyr::transmute(
      Patient = .data$Patient,
      Time = .data$Day,
      Score = suppressWarnings(as.numeric(.data[[score_col]]))
    ) %>%
    tidyr::drop_na(.data$Patient, .data$Time, .data$Score) %>%
    dplyr::mutate(
      Iteration = get_fc_iteration(.data$Time, horizon_for_counts)
    )

  if (!nrow(df)) return(tibble())

  # get_fc_training_iteration() returns the SET of valid forecasting
  # iterations. It is not a row-wise transformation, so call it once.
  train_it <- get_fc_training_iteration(df$Iteration)

  purrr::map_dfr(train_it, function(it) {
    split <- split_fc_dataset(df, it)
    test <- split$Testing

    tibble::tibble(
      Dataset = dataset,
      score = score_name,
      iter_plot = as.integer(it),
      n_active_patients = dplyr::n_distinct(test$Patient)
    )
  })
}


rs_get_active_patients <- function(
    datasets = rs_dataset_levels(),
    items = rs_default_items(),
    horizon_for_counts = 1L
) {
  purrr::map_dfr(datasets, function(dataset) {
    purrr::map_dfr(items, function(item) {
      rs_get_active_patients_per_iter_score(
        dataset = dataset,
        score_name = item,
        horizon_for_counts = horizon_for_counts
      )
    })
  })
}

# -----------------------------------------------------------------------
# Actual accumulated training-observation counts
# -----------------------------------------------------------------------
# A dataset-level observation is one observed patient-day. A patient-day is
# counted once if at least one of the nine PO-SCORAD items is observed.
# This gives a common training-size axis for all nine item models.
# -----------------------------------------------------------------------

rs_get_training_observation_lookup <- function(
    dataset,
    iter_values,
    items = rs_default_items(),
    horizon_for_counts = 1L
) {
  poscorad <- rs_load_dataset_checked(dataset)

  if (!all(c("Patient", "Day") %in% names(poscorad))) {
    stop(
      "Dataset ", dataset,
      " must contain Patient and Day columns to calculate accumulated observations."
    )
  }

  score_cols <- unique(
    vapply(
      items,
      function(item) {
        tryCatch(
          rs_get_item_score_column(item),
          error = function(e) NA_character_
        )
      },
      character(1)
    )
  )

  score_cols <- score_cols[
    !is.na(score_cols) & score_cols %in% names(poscorad)
  ]

  if (!length(score_cols)) {
    stop(
      "Could not resolve any PO-SCORAD item columns in dataset ", dataset,
      " when calculating accumulated training observations."
    )
  }

  visits <- poscorad %>%
    dplyr::filter(
      !is.na(.data$Patient),
      !is.na(.data$Day),
      is.finite(suppressWarnings(as.numeric(.data$Day)))
    ) %>%
    dplyr::filter(
      dplyr::if_any(
        dplyr::all_of(score_cols),
        ~ !is.na(.x)
      )
    ) %>%
    dplyr::transmute(
      Patient = as.character(.data$Patient),
      Time = suppressWarnings(as.numeric(.data$Day))
    ) %>%
    dplyr::distinct(.data$Patient, .data$Time)

  iters <- sort(
    unique(
      as.integer(
        iter_values[!is.na(iter_values) & is.finite(iter_values)]
      )
    )
  )

  purrr::map_dfr(iters, function(it) {
    # Existing runtime convention: iter_plot = 0 corresponds to training day 1.
    training_day <- as.integer(it) + 1L

    tibble::tibble(
      Dataset = dataset,
      iter_plot = as.integer(it),
      Training_day = training_day,
      N_training_observations = sum(visits$Time <= training_day, na.rm = TRUE)
    )
  })
}

rs_attach_training_observation_counts <- function(
    speed_df,
    datasets = rs_dataset_levels(),
    items = rs_default_items(),
    horizon_for_counts = 1L
) {
  lookup <- purrr::map_dfr(datasets, function(dataset) {
    iter_values <- speed_df %>%
      dplyr::filter(as.character(.data$Dataset) == dataset) %>%
      dplyr::pull(.data$iter_plot)

    if (!length(iter_values)) return(tibble::tibble())

    rs_get_training_observation_lookup(
      dataset = dataset,
      iter_values = iter_values,
      items = items,
      horizon_for_counts = horizon_for_counts
    )
  })

  if (!nrow(lookup)) {
    stop("Could not calculate accumulated training-observation counts.")
  }

  speed_df %>%
    dplyr::mutate(
      Dataset = as.character(.data$Dataset),
      iter_plot = as.integer(.data$iter_plot)
    ) %>%
    dplyr::left_join(
      lookup,
      by = c("Dataset", "iter_plot")
    )
}


rs_observation_x_scale <- function(
    df,
    show_training_day_axis = TRUE,
    t_horizon = 4L
) {

  if (isTRUE(show_training_day_axis)) {

    ds <- unique(as.character(df$Dataset))
    ds <- ds[!is.na(ds)]

    if (length(ds) != 1L) {
      stop(
        "The Training days secondary axis requires exactly one dataset."
      )
    }

    training_axis <- make_training_axis_breaks(
      dataset = ds[[1]],
      t_horizon = t_horizon,
      n_breaks = 10
    )

    return(
      ggplot2::scale_x_continuous(
        sec.axis = ggplot2::dup_axis(
          breaks = training_axis$N,
          labels = training_axis$LastTime,
          name = "Training days"
        )
      )
    )
  }

  ggplot2::scale_x_continuous()
}

# -----------------------------------------------------------------------
# Build standard speed-up tables
# -----------------------------------------------------------------------

rs_make_direct_speedup_data <- function(
    datasets = rs_dataset_levels(),
    items = rs_default_items(),
    result_root = "results"
) {
  settings <- rs_speedup_settings("xemapred_vs_eczemapred")

  runtime <- bind_rows(
    purrr::map_dfr(datasets, ~ rs_load_eczemapred_runtime(.x, result_root = result_root, items = items)),
    purrr::map_dfr(datasets, ~ rs_load_xemapred_runtime(.x, result_root = result_root, items = items))
  ) %>%
    filter(
      .data$Model %in% c("EczemaPred", "XemaPred"),
      .data$score %in% items,
      is.finite(.data$runtime_sec),
      .data$runtime_sec > 0,
      !is.na(.data$iter_plot)
    )

  if (!nrow(runtime)) stop("No direct EczemaPred/XemaPred runtime data loaded.")

  wide <- runtime %>%
    select(Dataset, score, iter, iter_plot, Model, runtime_sec) %>%
    tidyr::pivot_wider(names_from = Model, values_from = runtime_sec) %>%
    filter(
      is.finite(.data$EczemaPred),
      is.finite(.data$XemaPred),
      .data$EczemaPred > 0,
      .data$XemaPred > 0
    )

  wide %>%
    transmute(
      comparison = settings$comparison,
      Dataset = as.character(.data$Dataset),
      Dataset_label = rs_factor_dataset_label(.data$Dataset),
      score = factor(as.character(.data$score), levels = items),
      score_label = rs_score_label(.data$score, items = items),
      iter = .data$iter,
      iter_plot = as.integer(.data$iter_plot),
      model_a = settings$model_a,
      model_b = settings$model_b,
      runtime_a_sec = .data$EczemaPred,
      runtime_b_sec = .data$XemaPred,
      EczemaPred_runtime_sec = .data$EczemaPred,
      XemaPred_runtime_sec = .data$XemaPred,
      speedup = .data$EczemaPred / .data$XemaPred,
      runtime_reduction_pct = 100 * (1 - .data$XemaPred / .data$EczemaPred),
      n_active_patients = NA_integer_,
      n_patient_specific_patients = NA_integer_
    ) %>%
    filter(is.finite(.data$speedup), .data$speedup > 0)
}

rs_make_patientwise_speedup_data <- function(
    datasets = rs_dataset_levels(),
    items = rs_default_items(),
    result_root = "results",
    horizon_for_counts = 1L
) {
  settings <- rs_speedup_settings(
    "xemapred_vs_patient_specific"
  )

  # ------------------------------------------------------------
  # 1. Population-level XemaPred cohort runtime
  # ------------------------------------------------------------

  xemapred <- purrr::map_dfr(
    datasets,
    ~ rs_load_xemapred_runtime(
      dataset = .x,
      result_root = result_root,
      items = items
    )
  ) %>%
    dplyr::filter(
      .data$score %in% items,
      is.finite(.data$runtime_sec),
      .data$runtime_sec > 0,
      !is.na(.data$iter_plot)
    )

  # ------------------------------------------------------------
  # 2. Actual patient-specific runtimes
  # ------------------------------------------------------------

  patient_specific_raw <- purrr::map_dfr(
    datasets,
    ~ rs_load_patient_specific_runtime(
      dataset = .x,
      result_root = result_root,
      items = items
    )
  ) %>%
    dplyr::filter(
      .data$score %in% items,
      !is.na(.data$patient),
      is.finite(.data$runtime_sec),
      .data$runtime_sec > 0,
      !is.na(.data$iter_plot)
    )

  if (!nrow(xemapred)) {
    stop(
      "No XemaPred runtime data loaded for patient-wise comparison."
    )
  }

  if (!nrow(patient_specific_raw)) {
    stop(
      "No patient-specific XemaPred runtime data loaded."
    )
  }

  # ------------------------------------------------------------
  # 3. Summarise patient-specific runtime by item and update
  # ------------------------------------------------------------

  patient_specific <- patient_specific_raw %>%
    dplyr::group_by(
      .data$Dataset,
      .data$score,
      .data$iter_plot
    ) %>%
    dplyr::summarise(
      PatientSpecific_mean_runtime_sec =
        mean(.data$runtime_sec, na.rm = TRUE),

      PatientSpecific_median_runtime_sec =
        median(.data$runtime_sec, na.rm = TRUE),

      PatientSpecific_q25_runtime_sec =
        stats::quantile(
          .data$runtime_sec,
          0.25,
          na.rm = TRUE,
          names = FALSE
        ),

      PatientSpecific_q75_runtime_sec =
        stats::quantile(
          .data$runtime_sec,
          0.75,
          na.rm = TRUE,
          names = FALSE
        ),

      PatientSpecific_q95_runtime_sec =
        stats::quantile(
          .data$runtime_sec,
          0.95,
          na.rm = TRUE,
          names = FALSE
        ),

      PatientSpecific_max_runtime_sec =
        max(.data$runtime_sec, na.rm = TRUE),

      n_patient_specific_patients =
        dplyr::n_distinct(.data$patient),

      .groups = "drop"
    ) %>%
    dplyr::filter(
      is.finite(.data$PatientSpecific_mean_runtime_sec),
      is.finite(.data$PatientSpecific_median_runtime_sec),
      .data$PatientSpecific_mean_runtime_sec > 0,
      .data$PatientSpecific_median_runtime_sec > 0,
      .data$n_patient_specific_patients > 0
    )

  # ------------------------------------------------------------
  # 4. Convert population XemaPred runtime to per-patient runtime
  # ------------------------------------------------------------

  xemapred_patientwise <- xemapred %>%
    dplyr::transmute(
      Dataset = as.character(.data$Dataset),
      score = as.character(.data$score),
      iter = .data$iter,
      iter_plot = as.integer(.data$iter_plot),
      cohort_runtime_sec = .data$runtime_sec
    ) %>%
    dplyr::inner_join(
      patient_specific,
      by = c("Dataset", "score", "iter_plot")
    ) %>%
    dplyr::mutate(
      XemaPred_per_patient_sec =
        .data$cohort_runtime_sec /
        .data$n_patient_specific_patients
    )

  # ------------------------------------------------------------
  # 5. Calculate patient-wise speed-up
  # ------------------------------------------------------------

  xemapred_patientwise %>%
    dplyr::filter(
      is.finite(.data$XemaPred_per_patient_sec),
      is.finite(.data$PatientSpecific_mean_runtime_sec),
      is.finite(.data$PatientSpecific_median_runtime_sec),
      .data$XemaPred_per_patient_sec > 0,
      .data$PatientSpecific_mean_runtime_sec > 0,
      .data$PatientSpecific_median_runtime_sec > 0
    ) %>%
    dplyr::mutate(
      speedup_mean_patient_specific =
        .data$XemaPred_per_patient_sec /
        .data$PatientSpecific_mean_runtime_sec,

      speedup_median_patient_specific =
        .data$XemaPred_per_patient_sec /
        .data$PatientSpecific_median_runtime_sec,

      runtime_reduction_pct_mean =
        100 * (
          1 -
          .data$PatientSpecific_mean_runtime_sec /
            .data$XemaPred_per_patient_sec
        ),

      runtime_reduction_pct_median =
        100 * (
          1 -
          .data$PatientSpecific_median_runtime_sec /
            .data$XemaPred_per_patient_sec
        )
    ) %>%
    dplyr::transmute(
      comparison = settings$comparison,

      Dataset = as.character(.data$Dataset),
      Dataset_label =
        rs_factor_dataset_label(.data$Dataset),

      score =
        factor(as.character(.data$score), levels = items),

      score_label =
        rs_score_label(.data$score, items = items),

      iter = .data$iter,
      iter_plot = as.integer(.data$iter_plot),

      model_a = settings$model_a,
      model_b = settings$model_b,

      runtime_a_sec =
        .data$XemaPred_per_patient_sec,

      runtime_b_sec =
        .data$PatientSpecific_mean_runtime_sec,

      cohort_runtime_sec =
        .data$cohort_runtime_sec,

      XemaPred_per_patient_sec =
        .data$XemaPred_per_patient_sec,

      PatientSpecific_mean_runtime_sec =
        .data$PatientSpecific_mean_runtime_sec,

      PatientSpecific_median_runtime_sec =
        .data$PatientSpecific_median_runtime_sec,

      PatientSpecific_q25_runtime_sec =
        .data$PatientSpecific_q25_runtime_sec,

      PatientSpecific_q75_runtime_sec =
        .data$PatientSpecific_q75_runtime_sec,

      PatientSpecific_q95_runtime_sec =
        .data$PatientSpecific_q95_runtime_sec,

      PatientSpecific_max_runtime_sec =
        .data$PatientSpecific_max_runtime_sec,

      speedup =
        .data$speedup_mean_patient_specific,

      speedup_mean_patient_specific =
        .data$speedup_mean_patient_specific,

      speedup_median_patient_specific =
        .data$speedup_median_patient_specific,

      runtime_reduction_pct =
        .data$runtime_reduction_pct_mean,

      runtime_reduction_pct_mean =
        .data$runtime_reduction_pct_mean,

      runtime_reduction_pct_median =
        .data$runtime_reduction_pct_median,

      n_active_patients =
        .data$n_patient_specific_patients,

      n_patient_specific_patients =
        .data$n_patient_specific_patients
    ) %>%
    dplyr::filter(
      is.finite(.data$speedup),
      .data$speedup > 0
    )
}

rs_make_eczemapred_vs_patient_specific_data <- function(
    datasets = rs_dataset_levels(),
    items = rs_default_items(),
    result_root = "results",
    horizon_for_counts = 1L
) {
  settings <- rs_speedup_settings(
    "eczemapred_vs_patient_specific"
  )

  # ------------------------------------------------------------
  # 1. Load cohort-level EczemaPred runtime
  # ------------------------------------------------------------

  eczemapred <- purrr::map_dfr(
    datasets,
    ~ rs_load_eczemapred_runtime(
      dataset = .x,
      result_root = result_root,
      items = items
    )
  ) %>%
    dplyr::filter(
      .data$score %in% items,
      is.finite(.data$runtime_sec),
      .data$runtime_sec > 0,
      !is.na(.data$iter_plot)
    ) %>%
    dplyr::transmute(
      Dataset = as.character(.data$Dataset),
      score = as.character(.data$score),
      iter = .data$iter,
      iter_plot = as.integer(.data$iter_plot),
      EczemaPred_cohort_runtime_sec = .data$runtime_sec
    )

  # ------------------------------------------------------------
  # 2. Load measured patient-specific XemaPred runtimes
  # ------------------------------------------------------------

  patient_specific_raw <- purrr::map_dfr(
    datasets,
    ~ rs_load_patient_specific_runtime(
      dataset = .x,
      result_root = result_root,
      items = items
    )
  ) %>%
    dplyr::filter(
      .data$score %in% items,
      !is.na(.data$patient),
      is.finite(.data$runtime_sec),
      .data$runtime_sec > 0,
      !is.na(.data$iter_plot)
    )

  if (!nrow(eczemapred)) {
    stop("No EczemaPred runtime data loaded.")
  }

  if (!nrow(patient_specific_raw)) {
    stop("No patient-specific XemaPred runtime data loaded.")
  }

  # ------------------------------------------------------------
  # 3. Summarise patient-specific runtimes and patient counts
  # ------------------------------------------------------------

  patient_specific <- patient_specific_raw %>%
    dplyr::group_by(
      .data$Dataset,
      .data$score,
      .data$iter_plot
    ) %>%
    dplyr::summarise(
      PatientSpecific_mean_runtime_sec =
        mean(.data$runtime_sec, na.rm = TRUE),

      PatientSpecific_median_runtime_sec =
        median(.data$runtime_sec, na.rm = TRUE),

      PatientSpecific_q25_runtime_sec =
        stats::quantile(
          .data$runtime_sec,
          0.25,
          na.rm = TRUE,
          names = FALSE
        ),

      PatientSpecific_q75_runtime_sec =
        stats::quantile(
          .data$runtime_sec,
          0.75,
          na.rm = TRUE,
          names = FALSE
        ),

      PatientSpecific_q95_runtime_sec =
        stats::quantile(
          .data$runtime_sec,
          0.95,
          na.rm = TRUE,
          names = FALSE
        ),

      PatientSpecific_max_runtime_sec =
        max(.data$runtime_sec, na.rm = TRUE),

      n_patient_specific_patients =
        dplyr::n_distinct(.data$patient),

      .groups = "drop"
    ) %>%
    dplyr::filter(
      is.finite(.data$PatientSpecific_mean_runtime_sec),
      is.finite(.data$PatientSpecific_median_runtime_sec),
      .data$PatientSpecific_mean_runtime_sec > 0,
      .data$PatientSpecific_median_runtime_sec > 0,
      .data$n_patient_specific_patients > 0
    )

  # ------------------------------------------------------------
  # 4. Match EczemaPred with patient-specific XemaPred
  #    and amortise EczemaPred over matched patients
  # ------------------------------------------------------------

  eczemapred %>%
    dplyr::inner_join(
      patient_specific,
      by = c("Dataset", "score", "iter_plot")
    ) %>%
    dplyr::mutate(
      EczemaPred_per_patient_sec =
        .data$EczemaPred_cohort_runtime_sec /
        .data$n_patient_specific_patients
    ) %>%
    dplyr::filter(
      is.finite(.data$EczemaPred_per_patient_sec),
      is.finite(.data$PatientSpecific_mean_runtime_sec),
      is.finite(.data$PatientSpecific_median_runtime_sec),
      .data$EczemaPred_per_patient_sec > 0,
      .data$PatientSpecific_mean_runtime_sec > 0,
      .data$PatientSpecific_median_runtime_sec > 0
    ) %>%
    dplyr::mutate(
      speedup_mean_patient_specific =
        .data$EczemaPred_per_patient_sec /
        .data$PatientSpecific_mean_runtime_sec,

      speedup_median_patient_specific =
        .data$EczemaPred_per_patient_sec /
        .data$PatientSpecific_median_runtime_sec,

      runtime_reduction_pct_mean =
        100 * (
          1 -
            .data$PatientSpecific_mean_runtime_sec /
            .data$EczemaPred_per_patient_sec
        ),

      runtime_reduction_pct_median =
        100 * (
          1 -
            .data$PatientSpecific_median_runtime_sec /
            .data$EczemaPred_per_patient_sec
        )
    ) %>%
    dplyr::transmute(
      comparison = settings$comparison,

      Dataset = as.character(.data$Dataset),
      Dataset_label =
        rs_factor_dataset_label(.data$Dataset),

      score =
        factor(as.character(.data$score), levels = items),

      score_label =
        rs_score_label(.data$score, items = items),

      iter = .data$iter,
      iter_plot = as.integer(.data$iter_plot),

      model_a = settings$model_a,
      model_b = settings$model_b,

      runtime_a_sec =
        .data$EczemaPred_per_patient_sec,

      runtime_b_sec =
        .data$PatientSpecific_mean_runtime_sec,

      EczemaPred_cohort_runtime_sec =
        .data$EczemaPred_cohort_runtime_sec,

      EczemaPred_per_patient_sec =
        .data$EczemaPred_per_patient_sec,

      PatientSpecific_mean_runtime_sec =
        .data$PatientSpecific_mean_runtime_sec,

      PatientSpecific_median_runtime_sec =
        .data$PatientSpecific_median_runtime_sec,

      PatientSpecific_q25_runtime_sec =
        .data$PatientSpecific_q25_runtime_sec,

      PatientSpecific_q75_runtime_sec =
        .data$PatientSpecific_q75_runtime_sec,

      PatientSpecific_q95_runtime_sec =
        .data$PatientSpecific_q95_runtime_sec,

      PatientSpecific_max_runtime_sec =
        .data$PatientSpecific_max_runtime_sec,

      speedup =
        .data$speedup_mean_patient_specific,

      speedup_mean_patient_specific =
        .data$speedup_mean_patient_specific,

      speedup_median_patient_specific =
        .data$speedup_median_patient_specific,

      runtime_reduction_pct =
        .data$runtime_reduction_pct_mean,

      runtime_reduction_pct_mean =
        .data$runtime_reduction_pct_mean,

      runtime_reduction_pct_median =
        .data$runtime_reduction_pct_median,

      n_active_patients =
        .data$n_patient_specific_patients,

      n_patient_specific_patients =
        .data$n_patient_specific_patients
    ) %>%
    dplyr::filter(
      is.finite(.data$speedup),
      .data$speedup > 0
    )
}

rs_make_speedup_data <- function(
    comparison = c(
      "xemapred_vs_eczemapred",
      "xemapred_vs_patient_specific",
      "eczemapred_vs_patient_specific"
    ),
    datasets = rs_dataset_levels(),
    items = rs_default_items(),
    result_root = "results",
    horizon_for_counts = 1L
) {
  comparison <- match.arg(comparison)

  if (comparison == "xemapred_vs_eczemapred") {
    return(
      rs_make_direct_speedup_data(
        datasets = datasets,
        items = items,
        result_root = result_root
      )
    )
  }

  if (comparison == "xemapred_vs_patient_specific") {
    return(
      rs_make_patientwise_speedup_data(
        datasets = datasets,
        items = items,
        result_root = result_root,
        horizon_for_counts = horizon_for_counts
      )
    )
  }

  rs_make_eczemapred_vs_patient_specific_data(
    datasets = datasets,
    items = items,
    result_root = result_root,
    horizon_for_counts = horizon_for_counts
  )
}

# -----------------------------------------------------------------------
# Summaries
# -----------------------------------------------------------------------

rs_speedup_summary <- function(speed_df) {
  speed_df %>%
    group_by(.data$Dataset, .data$Dataset_label, .data$model_a, .data$model_b) %>%
    summarise(
      n_points = n(),
      n_iters = n_distinct(.data$iter_plot),
      n_items = n_distinct(.data$score),
      mean_model_a_runtime_sec = mean(.data$runtime_a_sec, na.rm = TRUE),
      median_model_a_runtime_sec = median(.data$runtime_a_sec, na.rm = TRUE),
      mean_model_b_runtime_sec = mean(.data$runtime_b_sec, na.rm = TRUE),
      median_model_b_runtime_sec = median(.data$runtime_b_sec, na.rm = TRUE),
      mean_runtime_saved_sec = mean(.data$runtime_a_sec - .data$runtime_b_sec, na.rm = TRUE),
      median_runtime_saved_sec = median(.data$runtime_a_sec - .data$runtime_b_sec, na.rm = TRUE),
      median_speedup = median(.data$speedup, na.rm = TRUE),
      p25_speedup = quantile(.data$speedup, 0.25, na.rm = TRUE, names = FALSE),
      p75_speedup = quantile(.data$speedup, 0.75, na.rm = TRUE, names = FALSE),
      mean_speedup = mean(.data$speedup, na.rm = TRUE),
      max_speedup = max(.data$speedup, na.rm = TRUE),
      median_runtime_reduction_pct = median(.data$runtime_reduction_pct, na.rm = TRUE),
      final_iter = max(iter_plot, na.rm = TRUE),
      final_speedup = median(speedup[iter_plot == max(iter_plot, na.rm = TRUE)], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      median_speedup_label = paste0(round(.data$median_speedup, 2), "x"),
      final_speedup_label = paste0(round(.data$final_speedup, 2), "x"),
      median_reduction_label = paste0(round(.data$median_runtime_reduction_pct, 1), "%")
    )
}

rs_speedup_summary_by_score <- function(speed_df) {
  speed_df %>%
    group_by(.data$Dataset, .data$Dataset_label, .data$score, .data$score_label) %>%
    summarise(
      n_points = n(),
      median_speedup = median(.data$speedup, na.rm = TRUE),
      p25_speedup = quantile(.data$speedup, 0.25, na.rm = TRUE, names = FALSE),
      p75_speedup = quantile(.data$speedup, 0.75, na.rm = TRUE, names = FALSE),
      mean_speedup = mean(.data$speedup, na.rm = TRUE),
      median_runtime_reduction_pct = median(.data$runtime_reduction_pct, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(median_speedup_label = paste0(round(.data$median_speedup, 2), "x"))
}

# -----------------------------------------------------------------------
# Cohort-level runtime reduction based on the slowest item
#
# The nine PO-SCORAD item models are run in parallel, so the wall-clock
# latency of one forecasting update is determined by the slowest item.
#
# Therefore:
#
#   cohort-level reduction =
#     max(EczemaPred item runtime) / max(XemaPred item runtime)
#
# This is the preferred headline runtime-reduction quantity for the
# EczemaPred vs XemaPred comparison.
# -----------------------------------------------------------------------

rs_max_runtime_speedup_by_iteration <- function(speed_df) {
  required_cols <- c(
    "Dataset",
    "Dataset_label",
    "iter_plot",
    "Training_day",
    "N_training_observations",
    "score",
    "runtime_a_sec",
    "runtime_b_sec"
  )

  missing_cols <- setdiff(required_cols, names(speed_df))
  if (length(missing_cols)) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  speed_df %>%
    dplyr::filter(
      is.finite(.data$runtime_a_sec),
      is.finite(.data$runtime_b_sec),
      .data$runtime_a_sec > 0,
      .data$runtime_b_sec > 0,
      is.finite(.data$N_training_observations)
    ) %>%
    dplyr::group_by(
      .data$Dataset,
      .data$Dataset_label,
      .data$iter_plot,
      .data$Training_day,
      .data$N_training_observations
    ) %>%
    dplyr::summarise(
      EczemaPred_max_runtime_sec = max(.data$runtime_a_sec, na.rm = TRUE),
      XemaPred_max_runtime_sec = max(.data$runtime_b_sec, na.rm = TRUE),
      N_items = dplyr::n_distinct(.data$score),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      max_runtime_speedup =
        .data$EczemaPred_max_runtime_sec / .data$XemaPred_max_runtime_sec,
      runtime_reduction_pct =
        100 * (1 - .data$XemaPred_max_runtime_sec / .data$EczemaPred_max_runtime_sec)
    ) %>%
    dplyr::arrange(.data$Dataset, .data$N_training_observations)
}

rs_max_runtime_patient_specific_by_iteration <- function(
    speed_df,
    comparison = c(
      "xemapred_vs_patient_specific",
      "eczemapred_vs_patient_specific"
    )
) {
  comparison <- match.arg(comparison)

  common_required <- c(
    "Dataset",
    "Dataset_label",
    "iter_plot",
    "Training_day",
    "N_training_observations",
    "score",
    "PatientSpecific_max_runtime_sec"
  )

  if (comparison == "xemapred_vs_patient_specific") {

    required_cols <- c(
      common_required,
      "cohort_runtime_sec"
    )

    missing_cols <- setdiff(required_cols, names(speed_df))

    if (length(missing_cols)) {
      stop(
        "Missing required columns: ",
        paste(missing_cols, collapse = ", ")
      )
    }

    return(
      speed_df %>%
        dplyr::filter(
          is.finite(.data$cohort_runtime_sec),
          is.finite(.data$PatientSpecific_max_runtime_sec),
          is.finite(.data$N_training_observations),
          .data$cohort_runtime_sec > 0,
          .data$PatientSpecific_max_runtime_sec > 0
        ) %>%
        dplyr::group_by(
          .data$Dataset,
          .data$Dataset_label,
          .data$iter_plot,
          .data$Training_day,
          .data$N_training_observations
        ) %>%
        dplyr::summarise(
          XemaPred_max_runtime_sec =
            max(.data$cohort_runtime_sec, na.rm = TRUE),

          PatientSpecific_max_runtime_sec =
            max(.data$PatientSpecific_max_runtime_sec, na.rm = TRUE),

          N_items =
            dplyr::n_distinct(.data$score),

          .groups = "drop"
        ) %>%
        dplyr::mutate(
          max_runtime_speedup =
            .data$XemaPred_max_runtime_sec /
            .data$PatientSpecific_max_runtime_sec,

          runtime_reduction_pct =
            100 * (
              1 -
              .data$PatientSpecific_max_runtime_sec /
                .data$XemaPred_max_runtime_sec
            )
        ) %>%
        dplyr::arrange(
          .data$Dataset,
          .data$N_training_observations
        )
    )
  }

  # ------------------------------------------------------------
  # EczemaPred vs patient-specific XemaPred
  # ------------------------------------------------------------

  required_cols <- c(
    common_required,
    "EczemaPred_cohort_runtime_sec"
  )

  missing_cols <- setdiff(required_cols, names(speed_df))

  if (length(missing_cols)) {
    stop(
      "Missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  speed_df %>%
    dplyr::filter(
      is.finite(.data$EczemaPred_cohort_runtime_sec),
      is.finite(.data$PatientSpecific_max_runtime_sec),
      is.finite(.data$N_training_observations),
      .data$EczemaPred_cohort_runtime_sec > 0,
      .data$PatientSpecific_max_runtime_sec > 0
    ) %>%
    dplyr::group_by(
      .data$Dataset,
      .data$Dataset_label,
      .data$iter_plot,
      .data$Training_day,
      .data$N_training_observations
    ) %>%
    dplyr::summarise(
      EczemaPred_max_runtime_sec =
        max(.data$EczemaPred_cohort_runtime_sec, na.rm = TRUE),

      PatientSpecific_max_runtime_sec =
        max(.data$PatientSpecific_max_runtime_sec, na.rm = TRUE),

      N_items =
        dplyr::n_distinct(.data$score),

      .groups = "drop"
    ) %>%
    dplyr::mutate(
      max_runtime_speedup =
        .data$EczemaPred_max_runtime_sec /
        .data$PatientSpecific_max_runtime_sec,

      runtime_reduction_pct =
        100 * (
          1 -
          .data$PatientSpecific_max_runtime_sec /
            .data$EczemaPred_max_runtime_sec
        )
    ) %>%
    dplyr::arrange(
      .data$Dataset,
      .data$N_training_observations
    )
}


rs_plot_max_runtime_speedup <- function(
    max_df,
    settings
) {
  cols <- rs_dataset_cols()

  ggplot2::ggplot(
    max_df,
    ggplot2::aes(
      x = .data$iter_plot,
      y = .data$max_runtime_speedup,
      colour = .data$Dataset_label
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 1,
      linetype = "dashed",
      colour = "grey50",
      linewidth = 0.5
    ) +
    ggplot2::geom_line(
      linewidth = 1.2
    ) +
    ggplot2::geom_point(
      size = 1.3,
      alpha = 0.65
    ) +
    ggplot2::scale_colour_manual(
      values = cols,
      drop = FALSE
    ) +
    ggplot2::scale_y_log10(
      breaks = c(1, 3, 10, 30, 100, 300, 1000),
      labels = rs_runtime_axis_labels
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(0, 120, by = 20),
      expand = ggplot2::expansion(mult = c(0.01, 0.01))
    ) +
    ggplot2::labs(
      x = "Training days",
      y = if (settings$comparison == "xemapred_vs_eczemapred") {
        "Full-update runtime ratio"
      } else {
        "Runtime ratio"
      },
      colour = NULL
    ) +
    rs_theme_runtime(base_size = 15)
}

# -----------------------------------------------------------------------
# Predetermined checkpoints for the maximum-runtime comparison
# -----------------------------------------------------------------------

rs_max_runtime_checkpoints <- function(max_speedup_df) {
  max_speedup_df %>%
    dplyr::group_by(.data$Dataset) %>%
    dplyr::group_modify(~ {
      its <- sort(unique(.x$iter_plot))
      n_it <- length(its)

      tibble::tibble(
        iter_plot = c(
          its[[1]],
          its[[ceiling(n_it / 2)]],
          its[[n_it]]
        ),
        Stage = c("Early", "Middle", "Final")
      )
    }) %>%
    dplyr::ungroup() %>%
    dplyr::inner_join(
      max_speedup_df,
      by = c("Dataset", "iter_plot")
    ) %>%
    dplyr::select(
      .data$Dataset,
      .data$Dataset_label,
      .data$Stage,
      .data$iter_plot,
      .data$Training_day,
      .data$N_training_observations,
      .data$EczemaPred_max_runtime_sec,
      .data$XemaPred_max_runtime_sec,
      .data$max_runtime_speedup,
      .data$runtime_reduction_pct,
      .data$N_items
    ) %>%
    dplyr::arrange(
      .data$Dataset,
      factor(.data$Stage, levels = c("Early", "Middle", "Final"))
    )
}


rs_print_max_runtime_checkpoints <- function(df) {
  cat(sprintf(
    "%-10s %-7s %5s %9s %12s %12s %12s %7s\n",
    "Dataset", "Stage", "Day", "N obs", "EczemaPred", "XemaPred", "Reduction", "Items"
  ))

  cat(paste(rep("-", 92), collapse = ""), "\n", sep = "")

  for (i in seq_len(nrow(df))) {
    cat(sprintf(
      "%-10s %-7s %5d %9s %12s %12s %11.1fx %7s\n",
      as.character(df$Dataset_label[[i]]),
      df$Stage[[i]],
      as.integer(df$Training_day[[i]]),
      scales::comma(df$N_training_observations[[i]], accuracy = 1),
      rs_sec_label(df$EczemaPred_max_runtime_sec[[i]]),
      rs_sec_label(df$XemaPred_max_runtime_sec[[i]]),
      df$max_runtime_speedup[[i]],
      paste0(df$N_items[[i]], "/9")
    ))
  }

  invisible(df)
}


rs_speedup_overlay_data <- function(speed_df, smooth_k = 6L) {
  speed_df %>%
    dplyr::filter(
      is.finite(.data$speedup),
      .data$speedup > 0,
      is.finite(.data$N_training_observations),
      is.finite(.data$Training_day)
    ) %>%
    dplyr::group_by(
      .data$Dataset,
      .data$Dataset_label,
      .data$iter_plot,
      .data$Training_day,
      .data$N_training_observations
    ) %>%
    dplyr::summarise(
      n_items = dplyr::n_distinct(.data$score),
      median_speedup = median(.data$speedup, na.rm = TRUE),
      p25_speedup = stats::quantile(.data$speedup, 0.25, na.rm = TRUE, names = FALSE),
      p75_speedup = stats::quantile(.data$speedup, 0.75, na.rm = TRUE, names = FALSE),
      mean_speedup = mean(.data$speedup, na.rm = TRUE),
      sd_speedup = stats::sd(.data$speedup, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      mean_minus_sd = pmax(.data$mean_speedup - .data$sd_speedup, .Machine$double.eps),
      mean_plus_sd = .data$mean_speedup + .data$sd_speedup
    ) %>%
    dplyr::arrange(.data$Dataset, .data$N_training_observations) %>%
    dplyr::group_by(.data$Dataset) %>%
    dplyr::mutate(
      median_speedup_smooth = rs_rolling_median(.data$median_speedup, k = smooth_k),
      ymin_smooth = rs_rolling_median(.data$p25_speedup, k = smooth_k),
      ymax_smooth = rs_rolling_median(.data$p75_speedup, k = smooth_k),
      mean_speedup_smooth = rs_rolling_mean(.data$mean_speedup, k = smooth_k),
      mean_minus_sd_smooth = rs_rolling_mean(.data$mean_minus_sd, k = smooth_k),
      mean_plus_sd_smooth = rs_rolling_mean(.data$mean_plus_sd, k = smooth_k)
    ) %>%
    dplyr::ungroup()
}


rs_panel_summary <- function(speedup_overlay_df) {
  speedup_overlay_df %>%
    dplyr::group_by(.data$Dataset, .data$Dataset_label) %>%
    dplyr::summarise(
      n_training_points = dplyr::n_distinct(.data$iter_plot),
      first_day = min(.data$Training_day, na.rm = TRUE),
      final_day = max(.data$Training_day, na.rm = TRUE),
      first_n_observations = .data$N_training_observations[which.min(.data$iter_plot)],
      final_n_observations = .data$N_training_observations[which.max(.data$iter_plot)],
      first_median_speedup = .data$median_speedup[which.min(.data$iter_plot)],
      final_median_speedup = .data$median_speedup[which.max(.data$iter_plot)],
      first_mean_speedup = .data$mean_speedup[which.min(.data$iter_plot)],
      final_mean_speedup = .data$mean_speedup[which.max(.data$iter_plot)],
      final_p25_speedup = .data$p25_speedup[which.max(.data$iter_plot)],
      final_p75_speedup = .data$p75_speedup[which.max(.data$iter_plot)],
      final_sd_speedup = .data$sd_speedup[which.max(.data$iter_plot)],
      max_smoothed_speedup = max(.data$median_speedup_smooth, na.rm = TRUE),
      observations_at_max_smoothed_speedup =
        .data$N_training_observations[which.max(.data$median_speedup_smooth)],
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      final_median_label = paste0(
        round(.data$final_median_speedup, 1), "x [IQR ",
        round(.data$final_p25_speedup, 1), "–",
        round(.data$final_p75_speedup, 1), "]"
      ),
      final_mean_sd_label = paste0(
        round(.data$final_mean_speedup, 1), "x ± ",
        round(.data$final_sd_speedup, 1)
      )
    )
}


# -----------------------------------------------------------------------
# Plot themes
# -----------------------------------------------------------------------

rs_theme_runtime <- function(base_size = 15) {
  theme_classic(base_size = base_size) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 12, colour = "grey20"),
      legend.key.width = grid::unit(0.9, "cm"),
      legend.spacing.x = grid::unit(4, "pt"),
      legend.margin = margin(t = 4, b = 4),
      axis.title.y = element_text(face = "bold", size = 13, margin = margin(r = 8)),
      axis.title.x = element_text(size = 12, margin = margin(t = 8)),
      axis.title.x.top = element_text(size = 11, face = "bold", margin = margin(b = 6)),
      axis.text = element_text(colour = "grey20"),
      axis.text.x.top = element_text(size = 10, colour = "grey20"),
      axis.line = element_line(colour = "grey20", linewidth = 0.6),
      axis.ticks = element_line(colour = "grey20", linewidth = 0.5),
      panel.grid.major.y = element_line(
        colour = "grey86", linetype = "dashed", linewidth = 0.35
      ),
      panel.grid.major.x = element_line(
        colour = "grey90", linetype = "dashed", linewidth = 0.30
      ),
      panel.grid.minor = element_blank(),
      plot.margin = margin(8, 14, 10, 12)
    )
}


# -----------------------------------------------------------------------
# Plots
# -----------------------------------------------------------------------

rs_plot_speedup_overlay <- function(
    speed_df,
    settings,
    smooth_k = 6L,
    show_training_day_axis = FALSE
) {
  overlay_df <- rs_speedup_overlay_data(speed_df, smooth_k = smooth_k)
  cols <- rs_dataset_cols()

  ggplot(
    overlay_df,
    aes(
      x = .data$N_training_observations,
      colour = .data$Dataset_label,
      fill = .data$Dataset_label,
      group = .data$Dataset_label
    )
  ) +
    geom_hline(
      yintercept = 1,
      linetype = "dashed",
      colour = "grey50",
      linewidth = 0.5
    ) +
    geom_ribbon(
      aes(ymin = .data$ymin_smooth, ymax = .data$ymax_smooth),
      alpha = 0.10,
      colour = NA
    ) +
    geom_line(
      aes(y = .data$median_speedup_smooth),
      linewidth = 1.25,
      lineend = "round"
    ) +
    scale_colour_manual(values = cols, drop = FALSE) +
    scale_fill_manual(values = cols, drop = FALSE) +
    scale_y_log10(
      breaks = if (!is.null(rs_speedup_overlay_axis_breaks(settings))) {
        rs_speedup_overlay_axis_breaks(settings)
      } else {
        rs_runtime_axis_breaks(speed_df$speedup)
      },
      labels = rs_runtime_axis_labels,
      expand = expansion(mult = c(0.04, 0.08))
    ) +
    rs_observation_x_scale(
      overlay_df,
      show_training_day_axis = show_training_day_axis
    ) +
    labs(
      x = "Number of training observations",
      y = if (settings$comparison == "xemapred_vs_eczemapred") {
        "Median item-level runtime ratio"
      } else {
        settings$ratio_label
      },
      colour = NULL,
      fill = NULL
    ) +
    rs_theme_runtime(base_size = 15)
}


rs_plot_speedup_raw_median_and_rolling <- function(
    speed_df,
    settings,
    smooth_k = 6L,
    show_training_day_axis = FALSE
) {
  overlay_df <- rs_speedup_overlay_data(speed_df, smooth_k = smooth_k)
  cols <- rs_dataset_cols()

  ggplot(
    overlay_df,
    aes(
      x = .data$N_training_observations,
      colour = .data$Dataset_label,
      fill = .data$Dataset_label,
      group = .data$Dataset_label
    )
  ) +
    geom_hline(
      yintercept = 1,
      linetype = "dashed",
      colour = "grey50",
      linewidth = 0.5
    ) +
    geom_ribbon(
      aes(ymin = .data$ymin_smooth, ymax = .data$ymax_smooth),
      alpha = 0.08,
      colour = NA
    ) +
    geom_line(
      aes(y = .data$median_speedup),
      linewidth = 0.35,
      alpha = 0.35,
      lineend = "round"
    ) +
    geom_point(
      aes(y = .data$median_speedup),
      size = 0.85,
      alpha = 0.35
    ) +
    geom_line(
      aes(y = .data$median_speedup_smooth),
      linewidth = 1.55,
      alpha = 1,
      lineend = "round"
    ) +
    scale_colour_manual(values = cols, drop = FALSE) +
    scale_fill_manual(values = cols, drop = FALSE) +
    scale_y_log10(
      breaks = if (!is.null(rs_speedup_overlay_axis_breaks(settings))) {
        rs_speedup_overlay_axis_breaks(settings)
      } else {
        rs_runtime_axis_breaks(speed_df$speedup)
      },
      labels = rs_runtime_axis_labels,
      expand = expansion(mult = c(0.04, 0.08))
    ) +
    rs_observation_x_scale(
      overlay_df,
      show_training_day_axis = show_training_day_axis
    ) +
    labs(
      x = "Number of training observations",
      y = if (settings$comparison == "xemapred_vs_eczemapred") {
        "Median item-level runtime ratio"
      } else {
        settings$ratio_label
      },
      colour = NULL,
      fill = NULL
    ) +
    rs_theme_runtime(base_size = 15)
}


rs_plot_speedup_mean_sd <- function(
    speed_df,
    settings,
    smooth_k = 6L,
    show_training_day_axis = FALSE
) {
  overlay_df <- rs_speedup_overlay_data(speed_df, smooth_k = smooth_k)
  cols <- rs_dataset_cols()

  ggplot(
    overlay_df,
    aes(
      x = .data$N_training_observations,
      colour = .data$Dataset_label,
      fill = .data$Dataset_label,
      group = .data$Dataset_label
    )
  ) +
    geom_hline(
      yintercept = 1,
      linetype = "dashed",
      colour = "grey50",
      linewidth = 0.5
    ) +
    geom_ribbon(
      aes(
        ymin = .data$mean_minus_sd_smooth,
        ymax = .data$mean_plus_sd_smooth
      ),
      alpha = 0.09,
      colour = NA
    ) +
    geom_line(
      aes(y = .data$mean_speedup),
      linewidth = 0.30,
      alpha = 0.25
    ) +
    geom_point(
      aes(y = .data$mean_speedup),
      size = 0.75,
      alpha = 0.25
    ) +
    geom_line(
      aes(y = .data$mean_speedup_smooth),
      linewidth = 1.35,
      lineend = "round"
    ) +
    scale_colour_manual(values = cols, drop = FALSE) +
    scale_fill_manual(values = cols, drop = FALSE) +
    scale_y_log10(
      breaks = if (!is.null(rs_speedup_overlay_axis_breaks(settings))) {
        rs_speedup_overlay_axis_breaks(settings)
      } else {
        rs_runtime_axis_breaks(speed_df$speedup)
      },
      labels = rs_runtime_axis_labels,
      expand = expansion(mult = c(0.04, 0.08))
    ) +
    rs_observation_x_scale(
      overlay_df,
      show_training_day_axis = show_training_day_axis
    ) +
    labs(
      x = "Number of training observations",
      y = settings$ratio_label,
      colour = NULL,
      fill = NULL
    ) +
    rs_theme_runtime(base_size = 15)
}

rs_plot_max_runtime_reduction <- function(
    max_speedup_df,
    smooth_k = 6L,
    show_training_day_axis = FALSE,
    y_label = "Full-update runtime ratio"
) {
  cols <- rs_dataset_cols()

  df <- max_speedup_df %>%
    dplyr::filter(
      is.finite(.data$max_runtime_speedup),
      .data$max_runtime_speedup > 0,
      is.finite(.data$N_training_observations)
    ) %>%
    dplyr::arrange(.data$Dataset, .data$N_training_observations) %>%
    dplyr::group_by(.data$Dataset) %>%
    dplyr::mutate(
      max_runtime_speedup_smooth =
        rs_rolling_median(.data$max_runtime_speedup, k = smooth_k)
    ) %>%
    dplyr::ungroup()

  if (!nrow(df)) return(NULL)

  ggplot(
    df,
    aes(
      x = .data$N_training_observations,
      colour = .data$Dataset_label,
      group = .data$Dataset_label
    )
  ) +
    geom_hline(
      yintercept = 1,
      linetype = "dashed",
      colour = "grey50",
      linewidth = 0.5
    ) +
    geom_line(
      aes(y = .data$max_runtime_speedup),
      linewidth = 0.35,
      alpha = 0.30
    ) +
    geom_point(
      aes(y = .data$max_runtime_speedup),
      size = 0.85,
      alpha = 0.35
    ) +
    geom_line(
      aes(y = .data$max_runtime_speedup_smooth),
      linewidth = 1.45,
      lineend = "round"
    ) +
    scale_colour_manual(values = cols, drop = FALSE) +
    scale_y_log10(
      breaks = c(0.1, 1, 10, 100, 1000, 10000),
      labels = rs_runtime_axis_labels,
      expand = expansion(mult = c(0.04, 0.08))
    ) +
    rs_observation_x_scale(
      df,
      show_training_day_axis = show_training_day_axis
    ) +
    labs(
      x = "Number of training observations",
      y = y_label,
      colour = NULL
    ) +
    rs_theme_runtime(base_size = 15)
}

rs_plot_speedup_boxplot <- function(speed_df, settings) {
  cols <- rs_dataset_cols()

  speed_df %>%
    mutate(Dataset_label = factor(.data$Dataset_label, levels = unname(rs_dataset_labels()))) %>%
    ggplot(aes(x = .data$Dataset_label, y = .data$speedup, colour = .data$Dataset_label, fill = .data$Dataset_label)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey45", linewidth = 0.55) +
    geom_jitter(width = 0.17, alpha = 0.22, size = 0.8, show.legend = FALSE) +
    geom_boxplot(
      width = 0.45,
      alpha = 0.68,
      outlier.shape = NA,
      colour = "grey25",
      linewidth = 0.55,
      show.legend = FALSE
    ) +
    scale_colour_manual(values = cols, drop = FALSE) +
    scale_fill_manual(values = cols, drop = FALSE) +
    scale_y_log10(
      breaks = if (!is.null(rs_speedup_distribution_axis_breaks(settings))) {
        rs_speedup_distribution_axis_breaks(settings)
      } else {
        rs_runtime_axis_breaks(speed_df$speedup)
      },
      labels = rs_runtime_axis_labels,
      name = settings$ratio_label,
      expand = expansion(mult = c(0.04, 0.10))
    ) +
      {
        if (settings$comparison == "xemapred_vs_eczemapred") {
          coord_cartesian(ylim = c(0.3, 1000))
        } else if (settings$comparison == "xemapred_vs_patient_specific") {
          coord_cartesian(ylim = c(0.2, 1000))
        } else {
          NULL
        }
      } +
    labs(x = NULL) +
    rs_theme_runtime(base_size = 16) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(face = "bold", size = 12, colour = "grey20")
    )
}

rs_plot_speedup_per_sign <- function(speed_df, dataset, settings) {
  df <- speed_df %>% filter(.data$Dataset == dataset)
  if (!nrow(df)) return(NULL)

  ggplot(df, aes(x = .data$score_label, y = .data$speedup)) +
    geom_jitter(width = 0.16, alpha = 0.18, size = 0.55, colour = "grey35") +
    geom_boxplot(
      width = 0.5,
      alpha = 0.7,
      fill = settings$high_colour,
      colour = "grey25",
      outlier.shape = NA,
      linewidth = 0.35
    ) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50", linewidth = 0.45) +
    scale_y_log10(
      breaks = if (!is.null(rs_speedup_distribution_axis_breaks(settings))) {
        rs_speedup_distribution_axis_breaks(settings)
      } else {
        rs_runtime_axis_breaks(df$speedup)
      },
      labels = rs_runtime_axis_labels,
      name = settings$ratio_label
    ) +
    labs(
      title = paste0(settings$title, " - ", rs_dataset_labels()[[dataset]]),
      subtitle = "Each point = one PO-SCORAD item at one matched forecasting iteration",
      x = NULL
    ) +
    rs_theme_runtime(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
      legend.position = "none"
    )
}

rs_plot_score_heatmap <- function(speed_summary_by_score, settings, items = rs_default_items()) {
  hm <- speed_summary_by_score %>%
    mutate(
      Dataset_label = factor(.data$Dataset_label, levels = unname(rs_dataset_labels())),
      score_label = rs_score_label(.data$score, items = items)
    )

  ggplot(hm, aes(x = .data$score_label, y = .data$Dataset_label, fill = .data$median_speedup)) +
    geom_tile(colour = "white", linewidth = 0.8) +
    geom_text(aes(label = paste0(round(.data$median_speedup, 1), "x")), fontface = "bold", size = 3.5) +
    scale_fill_gradient(low = "white", high = settings$high_colour, trans = "log10", name = "Median\nreduction") +
    labs(
      title = "Median runtime reduction by severity component",
      subtitle = paste("Values show", settings$heatmap_label),
      x = NULL,
      y = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 15, hjust = 0),
      plot.subtitle = element_text(size = 10, colour = "grey40", hjust = 0),
      axis.text.x = element_text(angle = 30, hjust = 1, face = "bold"),
      axis.text.y = element_text(face = "bold"),
      legend.position = "right"
    )
}

rs_plot_runtime_reduction_heatmap_dataset <- function(speed_df, dataset, settings, items = rs_default_items()) {
  fill_limits <- range(log10(speed_df$speedup), na.rm = TRUE)
  fill_limits <- fill_limits[is.finite(fill_limits)]

  if (length(fill_limits) != 2 || any(!is.finite(fill_limits))) {
    fill_limits <- c(-1, 1)
  }

  df <- speed_df %>%
    filter(.data$Dataset == dataset) %>%
    mutate(
      score_label = factor(as.character(.data$score_label), levels = rev(unname(rs_item_labels()[items]))),
      log10_speedup = log10(.data$speedup)
    ) %>%
    filter(is.finite(.data$log10_speedup))

  if (!nrow(df)) return(NULL)

  ggplot(df, aes(x = .data$iter_plot, y = .data$score_label, fill = .data$log10_speedup)) +
    geom_tile(colour = "white", linewidth = 0.08, height = 0.92) +
    scale_fill_gradient2(
      low = "#3B4CC0",
      mid = "white",
      high = settings$high_colour,
      midpoint = 0,
      name = "Runtime\nreduction",
      breaks = log10(c(0.1, 1, 10, 100, 1000)),
      labels = c("0.1x", "1x", "10x", "100x", "1000x"),
      limits = fill_limits,
      oob = scales::squish
    ) +
    scale_x_continuous(breaks = seq(0, 120, by = 20), expand = expansion(mult = c(0.01, 0.01))) +
    labs(title = rs_dataset_labels()[[dataset]], x = "Training days", y = NULL) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 15, hjust = 0.5, margin = margin(b = 8)),
      axis.title.x = element_text(face = "bold", margin = margin(t = 8)),
      axis.text.x = element_text(size = 10, colour = "grey20"),
      axis.text.y = element_text(face = "bold", size = 10, colour = "grey20"),
      axis.line = element_line(colour = "grey20", linewidth = 0.5),
      legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 9),
      plot.margin = margin(8, 12, 8, 12)
    )
}

rs_plot_runtime_reduction_heatmap_combined <- function(speed_df, settings, items = rs_default_items()) {
  p_derexyl <- rs_plot_runtime_reduction_heatmap_dataset(speed_df, "Derexyl", settings, items = items)
  p_pfdc <- rs_plot_runtime_reduction_heatmap_dataset(speed_df, "PFDC", settings, items = items)

  if (is.null(p_derexyl) || is.null(p_pfdc)) return(NULL)

  (p_derexyl / p_pfdc) +
    plot_layout(ncol = 1, guides = "collect") &
    theme(legend.position = "right")
}

# -----------------------------------------------------------------------
# Saving
# -----------------------------------------------------------------------

rs_save_plot <- function(
    plot,
    out_dir,
    filename_stem,
    width,
    height,
    formats = c("png", "pdf"),
    dpi = 320
) {
  rs_dir_create(out_dir)

  out_files <- character(0)

  for (ext in formats) {
    out_file <- file.path(out_dir, paste0(filename_stem, ".", ext))

    if (ext == "pdf") {
      ggplot2::ggsave(out_file, plot, width = width, height = height, device = cairo_pdf)
    } else {
      ggplot2::ggsave(out_file, plot, width = width, height = height, dpi = dpi)
    }

    out_files <- c(out_files, out_file)
  }

  cat("[SAVED] ", filename_stem, "\n", sep = "")

  invisible(out_files)
}
rs_sec_label <- function(x) {
  dplyr::case_when(
    is.na(x) ~ NA_character_,
    x >= 60 ~ paste0(round(x / 60, 2), " min"),
    x >= 1 ~ paste0(round(x, 2), " s"),
    TRUE ~ paste0(round(1000 * x, 1), " ms")
  )
}

rs_x_label <- function(x, digits = 2) {
  paste0(round(x, digits), "x")
}

rs_pct_label <- function(x, digits = 1) {
  paste0(round(x, digits), "%")
}

rs_print_speedup_summary <- function(speed_summary, settings) {
  out <- speed_summary %>%
    dplyr::transmute(
      Dataset = as.character(.data$Dataset_label),
      Days = .data$n_iters,
      Items = .data$n_items,
      MedianA = rs_sec_label(.data$median_model_a_runtime_sec),
      MedianB = rs_sec_label(.data$median_model_b_runtime_sec),
      Saved = rs_sec_label(.data$median_runtime_saved_sec),
      Speedup = rs_x_label(.data$median_speedup),
      IQR = paste0(rs_x_label(.data$p25_speedup), "-", rs_x_label(.data$p75_speedup)),
      Final = rs_x_label(.data$final_speedup),
      Reduction = rs_pct_label(.data$median_runtime_reduction_pct)
    )

  cat("\n")
  cat(sprintf(
    "%-10s %5s %5s %12s %12s %12s %10s %17s %10s %11s\n",
    "Dataset", "Days", "Items", "Median A", "Median B", "Saved", "Speed-up", "IQR", "Final", "Reduction"
  ))

  apply(out, 1, function(r) {
    cat(sprintf(
      "%-10s %5s %5s %12s %12s %12s %10s %17s %10s %11s\n",
      r[["Dataset"]],
      r[["Days"]],
      r[["Items"]],
      r[["MedianA"]],
      r[["MedianB"]],
      r[["Saved"]],
      r[["Speedup"]],
      r[["IQR"]],
      r[["Final"]],
      r[["Reduction"]]
    ))
  })

  invisible(out)
}

rs_print_per_sign_summary <- function(speed_summary_by_score) {
  out <- speed_summary_by_score %>%
    dplyr::arrange(.data$Dataset_label, dplyr::desc(.data$median_speedup)) %>%
    dplyr::transmute(
      Dataset = as.character(.data$Dataset_label),
      Component = as.character(.data$score_label),
      Speedup = rs_x_label(.data$median_speedup),
      IQR = paste0(rs_x_label(.data$p25_speedup), "-", rs_x_label(.data$p75_speedup)),
      Reduction = rs_pct_label(.data$median_runtime_reduction_pct)
    )

  cat("\n")
  cat(sprintf(
    "%-10s %-12s %10s %17s %11s\n",
    "Dataset", "Component", "Speed-up", "IQR", "Reduction"
  ))

  apply(out, 1, function(r) {
    cat(sprintf(
      "%-10s %-12s %10s %17s %11s\n",
      r[["Dataset"]],
      r[["Component"]],
      r[["Speedup"]],
      r[["IQR"]],
      r[["Reduction"]]
    ))
  })

  invisible(out)
}

# -----------------------------------------------------------------------
# Main runner
# -----------------------------------------------------------------------

run_runtime_speedup_analysis <- function(
    comparison = c(
      "xemapred_vs_eczemapred",
      "xemapred_vs_patient_specific",
      "eczemapred_vs_patient_specific"
    ),
    datasets = rs_dataset_levels(),
    items = rs_default_items(),
    result_root = "results",
    plot_root = "plots",
    horizon_for_counts = 1L,
    smooth_k = 6L,
    formats = c("png", "pdf")
) {
  comparison <- match.arg(comparison)
  settings <- rs_speedup_settings(comparison)

  out_dir <- rs_plot_dir(
    "runtime_speedup",
    settings$out_stub,
    plot_root = plot_root
  )

  csv_dir <- file.path(out_dir, "csv")
  combined_dir <- file.path(out_dir, "combined")
  dataset_dirs <- stats::setNames(file.path(out_dir, datasets), datasets)

  # Remove stale directories created by the previous manuscript/supplementary
  # layout. The renewed analysis writes only to combined/, per-dataset/, csv/.
  legacy_dirs <- file.path(out_dir, c("main", "supplementary"))
  purrr::walk(legacy_dirs, function(path) {
    if (dir.exists(path)) unlink(path, recursive = TRUE, force = TRUE)
  })

  purrr::walk(
    c(csv_dir, combined_dir, unname(dataset_dirs)),
    rs_dir_create
  )

  cat("============================================================\n")
  cat("[RUNTIME SPEED-UP] ", settings$title, "\n", sep = "")
  cat("============================================================\n")

  speed_df <- rs_make_speedup_data(
    comparison = comparison,
    datasets = datasets,
    items = items,
    result_root = result_root,
    horizon_for_counts = horizon_for_counts
  )

  if (!nrow(speed_df)) {
    stop("No matched runtime speed-up rows were produced for comparison: ", comparison)
  }

  # Attach the true accumulated dataset-level training size to every matched
  # runtime row. This changes the trajectory x-axis from training day to the
  # actual number of observed patient-days available for training.
  speed_df <- rs_attach_training_observation_counts(
    speed_df = speed_df,
    datasets = datasets,
    items = items,
    horizon_for_counts = horizon_for_counts
  )

  if (any(!is.finite(speed_df$N_training_observations))) {
    warning("Some runtime rows could not be matched to a training-observation count.")
  }

  speed_summary <- rs_speedup_summary(speed_df)
  speed_summary_by_score <- rs_speedup_summary_by_score(speed_df)
  overlay_df <- rs_speedup_overlay_data(speed_df, smooth_k = smooth_k)
  panel_summary <- rs_panel_summary(overlay_df)

  max_runtime_speedup <- NULL
  max_runtime_checkpoints <- NULL

  if (comparison == "xemapred_vs_eczemapred") {

    max_runtime_speedup <-
      rs_max_runtime_speedup_by_iteration(
        speed_df
      )

    max_runtime_checkpoints <-
      rs_max_runtime_checkpoints(
        max_runtime_speedup
      )

  } else if (
    comparison %in% c(
      "xemapred_vs_patient_specific",
      "eczemapred_vs_patient_specific"
    )
  ) {

    max_runtime_speedup <-
      rs_max_runtime_patient_specific_by_iteration(
        speed_df = speed_df,
        comparison = comparison
      )
  }

  # -------------------------------------------------------------------
  # CSV outputs
  # -------------------------------------------------------------------

  readr::write_csv(
    speed_df,
    file.path(
      csv_dir,
      paste0("runtime_speedup_per_item_iteration_", settings$out_stub, ".csv")
    )
  )

  readr::write_csv(
    speed_summary,
    file.path(
      csv_dir,
      paste0("runtime_speedup_summary_by_dataset_", settings$out_stub, ".csv")
    )
  )

  readr::write_csv(
    speed_summary_by_score,
    file.path(
      csv_dir,
      paste0("runtime_speedup_summary_by_score_", settings$out_stub, ".csv")
    )
  )

  readr::write_csv(
    overlay_df,
    file.path(
      csv_dir,
      paste0("runtime_speedup_overlay_by_training_observations_", settings$out_stub, ".csv")
    )
  )

  readr::write_csv(
    panel_summary,
    file.path(
      csv_dir,
      paste0("runtime_speedup_panel_summary_", settings$out_stub, ".csv")
    )
  )

  heatmap_data <- speed_df %>%
    dplyr::mutate(
      Model_comparison = settings$heatmap_label,
      runtime_reduction_factor = .data$speedup,
      log10_runtime_reduction_factor = log10(.data$speedup)
    )

  readr::write_csv(
    heatmap_data,
    file.path(
      csv_dir,
      paste0("runtime_reduction_heatmap_data_", settings$out_stub, ".csv")
    )
  )

  if (!is.null(max_runtime_speedup)) {

    readr::write_csv(
      max_runtime_speedup,
      file.path(
        csv_dir,
        paste0(
          "runtime_reduction_max_runtime_by_iteration_",
          settings$out_stub,
          ".csv"
        )
      )
    )
  }

  if (!is.null(max_runtime_checkpoints)) {

    readr::write_csv(
      max_runtime_checkpoints,
      file.path(
        csv_dir,
        paste0(
          "runtime_reduction_max_runtime_checkpoints_",
          settings$out_stub,
          ".csv"
        )
      )
    )
  }

  # -------------------------------------------------------------------
  # Combined plots
  # -------------------------------------------------------------------
  # No training-day secondary axis is shown in the combined two-dataset
  # plots because the mapping from observation count to day differs by
  # dataset. Per-dataset plots show training day on the upper axis.

  p_overlay <- rs_plot_speedup_overlay(
    speed_df,
    settings = settings,
    smooth_k = smooth_k,
    show_training_day_axis = FALSE
  )

  rs_save_plot(
    p_overlay,
    combined_dir,
    paste0("runtime_speedup_overlay_", settings$out_stub),
    width = 8.5,
    height = 5.2,
    formats = formats
  )

  p_overlay_raw <- rs_plot_speedup_raw_median_and_rolling(
    speed_df,
    settings = settings,
    smooth_k = smooth_k,
    show_training_day_axis = FALSE
  )

  rs_save_plot(
    p_overlay_raw,
    combined_dir,
    paste0("runtime_speedup_raw_median_and_rolling_", settings$out_stub),
    width = 8.5,
    height = 5.2,
    formats = formats
  )

  if (!is.null(max_runtime_speedup)) {

    p_max_runtime <- rs_plot_max_runtime_speedup(
      max_runtime_speedup,
      settings = settings
    )

    rs_save_plot(
      p_max_runtime,
      combined_dir,
      paste0(
        "runtime_speedup_max_runtime_",
        settings$out_stub
      ),
      width = 8.5,
      height = 5.2,
      formats = formats
    )
  }


  # Optional mean ± SD diagnostic. For runtime ratios the median + IQR plot
  # above is generally more robust and is the preferred manuscript summary.
  p_mean_sd <- rs_plot_speedup_mean_sd(
    speed_df,
    settings = settings,
    smooth_k = smooth_k,
    show_training_day_axis = FALSE
  )

  rs_save_plot(
    p_mean_sd,
    combined_dir,
    paste0("runtime_speedup_mean_sd_", settings$out_stub),
    width = 8.5,
    height = 5.2,
    formats = formats
  )

  # Retained for backward compatibility/diagnostics. This pools iterations
  # with different training sizes, so it should not be the headline figure.
  p_box <- rs_plot_speedup_boxplot(speed_df, settings = settings)

  rs_save_plot(
    p_box,
    combined_dir,
    paste0("runtime_speedup_distribution_boxplot_", settings$out_stub),
    width = 6.8,
    height = 4.4,
    formats = formats
  )

  p_score_heatmap <- rs_plot_score_heatmap(
    speed_summary_by_score,
    settings = settings,
    items = items
  )

  rs_save_plot(
    p_score_heatmap,
    combined_dir,
    paste0("runtime_speedup_score_heatmap_", settings$out_stub),
    width = 8.5,
    height = 3.5,
    formats = formats
  )

  p_heat_combined <- rs_plot_runtime_reduction_heatmap_combined(
    speed_df,
    settings = settings,
    items = items
  )

  if (!is.null(p_heat_combined)) {
    rs_save_plot(
      p_heat_combined,
      combined_dir,
      paste0("runtime_reduction_heatmap_combined_", settings$out_stub),
      width = 10,
      height = 8.5,
      formats = formats
    )
  }

  if (!is.null(max_runtime_speedup)) {
    p_max_combined <- rs_plot_max_runtime_reduction(
      max_runtime_speedup,
      smooth_k = smooth_k,
      show_training_day_axis = FALSE,
      y_label = if (comparison == "xemapred_vs_eczemapred") {
        "Full-update runtime ratio"
      } else {
        "Runtime ratio"
      }
    )

    if (!is.null(p_max_combined)) {
      rs_save_plot(
        p_max_combined,
        combined_dir,
        paste0("runtime_speedup_max_runtime_", settings$out_stub),
        width = 8.5,
        height = 5.2,
        formats = formats
      )
    }
  }

  # -------------------------------------------------------------------
  # Per-dataset plots
  # -------------------------------------------------------------------

  for (dataset in datasets) {
    ds_out_dir <- dataset_dirs[[dataset]]
    ds_speed <- speed_df %>%
      dplyr::filter(as.character(.data$Dataset) == dataset)

    if (!nrow(ds_speed)) next

    p_overlay_ds <- rs_plot_speedup_overlay(
      ds_speed,
      settings = settings,
      smooth_k = smooth_k,
      show_training_day_axis = TRUE
    ) +
      ggplot2::theme(legend.position = "none")

    rs_save_plot(
      p_overlay_ds,
      ds_out_dir,
      paste0("runtime_speedup_overlay_", settings$out_stub, "_", dataset),
      width = 8.5,
      height = 5.2,
      formats = formats
    )

    p_overlay_raw_ds <- rs_plot_speedup_raw_median_and_rolling(
      ds_speed,
      settings = settings,
      smooth_k = smooth_k,
      show_training_day_axis = TRUE
    ) +
      ggplot2::theme(legend.position = "none")

    rs_save_plot(
      p_overlay_raw_ds,
      ds_out_dir,
      paste0("runtime_speedup_raw_median_and_rolling_", settings$out_stub, "_", dataset),
      width = 8.5,
      height = 5.2,
      formats = formats
    )

    p_mean_sd_ds <- rs_plot_speedup_mean_sd(
      ds_speed,
      settings = settings,
      smooth_k = smooth_k,
      show_training_day_axis = TRUE
    ) +
      ggplot2::theme(legend.position = "none")

    rs_save_plot(
      p_mean_sd_ds,
      ds_out_dir,
      paste0("runtime_speedup_mean_sd_", settings$out_stub, "_", dataset),
      width = 8.5,
      height = 5.2,
      formats = formats
    )

    p_sign <- rs_plot_speedup_per_sign(
      speed_df,
      dataset = dataset,
      settings = settings
    )

    if (!is.null(p_sign)) {
      rs_save_plot(
        p_sign,
        ds_out_dir,
        paste0("runtime_speedup_per_sign_", settings$out_stub, "_", dataset),
        width = 8.5,
        height = 4.8,
        formats = formats
      )
    }

    p_heat_ds <- rs_plot_runtime_reduction_heatmap_dataset(
      speed_df,
      dataset = dataset,
      settings = settings,
      items = items
    )

    if (!is.null(p_heat_ds)) {
      rs_save_plot(
        p_heat_ds,
        ds_out_dir,
        paste0("runtime_reduction_heatmap_", settings$out_stub, "_", dataset),
        width = 9,
        height = 4.8,
        formats = formats
      )
    }

    if (!is.null(max_runtime_speedup)) {
      ds_max <- max_runtime_speedup %>%
        dplyr::filter(as.character(.data$Dataset) == dataset)

      p_max_ds <- rs_plot_max_runtime_reduction(
        ds_max,
        smooth_k = smooth_k,
        show_training_day_axis = TRUE,
        y_label = if (comparison == "xemapred_vs_eczemapred") {
          "Full-update runtime ratio"
        } else {
          "Runtime ratio"
        }
      ) +
        ggplot2::theme(legend.position = "none")

      if (!is.null(p_max_ds)) {
        rs_save_plot(
          p_max_ds,
          ds_out_dir,
          paste0("runtime_speedup_max_runtime_", settings$out_stub, "_", dataset),
          width = 8.5,
          height = 5.2,
          formats = formats
        )
      }
    }
  }

  # -------------------------------------------------------------------
  # Console summaries
  # -------------------------------------------------------------------

  cat("\n============================================================\n")
  cat("[SUMMARY] ", settings$title, "\n", sep = "")
  cat("============================================================\n")
  rs_print_speedup_summary(speed_summary, settings)

  cat("\n============================================================\n")
  cat("[PER-COMPONENT SUMMARY]\n")
  cat("============================================================\n")
  rs_print_per_sign_summary(speed_summary_by_score)

  if (!is.null(max_runtime_checkpoints)) {
    cat("\n============================================================\n")
    cat("[MAXIMUM ITEM RUNTIME REDUCTION — ADDITIONAL ANALYSIS]\n")
    cat("============================================================\n\n")
    rs_print_max_runtime_checkpoints(max_runtime_checkpoints)
  }

  cat("\n[DONE]\n")
  cat("[OUT] ", out_dir, "\n", sep = "")

  invisible(list(
    settings = settings,
    speed_df = speed_df,
    speed_summary = speed_summary,
    speed_summary_by_score = speed_summary_by_score,
    overlay_df = overlay_df,
    panel_summary = panel_summary,
    max_runtime_speedup = max_runtime_speedup,
    max_runtime_checkpoints = max_runtime_checkpoints,
    out_dir = out_dir,
    combined_dir = combined_dir,
    dataset_dirs = dataset_dirs
  ))
}

