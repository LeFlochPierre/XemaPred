#!/usr/bin/env Rscript

# -----------------------------------------------------------------------
# Forecasting-performance learning curves
# -----------------------------------------------------------------------

source(here::here("scripts", "03_plotting", "shared", "00_plot_setup.R"))
source(here::here("scripts", "03_plotting", "shared", "plot_config.R"))
source(here::here("scripts", "03_plotting", "shared", "score_metadata.R"))
source(here::here("scripts", "03_plotting", "shared", "plot_paths.R"))
source(here::here("scripts", "03_plotting", "shared", "plot_scales.R"))
source(here::here("scripts", "03_plotting", "shared", "plot_themes.R"))
source(here::here("scripts", "03_plotting", "shared", "save_plots.R"))
source(here::here("scripts", "03_plotting", "shared", "observed_data.R"))
source(here::here("scripts", "03_plotting", "shared", "learning_curve_data.R"))
source(here::here("scripts", "03_plotting", "shared", "learning_curve_plot_helpers.R"))

# -----------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------

DATASETS <- DATASETS_DEFAULT

PRED_HORIZON <- 4
TRAINING_HORIZON <- 4

METRICS <- c(
  "lpd",
  "crps",
  "accuracy_median",
  "accuracy_map",
  "accuracy_prob"
)

# Main comparison.
#
# Use only clean model names:
#   EczemaPred
#   XemaPred
#   PatientSpecificXemaPred
#
# Examples:
#   XemaPred vs EczemaPred:
#     MODELS_TO_PLOT <- c("EczemaPred", "XemaPred")
#
#   Population-level vs patient-specific XemaPred:
#     MODELS_TO_PLOT <- c("XemaPred", "PatientSpecificXemaPred")
#
#   All three:
#     MODELS_TO_PLOT <- c("EczemaPred", "XemaPred", "PatientSpecificXemaPred")

# MODELS_TO_PLOT <- c("EczemaPred", "XemaPred")#, "PatientSpecificXemaPred")
MODELS_TO_PLOT <- c(
  "EczemaPred",
  "XemaPred",
  "PatientSpecificXemaPredCanonical"
)

ADD_REF_MODELS <- TRUE          # cleaner for a before/after comparison
REF_MODELS <- REFERENCE_MODEL_DEFAULTS

# TRUE adds any reference-like models that are actually present in the
# performance tables, for example RW, MC, historical, uniform, AR1,
# MixedAR1, or Smoothing.
ADD_ALL_REF_LIKE_MODELS <- TRUE

RUN_ITEM_LEVEL <- TRUE
RUN_POSCORAD <- TRUE

ITEMS_TO_PLOT <- c(
  "dryness",
  "extent",
  "itching",
  "redness",
  "oozing",
  "sleep",
  "swelling",
  "thickening",
  "scratching"
)

Y_MODE <- "fixed"

# The old item-level plots used log-probability tick labels such as log(0.1).
# The old PO-SCORAD plots often used numeric LPD ticks.
ITEM_LPD_STYLE <- "log_label"
POSCORAD_LPD_STYLE <- "numeric"

LPD_BREAKS <- c(0.01, 0.05, 0.1, 0.25, 0.5, 1)

RESULT_ROOT <- "results"
PLOT_ROOT <- "plots"

PLOT_FORMATS <- c("png")

# -----------------------------------------------------------------------
# Console helpers
# -----------------------------------------------------------------------

hr <- function(char = "=") {
  cat(paste0("\n", paste(rep(char, 70), collapse = ""), "\n"))
}

ok <- function(msg) {
  cat(paste0("[OK]   ", msg, "\n"))
}

info <- function(msg) {
  cat(paste0("[INFO] ", msg, "\n"))
}

warn <- function(msg) {
  cat(paste0("[WARN] ", msg, "\n"))
}

#
# -----------------------------------------------------------------------
# Runner
# -----------------------------------------------------------------------

run_learning_curve_scope <- function(
    item_scope = c("items", "poscorad"),
    datasets = DATASETS,
    metrics = METRICS,
    pred_horizon = PRED_HORIZON,
    training_horizon = TRAINING_HORIZON,
    models_to_plot = MODELS_TO_PLOT,
    add_ref_models = ADD_REF_MODELS,
    ref_models = REF_MODELS,
    add_all_ref_like_models = ADD_ALL_REF_LIKE_MODELS,
    items_to_plot = ITEMS_TO_PLOT,
    y_mode = Y_MODE,
    lpd_style = "log_label",
    lpd_breaks = LPD_BREAKS,
    result_root = RESULT_ROOT,
    plot_root = PLOT_ROOT,
    formats = PLOT_FORMATS,
    save = TRUE
) {
  item_scope <- match.arg(item_scope)

  scope_label <- dplyr::case_when(
    item_scope == "items"    ~ "item_level",
    item_scope == "poscorad" ~ "poscorad",
    TRUE                     ~ item_scope
  )

  hr("=")
  info(paste("Reading", scope_label, "forecasting-performance learning curves"))

  perf <- read_learning_curve_perf(
    datasets = datasets,
    item_scope = item_scope,
    result_root = result_root
  )

  # ---------------------------------------------------------------------
  # Canonical patient-specific XemaPred
  #
  # PFDC    -> cross-cohort prior
  # Derexyl -> population/half-cohort prior
  # ---------------------------------------------------------------------

  perf <- perf %>%
    dplyr::mutate(
      Model_type = dplyr::case_when(

        .data$Dataset == "PFDC" &
          .data$Model_type == "PatientSpecificXemaPredCrossCohort" ~
          "PatientSpecificXemaPredCanonical",

        .data$Dataset == "Derexyl" &
          .data$Model_type == "PatientSpecificXemaPredPopPrior" ~
          "PatientSpecificXemaPredCanonical",

        TRUE ~ .data$Model_type
      )
    )

  if (item_scope == "items") {
    perf <- perf %>%
      filter(.data$Item %in% items_to_plot)
  }

  if (!nrow(perf)) {
    stop("No rows left after applying item-scope filters.")
  }

  model_info <- resolve_learning_curve_models(
    perf = perf,
    requested_models = models_to_plot,
    add_ref_models = add_ref_models,
    ref_models = ref_models,
    add_all_ref_like_models = add_all_ref_like_models
  )

  if (length(model_info$missing_models)) {
    warn(
      paste(
        "Requested models not found and skipped:",
        paste(model_info$missing_models, collapse = ", ")
      )
    )
  }

  perf <- perf %>%
    filter(.data$Model_type %in% model_info$model_order_raw)

  experiment_name <- make_learning_curve_experiment_name(
    main_models = model_info$requested_present,
    include_references = model_info$has_references
  )

  info(
    paste(
      "Models plotted:",
      paste(model_info$model_order_raw, collapse = ", ")
    )
  )

  info(
    paste(
      "Legend labels:",
      paste(learning_curve_model_label(model_info$model_order_raw), collapse = ", ")
    )
  )

  info(
    paste(
      "Output folder:",
      file.path(
        plot_root,
        "learning_curves",
        "forecasting_performance",
        experiment_name,
        scope_label
      )
    )
  )

  for (metric in normalize_learning_curve_metric(metrics)) {
    info(paste("Metric:", metric))

    faceted <- item_scope == "items"
    item_order <- if (faceted) items_to_plot else NULL

    plots <- setNames(
      lapply(datasets, function(dataset) {
        plot_learning_curve(
          perf = perf,
          dataset = dataset,
          pred_horizon = pred_horizon,
          t_horizon = training_horizon,
          metric = metric,
          model_order_raw = model_info$model_order_raw,
          item_order = item_order,
          faceted = faceted,
          y_mode = y_mode,
          lpd_style = lpd_style,
          lpd_breaks = lpd_breaks,
          base_size = if (faceted) 14 else 15
        )
      }),
      datasets
    )

    plots <- plots[!vapply(plots, is.null, logical(1))]

    if (!length(plots)) {
      warn(paste("No plots generated for metric:", metric))
      next
    }

    if (!isTRUE(save)) {
      next
    }

    for (dataset in names(plots)) {
      out_dir <- get_plot_dir(
        "learning_curves",
        "forecasting_performance",
        experiment_name,
        scope_label,
        metric,
        dataset,
        plot_root = plot_root
      )

      filename_stem <- paste0(
        scope_label,
        "_learning_curve_",
        dataset,
        "_",
        metric,
        "_h",
        pred_horizon
      )

      out_file <- file.path(out_dir, paste0(filename_stem, ".", formats[1]))

      save_plot(
        plot = plots[[dataset]],
        out_dir = out_dir,
        filename_stem = filename_stem,
        width = if (faceted) 12 else 12,
        height = if (faceted) 8 else 7,
        dpi = 300,
        formats = formats,
        verbose = FALSE
      )

      ok(paste("Saved:", basename(out_file)))
    }

    if (all(c("Derexyl", "PFDC") %in% names(plots))) {
      combined <- combine_learning_curve_plots(
        plots = plots,
        orientation = if (faceted) "vertical" else "horizontal"
      )

      combined_dir <- get_plot_dir(
        "learning_curves",
        "forecasting_performance",
        experiment_name,
        scope_label,
        metric,
        "COMBINED",
        plot_root = plot_root
      )

      combined_stem <- paste0(
        scope_label,
        "_learning_curve_COMBINED_",
        metric,
        "_h",
        pred_horizon
      )

      combined_file <- file.path(combined_dir, paste0(combined_stem, ".", formats[1]))

      save_plot(
        plot = combined,
        out_dir = combined_dir,
        filename_stem = combined_stem,
        width = if (faceted) 12 else 18,
        height = if (faceted) 16 else 7,
        dpi = 300,
        formats = formats,
        verbose = FALSE
      )

      ok(paste("Saved:", basename(combined_file)))
    }
  }

  invisible(TRUE)
}

main <- function() {
  hr("=")
  cat("FORECASTING-PERFORMANCE LEARNING-CURVE PLOTS\n")
  hr("=")

  if (isTRUE(RUN_ITEM_LEVEL)) {
    run_learning_curve_scope(
      item_scope = "items",
      lpd_style = ITEM_LPD_STYLE
    )
  }

  if (isTRUE(RUN_POSCORAD)) {
    run_learning_curve_scope(
      item_scope = "poscorad",
      lpd_style = POSCORAD_LPD_STYLE
    )
  }

  hr("=")
  ok("Finished forecasting-performance learning-curve plots.")
  hr("=")

  invisible(TRUE)
}

main()