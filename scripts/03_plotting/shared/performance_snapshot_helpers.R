# -----------------------------------------------------------------------
# Shared forecasting-performance snapshot helpers
# -----------------------------------------------------------------------

fps_accuracy_metrics <- function() {
  c("accuracy", "accuracy_median", "accuracy_map", "accuracy_prob")
}

fps_is_accuracy_metric <- function(metric) {
  normalize_learning_curve_metric(metric) %in% fps_accuracy_metrics()
}

fps_metric_label_short <- function(metric) {
  metric <- normalize_learning_curve_metric(metric)

  dplyr::case_when(
    metric == "lpd" ~ "LPD",
    metric == "crps" ~ "CRPS",
    metric == "accuracy_median" ~ "Median accuracy",
    metric == "accuracy_map" ~ "MAP accuracy",
    metric == "accuracy_prob" ~ "Probabilistic accuracy",
    metric == "accuracy" ~ "Point-forecast accuracy",
    TRUE ~ metric
  )
}

fps_metric_axis_label <- function(metric, percent_accuracy = TRUE) {
  metric <- normalize_learning_curve_metric(metric)
  lab <- fps_metric_label_short(metric)

  if (isTRUE(percent_accuracy)) {
    lab <- ifelse(
      fps_is_accuracy_metric(metric),
      paste0(lab, " (%)"),
      lab
    )
  }

  lab
}

fps_model_palette <- function(model_labels) {
  model_labels <- unique(as.character(model_labels))

  palette <- LEARNING_CURVE_MODEL_COLOURS

  missing_labels <- setdiff(model_labels, names(palette))

  if (length(missing_labels)) {
    extra_cols <- grDevices::hcl.colors(length(missing_labels), palette = "Set 2")
    names(extra_cols) <- missing_labels
    palette <- c(palette, extra_cols)
  }

  palette[model_labels]
}

fps_get_snapshot_N <- function(dataset, target_day, horizon = 4) {
  df_dataset <- load_dataset_checked(dataset) %>%
    dplyr::rename(Time = Day)

  fc <- detail_fc_training(df_dataset, horizon)

  id <- which.min(abs(fc$LastTime - target_day))

  tibble::tibble(
    Dataset = dataset,
    target_day = target_day,
    matched_day = fc$LastTime[id],
    matched_N = fc$N[id]
  )
}

fps_get_snapshot_table <- function(datasets, snapshot_day, horizon = 4) {
  purrr::map_dfr(datasets, function(dataset) {
    fps_get_snapshot_N(
      dataset = dataset,
      target_day = snapshot_day[[dataset]],
      horizon = horizon
    )
  })
}

fps_prepare_poscorad_summary_data <- function(
    perf,
    snapshot_tbl,
    metrics = c("lpd", "accuracy_prob", "accuracy_map"),
    model_types = c("PatientSpecificXemaPred", "XemaPred", "EczemaPred"),
    model_order = model_types,
    horizon = 4,
    metric_order = c("lpd", "accuracy_prob", "accuracy_map")
) {
  metrics <- normalize_learning_curve_metric(metrics)
  metric_order <- normalize_learning_curve_metric(metric_order)

  model_label_levels <- learning_curve_model_label(model_order)

  perf_use <- perf %>%
    dplyr::mutate(
      metric = normalize_learning_curve_metric(.data$metric),
      Variable = tolower(trimws(as.character(.data$Variable)))
    ) %>%
    dplyr::filter(
      .data$Horizon == .env$horizon,
      .data$Variable == "fit",
      .data$metric %in% .env$metrics,
      .data$Model_type %in% .env$model_types
    ) %>%
    dplyr::left_join(
      snapshot_tbl %>%
        dplyr::select(Dataset, target_day, matched_day, matched_N),
      by = "Dataset"
    )

  if (!nrow(perf_use)) {
    stop("No POSCORAD performance rows found after filtering.")
  }

  snapshot_values <- perf_use %>%
    dplyr::mutate(
      N_dist = abs(.data$N - .data$matched_N)
    ) %>%
    dplyr::group_by(.data$Dataset, .data$metric, .data$Model_type) %>%
    dplyr::slice_min(order_by = .data$N_dist, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(EvalPoint = "80% training data")

  final_values <- perf_use %>%
    dplyr::group_by(.data$Dataset, .data$metric, .data$Model_type) %>%
    dplyr::slice_max(order_by = .data$N, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      EvalPoint = "100% training data",
      N_dist = abs(.data$N - .data$matched_N)
    )

  out <- dplyr::bind_rows(snapshot_values, final_values) %>%
    dplyr::mutate(
      Dataset_label = dataset_label(.data$Dataset),
      Metric_label = fps_metric_label_short(.data$metric),
      Metric_axis = fps_metric_axis_label(.data$metric, percent_accuracy = TRUE),
      Model_label = learning_curve_model_label(.data$Model_type),

      EvalPoint = factor(
        .data$EvalPoint,
        levels = c("80% training data", "100% training data")
      ),

      Dataset_label = factor(
        .data$Dataset_label,
        levels = unname(DATASET_LABELS)
      ),

      Metric_label = factor(
        .data$Metric_label,
        levels = fps_metric_label_short(metric_order)
      ),

      Metric_axis = factor(
        .data$Metric_axis,
        levels = fps_metric_axis_label(metric_order, percent_accuracy = TRUE)
      ),

      Model_label = factor(
        .data$Model_label,
        levels = model_label_levels
      ),

      is_accuracy = fps_is_accuracy_metric(.data$metric),

      Mean_plot = ifelse(.data$is_accuracy, 100 * .data$Mean, .data$Mean),
      SE_plot = ifelse(.data$is_accuracy, 100 * .data$SE, .data$SE),

      xmin = ifelse(
        .data$is_accuracy,
        pmax(0, 100 * (.data$Mean - .data$SE)),
        .data$Mean - .data$SE
      ),

      xmax = ifelse(
        .data$is_accuracy,
        pmin(100, 100 * (.data$Mean + .data$SE)),
        .data$Mean + .data$SE
      )
    )

  model_index_lookup <- seq_along(model_label_levels)
  names(model_index_lookup) <- model_label_levels

  dataset_base_y <- seq_along(rev(levels(out$Dataset_label)))
  names(dataset_base_y) <- rev(levels(out$Dataset_label))

  out %>%
    dplyr::mutate(
      model_index = unname(model_index_lookup[as.character(.data$Model_label)]),
      n_models = length(model_label_levels),

      # First model in model_order appears highest.
      y_offset = ((.data$n_models + 1) / 2 - .data$model_index) * 0.20,
      y_base = unname(dataset_base_y[as.character(.data$Dataset_label)]),
      y_pos = .data$y_base + .data$y_offset
    )
}

theme_poscorad_summary <- function(base_size = 13) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(
        face = "bold",
        size = 10.8,
        colour = "grey15",
        margin = ggplot2::margin(t = 3, b = 5)
      ),

      legend.position = "top",
      legend.justification = "center",
      legend.box = "horizontal",
      legend.text = ggplot2::element_text(size = 10),
      legend.key.width = grid::unit(0.55, "cm"),
      legend.key.height = grid::unit(0.35, "cm"),
      legend.margin = ggplot2::margin(b = 2),

      axis.text.y = ggplot2::element_text(
        size = 10.5,
        colour = "grey20",
        margin = ggplot2::margin(r = 5)
      ),
      axis.text.x = ggplot2::element_text(
        size = 9.2,
        colour = "grey20"
      ),
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_line(
        colour = "grey35",
        linewidth = 0.30
      ),

      panel.grid.major.x = ggplot2::element_line(
        colour = "grey90",
        linewidth = 0.35
      ),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),

      panel.border = ggplot2::element_rect(
        colour = "grey45",
        fill = NA,
        linewidth = 0.35
      ),

      panel.background = ggplot2::element_rect(
        fill = "white",
        colour = NA
      ),
      plot.background = ggplot2::element_rect(
        fill = "white",
        colour = NA
      ),

      panel.spacing.x = grid::unit(1.0, "lines"),
      panel.spacing.y = grid::unit(0.45, "lines"),

      plot.margin = ggplot2::margin(6, 12, 6, 8)
    )
}

fps_plot_poscorad_summary <- function(
  plot_df,
  base_size = 13,
  dataset_labels = c(
    Derexyl = "Dataset 1",
    PFDC = "Dataset 2"
  )
) {
  model_levels <- levels(plot_df$Model_label)
  plot_cols <- fps_model_palette(model_levels)

  y_breaks <- sort(unique(plot_df$y_base))

  # Original dataset labels in the same order as y_breaks
  y_labels_original <- vapply(
    y_breaks,
    function(y) {
      as.character(plot_df$Dataset_label[match(y, plot_df$y_base)])
    },
    character(1)
  )

  # Replace Derexyl/PFDC with Dataset 1/Dataset 2
  y_labels <- dplyr::recode(
    y_labels_original,
    !!!dataset_labels,
    .default = y_labels_original
  )

  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = .data$Mean_plot,
      y = .data$y_pos,
      colour = .data$Model_label
    )
  ) +
    ggplot2::geom_hline(
      yintercept = y_breaks,
      colour = "grey90",
      linewidth = 0.35
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = .data$xmin,
        xend = .data$xmax,
        y = .data$y_pos,
        yend = .data$y_pos
      ),
      linewidth = 0.75,
      alpha = 0.90,
      lineend = "round"
    ) +
    ggplot2::geom_point(
      size = 3.0,
      shape = 16,
      alpha = 0.98
    ) +
    ggplot2::facet_grid(
      .data$EvalPoint ~ .data$Metric_axis,
      scales = "free_x"
    ) +
    ggplot2::scale_colour_manual(
      values = plot_cols,
      limits = model_levels,
      breaks = model_levels,
      drop = FALSE
    ) +
    ggplot2::scale_y_continuous(
      breaks = y_breaks,
      labels = y_labels,
      limits = range(y_breaks) + c(-0.40, 0.40),
      expand = ggplot2::expansion(mult = c(0.01, 0.01))
    ) +
    ggplot2::labs(
      x = NULL,
      y = NULL,
      colour = NULL
    ) +
    theme_poscorad_summary(base_size = base_size)
}

# -----------------------------------------------------------------------
# Item-level snapshot helpers
# -----------------------------------------------------------------------

fps_item_group <- function(item) {
  dplyr::case_when(
    item == "extent" ~ "Extent",
    item %in% INTENSITY_ITEMS ~ "Intensity signs",
    item %in% c("itching", "sleep") ~ "Subjective symptoms",
    TRUE ~ "Other"
  )
}

fps_model_order_for_item <- function(
    item,
    main_models,
    ref_extent_subjective = c("historical", "RW", "uniform"),
    ref_intensity = c("MC", "historical", "uniform")
) {
  item <- as.character(item)

  if (item %in% c("extent", "itching", "sleep")) {
    return(c(main_models, ref_extent_subjective))
  }

  if (item %in% INTENSITY_ITEMS) {
    return(c(main_models, ref_intensity))
  }

  main_models
}

fps_model_rank_for_item <- function(
    item,
    model_type,
    main_models,
    ref_extent_subjective = c("historical", "RW", "uniform"),
    ref_intensity = c("MC", "historical", "uniform")
) {
  purrr::map2_int(
    as.character(item),
    as.character(model_type),
    function(item_i, model_i) {
      ord <- fps_model_order_for_item(
        item = item_i,
        main_models = main_models,
        ref_extent_subjective = ref_extent_subjective,
        ref_intensity = ref_intensity
      )

      idx <- match(model_i, ord)

      if (is.na(idx)) {
        length(ord) + 999L
      } else {
        idx
      }
    }
  )
}

fps_item_model_levels <- function(
    model_types,
    main_models,
    ref_models = c("MC", "historical", "RW", "uniform")
) {
  fixed_order <- unique(c(main_models, ref_models))
  model_types <- unique(as.character(model_types))

  ordered <- fixed_order[fixed_order %in% model_types]
  extra <- sort(setdiff(model_types, ordered))

  learning_curve_model_label(c(ordered, extra))
}

fps_snapshot_x_scale <- function(
    metric,
    lower,
    upper,
    x_label = NULL,
    lpd_breaks = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1)
) {
  metric <- normalize_learning_curve_metric(metric)

  if (is.null(x_label)) {
    x_label <- fps_metric_label_short(metric)
  }

  if (metric == "lpd") {
    x_range <- range(c(lower, upper, log(lpd_breaks)), na.rm = TRUE)
    x_pad <- diff(x_range) * 0.03

    if (!is.finite(x_pad) || x_pad == 0) {
      x_pad <- 0.1
    }

    return(
      ggplot2::scale_x_continuous(
        breaks = log(lpd_breaks),
        labels = paste0("log(", lpd_breaks, ")"),
        name = x_label,
        limits = c(x_range[1] - x_pad, x_range[2] + x_pad),
        expand = ggplot2::expansion(mult = c(0.01, 0.03))
      )
    )
  }

  if (fps_is_accuracy_metric(metric)) {
    return(
      ggplot2::scale_x_continuous(
        breaks = seq(0, 1, by = 0.2),
        labels = scales::percent_format(accuracy = 1),
        limits = c(0, 1),
        name = x_label,
        expand = ggplot2::expansion(mult = c(0.01, 0.03))
      )
    )
  }

  ggplot2::scale_x_continuous(
    name = x_label,
    expand = ggplot2::expansion(mult = c(0.03, 0.05))
  )
}

fps_prepare_item_snapshot_data <- function(
    perf,
    snapshot_tbl,
    metric,
    model_types,
    main_models,
    items_to_plot,
    horizon = 4,
    ref_extent_subjective = c("historical", "RW", "uniform"),
    ref_intensity = c("MC", "historical", "uniform")
) {
  metric <- normalize_learning_curve_metric(metric)
  is_accuracy <- fps_is_accuracy_metric(metric)

  dat <- perf %>%
    dplyr::mutate(
      metric = normalize_learning_curve_metric(.data$metric),
      Variable = tolower(trimws(as.character(.data$Variable)))
    ) %>%
    dplyr::filter(
      .data$Horizon == .env$horizon,
      .data$metric == .env$metric,
      .data$Variable == "fit",
      .data$Item %in% .env$items_to_plot,
      .data$Model_type %in% .env$model_types
    ) %>%
    dplyr::left_join(
      snapshot_tbl %>%
        dplyr::select(Dataset, target_day, matched_day, matched_N),
      by = "Dataset"
    )

  if (!nrow(dat)) {
    return(list(data = dat, model_levels = character(0)))
  }

  dat <- dat %>%
    dplyr::mutate(
      N_dist = abs(.data$N - .data$matched_N)
    ) %>%
    dplyr::group_by(.data$Dataset, .data$Model_type, .data$Item) %>%
    dplyr::slice_min(order_by = .data$N_dist, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      Dataset_label = dataset_label(.data$Dataset),
      Model_label = learning_curve_model_label(.data$Model_type),
      Item_label = item_label(.data$Item),
      GroupBlock = fps_item_group(.data$Item),

      Lower = if (is_accuracy) {
        pmax(0, .data$Mean - .data$SE)
      } else {
        .data$Mean - .data$SE
      },

      Upper = if (is_accuracy) {
        pmin(1, .data$Mean + .data$SE)
      } else {
        .data$Mean + .data$SE
      },

      is_ref = is_reference_model(.data$Model_type)
    )

  model_levels <- fps_item_model_levels(
    model_types = dat$Model_type,
    main_models = main_models
  )

  dat <- dat %>%
    dplyr::mutate(
      Dataset_label = factor(.data$Dataset_label, levels = unname(DATASET_LABELS)),
      Model_label = factor(.data$Model_label, levels = model_levels),
      GroupBlock = factor(
        .data$GroupBlock,
        levels = c("Extent", "Intensity signs", "Subjective symptoms")
      )
    )

  list(data = dat, model_levels = model_levels)
}

fps_make_item_axis <- function(items_to_plot) {
  default_axis <- tibble::tibble(
    Item = c(
      "extent",
      "thickening", "swelling", "scratching", "redness", "oozing", "dryness",
      "sleep", "itching"
    ),
    Item_y_base = c(
      10,
      8, 7, 6, 5, 4, 3,
      1.5, 0.5
    )
  )

  default_axis %>%
    dplyr::filter(.data$Item %in% items_to_plot) %>%
    dplyr::mutate(
      GroupBlock = factor(
        fps_item_group(.data$Item),
        levels = c("Extent", "Intensity signs", "Subjective symptoms")
      ),
      Item_lab = item_label(.data$Item)
    )
}

fps_plot_item_snapshot <- function(
    snapshot_obj,
    metric,
    items_to_plot,
    main_models,
    ref_extent_subjective = c("historical", "RW", "uniform"),
    ref_intensity = c("MC", "historical", "uniform"),
    base_size = 13
) {
  dat <- snapshot_obj$data
  model_levels <- snapshot_obj$model_levels

  if (!nrow(dat)) {
    return(NULL)
  }

  item_axis <- fps_make_item_axis(items_to_plot)

  dat <- dat %>%
    dplyr::mutate(
      Item = as.character(.data$Item),
      Model_type_chr = as.character(.data$Model_type),
      model_rank_top = fps_model_rank_for_item(
        item = .data$Item,
        model_type = .data$Model_type_chr,
        main_models = main_models,
        ref_extent_subjective = ref_extent_subjective,
        ref_intensity = ref_intensity
      )
    ) %>%
    dplyr::left_join(
      item_axis %>%
        dplyr::select(Item, GroupBlock, Item_y_base, Item_lab),
      by = c("Item", "GroupBlock")
    ) %>%
    dplyr::group_by(.data$Dataset_label, .data$GroupBlock, .data$Item) %>%
    dplyr::arrange(.data$model_rank_top, .data$Model_type_chr, .by_group = TRUE) %>%
    dplyr::mutate(
      model_index = dplyr::row_number(),

      # First model in the requested order appears highest.
      offset = (mean(.data$model_index) - .data$model_index) * 0.13,
      Item_y = .data$Item_y_base + .data$offset
    ) %>%
    dplyr::ungroup()

  row_space_df <- tibble::tribble(
    ~GroupBlock,              ~ymin, ~ymax,
    "Extent",                  9.15, 10.85,
    "Intensity signs",         2.45,  8.55,
    "Subjective symptoms",     0.00,  2.00
  ) %>%
    dplyr::mutate(
      GroupBlock = factor(
        .data$GroupBlock,
        levels = c("Extent", "Intensity signs", "Subjective symptoms")
      )
    ) %>%
    tidyr::expand_grid(
      Dataset_label = unique(dat$Dataset_label)
    ) %>%
    tidyr::pivot_longer(
      cols = c("ymin", "ymax"),
      names_to = "bound",
      values_to = "Item_y"
    ) %>%
    dplyr::mutate(
      x_dummy = min(dat$Lower, na.rm = TRUE)
    )

  plot_cols <- fps_model_palette(model_levels)
  main_labels <- learning_curve_model_label(main_models)

  ggplot2::ggplot(
    dat,
    ggplot2::aes(
      y = .data$Item_y,
      colour = .data$Model_label
    )
  ) +
    ggplot2::geom_blank(
      data = row_space_df,
      ggplot2::aes(x = .data$x_dummy, y = .data$Item_y),
      inherit.aes = FALSE
    ) +
    ggplot2::geom_hline(
    data = item_axis,
    ggplot2::aes(yintercept = .data$Item_y_base),
    colour = "grey90",
    linewidth = 0.35
    ) +
    ggplot2::geom_segment(
    ggplot2::aes(
        x = .data$Lower,
        xend = .data$Upper,
        y = .data$Item_y,
        yend = .data$Item_y
    ),
    linewidth = 0.70,
    alpha = 0.85,
    lineend = "round",
    show.legend = FALSE,
    na.rm = TRUE
    ) +
    ggplot2::geom_point(
    data = dat %>% dplyr::filter(!.data$is_ref),
    ggplot2::aes(x = .data$Mean),
    size = 3.1,
    shape = 16,
    alpha = 0.98,
    na.rm = TRUE
    ) +
    ggplot2::geom_point(
    data = dat %>% dplyr::filter(.data$is_ref),
    ggplot2::aes(x = .data$Mean),
    size = 2.9,
    shape = 21,
    fill = "white",
    stroke = 1.05,
    alpha = 0.95,
    na.rm = TRUE
    ) +
    ggplot2::scale_colour_manual(
      values = plot_cols,
      limits = model_levels,
      breaks = model_levels,
      drop = FALSE
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(
          shape = ifelse(model_levels %in% main_labels, 16, 21),
          fill = ifelse(model_levels %in% main_labels, plot_cols[model_levels], "white"),
          colour = plot_cols[model_levels],
          size = 3,
          alpha = 1,
          stroke = 1.05
        )
      )
    ) +
    fps_snapshot_x_scale(
      metric = metric,
      lower = dat$Lower,
      upper = dat$Upper,
      x_label = fps_metric_label_short(metric)
    ) +
    ggplot2::scale_y_continuous(
      breaks = item_axis$Item_y_base,
      labels = item_axis$Item_lab,
      expand = ggplot2::expansion(mult = c(0.03, 0.03))
    ) +
    ggplot2::facet_grid(
      GroupBlock ~ Dataset_label,
      scales = "free_y",
      space = "free_y"
    ) +
    ggplot2::labs(
      y = NULL,
      colour = NULL
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(colour = "grey88", linewidth = 0.35),

      strip.background = ggplot2::element_rect(fill = "grey96", colour = "grey70", linewidth = 0.35),
      strip.text.x = ggplot2::element_text(face = "bold", size = 11, margin = ggplot2::margin(t = 3, b = 3)),
      strip.text.y.right = ggplot2::element_text(face = "bold", size = 8.5, angle = -90),

      axis.title.x = ggplot2::element_text(face = "bold", size = 12, margin = ggplot2::margin(t = 6)),
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1, size = 8.8, colour = "grey20"),
      axis.text.y = ggplot2::element_text(size = 10, colour = "grey20", margin = ggplot2::margin(r = 4)),
      axis.ticks.y = ggplot2::element_blank(),

      legend.position = "top",
      legend.justification = "center",
      legend.box.just = "center",
      legend.box = "horizontal",
      legend.text = ggplot2::element_text(size = 10),
      legend.key.width = grid::unit(0.55, "cm"),
      legend.key.height = grid::unit(0.35, "cm"),
      legend.spacing.x = grid::unit(3, "pt"),
      legend.margin = ggplot2::margin(b = 2),

      panel.spacing.x = grid::unit(0.25, "lines"),
      panel.spacing.y = grid::unit(0.55, "lines"),
      plot.margin = ggplot2::margin(8, 12, 8, 8),

      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank()
    )
}

fps_save_item_snapshot_numbers <- function(snapshot_obj, out_dir, filename_stem, metric) {
  dat <- snapshot_obj$data

  if (!nrow(dat)) {
    return(invisible(FALSE))
  }

  metric <- normalize_learning_curve_metric(metric)
  value_name <- fps_metric_label_short(metric)
  is_accuracy_metric <- fps_is_accuracy_metric(metric)

  numbers <- dat %>%
    dplyr::mutate(
      Value = .data$Mean,
      Value_percent = if (is_accuracy_metric) 100 * .data$Mean else NA_real_,
      SE_percent = if (is_accuracy_metric) 100 * .data$SE else NA_real_,
      Lower_percent = if (is_accuracy_metric) 100 * .data$Lower else NA_real_,
      Upper_percent = if (is_accuracy_metric) 100 * .data$Upper else NA_real_
    ) %>%
    dplyr::select(
      Dataset,
      Dataset_label,
      Model_type,
      Model_label,
      Item,
      Item_label,
      Horizon,
      !!value_name := Value,
      SE,
      Lower,
      Upper,
      Value_percent,
      SE_percent,
      Lower_percent,
      Upper_percent,
      N,
      LastTime,
      matched_N,
      matched_day,
      target_day,
      N_dist
    ) %>%
    dplyr::arrange(.data$Dataset, .data$Item, .data$Model_type)

  out_file <- file.path(out_dir, paste0(filename_stem, ".csv"))

  readr::write_csv(numbers, out_file)

  invisible(out_file)
}

# -----------------------------------------------------------------------
# PO-SCORAD forecasting-performance comparison tables
# -----------------------------------------------------------------------

fps_comparison_settings <- function(comparison) {
  if (comparison == "xemapred_vs_eczemapred") {
    return(list(
      model_a = "EczemaPred",
      model_b = "XemaPred",
      model_a_label = "EczemaPred",
      model_b_label = "XemaPred",
      comparison_label = "XemaPred vs EczemaPred",
      out_prefix = "xemapred_vs_eczemapred"
    ))
  }

  if (comparison == "patient_specific_vs_xemapred") {
    return(list(
      model_a = "XemaPred",
      model_b = "PatientSpecificXemaPred",
      model_a_label = "Population-level XemaPred",
      model_b_label = "Patient-specific XemaPred",
      comparison_label = "Patient-specific XemaPred vs population-level XemaPred",
      out_prefix = "patient_specific_xemapred_vs_population_xemapred"
    ))
  }

  stop(
    "Unknown comparison: ", comparison,
    "\nUse 'xemapred_vs_eczemapred' or 'patient_specific_vs_xemapred'."
  )
}

fps_higher_is_better <- function(metric) {
  metric <- normalize_learning_curve_metric(metric)
  !metric %in% c("crps", "mae", "mae_median", "mae_map")
}

fps_fmt_mean <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))
}

fps_fmt_acc_pct <- function(x, digits = 1) {
  ifelse(is.na(x), "NA", paste0(formatC(100 * x, format = "f", digits = digits), "%"))
}

fps_fmt_metric_value <- function(metric, x, digits_lpd = 3, digits_acc = 1) {
  metric <- normalize_learning_curve_metric(metric)

  ifelse(
    fps_is_accuracy_metric(metric),
    fps_fmt_acc_pct(x, digits = digits_acc),
    fps_fmt_mean(x, digits = digits_lpd)
  )
}

fps_fmt_metric_difference <- function(metric, x, digits_lpd = 3, digits_acc = 1) {
  metric <- normalize_learning_curve_metric(metric)

  ifelse(
    fps_is_accuracy_metric(metric),
    paste0(formatC(100 * x, format = "f", digits = digits_acc), " percentage points"),
    formatC(x, format = "f", digits = digits_lpd)
  )
}

fps_prepare_performance_values <- function(
    perf,
    datasets,
    metrics,
    model_types,
    ref_models = character(0),
    horizon = 4
) {
  metrics <- normalize_learning_curve_metric(metrics)

  out <- perf %>%
    dplyr::mutate(
      metric = normalize_learning_curve_metric(.data$metric),
      Variable = tolower(trimws(as.character(.data$Variable)))
    ) %>%
    dplyr::filter(
      .data$Dataset %in% .env$datasets,
      .data$Horizon == .env$horizon,
      .data$Variable == "fit",
      .data$metric %in% .env$metrics,
      .data$Model_type %in% c(.env$model_types, .env$ref_models)
    ) %>%
    dplyr::mutate(
      Dataset_label = dataset_label(.data$Dataset),
      Metric_label = fps_metric_label_short(.data$metric),
      Model_label = learning_curve_model_label(.data$Model_type)
    )

  if (!nrow(out)) {
    stop("No PO-SCORAD performance rows found after filtering.")
  }

  out
}

fps_final_values <- function(perf_use) {
  perf_use %>%
    dplyr::group_by(
      .data$Dataset,
      .data$Dataset_label,
      .data$metric,
      .data$Metric_label,
      .data$Model_type,
      .data$Model_label
    ) %>%
    dplyr::slice_max(order_by = .data$N, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(
      Dataset,
      Dataset_label,
      metric,
      Metric_label,
      Model_type,
      Model_label,
      N,
      LastTime,
      Mean,
      SE
    ) %>%
    dplyr::arrange(.data$Dataset, .data$metric, dplyr::desc(.data$Mean))
}

fps_snapshot_values <- function(perf_use, snapshot_tbl) {
  perf_use %>%
    dplyr::left_join(
      snapshot_tbl %>%
        dplyr::select(Dataset, target_day, matched_day, matched_N),
      by = "Dataset"
    ) %>%
    dplyr::mutate(
      N_dist = abs(.data$N - .data$matched_N)
    ) %>%
    dplyr::group_by(.data$Dataset, .data$metric, .data$Model_type) %>%
    dplyr::slice_min(order_by = .data$N_dist, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(
      Dataset,
      Dataset_label,
      metric,
      Metric_label,
      Model_type,
      Model_label,
      target_day,
      matched_day,
      matched_N,
      N,
      LastTime,
      Mean,
      SE,
      N_dist
    ) %>%
    dplyr::arrange(.data$Dataset, .data$metric, dplyr::desc(.data$Mean))
}

fps_make_model_comparison <- function(values_df, settings) {
  model_a <- settings$model_a
  model_b <- settings$model_b

  out <- values_df %>%
    dplyr::filter(.data$Model_type %in% c(model_a, model_b)) %>%
    dplyr::select(
      Dataset,
      Dataset_label,
      metric,
      Metric_label,
      Model_type,
      Mean,
      SE,
      N,
      LastTime,
      dplyr::any_of(c("target_day", "matched_day", "matched_N", "N_dist"))
    ) %>%
    tidyr::pivot_wider(
      names_from = Model_type,
      values_from = c(Mean, SE, N, LastTime),
      names_sep = "_"
    )

  mean_a <- paste0("Mean_", model_a)
  mean_b <- paste0("Mean_", model_b)

  if (!all(c(mean_a, mean_b) %in% names(out))) {
    stop(
      "Could not find both models in comparison table.\n",
      "Missing columns from: ", paste(c(mean_a, mean_b), collapse = ", "), "\n",
      "Available columns: ", paste(names(out), collapse = ", ")
    )
  }

  out %>%
    dplyr::mutate(
      difference = .data[[mean_b]] - .data[[mean_a]],
      higher_is_better = fps_higher_is_better(.data$metric),
      better_model = dplyr::if_else(
        (.data$higher_is_better & .data$difference >= 0) |
          (!.data$higher_is_better & .data$difference <= 0),
        settings$model_b_label,
        settings$model_a_label
      )
    ) %>%
    dplyr::arrange(.data$Dataset, .data$metric)
}

fps_early_comparison <- function(perf_use, settings, n_early = 3) {
  model_a <- settings$model_a
  model_b <- settings$model_b

  perf_use %>%
    dplyr::filter(.data$Model_type %in% c(model_a, model_b)) %>%
    dplyr::group_by(.data$Dataset, .data$Dataset_label, .data$metric, .data$Metric_label) %>%
    dplyr::arrange(.data$N, .by_group = TRUE) %>%
    dplyr::mutate(iter_rank = dplyr::dense_rank(.data$N)) %>%
    dplyr::filter(.data$iter_rank <= n_early) %>%
    dplyr::ungroup() %>%
    dplyr::select(
      Dataset,
      Dataset_label,
      metric,
      Metric_label,
      iter_rank,
      N,
      LastTime,
      Model_type,
      Mean,
      SE
    ) %>%
    tidyr::pivot_wider(
      names_from = Model_type,
      values_from = c(Mean, SE),
      names_sep = "_"
    ) %>%
    dplyr::filter(
      !is.na(.data[[paste0("Mean_", model_a)]]),
      !is.na(.data[[paste0("Mean_", model_b)]])
    ) %>%
    dplyr::mutate(
      difference = .data[[paste0("Mean_", model_b)]] - .data[[paste0("Mean_", model_a)]],
      higher_is_better = fps_higher_is_better(.data$metric),
      model_b_better = (.data$higher_is_better & .data$difference >= 0) |
        (!.data$higher_is_better & .data$difference <= 0)
    )
}

fps_early_summary <- function(early_compare) {
  early_compare %>%
    dplyr::group_by(.data$Dataset, .data$Dataset_label, .data$metric, .data$Metric_label) %>%
    dplyr::summarise(
      n_early_iterations = dplyr::n(),
      n_model_b_better = sum(.data$model_b_better, na.rm = TRUE),
      mean_difference = mean(.data$difference, na.rm = TRUE),
      median_difference = median(.data$difference, na.rm = TRUE),
      first_iteration_difference = .data$difference[which.min(.data$iter_rank)],
      last_early_iteration_difference = .data$difference[which.max(.data$iter_rank)],
      .groups = "drop"
    )
}

fps_early_days_summary <- function(early_compare) {
  early_compare %>%
    dplyr::distinct(.data$Dataset, .data$Dataset_label, .data$iter_rank, .data$N, .data$LastTime) %>%
    dplyr::arrange(.data$Dataset_label, .data$iter_rank) %>%
    dplyr::group_by(.data$Dataset, .data$Dataset_label) %>%
    dplyr::summarise(
      n_iterations = dplyr::n_distinct(.data$iter_rank),
      first_training_day = min(.data$LastTime, na.rm = TRUE),
      last_training_day = max(.data$LastTime, na.rm = TRUE),
      first_training_observations = min(.data$N, na.rm = TRUE),
      last_training_observations = max(.data$N, na.rm = TRUE),
      .groups = "drop"
    )
}

fps_best_reference_comparison <- function(final_values, settings, ref_models) {
  model_b <- settings$model_b

  best_reference <- final_values %>%
    dplyr::filter(.data$Model_type %in% ref_models) %>%
    dplyr::group_by(.data$Dataset, .data$Dataset_label, .data$metric, .data$Metric_label) %>%
    dplyr::arrange(
      dplyr::if_else(
        fps_higher_is_better(.data$metric),
        -.data$Mean,
        .data$Mean
      ),
      .by_group = TRUE
    ) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::rename(
      Best_reference = Model_label,
      Best_reference_mean = Mean,
      Best_reference_SE = SE
    ) %>%
    dplyr::select(
      Dataset,
      Dataset_label,
      metric,
      Metric_label,
      Best_reference,
      Best_reference_mean,
      Best_reference_SE
    )

  final_values %>%
    dplyr::filter(.data$Model_type == model_b) %>%
    dplyr::select(
      Dataset,
      Dataset_label,
      metric,
      Metric_label,
      Model_B_mean = Mean,
      Model_B_SE = SE
    ) %>%
    dplyr::left_join(
      best_reference,
      by = c("Dataset", "Dataset_label", "metric", "Metric_label")
    ) %>%
    dplyr::mutate(
      delta_model_b_minus_best_reference = .data$Model_B_mean - .data$Best_reference_mean,
      higher_is_better = fps_higher_is_better(.data$metric),
      model_b_better_than_best_reference = (.data$higher_is_better &
        .data$delta_model_b_minus_best_reference >= 0) |
        (!.data$higher_is_better & .data$delta_model_b_minus_best_reference <= 0)
    )
}

fps_print_model_comparison <- function(compare_df, settings, title, snapshot = FALSE) {
  model_a <- settings$model_a
  model_b <- settings$model_b

  cat("\n============================================================\n")
  cat(title, "\n", sep = "")
  cat("============================================================\n")

  compare_df %>%
    dplyr::mutate(
      value_b = mapply(fps_fmt_metric_value, .data$metric, .data[[paste0("Mean_", model_b)]]),
      se_b_txt = mapply(fps_fmt_metric_value, .data$metric, .data[[paste0("SE_", model_b)]]),
      value_a = mapply(fps_fmt_metric_value, .data$metric, .data[[paste0("Mean_", model_a)]]),
      se_a_txt = mapply(fps_fmt_metric_value, .data$metric, .data[[paste0("SE_", model_a)]]),
      diff_txt = mapply(fps_fmt_metric_difference, .data$metric, .data$difference),
      eval_txt = if (snapshot && "target_day" %in% names(.)) {
        paste0(" at day ", .data$target_day)
      } else {
        ""
      },
      sentence = paste0(
        .data$Dataset_label, " - ", .data$Metric_label, .data$eval_txt, ": ",
        settings$model_b_label, " = ", .data$value_b, " ± ", .data$se_b_txt,
        ", ", settings$model_a_label, " = ", .data$value_a, " ± ", .data$se_a_txt,
        ", difference = ", .data$diff_txt,
        " (", .data$better_model, " better)."
      )
    ) %>%
    dplyr::pull(.data$sentence) %>%
    cat(sep = "\n")
}

fps_print_early_summary <- function(early_summary, settings, n_early) {
  cat("\n\n============================================================\n")
  cat(settings$comparison_label, ", first ", n_early, " forecasting iterations\n", sep = "")
  cat("============================================================\n")

  early_summary %>%
    dplyr::mutate(
      diff_txt = mapply(fps_fmt_metric_difference, .data$metric, .data$median_difference),
      sentence = paste0(
        .data$Dataset_label, " - ", .data$Metric_label, ": ",
        .data$n_model_b_better, "/", .data$n_early_iterations,
        " early iterations were better with ", settings$model_b_label,
        "; median difference = ", .data$diff_txt, "."
      )
    ) %>%
    dplyr::pull(.data$sentence) %>%
    cat(sep = "\n")
}

fps_print_best_reference <- function(ref_compare, settings) {
  cat("\n\n============================================================\n")
  cat(settings$model_b_label, " vs best reference model, final learning-curve point\n", sep = "")
  cat("============================================================\n")

  ref_compare %>%
    dplyr::mutate(
      value_b = mapply(fps_fmt_metric_value, .data$metric, .data$Model_B_mean),
      se_b = mapply(fps_fmt_metric_value, .data$metric, .data$Model_B_SE),
      value_ref = mapply(fps_fmt_metric_value, .data$metric, .data$Best_reference_mean),
      se_ref = mapply(fps_fmt_metric_value, .data$metric, .data$Best_reference_SE),
      diff_txt = mapply(
        fps_fmt_metric_difference,
        .data$metric,
        .data$delta_model_b_minus_best_reference
      ),
      sentence = paste0(
        .data$Dataset_label, " - ", .data$Metric_label, ": ",
        settings$model_b_label, " = ", .data$value_b, " ± ", .data$se_b,
        "; best reference = ", .data$Best_reference, " (",
        .data$value_ref, " ± ", .data$se_ref, ")",
        "; difference = ", .data$diff_txt, "."
      )
    ) %>%
    dplyr::pull(.data$sentence) %>%
    cat(sep = "\n")
}

fps_fmt_table_cell <- function(metric, mean, se, digits_lpd = 3, digits_acc = 1) {
  metric <- normalize_learning_curve_metric(metric)

  if (fps_is_accuracy_metric(metric)) {
    return(paste0(
      formatC(100 * mean, format = "f", digits = digits_acc),
      " ± ",
      formatC(100 * se, format = "f", digits = digits_acc)
    ))
  }

  paste0(
    formatC(mean, format = "f", digits = digits_lpd),
    " ± ",
    formatC(se, format = "f", digits = digits_lpd)
  )
}

fps_three_model_80_100_table <- function(
    perf_use,
    snapshot_tbl,
    metrics = c("lpd", "accuracy_map", "accuracy_prob"),
    models = c("EczemaPred", "XemaPred", "PatientSpecificXemaPred")
) {
  metrics <- normalize_learning_curve_metric(metrics)

  final_values <- perf_use %>%
    dplyr::filter(
      .data$metric %in% .env$metrics,
      .data$Model_type %in% .env$models
    ) %>%
    dplyr::group_by(.data$Dataset, .data$metric, .data$Model_type) %>%
    dplyr::slice_max(order_by = .data$N, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(`Training data (%)` = "100")

  snapshot_values <- perf_use %>%
    dplyr::filter(
      .data$metric %in% .env$metrics,
      .data$Model_type %in% .env$models
    ) %>%
    dplyr::left_join(
      snapshot_tbl %>%
        dplyr::select(Dataset, target_day, matched_day, matched_N),
      by = "Dataset"
    ) %>%
    dplyr::mutate(
      N_dist = abs(.data$N - .data$matched_N)
    ) %>%
    dplyr::group_by(.data$Dataset, .data$metric, .data$Model_type) %>%
    dplyr::slice_min(order_by = .data$N_dist, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(`Training data (%)` = "80")

  out <- dplyr::bind_rows(final_values, snapshot_values) %>%
    dplyr::mutate(
      Dataset = dplyr::recode(
        as.character(.data$Dataset),
        Derexyl = "1",
        PFDC = "2",
        .default = as.character(.data$Dataset)
      ),
      Metric = fps_metric_label_short(.data$metric),
      Metric = dplyr::recode(
        .data$Metric,
        "MAP accuracy" = "MAP accuracy",
        "Probabilistic accuracy" = "Probabilistic accuracy",
        .default = .data$Metric
      ),
      Model_column = dplyr::case_when(
        .data$Model_type == "EczemaPred" ~ "EczemaPred",
        .data$Model_type == "XemaPred" ~ "XemaPred, population-level",
        .data$Model_type == "PatientSpecificXemaPred" ~ "XemaPred, Patient-specific",
        TRUE ~ learning_curve_model_label(.data$Model_type)
      ),
      Value = mapply(
        fps_fmt_table_cell,
        .data$metric,
        .data$Mean,
        .data$SE
      )
    ) %>%
    dplyr::select(
      Dataset,
      `Training data (%)`,
      Metric,
      Model_column,
      Value
    ) %>%
    tidyr::pivot_wider(
      names_from = Model_column,
      values_from = Value
    ) %>%
    dplyr::mutate(
      Dataset = factor(.data$Dataset, levels = c("1", "2")),
      `Training data (%)` = factor(.data$`Training data (%)`, levels = c("100", "80")),
      Metric = factor(
        .data$Metric,
        levels = c("LPD", "MAP accuracy", "Probabilistic accuracy")
      )
    ) %>%
    dplyr::arrange(
      .data$`Training data (%)`,
      .data$Dataset,
      .data$Metric
    ) %>%
    dplyr::mutate(
      Dataset = as.character(.data$Dataset),
      `Training data (%)` = as.character(.data$`Training data (%)`),
      Metric = as.character(.data$Metric)
    )

  expected_cols <- c(
    "Dataset",
    "Training data (%)",
    "Metric",
    "EczemaPred",
    "XemaPred, population-level",
    "XemaPred, Patient-specific"
  )

  missing_cols <- setdiff(expected_cols, names(out))

  for (col in missing_cols) {
    out[[col]] <- NA_character_
  }

  out %>%
    dplyr::select(dplyr::all_of(expected_cols))
}