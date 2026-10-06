# -----------------------------------------------------------------------
# Shared learning-curve plotting helpers
# -----------------------------------------------------------------------

LEARNING_CURVE_MODEL_COLOURS <- c(
  "EczemaPred"                       = "#0072B2",
  "XemaPred"                         = "#d66309",
  "Patient-specific XemaPred"        = "#556B2F",
  "Patient-specific XemaPred (non-power prior)" = "#076875",
  "Random walk"  = "#56B4E9",
  "Markov chain" = "#009E73",
  "Historical"   = "#808080",
  "Uniform"      = "#000000",
  "AR(1)"        = "#CC79A7",
  "Mixed AR(1)"  = "#999933",
  "Smoothing"    = "#8C564B"
)

learning_curve_model_scales <- function(model_labels) {
  model_labels <- unique(as.character(model_labels))

  palette <- LEARNING_CURVE_MODEL_COLOURS

  missing_labels <- setdiff(model_labels, names(palette))

  if (length(missing_labels)) {
    extra_cols <- grDevices::hcl.colors(length(missing_labels), palette = "Set 2")
    names(extra_cols) <- missing_labels
    palette <- c(palette, extra_cols)
  }

  values <- palette[model_labels]

  list(
    scale_colour_manual(values = values, drop = FALSE),
    scale_fill_manual(values = values, drop = FALSE)
  )
}

learning_curve_model_group <- function(model_label) {
  ifelse(
    grepl("XemaPred|EczemaPred", as.character(model_label)),
    "main",
    "reference"
  )
}

learning_curve_y_scale <- function(
    metric,
    y_mode = c("fixed", "free"),
    lpd_style = c("log_label", "numeric"),
    lpd_breaks = c(0.01, 0.05, 0.1, 0.25, 0.5, 1)
) {
  y_mode <- match.arg(y_mode)
  lpd_style <- match.arg(lpd_style)
  metric <- normalize_learning_curve_metric(metric)

  y_label <- learning_curve_metric_label(metric)

  if (metric == "lpd") {
    if (lpd_style == "log_label") {
      return(
        scale_y_continuous(
          breaks = log(lpd_breaks),
          labels = paste0("log(", lpd_breaks, ")"),
          limits = if (y_mode == "fixed") c(NA, 0.2) else NULL,
          expand = ggplot2::expansion(mult = c(0.03, 0.08)),
          name = y_label
        )
      )
    }

    return(
      scale_y_continuous(
        breaks = -1:-5,
        labels = as.character(-1:-5),
        limits = if (y_mode == "fixed") c(-5.5, -0.5) else NULL,
        expand = ggplot2::expansion(mult = c(0.03, 0.08)),
        name = y_label
      )
    )
  }

  if (metric == "crps") {
    return(
      scale_y_continuous(
        name = y_label,
        expand = ggplot2::expansion(mult = c(0.03, 0.08))
      )
    )
  }

  if (metric %in% c("accuracy", "accuracy_median", "accuracy_map", "accuracy_prob")) {
    return(
      scale_y_continuous(
        limits = if (y_mode == "fixed") c(0, 1) else NULL,
        breaks = c(0, 0.25, 0.5, 0.75, 1),
        labels = c("0", "0.25", "0.50", "0.75", "1.0"),
        expand = ggplot2::expansion(mult = c(0.02, 0.04)),
        name = y_label
      )
    )
  }

  scale_y_continuous(name = y_label)
}

make_training_axis_breaks <- function(dataset, t_horizon = 4, n_breaks = 10) {
  if (exists("detail_fc_training", mode = "function")) {
    df_dataset <- load_dataset_checked(dataset) %>%
      rename(Time = .data$Day)

    fc_it <- detail_fc_training(df_dataset, t_horizon)

    required_cols <- c("N", "LastTime", "Proportion")
    missing_cols <- setdiff(required_cols, names(fc_it))

    if (length(missing_cols)) {
      stop(
        "detail_fc_training() output is missing columns: ",
        paste(missing_cols, collapse = ", ")
      )
    }

    ids <- vapply(
      seq(0, 1, length.out = n_breaks),
      function(x) which.min((x - fc_it$Proportion)^2),
      numeric(1)
    )

    return(
      tibble::tibble(
        N = fc_it$N[ids],
        LastTime = fc_it$LastTime[ids]
      ) %>%
        distinct(.data$N, .data$LastTime)
    )
  }

  warning(
    "detail_fc_training() was not found. ",
    "Using an approximate secondary training-day axis."
  )

  df_dataset <- load_dataset_checked(dataset)

  day_counts <- df_dataset %>%
    filter(!is.na(.data$Day)) %>%
    group_by(.data$Day) %>%
    summarise(n_day = n(), .groups = "drop") %>%
    arrange(.data$Day) %>%
    mutate(
      N = cumsum(.data$n_day),
      LastTime = .data$Day,
      Proportion = (.data$N - min(.data$N)) / (max(.data$N) - min(.data$N))
    )

  ids <- vapply(
    seq(0, 1, length.out = n_breaks),
    function(x) which.min((x - day_counts$Proportion)^2),
    numeric(1)
  )

  day_counts %>%
    slice(ids) %>%
    transmute(
      N = .data$N,
      LastTime = .data$LastTime
    ) %>%
    distinct(.data$N, .data$LastTime)
}

prepare_learning_curve_plot_data <- function(
    perf,
    dataset,
    pred_horizon,
    metric,
    model_order_raw,
    item_order = NULL
) {
  metric <- normalize_learning_curve_metric(metric)

  dat <- perf %>%
    filter(
      .data$Dataset == .env$dataset,
      .data$Horizon == .env$pred_horizon,
      .data$metric == .env$metric
    )

  if ("Variable" %in% names(dat)) {
    dat <- dat %>%
      filter(tolower(trimws(as.character(.data$Variable))) == "fit")
  }

  if (!is.null(item_order)) {
    dat <- dat %>%
      filter(.data$Item %in% item_order)
  }

  if (!nrow(dat)) {
    return(dat)
  }

  dat <- dat %>%
    mutate(
      Dataset_label = dataset_label(.data$Dataset),
      Model_label = learning_curve_model_label(.data$Model_type)
    )

  requested_label_order <- unique(learning_curve_model_label(model_order_raw))
  present_labels <- unique(as.character(dat$Model_label))

  model_label_order <- requested_label_order[requested_label_order %in% present_labels]
  model_label_order <- c(model_label_order, setdiff(present_labels, model_label_order))

  dat <- dat %>%
    filter(.data$Model_label %in% model_label_order) %>%
    mutate(
      Model_label = factor(.data$Model_label, levels = model_label_order),
      Model_group = learning_curve_model_group(.data$Model_label)
    )

  if (!is.null(item_order)) {
    item_levels <- unname(ITEM_LABELS[item_order])

    dat <- dat %>%
      mutate(
        Item_label = item_label(.data$Item),
        Item_label = factor(.data$Item_label, levels = item_levels)
      )
  } else {
    dat <- dat %>%
      mutate(
        Item_label = item_label(.data$Item),
        Item_label = factor(.data$Item_label, levels = unique(.data$Item_label))
      )
  }

  dat
}

plot_learning_curve <- function(
    perf,
    dataset,
    pred_horizon = 4,
    t_horizon = 4,
    metric = "lpd",
    model_order_raw,
    item_order = NULL,
    faceted = TRUE,
    y_mode = c("fixed", "free"),
    lpd_style = c("log_label", "numeric"),
    lpd_breaks = c(0.01, 0.05, 0.1, 0.25, 0.5, 1),
    base_size = 14
) {
  y_mode <- match.arg(y_mode)
  lpd_style <- match.arg(lpd_style)

  dat <- prepare_learning_curve_plot_data(
    perf = perf,
    dataset = dataset,
    pred_horizon = pred_horizon,
    metric = metric,
    model_order_raw = model_order_raw,
    item_order = item_order
  )

  if (!nrow(dat)) {
    return(NULL)
  }

  training_axis <- make_training_axis_breaks(
    dataset = dataset,
    t_horizon = t_horizon
  )

  model_levels <- levels(dat$Model_label)

  p <- ggplot(
    dat,
    aes(
      x = .data$N,
      y = .data$Mean,
      ymin = .data$Mean - .data$SE,
      ymax = .data$Mean + .data$SE,
      colour = .data$Model_label,
      fill = .data$Model_label,
      group = .data$Model_label
    )
  ) +
    geom_ribbon(
      data = subset(dat, Model_group == "reference"),
      alpha = 0.08,
      colour = NA
    ) +
    geom_line(
      data = subset(dat, Model_group == "reference"),
      linewidth = 0.55,
      alpha = 0.65
    ) +
    geom_point(
      data = subset(dat, Model_group == "reference"),
      size = 1.15,
      shape = 17,
      alpha = 0.65
    ) +
    geom_ribbon(
      data = subset(dat, Model_group == "main"),
      alpha = 0.35,
      colour = NA
    ) +
    geom_line(
      data = subset(dat, Model_group == "main"),
      linewidth = 1.20,
      alpha = 1.00
    ) +
    geom_point(
      data = subset(dat, Model_group == "main"),
      size = 1.80,
      shape = 16,
      alpha = 1.00
    ) +
    {
      sc <- learning_curve_model_scales(model_levels)
      list(sc[[1]], sc[[2]])
    } +
    scale_x_continuous(
      sec.axis = dup_axis(
        breaks = training_axis$N,
        labels = training_axis$LastTime,
        name = "Training days"
      )
    ) +
    learning_curve_y_scale(
      metric = metric,
      y_mode = y_mode,
      lpd_style = lpd_style,
      lpd_breaks = lpd_breaks
    ) +
    labs(
      title = dataset_label(dataset),
      x = "Number of training observations",
      colour = NULL,
      fill = NULL
    ) +
    theme_learning_curve(base_size = base_size)

  if (isTRUE(faceted)) {
    p <- p +
      facet_wrap(
        ~Item_label,
        scales = "free_y",
        strip.position = "top"
      )
  }

  p
}

collect_learning_curve_labels <- function(plots) {
  labels <- unlist(lapply(plots, function(plot_obj) {
    dat <- plot_obj$data

    if (!"Model_label" %in% names(dat)) {
      return(character(0))
    }

    if (is.factor(dat$Model_label)) {
      return(levels(droplevels(dat$Model_label)))
    }

    unique(as.character(dat$Model_label))
  }))

  unique(as.character(labels))
}

build_learning_curve_legend_plot <- function(model_labels) {
  model_labels <- unique(as.character(model_labels))

  line_df <- tidyr::expand_grid(
    Model_label = factor(model_labels, levels = model_labels),
    x = c(1, 2)
  ) %>%
    mutate(
      y = 1,
      Model_group = learning_curve_model_group(.data$Model_label)
    )

  point_df <- tibble::tibble(
    Model_label = factor(model_labels, levels = model_labels),
    x = 1.5,
    y = 1,
    Model_group = learning_curve_model_group(model_labels)
  )

  ggplot() +
    geom_line(
      data = subset(line_df, Model_group == "reference"),
      aes(
        x = .data$x,
        y = .data$y,
        colour = .data$Model_label,
        group = .data$Model_label
      ),
      linewidth = 0.55,
      alpha = 0.65
    ) +
    geom_point(
      data = subset(point_df, Model_group == "reference"),
      aes(
        x = .data$x,
        y = .data$y,
        colour = .data$Model_label
      ),
      size = 1.15,
      shape = 17,
      alpha = 0.65
    ) +
    geom_line(
      data = subset(line_df, Model_group == "main"),
      aes(
        x = .data$x,
        y = .data$y,
        colour = .data$Model_label,
        group = .data$Model_label
      ),
      linewidth = 1.20,
      alpha = 1.00
    ) +
    geom_point(
      data = subset(point_df, Model_group == "main"),
      aes(
        x = .data$x,
        y = .data$y,
        colour = .data$Model_label
      ),
      size = 1.65,
      shape = 16,
      alpha = 1.00
    ) +
    {
      scales <- learning_curve_model_scales(model_labels)
      scales[[1]]
    } +
    guides(
      colour = guide_legend(
        nrow = 1,
        byrow = TRUE
      )
    ) +
    theme_void(base_size = 14) +
    theme(
      legend.position = "top",
      legend.direction = "horizontal",
      legend.key = element_rect(
        fill = "white",
        colour = "grey65",
        linewidth = 0.35
      ),
      legend.key.width = grid::unit(0.95, "cm"),
      legend.key.height = grid::unit(0.45, "cm"),
      legend.margin = margin(b = 2),
      legend.box.margin = margin(b = 2),
      legend.spacing.x = grid::unit(0.15, "cm")
    )
}

get_present_model_labels_from_plot <- function(plot_obj) {
  dat <- plot_obj$data

  if (!"Model_label" %in% names(dat)) {
    return(character(0))
  }

  present <- unique(as.character(dat$Model_label))
  present <- present[!is.na(present)]

  if (is.factor(dat$Model_label)) {
    level_order <- levels(dat$Model_label)
    return(level_order[level_order %in% present])
  }

  present
}

combine_learning_curve_plots <- function(
    plots,
    orientation = c("vertical", "horizontal")
) {
  orientation <- match.arg(orientation)

  if (!length(plots)) {
    return(NULL)
  }

  present_labels_by_plot <- lapply(plots, get_present_model_labels_from_plot)
  n_labels_by_plot <- vapply(present_labels_by_plot, length, integer(1))

  legend_source_id <- which.max(n_labels_by_plot)
  legend_source <- plots[[legend_source_id]]

  legend_plot <- legend_source +
    theme(legend.position = "top")

  legend <- cowplot::get_plot_component(
    legend_plot,
    "guide-box-top",
    return_all = TRUE
  )

  if (is.list(legend) && length(legend) == 1) {
    legend <- legend[[1]]
  }

  if (is.null(legend) || inherits(legend, "zeroGrob")) {
    legend <- cowplot::get_legend(legend_plot)
  }

  plots_no_legend <- lapply(plots, function(p) {
    p +
      theme(
        legend.position = "none",
        plot.margin = margin(6, 16, 6, 6)
      )
  })

  ncol <- if (orientation == "horizontal") length(plots_no_legend) else 1

  # Important:
  # Use patchwork for the plot body.
  # cowplot::plot_grid(..., align = "hv", axis = "tbl") can segfault with
  # secondary top axes and many facets.
  body <- patchwork::wrap_plots(
    plots_no_legend,
    ncol = ncol
  )

  legend_height <- if (orientation == "vertical") 0.07 else 0.16

  cowplot::plot_grid(
    legend,
    body,
    ncol = 1,
    rel_heights = c(legend_height, 1 - legend_height)
  )
}