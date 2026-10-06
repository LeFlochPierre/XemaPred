#!/usr/bin/env Rscript
# plot_patient_forecasted_trajectories.R

# -------------------------------------------------------------------------
# Patient-level all-item forecast plots
#
# Purpose:
# - Plot observed scores and predictive distributions for all PO-SCORAD items.
# - Optionally compare multiple models side by side.
#
# Output:
# plots/patient_examples/forecast_all_items/<dataset>/
#   patient_<id>_all_items_forecasts_<models>_H<horizon>.png
#   patient_<id>_all_items_forecasts_<models>_H<horizon>.pdf
#   patient_<id>_all_items_forecasts_<models>_H<horizon>.svg
# -------------------------------------------------------------------------

source(here::here("scripts", "03_plotting", "shared", "00_plot_setup.R"))
source(here::here("scripts", "03_plotting", "shared", "plot_config.R"))
source(here::here("scripts", "03_plotting", "shared", "score_metadata.R"))
source(here::here("scripts", "03_plotting", "shared", "plot_paths.R"))
source(here::here("scripts", "03_plotting", "shared", "plot_scales.R"))
source(here::here("scripts", "03_plotting", "shared", "plot_themes.R"))
source(here::here("scripts", "03_plotting", "shared", "save_plots.R"))
source(here::here("scripts", "03_plotting", "shared", "observed_data.R"))
source(here::here("scripts", "03_plotting", "shared", "prediction_data.R"))


# -------------------------------------------------------------------------
# Helper: print full nested errors from purrr/rlang
# -------------------------------------------------------------------------

print_error_chain <- function(e) {
  message(conditionMessage(e))

  parent <- e$parent
  level <- 1L

  while (!is.null(parent)) {
    message(
      paste0(
        "  Caused by [",
        level,
        "]: ",
        conditionMessage(parent)
      )
    )

    parent <- parent$parent
    level <- level + 1L
  }

  invisible(NULL)
}


# -------------------------------------------------------------------------
# Helper: clean graphics devices and memory between patients
# -------------------------------------------------------------------------

clean_plot_memory <- function() {
  try(grDevices::graphics.off(), silent = TRUE)
  invisible(gc(verbose = FALSE))
}


# -------------------------------------------------------------------------
# Display labels for model names
# -------------------------------------------------------------------------

model_display_label <- function(model_type) {
  get_model_display_label(model_type)
}


# -------------------------------------------------------------------------
# Plot one item panel
# -------------------------------------------------------------------------

plot_one_forecast_panel <- function(df_item, dataset, show_x = FALSE) {
  item <- unique(as.character(df_item$Item))[1]

  if (is.na(item) || length(item) == 0) {
    stop("Cannot plot empty item panel.")
  }

  ylims <- get_item_ylim(item, padded = FALSE)

  x_interval <- get_x_break_interval(
    dataset,
    plot_type = "patient_forecast"
  )

  x_min <- min(df_item$Time, na.rm = TRUE)
  x_max <- max(df_item$Time, na.rm = TRUE)

  if (!is.finite(x_min) || !is.finite(x_max)) {
    stop("Invalid x-axis range for item: ", item)
  }

  x_breaks <- seq(
    from = floor(x_min / x_interval) * x_interval,
    to = ceiling(x_max / x_interval) * x_interval,
    by = x_interval
  )

  p <- ggplot(df_item, aes(x = Time)) +
    geom_ribbon(
      aes(
        ymin = q05,
        ymax = q95,
        fill = "90% predictive interval"
      ),
      alpha = 0.90
    ) +
    geom_ribbon(
      aes(
        ymin = q10,
        ymax = q90,
        fill = "80% predictive interval"
      ),
      alpha = 0.80
    ) +
    geom_ribbon(
      aes(
        ymin = q25,
        ymax = q75,
        fill = "50% predictive interval"
      ),
      alpha = 0.55
    ) +
    geom_line(
      aes(y = pred_median, colour = "Median forecast"),
      linewidth = 0.75
    ) +
    geom_line(
      aes(y = Score, colour = "Observed score"),
      linewidth = 0.60
    ) +
    geom_point(
      aes(y = Score, colour = "Observed score"),
      size = 1.10
    ) +
    scale_fill_manual(
      values = c(
        "90% predictive interval" = "#DCEEFF",
        "80% predictive interval" = "#A9D6FF",
        "50% predictive interval" = "#5DADE2"
      )
    ) +
    scale_colour_manual(
      values = c(
        "Median forecast" = "#1565C0",
        "Observed score" = "black"
      )
    ) +
    scale_y_continuous(
      limits = ylims,
      name = unique(df_item$Item_label)
    ) +
    scale_x_continuous(
      breaks = x_breaks,
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    labs(
      x = if (show_x) "Day" else NULL,
      colour = NULL,
      fill = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(
      text = element_text(family = "sans"),
      axis.title.y = element_text(
        size = 10,
        face = "bold",
        angle = 90,
        margin = margin(r = 8)
      ),
      axis.text.y = element_text(size = 8, colour = "black"),
      axis.text.x = element_text(size = 8, colour = "black"),
      axis.title.x = element_text(
        size = 10,
        colour = "black",
        margin = margin(t = 6)
      ),
      axis.line.x = element_line(colour = "black", linewidth = 0.35),
      axis.line.y = element_line(colour = "black", linewidth = 0.35),
      axis.ticks.x = element_line(colour = "black", linewidth = 0.35),
      axis.ticks.y = element_line(colour = "black", linewidth = 0.35),
      axis.ticks.length = unit(2.5, "pt"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(colour = "grey88", linewidth = 0.25),
      panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.25),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.text = element_text(size = 8),
      legend.key.width = unit(1.2, "cm"),
      plot.margin = margin(2, 4, 2, 4)
    )

  if (!show_x) {
    p <- p +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank()
      )
  }

  p
}


# -------------------------------------------------------------------------
# Plot one model column for one patient
# -------------------------------------------------------------------------

plot_patient_model_all_items <- function(
    dataset,
    patient_id,
    model_type,
    horizon = 4,
    result_root = "results"
) {
  df_all <- load_patient_all_item_predictions(
    dataset = dataset,
    patient_id = patient_id,
    horizon = horizon,
    model_type = model_type,
    result_root = result_root
  )

  if (nrow(df_all) == 0) {
    stop(
      "No prediction rows available after filtering. ",
      "dataset=", dataset,
      " | patient=", patient_id,
      " | model=", model_type
    )
  }

  item_levels <- levels(df_all$Item_label)

  panels <- vector("list", length(item_levels))

  for (i in seq_along(item_levels)) {
    item_label_i <- item_levels[[i]]

    df_item <- df_all %>%
      filter(Item_label == item_label_i)

    if (nrow(df_item) == 0) {
      stop(
        "No rows available for item panel: ",
        item_label_i,
        " | dataset=",
        dataset,
        " | patient=",
        patient_id,
        " | model=",
        model_type
      )
    }

    panels[[i]] <- plot_one_forecast_panel(
      df_item = df_item,
      dataset = dataset,
      show_x = i == length(item_levels)
    )
  }

  panel_heights <- rep(0.8, length(item_levels))

  if (length(panel_heights) >= 1) {
    panel_heights[1] <- 1.1
    panel_heights[length(panel_heights)] <- 1.1
  }

  patchwork::wrap_plots(
    panels,
    ncol = 1,
    heights = panel_heights
  ) +
    patchwork::plot_annotation(
      title = model_display_label(model_type)
    ) &
    theme(
      text = element_text(family = "sans"),
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        size = 12
      )
    )
}


# -------------------------------------------------------------------------
# Save one patient forecast comparison
# -------------------------------------------------------------------------

save_patient_all_items_forecast <- function(
    dataset,
    patient_id,
    horizon = 4,
    model_types = c("EczemaPred", "XemaPred"),
    result_root = "results",
    plot_root = "plots"
) {
  clean_plot_memory()

  on.exit({
    clean_plot_memory()
  }, add = TRUE)

  out_dir <- get_plot_dir(
    "patient_examples",
    "forecast_all_items",
    dataset,
    plot_root = plot_root
  )

  plots <- vector("list", length(model_types))
  names(plots) <- model_types

  for (j in seq_along(model_types)) {
    model_type <- model_types[[j]]

    plots[[j]] <- plot_patient_model_all_items(
      dataset = dataset,
      patient_id = patient_id,
      model_type = model_type,
      horizon = horizon,
      result_root = result_root
    )
  }

  combined <- patchwork::wrap_plots(
    plots,
    ncol = length(plots),
    guides = "collect"
  ) +
    patchwork::plot_annotation(
      caption = paste(
        "Blue bands show nested equal-tailed predictive intervals.",
        "The blue line shows the median forecast.",
        "Black points and lines show observed scores."
      )
    ) &
    theme(
      text = element_text(family = "sans"),
      legend.position = "bottom",
      legend.text = element_text(size = 8),
      legend.key.width = unit(1.2, "cm"),
      plot.caption = element_text(size = 8, hjust = 0)
    )

  filename_stem <- paste0(
    "patient_",
    patient_id,
    "_all_items_forecasts_",
    paste(model_types, collapse = "_vs_"),
    "_H",
    horizon
  )

  save_plot(
    combined,
    out_dir = out_dir,
    filename_stem = filename_stem,
    width = max(5.5, 5.0 * length(model_types)),
    height = 12,
    dpi = 320,
    formats = c("png", "pdf", "svg"),
    verbose = TRUE
  )

  rm(plots, combined)
  clean_plot_memory()

  invisible(TRUE)
}


# -------------------------------------------------------------------------
# Save all selected patients
# -------------------------------------------------------------------------

save_dataset_all_items_forecasts <- function(
    dataset,
    patients = NULL,
    horizon = 4,
    model_types = c("EczemaPred", "XemaPred"),
    result_root = "results",
    plot_root = "plots"
) {
  if (is.null(patients)) {
    patients <- get_available_patients(dataset)
  }

  for (patient_id in patients) {
    cat("[INFO] ", dataset, " | patient ", patient_id, "\n", sep = "")

    tryCatch(
      {
        save_patient_all_items_forecast(
          dataset = dataset,
          patient_id = patient_id,
          horizon = horizon,
          model_types = model_types,
          result_root = result_root,
          plot_root = plot_root
        )

        cat("[OK] ", dataset, " | patient ", patient_id, "\n", sep = "")
      },
      error = function(e) {
        message(
          "[WARN] Could not save all-item forecast plot for ",
          dataset,
          " | patient ",
          patient_id,
          ". Full error:"
        )

        print_error_chain(e)
      }
    )

    clean_plot_memory()
  }

  invisible(TRUE)
}


# -------------------------------------------------------------------------
# Command-line interface
# -------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

datasets <- c("Derexyl")
patients <- NULL
horizon <- 4L
model_types <- c("EczemaPred", "XemaPred")
result_root <- "results"
plot_root <- "plots"

i <- 1

while (i <= length(args)) {
  key <- args[[i]]

  if (key %in% c("-h", "--help")) {
    cat(
      "
plot_patient_forecasted_trajectories.R
=======================================

Usage:
  Rscript scripts/03_plotting/patient_trajectories/plot_patient_forecasted_trajectories.R \\
    [--datasets PFDC,Derexyl] \\
    [--dataset PFDC] \\
    [--patients 1,2,3] \\
    [--patient 1] \\
    [--t_horizon 4] \\
    [--models EczemaPred,XemaPred] \\
    [--result_root results] \\
    [--plot_root plots]

Examples:
  Rscript scripts/03_plotting/patient_trajectories/plot_patient_forecasted_trajectories.R \\
    --dataset PFDC \\
    --patient 1 \\
    --models EczemaPred,XemaPred

  Rscript scripts/03_plotting/patient_trajectories/plot_patient_forecasted_trajectories.R \\
    --datasets PFDC,Derexyl \\
    --patients 1,2,3 \\
    --models EczemaPred,XemaPred

  Rscript scripts/03_plotting/patient_trajectories/plot_patient_forecasted_trajectories.R \\
    --dataset PFDC \\
    --patient 1 \\
    --models EczemaPred,XemaPred,PatientSpecificXemaPred
"
    )
    quit(status = 0)

  } else if (key %in% c("--dataset", "--datasets")) {
    datasets <- parse_csv_arg(args[[i + 1]])
    i <- i + 2

  } else if (key == "--patient") {
    patients <- as.integer(args[[i + 1]])
    i <- i + 2

  } else if (key == "--patients") {
    patients <- as.integer(parse_csv_arg(args[[i + 1]]))
    i <- i + 2

  } else if (key == "--t_horizon") {
    horizon <- as.integer(args[[i + 1]])
    i <- i + 2

  } else if (key == "--models") {
    model_types <- parse_csv_arg(args[[i + 1]])
    i <- i + 2

  } else if (key == "--result_root") {
    result_root <- args[[i + 1]]
    i <- i + 2

  } else if (key == "--plot_root") {
    plot_root <- args[[i + 1]]
    i <- i + 2

  } else {
    stop("Unknown argument: ", key)
  }
}


# -------------------------------------------------------------------------
# Validate arguments
# -------------------------------------------------------------------------

if (exists("DATASET_LABELS")) {
  valid_datasets <- names(DATASET_LABELS)
} else {
  valid_datasets <- c("Derexyl", "PFDC", "Fake")
}

bad_datasets <- setdiff(datasets, valid_datasets)

if (length(bad_datasets) > 0) {
  stop(
    "Unknown dataset(s): ",
    paste(bad_datasets, collapse = ", "),
    "\nAllowed datasets are: ",
    paste(valid_datasets, collapse = ", ")
  )
}

model_types <- vapply(
  model_types,
  check_model_type,
  character(1)
)


# -------------------------------------------------------------------------
# Run
# -------------------------------------------------------------------------

cat("============================================================\n")
cat("[INIT] Patient all-item forecast plots\n")
cat("============================================================\n")
cat("[INFO] datasets : ", paste(datasets, collapse = ","), "\n", sep = "")
cat(
  "[INFO] patients : ",
  if (is.null(patients)) "ALL" else paste(patients, collapse = ","),
  "\n",
  sep = ""
)
cat("[INFO] horizon : ", horizon, "\n", sep = "")
cat("[INFO] models : ", paste(model_types, collapse = ","), "\n", sep = "")
cat("[INFO] result_root : ", result_root, "\n", sep = "")
cat("[INFO] plot_root : ", plot_root, "\n", sep = "")

for (dataset in datasets) {
  cat("============================================================\n")
  cat("[DATASET] ", dataset, "\n", sep = "")
  cat("============================================================\n")

  save_dataset_all_items_forecasts(
    dataset = dataset,
    patients = patients,
    horizon = horizon,
    model_types = model_types,
    result_root = result_root,
    plot_root = plot_root
  )

  clean_plot_memory()
}

cat("============================================================\n")
cat("[OK] Done.\n")
cat("============================================================\n")