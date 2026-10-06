#!/usr/bin/env Rscript
# plot_patient_observed_trajectories.R

# -------------------------------------------------------------------------
# Patient-level observed PO-SCORAD trajectories
#
# Purpose:
# - Plot observed trajectories for all PO-SCORAD items and total PO-SCORAD.
# - Save one stacked figure per patient.
#
# Output:
# plots/patient_examples/observed_trajectories/<dataset>/
#   patient_<id>_observed_trajectories.png
#   patient_<id>_observed_trajectories.pdf
#   patient_<id>_observed_trajectories.svg
# -------------------------------------------------------------------------

source(here::here("scripts", "03_plotting", "shared", "00_plot_setup.R"))
source(here::here("scripts", "03_plotting", "shared", "plot_config.R"))
source(here::here("scripts", "03_plotting", "shared", "score_metadata.R"))
source(here::here("scripts", "03_plotting", "shared", "plot_paths.R"))
source(here::here("scripts", "03_plotting", "shared", "plot_scales.R"))
source(here::here("scripts", "03_plotting", "shared", "plot_themes.R"))
source(here::here("scripts", "03_plotting", "shared", "save_plots.R"))
source(here::here("scripts", "03_plotting", "shared", "observed_data.R"))

# -------------------------------------------------------------------------
# Helper: run each patient in a fresh R process
# -------------------------------------------------------------------------

get_current_script_path <- function() {
  args_full <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args_full, value = TRUE)

  if (length(file_arg) != 1L) {
    stop("Could not determine current script path from commandArgs().")
  }

  normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
}


run_patient_workers <- function(
    datasets,
    patients,
    horizon,
    model_types,
    result_root,
    plot_root
) {
  script_path <- get_current_script_path()
  rscript <- file.path(R.home("bin"), "Rscript")

  failures <- tibble::tibble(
    Dataset = character(),
    Patient = integer(),
    Status = integer()
  )

  for (dataset_i in datasets) {
    patient_ids <- patients

    if (is.null(patient_ids)) {
      patient_ids <- get_available_patients(dataset_i)
    }

    cat("============================================================\n")
    cat("[DATASET] ", dataset_i, "\n", sep = "")
    cat("[INFO] Running one patient per R process\n")
    cat("============================================================\n")

    for (patient_i in patient_ids) {
      cat("[RUN] ", dataset_i, " | patient ", patient_i, "\n", sep = "")

      status <- system2(
        command = rscript,
        args = c(
          script_path,
          "--worker",
          "--datasets", dataset_i,
          "--patient", as.character(patient_i),
          "--t_horizon", as.character(horizon),
          "--models", paste(model_types, collapse = ","),
          "--result_root", result_root,
          "--plot_root", plot_root
        )
      )

      if (!identical(status, 0L)) {
        failures <- dplyr::bind_rows(
          failures,
          tibble::tibble(
            Dataset = dataset_i,
            Patient = as.integer(patient_i),
            Status = as.integer(status)
          )
        )

        cat(
          "[WARN] Failed ",
          dataset_i,
          " | patient ",
          patient_i,
          " | status ",
          status,
          "\n",
          sep = ""
        )
      } else {
        cat("[OK] ", dataset_i, " | patient ", patient_i, "\n", sep = "")
      }

      invisible(gc())
    }
  }

  if (nrow(failures) > 0) {
    cat("============================================================\n")
    cat("[WARN] Some patient plots failed:\n")
    print(failures)
    cat("============================================================\n")
  }

  invisible(failures)
}

# -------------------------------------------------------------------------
# Helper: add invisible rows to force fixed y-axis limits per item
# -------------------------------------------------------------------------

make_axis_padding_rows <- function(df_long) {
  x_ref <- min(df_long$Time, na.rm = TRUE)

  purrr::map_dfr(unique(as.character(df_long$Item)), function(item) {
    ylims <- get_item_ylim(item, padded = TRUE)

    tibble::tibble(
      Time = x_ref,
      Item = item,
      Item_label = factor(
        ITEM_LABELS[[item]],
        levels = levels(df_long$Item_label)
      ),
      Score = ylims
    )
  })
}


# -------------------------------------------------------------------------
# Main plotting function
# -------------------------------------------------------------------------

plot_patient_observed_trajectories <- function(df_long, dataset) {
  axis_padding <- make_axis_padding_rows(df_long)

  x_interval <- get_x_break_interval(
    dataset,
    plot_type = "patient_observed"
  )

  x_min <- min(df_long$Time, na.rm = TRUE)
  x_max <- max(df_long$Time, na.rm = TRUE)

  x_breaks <- seq(
    from = ceiling(x_min / x_interval) * x_interval,
    to = floor(x_max / x_interval) * x_interval,
    by = x_interval
  )

  # Always include the first observed day.
  x_breaks <- sort(unique(c(x_min, x_breaks)))

  ggplot(df_long, aes(x = Time, y = Score)) +
    geom_blank(
      data = axis_padding,
      aes(x = Time, y = Score)
    ) +
    geom_line(
      colour = "grey15",
      linewidth = 0.48,
      lineend = "round"
    ) +
    geom_point(
      colour = "grey10",
      fill = "white",
      shape = 21,
      stroke = 0.35,
      size = 1.65
    ) +
    facet_grid(
      Item_label ~ .,
      scales = "free_y",
      switch = "y"
    ) +
    scale_x_continuous(
      breaks = x_breaks,
      limits = c(x_min, x_max),
      expand = expansion(mult = c(0.006, 0.006))
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0)),
      labels = function(x) ifelse(x < 0, "", x)
    ) +
    labs(
      x = "Day",
      y = NULL
    ) +
    theme_patient_facets(base_size = 11) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      axis.title.x = element_text(
        size = 10,
        colour = "grey10",
        margin = margin(t = 8)
      ),
      axis.text.x = element_text(size = 8, colour = "grey25"),
      axis.text.y = element_text(size = 7, colour = "grey30"),
      axis.ticks = element_line(colour = "grey45", linewidth = 0.25),
      axis.ticks.length = unit(2.2, "pt"),
      panel.spacing.y = unit(1.15, "lines"),
      plot.margin = margin(t = 12, r = 14, b = 10, l = 18)
    )
}


# -------------------------------------------------------------------------
# Save one patient
# -------------------------------------------------------------------------

save_patient_observed_trajectories <- function(
  dataset,
  patient_id,
  plot_root = "plots"
) {
  out_dir <- get_plot_dir(
    "patient_examples",
    "observed_trajectories",
    dataset,
    plot_root = plot_root
  )

  df_long <- load_patient_observed_long(
    dataset = dataset,
    patient_id = patient_id,
    items = ITEM_ORDER_WITH_POSCORAD
  )

  p <- plot_patient_observed_trajectories(
    df_long = df_long,
    dataset = dataset
  )

  filename_stem <- paste0("patient_", patient_id, "_observed_trajectories")

  save_plot(
    p,
    out_dir = out_dir,
    filename_stem = filename_stem,
    width = 7.6,
    height = 15.2,
    formats = c("png", "pdf", "svg")
  )

  invisible(p)
}


# -------------------------------------------------------------------------
# Save all selected patients from one dataset
# -------------------------------------------------------------------------

save_dataset_observed_trajectories <- function(
  dataset,
  patients = NULL,
  plot_root = "plots"
) {
  if (is.null(patients)) {
    patients <- get_available_patients(dataset)
  }

  purrr::walk(patients, function(patient_id) {
    cat("------------------------------------------------------------\n")
    cat("[INFO] Dataset=", dataset, " | Patient=", patient_id, "\n", sep = "")
    cat("------------------------------------------------------------\n")

    tryCatch(
      save_patient_observed_trajectories(
        dataset = dataset,
        patient_id = patient_id,
        plot_root = plot_root
      ),
      error = function(e) {
        message(
          "[WARN] Could not save observed trajectories for ",
          dataset,
          " | patient ",
          patient_id,
          " : ",
          e$message
        )
      }
    )
  })
}


# -------------------------------------------------------------------------
# Command-line interface
# -------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

datasets <- DATASETS_DEFAULT
patients <- NULL
plot_root <- "plots"

i <- 1

while (i <= length(args)) {
  key <- args[[i]]

  if (key %in% c("-h", "--help")) {
    cat(
      "
01_plot_patient_observed_trajectories.R
========================================

Usage:
  Rscript scripts/03_plotting/patient_examples/01_plot_patient_observed_trajectories.R \\
    [--datasets PFDC,Derexyl] \\
    [--patients 1,2,3] \\
    [--plot_root plots]

Examples:
  Rscript scripts/03_plotting/patient_examples/01_plot_patient_observed_trajectories.R \\
    --datasets PFDC \\
    --patients 1

  Rscript scripts/03_plotting/patient_examples/01_plot_patient_observed_trajectories.R \\
    --datasets Derexyl \\
    --patients 1,2,3

  Rscript scripts/03_plotting/patient_examples/01_plot_patient_observed_trajectories.R \\
    --datasets PFDC,Derexyl
"
    )
    quit(status = 0)

  } else if (key == "--datasets") {
    datasets <- parse_csv_arg(args[[i + 1]])
    i <- i + 2

  } else if (key == "--patients") {
    patients <- as.integer(parse_csv_arg(args[[i + 1]]))
    i <- i + 2

  } else if (key == "--plot_root") {
    plot_root <- args[[i + 1]]
    i <- i + 2

  } else {
    stop("Unknown argument: ", key)
  }
}


cat("============================================================\n")
cat("[INIT] Observed patient trajectory plots\n")
cat("============================================================\n")
cat("[INFO] datasets : ", paste(datasets, collapse = ","), "\n", sep = "")
cat(
  "[INFO] patients : ",
  if (is.null(patients)) "ALL" else paste(patients, collapse = ","),
  "\n",
  sep = ""
)
cat("[INFO] plot_root : ", plot_root, "\n", sep = "")

for (dataset in datasets) {
  save_dataset_observed_trajectories(
    dataset = dataset,
    patients = patients,
    plot_root = plot_root
  )
}

cat("============================================================\n")
cat("[OK] Done.\n")
cat("============================================================\n")