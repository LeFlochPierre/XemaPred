# -----------------------------------------------------------------------
# Shared plotting and prediction paths
# -----------------------------------------------------------------------

dir_create <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}


# -----------------------------------------------------------------------
# Roots
# -----------------------------------------------------------------------

get_results_root <- function(result_root = "results") {
  here::here(result_root)
}


get_plot_root <- function(plot_root = "plots") {
  out <- here::here(plot_root)
  dir_create(out)
  out
}


get_plot_dir <- function(..., plot_root = "plots") {
  out <- file.path(get_plot_root(plot_root), ...)
  dir_create(out)
  out
}


# -----------------------------------------------------------------------
# Clean model names
# -----------------------------------------------------------------------

canonical_model_type <- function(model_type) {
  key <- tolower(as.character(model_type))
  key <- gsub("[_ -]", "", key)

  dplyr::case_when(
    key == "xemapred" ~ "XemaPred",
    key == "patientspecificxemapred" ~ "PatientSpecificXemaPred",
    key == "eczemapred" ~ "EczemaPred",
    key == "uniform" ~ "uniform",
    key == "historical" ~ "historical",
    key == "mc" ~ "MC",
    TRUE ~ NA_character_
  )
}


check_model_type <- function(model_type) {
  model_type_clean <- canonical_model_type(model_type)

  if (is.na(model_type_clean)) {
    stop(
      "Unknown model_type: ", model_type, "\n",
      "Allowed model types are: ",
      "XemaPred, PatientSpecificXemaPred, EczemaPred, uniform, historical, MC"
    )
  }

  model_type_clean
}


get_model_display_label <- function(model_type) {
  model_type <- check_model_type(model_type)

  dplyr::case_when(
    model_type == "XemaPred" ~ "XemaPred",
    model_type == "PatientSpecificXemaPred" ~ "Patient-specific XemaPred",
    model_type == "EczemaPred" ~ "EczemaPred",
    model_type == "uniform" ~ "Uniform",
    model_type == "historical" ~ "Historical",
    model_type == "MC" ~ "Markov chain",
    TRUE ~ model_type
  )
}


# -----------------------------------------------------------------------
# Item-specific base model
# -----------------------------------------------------------------------

get_item_model_base <- function(item) {
  item <- as.character(item)

  if (item == "extent") {
    return("BinMC")
  }

  if (item %in% c("itching", "sleep")) {
    return("BinRW")
  }

  if (is_intensity_item(item)) {
    return("OrderedRW")
  }

  stop("Unknown item: ", item)
}


# -----------------------------------------------------------------------
# Model directory
# -----------------------------------------------------------------------

get_model_dir <- function(item, model_type) {
  model_type <- check_model_type(model_type)
  base_model <- get_item_model_base(item)

  dplyr::case_when(
    model_type == "XemaPred" ~ paste0(base_model, "-XemaPred"),
    model_type == "PatientSpecificXemaPred" ~ paste0(base_model, "-PatientSpecificXemaPred"),
    model_type == "EczemaPred" ~ base_model,
    model_type == "uniform" ~ "uniform",
    model_type == "historical" ~ "historical",
    model_type == "MC" ~ "MC",
    TRUE ~ NA_character_
  )
}


# -----------------------------------------------------------------------
# Prediction file name
# -----------------------------------------------------------------------

get_prediction_file_name <- function(model_type) {
  model_type <- check_model_type(model_type)

  if (model_type %in% c("XemaPred", "PatientSpecificXemaPred")) {
    return("predictions_ALL_PATIENTS.rds")
  }

  "predictions.rds"
}


# -----------------------------------------------------------------------
# Full prediction path
# -----------------------------------------------------------------------

get_prediction_path <- function(
    dataset,
    item,
    horizon,
    model_type,
    result_root = "results"
) {
  model_type <- check_model_type(model_type)

  path <- file.path(
    get_results_root(result_root),
    dataset,
    item,
    get_model_dir(item, model_type),
    paste0("H", horizon),
    "final",
    get_prediction_file_name(model_type)
  )

  if (!file.exists(path)) {
    stop(
      "Missing prediction file.\n",
      "dataset: ", dataset, "\n",
      "item: ", item, "\n",
      "horizon: ", horizon, "\n",
      "model_type: ", model_type, "\n",
      "expected path:\n  ", path
    )
  }

  path
}


# -----------------------------------------------------------------------
# Debug helper
# -----------------------------------------------------------------------

print_prediction_paths <- function(
    dataset,
    horizon,
    model_types,
    items = ITEM_ORDER_MAIN,
    result_root = "results"
) {
  for (model_type in model_types) {
    model_type <- check_model_type(model_type)

    cat("============================================================\n")
    cat("[MODEL] ", model_type, " -> ", get_model_display_label(model_type), "\n", sep = "")
    cat("============================================================\n")

    for (item in items) {
      path <- file.path(
        get_results_root(result_root),
        dataset,
        item,
        get_model_dir(item, model_type),
        paste0("H", horizon),
        "final",
        get_prediction_file_name(model_type)
      )

      cat("[ITEM] ", item, "\n", sep = "")

      if (file.exists(path)) {
        cat("  FOUND:   ", path, "\n", sep = "")
      } else {
        cat("  MISSING: ", path, "\n", sep = "")
      }
    }
  }

  invisible(TRUE)
}