# -----------------------------------------------------------------------
# Shared learning-curve data helpers
# -----------------------------------------------------------------------

LEARNING_CURVE_METRICS <- c(
  "lpd",
  "crps",
  "accuracy",
  "accuracy_median",
  "accuracy_map",
  "accuracy_prob"
)

LEARNING_CURVE_METRIC_PATTERN <- paste(
  c(
    "lpd", "LPD",
    "crps", "CRPS",
    "accuracy", "Accuracy",
    "accuracy_median", "Accuracy_median",
    "accuracy_map", "Accuracy_map",
    "accuracy_prob", "Accuracy_prob"
  ),
  collapse = "|"
)

REFERENCE_MODEL_DEFAULTS <- c("RW", "MC", "historical", "uniform")

REFERENCE_MODEL_PATTERNS <- c(
  "^RW$",
  "^MC$",
  "^historical$",
  "^uniform$",
  "^AR1$",
  "^MixedAR1$",
  "^Smoothing$"
)

lc_sanitize_token <- function(x) {
  x <- gsub("[/\\\\]", "_", x)
  x <- gsub("[()\\[\\]]", "", x)
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)

  if (!nzchar(x)) "experiment" else x
}

normalize_learning_curve_metric <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x[x == ""] <- NA_character_
  x
}

learning_curve_metric_label <- function(metric) {
  metric <- normalize_learning_curve_metric(metric)

  dplyr::case_when(
    metric == "lpd"             ~ "LPD",
    metric == "crps"            ~ "CRPS",
    metric == "accuracy"        ~ "Point-forecast accuracy",
    metric == "accuracy_median" ~ "Median point-forecast accuracy",
    metric == "accuracy_map"    ~ "MAP point-forecast accuracy",
    metric == "accuracy_prob"   ~ "Probabilistic accuracy",
    TRUE                        ~ metric
  )
}

learning_curve_model_label <- function(model_type) {
  model_type <- as.character(model_type)

  dplyr::case_when(
    model_type == "XemaPred" ~ "XemaPred",

    # Canonical patient-specific model:
    # PFDC    = CrossCohort
    # Derexyl = PopPrior
    model_type == "PatientSpecificXemaPredCanonical" ~
      "Patient-specific XemaPred",

    # Original / incorrect non-power-prior implementation
    model_type == "PatientSpecificXemaPred" ~
      "Patient-specific XemaPred (non-power prior)",

    model_type == "PatientSpecificXemaPredPopPrior" ~
      "Patient-specific XemaPred (PopPrior)",

    model_type == "PatientSpecificXemaPredCrossCohort" ~
      "Patient-specific XemaPred (cross-cohort)",

    model_type == "EczemaPred" ~ "EczemaPred",
    model_type == "RW" ~ "Random walk",
    model_type == "MC" ~ "Markov chain",
    model_type == "historical" ~ "Historical",
    model_type == "uniform" ~ "Uniform",
    model_type == "AR1" ~ "AR(1)",
    model_type == "MixedAR1" ~ "Mixed AR(1)",
    model_type == "Smoothing" ~ "Smoothing",
    TRUE ~ model_type
  )
}

is_reference_model <- function(model_type) {
  model_type <- as.character(model_type)

  purrr::map_lgl(model_type, function(x) {
    any(grepl(paste(REFERENCE_MODEL_PATTERNS, collapse = "|"), x))
  })
}

detect_reference_like_models <- function(present_models) {
  present_models <- as.character(present_models)

  unique(unlist(lapply(REFERENCE_MODEL_PATTERNS, function(pattern) {
    grep(pattern, present_models, value = TRUE)
  })))
}

load_rdata_first_object <- function(path, prefer = NULL) {
  if (!file.exists(path)) {
    stop("Missing file: ", path)
  }

  env <- new.env(parent = emptyenv())
  load(path, envir = env)

  if (!is.null(prefer) && exists(prefer, envir = env, inherits = FALSE)) {
    return(env[[prefer]])
  }

  objects <- ls(envir = env)

  if (!length(objects)) {
    stop("No objects found in RData file: ", path)
  }

  preferred_objects <- c(
    "perf_all",
    "perf_all_wref",
    "perf",
    "perf_items",
    "summary_lpd",
    "perf_model"
  )

  for (object_name in preferred_objects) {
    if (exists(object_name, envir = env, inherits = FALSE)) {
      return(env[[object_name]])
    }
  }

  env[[objects[[1]]]]
}

extract_model_from_perf_filename <- function(path, dataset = NULL) {
  filename <- basename(path)
  filename <- sub("^perf_", "", filename)
  filename <- sub("\\.RData$", "", filename, ignore.case = TRUE)

  if (!is.null(dataset)) {
    filename <- sub(
      paste0("_", dataset, "_(", LEARNING_CURVE_METRIC_PATTERN, ")$"),
      "",
      filename
    )
  } else {
    filename <- sub(
      paste0("_(Derexyl|PFDC)_(", LEARNING_CURVE_METRIC_PATTERN, ")$"),
      "",
      filename
    )
  }

  filename
}

infer_metric_from_perf_path <- function(path) {
  path_norm <- gsub("\\\\", "/", path)

  metric_from_folder <- stringr::str_match(
    path_norm,
    "/(?:ITEMS|POSCORAD)/([^/]+)/"
  )[, 2]

  metric_from_folder <- normalize_learning_curve_metric(metric_from_folder)

  if (!is.na(metric_from_folder) && metric_from_folder %in% LEARNING_CURVE_METRICS) {
    return(metric_from_folder)
  }

  metric_from_file <- stringr::str_match(
    basename(path),
    paste0("_(", LEARNING_CURVE_METRIC_PATTERN, ")\\.RData$")
  )[, 2]

  normalize_learning_curve_metric(metric_from_file)
}

normalize_learning_curve_columns <- function(df, item_scope = c("items", "poscorad")) {
  item_scope <- match.arg(item_scope)

  if (!is.data.frame(df)) {
    stop("Learning-curve object is not a data frame.")
  }

  if ("Metric" %in% names(df) && !"metric" %in% names(df)) {
    df <- dplyr::rename(df, metric = Metric)
  }

  if (!"metric" %in% names(df)) {
    df$metric <- NA_character_
  }

  if (!"Variable" %in% names(df)) {
    df$Variable <- "Fit"
  }

  if (!"Model_type" %in% names(df)) {
    if ("Model" %in% names(df)) {
      df$Model_type <- as.character(df$Model)
    } else {
      df$Model_type <- NA_character_
    }
  }

  if (!"Dataset" %in% names(df)) {
    df$Dataset <- NA_character_
  }

  if (!"Item" %in% names(df)) {
    df$Item <- if (item_scope == "poscorad") "SCORAD" else NA_character_
  }

  if (item_scope == "poscorad") {
    df$Item <- "SCORAD"
  }

  df %>%
    mutate(
      Dataset = as.character(.data$Dataset),
      Model_type = as.character(.data$Model_type),
      Item = if (item_scope == "poscorad") {
        "SCORAD"
      } else {
        tolower(as.character(.data$Item))
      },
      Variable = tolower(trimws(as.character(.data$Variable))),
      Horizon = suppressWarnings(as.integer(.data$Horizon)),
      N = suppressWarnings(as.integer(.data$N)),
      Mean = suppressWarnings(as.numeric(.data$Mean)),
      SE = suppressWarnings(as.numeric(.data$SE)),
      metric = normalize_learning_curve_metric(.data$metric)
    ) %>%
    filter(
      !is.na(.data$Dataset),
      !is.na(.data$Model_type),
      !is.na(.data$Item),
      !is.na(.data$Horizon),
      !is.na(.data$N),
      !is.na(.data$Mean),
      !is.na(.data$SE)
    )
}

assert_no_duplicate_learning_curve_rows <- function(df, item_scope = c("items", "poscorad")) {
  item_scope <- match.arg(item_scope)

  key_cols <- c("Dataset", "metric", "Model_type", "Item", "Horizon", "N", "Variable")
  missing_cols <- setdiff(key_cols, names(df))

  if (length(missing_cols)) {
    stop(
      "Learning-curve table is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  duplicates <- df %>%
    count(across(all_of(key_cols)), name = "n") %>%
    filter(.data$n > 1)

  if (nrow(duplicates)) {
    print(duplicates, n = 200)
    stop(
      "Duplicate learning-curve rows detected for the same ",
      "Dataset/metric/Model/Item/Horizon/N/Variable."
    )
  }

  invisible(TRUE)
}

read_learning_curve_folder <- function(folder, dataset, item_scope = c("items", "poscorad")) {
  item_scope <- match.arg(item_scope)

  files <- list.files(
    folder,
    pattern = paste0(
      "^perf_.+_",
      dataset,
      "_(",
      LEARNING_CURVE_METRIC_PATTERN,
      ")\\.RData$"
    ),
    full.names = TRUE,
    recursive = TRUE
  )

  if (!length(files)) {
    return(NULL)
  }

  purrr::map_dfr(files, function(path) {
    model_type <- extract_model_from_perf_filename(path, dataset = dataset)
    metric <- infer_metric_from_perf_path(path)

    df <- load_rdata_first_object(path)

    if (!is.data.frame(df)) {
      warning("Skipping non-data-frame object in: ", path)
      return(NULL)
    }

    df$Dataset <- dataset
    df$Model_type <- model_type

    df <- normalize_learning_curve_columns(df, item_scope = item_scope)

    if (!"metric" %in% names(df) || all(is.na(df$metric))) {
      df$metric <- metric
    } else {
      df$metric[is.na(df$metric)] <- metric
    }

    df
  })
}

read_learning_curve_perf <- function(
    datasets = DATASETS_DEFAULT,
    item_scope = c("items", "poscorad"),
    result_root = "results"
) {
  item_scope <- match.arg(item_scope)

  folder_name <- dplyr::case_when(
    item_scope == "items"    ~ "ITEMS",
    item_scope == "poscorad" ~ "POSCORAD",
    TRUE                     ~ NA_character_
  )

  out <- purrr::map_dfr(datasets, function(dataset) {
    folder <- file.path(
      get_results_root(result_root),
      dataset,
      "ALL_MODELS",
      folder_name
    )

    if (!dir.exists(folder)) {
      warning("Skipping missing learning-curve folder: ", folder)
      return(NULL)
    }

    read_learning_curve_folder(
      folder = folder,
      dataset = dataset,
      item_scope = item_scope
    )
  })

  if (is.null(out) || !nrow(out)) {
    stop(
      "No learning-curve performance files found for scope '",
      item_scope,
      "' under results/<dataset>/ALL_MODELS/",
      folder_name,
      "."
    )
  }

  out <- normalize_learning_curve_columns(out, item_scope = item_scope)
  assert_no_duplicate_learning_curve_rows(out, item_scope = item_scope)

  out
}

resolve_learning_curve_models <- function(
    perf,
    requested_models,
    add_ref_models = FALSE,
    ref_models = REFERENCE_MODEL_DEFAULTS,
    add_all_ref_like_models = FALSE
) {
  present_models <- sort(unique(as.character(perf$Model_type)))

  models_full <- unique(as.character(requested_models))

  if (isTRUE(add_ref_models)) {
    models_full <- unique(c(models_full, ref_models))
  }

  if (isTRUE(add_all_ref_like_models)) {
    models_full <- unique(c(models_full, detect_reference_like_models(present_models)))
  }

  missing_models <- setdiff(models_full, present_models)
  model_order_raw <- models_full[models_full %in% present_models]

  if (!length(model_order_raw)) {
    stop("None of the requested models were found in the learning-curve table.")
  }

  requested_present <- requested_models[requested_models %in% present_models]
  has_references <- any(is_reference_model(model_order_raw))

  list(
    present_models = present_models,
    missing_models = missing_models,
    model_order_raw = model_order_raw,
    requested_present = requested_present,
    has_references = has_references
  )
}

make_learning_curve_experiment_name <- function(main_models, include_references = FALSE) {
  if (!length(main_models)) {
    main_models <- "models"
  }

  out <- lc_sanitize_token(
    paste(learning_curve_model_label(main_models), collapse = "_vs_")
  )

  if (isTRUE(include_references)) {
    out <- paste0(out, "_with_references")
  }

  out
}