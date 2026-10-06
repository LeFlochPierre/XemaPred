# ======================================================================
# Shared helpers for EczemaPred / reference-model result processing
#
# Source this file after:
#   source(here::here("scripts", "00_Setup", "00_init.R"))
# and after loading dplyr/tidyr/purrr/glue/tibble/scoringRules.
# ======================================================================

# ----------------------------------------------------------------------
# Logging
# ----------------------------------------------------------------------

log_header <- function(title) {
  cat("\n", strrep("=", 72), "\n", sep = "")
  cat(title, "\n", sep = "")
  cat(strrep("=", 72), "\n", sep = "")
}

log_section <- function(title) {
  cat("\n", strrep("-", 72), "\n", sep = "")
  cat(title, "\n", sep = "")
  cat(strrep("-", 72), "\n", sep = "")
}

log_info <- function(...) {
  cat("[INFO] ", ..., "\n", sep = "")
}

log_warn <- function(...) {
  cat("[WARN] ", ..., "\n", sep = "")
}

log_save <- function(path) {
  cat("[SAVE] ", path, "\n", sep = "")
}

# ----------------------------------------------------------------------
# Small utilities
# ----------------------------------------------------------------------

dir_create <- function(p) {
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
}

sanitize_token <- function(x) {
  x <- gsub("[/\\\\]", "_", x)
  x <- gsub("[()\\[\\]]", "", x)
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)

  if (!nzchar(x)) "MODEL" else x
}

parse_csv_arg <- function(x) {
  x <- gsub("\\s+", "", x)
  if (!nzchar(x)) character(0) else strsplit(x, ",")[[1]]
}

bool_arg <- function(x) {
  as.integer(x) == 1L
}

get_max_train_day <- function(ds, dict_datasets) {
  dict_datasets %>%
    filter(.data$Dataset == ds) %>%
    pull(.data$Max_train_day) %>%
    .[[1]]
}

# ----------------------------------------------------------------------
# Output paths
# ----------------------------------------------------------------------

processed_root <- function(ds) {
  here::here("results", ds, "ALL_MODELS")
}

items_metric_dir <- function(ds, metric) {
  file.path(processed_root(ds), "ITEMS", metric)
}

poscorad_metric_dir <- function(ds, metric) {
  file.path(processed_root(ds), "POSCORAD", metric)
}

comp_time_dir <- function(ds) {
  file.path(processed_root(ds), "comp_time")
}

# ----------------------------------------------------------------------
# Model mappings
# ----------------------------------------------------------------------

is_intensity_item <- function(item) {
  item %in% detail_POSCORAD("Intensity signs")$Name
}

eczemapred_model_for_item <- function(item) {
  if (item == "extent") return("BinMC")
  if (item %in% c("itching", "sleep")) return("BinRW")
  if (is_intensity_item(item)) return("OrderedRW")

  stop("Unknown EczemaPred item: ", item)
}

eczemapred_items <- function() {
  c(
    "extent",
    detail_POSCORAD("Intensity signs")$Name,
    "itching",
    "sleep"
  )
}

reference_models_for_item <- function(item) {
  if (item == "SCORAD") {
    return(c("uniform", "historical", "RW", "AR1", "MixedAR1", "Smoothing"))
  }

  if (item == "extent" || item %in% c("itching", "sleep")) {
    return(c("uniform", "historical", "RW"))
  }

  if (is_intensity_item(item)) {
    return(c("uniform", "historical", "MC"))
  }

  stop("Unknown reference-model item: ", item)
}

# ----------------------------------------------------------------------
# Raw prediction / diagnostics paths
# ----------------------------------------------------------------------

compare_old_new_eczemapred_join <- function(ds, perf_horizon = 4, dict_datasets) {
  tmp <- build_poscorad_testing_frame(ds, perf_horizon, dict_datasets)

  pred_sc_old <- tmp$pred
  pred_sc_new <- tmp$pred

  for (severity_item in eczemapred_items()) {
    f <- eczemapred_prediction_file(ds, severity_item, perf_horizon)

    res_old <- readRDS(f) %>%
      dplyr::select(Patient, Time, Iteration, Samples) %>%
      dplyr::rename(!!paste0(severity_item, "_pred") := Samples)

    res_new <- readRDS(f) %>%
      dplyr::filter(.data$Horizon <= perf_horizon) %>%
      dplyr::select(Patient, Time, Horizon, Iteration, Samples) %>%
      dplyr::rename(!!paste0(severity_item, "_pred") := Samples)

    pred_sc_old <- dplyr::left_join(
      pred_sc_old,
      res_old,
      by = c("Patient", "Time", "Iteration")
    )

    pred_sc_new <- dplyr::left_join(
      pred_sc_new,
      res_new,
      by = c("Patient", "Time", "Horizon", "Iteration")
    )
  }

  pred_sc_old <- pred_sc_old %>% tidyr::drop_na()
  pred_sc_new <- pred_sc_new %>% tidyr::drop_na()

  res_old <- build_scorad_samples(pred_sc_old) %>%
    ensure_metric_column(
      item = "SCORAD",
      metric = "Accuracy_prob",
      allow_no_samples = FALSE
    )

  res_new <- build_scorad_samples(pred_sc_new) %>%
    ensure_metric_column(
      item = "SCORAD",
      metric = "Accuracy_prob",
      allow_no_samples = FALSE
    )

  fc_it_sc <- detail_fc_training(
    tmp$POSCORAD %>% dplyr::rename(Time = Day),
    perf_horizon
  )

  perf_old <- estimate_performance(
    "Accuracy_prob",
    res_old,
    fc_it_sc,
    adjust_horizon = TRUE
  ) %>%
    dplyr::filter(.data$Variable == "Fit") %>%
    dplyr::mutate(join_style = "old_no_Horizon")

  perf_new <- estimate_performance(
    "Accuracy_prob",
    res_new,
    fc_it_sc,
    adjust_horizon = TRUE
  ) %>%
    dplyr::filter(.data$Variable == "Fit") %>%
    dplyr::mutate(join_style = "new_with_Horizon")

  dplyr::bind_rows(perf_old, perf_new) %>%
    dplyr::mutate(
      Dataset = ds,
      Accuracy_percent = 100 * .data$Mean,
      SE_percent = 100 * .data$SE
    ) %>%
    dplyr::select(
      Dataset,
      join_style,
      N,
      LastTime,
      Accuracy_percent,
      SE_percent,
      Mean,
      SE
    ) %>%
    dplyr::arrange(.data$LastTime, .data$join_style)
}

eczemapred_prediction_file <- function(ds, item, horizon) {
  model <- eczemapred_model_for_item(item)

  file.path(
    "results",
    ds,
    item,
    model,
    paste0("H", horizon),
    "final",
    "predictions.rds"
  )
}

eczemapred_diag_file <- function(ds, item, horizon) {
  model <- eczemapred_model_for_item(item)

  file.path(
    "results",
    ds,
    item,
    model,
    paste0("H", horizon),
    "final",
    "diagnostics.rds"
  )
}

reference_prediction_file <- function(ds, item, model, horizon) {
  file.path(
    "results",
    ds,
    item,
    model,
    paste0("H", horizon),
    "final",
    "predictions.rds"
  )
}

# ----------------------------------------------------------------------
# Expected file grids
# ----------------------------------------------------------------------

build_eczemapred_item_grid <- function(ds, horizon) {
  tibble(
    Dataset = ds,
    Item = eczemapred_items()
  ) %>%
    mutate(
      Model = vapply(.data$Item, eczemapred_model_for_item, character(1)),
      Model_type = "EczemaPred",
      Component = case_when(
        .data$Item == "extent" ~ "Extent",
        .data$Item %in% c("itching", "sleep") ~ "Subjective symptoms",
        TRUE ~ "Intensity signs"
      ),
      Inference = "HMC",
      Label = "Main model",
      File = vapply(
        .data$Item,
        function(item) eczemapred_prediction_file(ds, item, horizon),
        character(1)
      )
    )
}

build_eczemapred_diag_grid <- function(ds, horizon) {
  build_eczemapred_item_grid(ds, horizon) %>%
    mutate(
      File = vapply(
        .data$Item,
        function(item) eczemapred_diag_file(ds, item, horizon),
        character(1)
      )
    )
}

reference_item_grid <- function(ds, horizon) {
  dplyr::bind_rows(
    tibble(
      Dataset = ds,
      Item = "extent",
      Model = reference_models_for_item("extent"),
      Component = "Extent",
      Inference = "Reference",
      Label = "Reference model"
    ),
    tidyr::expand_grid(
      Dataset = ds,
      Item = c("itching", "sleep"),
      Model = reference_models_for_item("itching"),
      Component = "Subjective symptoms",
      Inference = "Reference",
      Label = "Reference model"
    ),
    tidyr::expand_grid(
      Dataset = ds,
      Item = detail_POSCORAD("Intensity signs")$Name,
      Model = reference_models_for_item("dryness"),
      Component = "Intensity signs",
      Inference = "Reference",
      Label = "Reference model"
    )
  ) %>%
    mutate(
      Model_type = as.character(.data$Model),
      File = pmap_chr(
        list(.data$Dataset, .data$Item, .data$Model),
        function(Dataset, Item, Model) {
          reference_prediction_file(
            ds = Dataset,
            item = Item,
            model = Model,
            horizon = horizon
          )
        }
      )
    )
}

reference_poscorad_grid <- function(ds, horizon) {
  tibble(
    Dataset = ds,
    Item = "SCORAD",
    Model = reference_models_for_item("SCORAD"),
    Component = "SCORAD",
    Inference = "Reference",
    Label = "Reference model",
    Model_type = as.character(reference_models_for_item("SCORAD"))
  ) %>%
    mutate(
      File = pmap_chr(
        list(.data$Dataset, .data$Item, .data$Model),
        function(Dataset, Item, Model) {
          reference_prediction_file(
            ds = Dataset,
            item = Item,
            model = Model,
            horizon = horizon
          )
        }
      )
    )
}

# ----------------------------------------------------------------------
# File preflight
# ----------------------------------------------------------------------

preflight_files <- function(grids, strict = TRUE) {
  if (!length(grids)) {
    log_warn("No file grids selected for preflight.")
    return(invisible(TRUE))
  }

  grid <- bind_rows(grids) %>%
    mutate(exists = file.exists(.data$File))

  log_section("Input file availability")

  summary <- grid %>%
    count(Dataset, Section, exists, name = "n") %>%
    arrange(Dataset, Section, desc(exists))

  print(summary, n = Inf)

  missing <- grid %>%
    filter(!.data$exists) %>%
    select(Section, Dataset, Item, Model, File)

  if (nrow(missing) > 0) {
    log_warn(nrow(missing), " required input file(s) missing.")
    print(missing, n = Inf)

    if (isTRUE(strict)) {
      stop("Stopping because --strict_files is enabled and some input files are missing.")
    }
  } else {
    log_info("All selected input files are available.")
  }

  invisible(nrow(missing) == 0)
}

# ----------------------------------------------------------------------
# Accuracy / point prediction utilities
# ----------------------------------------------------------------------

acc_thr_item <- function(item) {
  if (item == "SCORAD") return(5)
  if (item == "extent") return(5)
  if (item %in% c("itching", "sleep")) return(1)

  # Intensity signs are ordinal grades; exact match is encoded using
  # a tiny threshold to avoid edge-case equality issues.
  1e-12
}

sample_matrix <- function(samples) {
  do.call(rbind, samples)
}

median_point_accuracy <- function(y, pred, ct) {
  med <- vapply(pred, function(x) stats::median(x, na.rm = TRUE), numeric(1))
  as.numeric(abs(y - med) <= ct)
}

median_abs_error <- function(y, pred) {
  med <- vapply(pred, function(x) stats::median(x, na.rm = TRUE), numeric(1))
  abs(y - med)
}

safe_mode <- function(x, digits_fallback = 2, unique_frac_cutoff = 0.9) {
  x <- suppressWarnings(as.numeric(unlist(x)))
  x <- x[is.finite(x)]

  if (!length(x)) return(NA_real_)

  unique_frac <- length(unique(x)) / length(x)

  if (unique_frac < unique_frac_cutoff) {
    ux <- sort(unique(x))
    return(ux[which.max(tabulate(match(x, ux)))])
  }

  xr <- round(x, digits = digits_fallback)
  ux <- sort(unique(xr))
  ux[which.max(tabulate(match(xr, ux)))]
}

map_point_accuracy <- function(y, pred, ct) {
  mp <- vapply(pred, safe_mode, numeric(1))
  as.numeric(abs(y - mp) <= ct)
}

map_abs_error <- function(y, pred) {
  mp <- vapply(pred, safe_mode, numeric(1))
  abs(y - mp)
}

prob_within_tol <- function(y, pred, ct) {
  vapply(seq_along(pred), function(i) {
    mean(abs(pred[[i]] - y[i]) <= ct, na.rm = TRUE)
  }, numeric(1))
}

clamp_scorad_lpd <- function(res) {
  if (!"lpd" %in% names(res)) return(res)

  max_score <- 103
  lpd_lb <- -log(max_score * 100)
  lpd_ub <- -log(0.1)

  res %>%
    mutate(
      lpd = replace(.data$lpd, .data$lpd < lpd_lb, lpd_lb),
      lpd = replace(.data$lpd, .data$lpd > lpd_ub, lpd_ub)
    )
}

ensure_metric_column <- function(res, item, metric, allow_no_samples = FALSE) {
  stopifnot(is.data.frame(res))
  stopifnot(all(c("Iteration", "Horizon") %in% names(res)))

  if (metric %in% names(res)) {
    if (item == "SCORAD") return(clamp_scorad_lpd(res))
    return(res)
  }

  has_samples <- all(c("Score", "Samples") %in% names(res))

  if (!has_samples) {
    if (isTRUE(allow_no_samples)) return(res)

    stop(
      "Missing Score/Samples needed to compute metric = ", metric,
      " for item = ", item,
      ". Columns are: ", paste(names(res), collapse = ", ")
    )
  }

  if (metric == "lpd") {
    res <- res %>%
      mutate(lpd = -scoringRules::logs_sample(.data$Score, sample_matrix(.data$Samples)))

  } else if (metric == "CRPS") {
    res <- res %>%
      mutate(CRPS = scoringRules::crps_sample(.data$Score, sample_matrix(.data$Samples)))

  } else if (metric == "Accuracy_median") {
    thr <- acc_thr_item(item)
    res <- res %>%
      mutate(
        Accuracy_median = median_point_accuracy(.data$Score, .data$Samples, thr),
        MAE_median = median_abs_error(.data$Score, .data$Samples)
      )

  } else if (metric == "Accuracy_map") {
    thr <- acc_thr_item(item)
    res <- res %>%
      mutate(
        Accuracy_map = map_point_accuracy(.data$Score, .data$Samples, thr),
        MAE_map = map_abs_error(.data$Score, .data$Samples)
      )

  } else if (metric == "Accuracy_prob") {
    thr <- acc_thr_item(item)
    res <- res %>%
      mutate(Accuracy_prob = prob_within_tol(.data$Score, .data$Samples, thr))

  } else if (metric == "MAE_median") {
    res <- res %>%
      mutate(MAE_median = median_abs_error(.data$Score, .data$Samples))

  } else if (metric == "MAE_map") {
    res <- res %>%
      mutate(MAE_map = map_abs_error(.data$Score, .data$Samples))

  } else {
    stop("Unknown metric: ", metric)
  }

  if (item == "SCORAD") clamp_scorad_lpd(res) else res
}

# ----------------------------------------------------------------------
# Generic performance computation
# ----------------------------------------------------------------------

compute_perf_from_grid <- function(
    grid,
    ds,
    metric,
    perf_horizon,
    dict_datasets,
    allow_no_samples,
    filter_rw_iteration_zero,
    adjust_horizon_fn
) {
  if (!nrow(grid)) return(tibble())

  df_ds <- load_dataset(ds)
  fc_it <- detail_fc_training(df_ds %>% rename(Time = Day), perf_horizon)
  max_day <- get_max_train_day(ds, dict_datasets)

  perf_list <- vector("list", nrow(grid))
  skipped <- tibble(Item = character(), Model = character(), Reason = character())

  for (i in seq_len(nrow(grid))) {
    item <- grid$Item[i]
    model <- grid$Model[i]
    f <- grid$File[i]

    if (!file.exists(f)) {
      stop("Missing input file after preflight: ", f)
    }

    res <- readRDS(f) %>%
      filter(
        .data$Horizon <= perf_horizon,
        (.data$Time - .data$Horizon) <= max_day
      )

    if (isTRUE(filter_rw_iteration_zero) && model == "RW") {
      res <- res %>% filter(.data$Iteration > 0)
    }

    res <- ensure_metric_column(
      res = res,
      item = item,
      metric = metric,
      allow_no_samples = allow_no_samples
    )

    if (!(metric %in% names(res))) {
      skipped <- bind_rows(
        skipped,
        tibble(
          Item = item,
          Model = model,
          Reason = "Metric column and Samples column are unavailable"
        )
      )
      next
    }

    perf_list[[i]] <- estimate_performance(
      metric,
      res,
      fc_it,
      adjust_horizon = adjust_horizon_fn(model)
    ) %>%
      bind_cols(grid[i, ])
  }

  if (nrow(skipped)) {
    expected_skipped <- skipped %>%
      filter(
        .data$Item != "SCORAD",
        .data$Model %in% c("uniform", "historical"),
        metric %in% c(
          "CRPS",
          "Accuracy_median",
          "Accuracy_map",
          "Accuracy_prob",
          "MAE_median",
          "MAE_map"
        )
      )

    unexpected_skipped <- skipped %>%
      anti_join(
        expected_skipped,
        by = c("Item", "Model", "Reason")
      )

    if (nrow(expected_skipped)) {
      log_info(
        nrow(expected_skipped),
        " expected item-level uniform/historical skip(s) for ",
        metric,
        ": we do not save Samples for these baselines."
      )
    }

    if (nrow(unexpected_skipped)) {
      log_warn(
        nrow(unexpected_skipped),
        " unexpected model(s) skipped for metric ",
        metric,
        " because required sample-level predictions were unavailable."
      )
      print(unexpected_skipped, n = Inf)
    }
  }

  perf_raw <- bind_rows(perf_list)

  if (!nrow(perf_raw) || !"Variable" %in% names(perf_raw)) {
    return(tibble())
  }

  perf_raw %>%
    filter(.data$Variable == "Fit") %>%
    mutate(
      Model_type = as.character(.data$Model_type),
      Metric = metric
    )
}

# ----------------------------------------------------------------------
# Save helpers
# ----------------------------------------------------------------------

write_split_perf_by_model <- function(out_root, perf_df, ds, metric, overwrite = TRUE) {
  stopifnot(is.data.frame(perf_df))

  if (!nrow(perf_df)) {
    log_warn("No performance rows to save in ", out_root)
    return(invisible(FALSE))
  }

  stopifnot("Model_type" %in% names(perf_df))

  dir_create(out_root)

  models <- sort(unique(as.character(perf_df$Model_type)))

  for (m in models) {
    df_m <- perf_df %>% filter(.data$Model_type == m)
    if (!nrow(df_m)) next

    out_path <- file.path(
      out_root,
      paste0("perf_", sanitize_token(m), "_", ds, "_", metric, ".RData")
    )

    if (!file.exists(out_path) || isTRUE(overwrite)) {
      perf_model <- df_m
      save(perf_model, file = out_path)
      log_save(out_path)
    } else {
      log_info("Existing output kept: ", out_path)
    }
  }

  invisible(TRUE)
}

# ----------------------------------------------------------------------
# PO-SCORAD reconstruction from item-level EczemaPred predictions
# ----------------------------------------------------------------------

build_poscorad_testing_frame <- function(ds, perf_horizon, dict_datasets) {
  POSCORAD <- load_dataset(ds)

  pred <- POSCORAD %>%
    rename(Time = Day, Score = SCORAD) %>%
    select(Patient, Time, Score) %>%
    drop_na() %>%
    mutate(Iteration = get_fc_iteration(.data$Time, perf_horizon))

  pred <- lapply(unique(pred[["Iteration"]]), function(i) {
    split_fc_dataset(pred, i)$Testing
  }) %>%
    bind_rows() %>%
    select(-LastScore, -LastTime) %>%
    filter(.data$Horizon <= perf_horizon)

  max_day <- get_max_train_day(ds, dict_datasets)

  pred <- pred %>%
    filter((.data$Time - .data$Horizon) <= max_day)

  list(
    pred = pred,
    POSCORAD = POSCORAD,
    max_day = max_day
  )
}

attach_eczemapred_item_predictions <- function(ds, pred_sc, perf_horizon) {
  for (severity_item in eczemapred_items()) {
    f <- eczemapred_prediction_file(ds, severity_item, perf_horizon)

    if (!file.exists(f)) {
      stop("Missing EczemaPred item prediction file: ", f)
    }

    res_item <- readRDS(f) %>%
      dplyr::filter(.data$Horizon <= perf_horizon) %>%
      dplyr::select(Patient, Time, Horizon, Iteration, Samples) %>%
      dplyr::rename(!!paste0(severity_item, "_pred") := Samples)

    dup <- res_item %>%
      dplyr::count(.data$Patient, .data$Time, .data$Horizon, .data$Iteration) %>%
      dplyr::filter(.data$n > 1)

    if (nrow(dup) > 0) {
      stop(
        "Duplicate EczemaPred prediction rows found for ",
        severity_item,
        " using Patient/Time/Horizon/Iteration keys."
      )
    }

    pred_sc <- dplyr::left_join(
      pred_sc,
      res_item,
      by = c("Patient", "Time", "Horizon", "Iteration")
    )
  }

  pred_sc %>% tidyr::drop_na()
}

build_scorad_samples <- function(pred_sc) {
  pred_sc <- pred_sc %>%
    mutate(
      B_pred = pmap(
        list(
          redness_pred,
          dryness_pred,
          swelling_pred,
          oozing_pred,
          scratching_pred,
          thickening_pred
        ),
        function(...) Reduce(`+`, list(...))
      ),
      C_pred = pmap(
        list(itching_pred, sleep_pred),
        function(...) Reduce(`+`, list(...))
      ),
      oSCORAD_pred = pmap(
        list(extent_pred, B_pred),
        function(x, y) 0.2 * x + 3.5 * y
      ),
      SCORAD_pred = pmap(
        list(oSCORAD_pred, C_pred),
        function(x, y) x + y
      )
    )

  pred_sc %>%
    rename(Samples = SCORAD_pred) %>%
    select(Patient, Time, Score, Horizon, Iteration, Samples)
}
