# -----------------------------------------------------------------------
# Shared scalability plotting helpers
# -----------------------------------------------------------------------

sc_default_items <- function() {
  c(
    "dryness", "extent", "itching", "redness", "oozing",
    "sleep", "swelling", "thickening", "scratching"
  )
}

sc_print_console_table <- function(x, width = 180) {
  old_width <- getOption("width")
  on.exit(options(width = old_width), add = TRUE)

  options(width = width)

  print(
    as.data.frame(x),
    row.names = FALSE,
    right = FALSE
  )
}

sc_dataset_labels <- function() {
  c(
    Derexyl = "Dataset 1",
    PFDC = "Dataset 2"
  )
}

sc_item_labels <- function() {
  c(
    dryness = "Dryness",
    extent = "Extent",
    itching = "Itch",
    redness = "Redness",
    oozing = "Oozing",
    sleep = "Sleep",
    swelling = "Swelling",
    thickening = "Thickening",
    scratching = "Scratching"
  )
}

sc_model_colours <- function() {
  c(
    "EczemaPred" = "#0072B2",
    "XemaPred" = "#D55E00",
    "Patient-specific XemaPred" = "#556B2F"
  )
}

sc_comparison_settings <- function(comparison) {
  if (comparison == "xemapred_vs_eczemapred") {
    return(list(
      comparison = comparison,
      model_a = "EczemaPred",
      model_b = "XemaPred",
      model_order = c("EczemaPred", "XemaPred"),
      out_stub = "xemapred_vs_eczemapred",
      title = "XemaPred vs EczemaPred",
      unit = "min",
      unit_divisor = 60,
      y_label = "Full-update runtime (min)"
    ))
  }

  if (comparison == "xemapred_vs_patient_specific") {
    return(list(
      comparison = comparison,
      model_a = "XemaPred",
      model_b = "Patient-specific XemaPred",
      model_order = c("XemaPred", "Patient-specific XemaPred"),
      out_stub = "xemapred_vs_patient_specific",
      title = "Patient-specific XemaPred vs population-level XemaPred",
      unit = "s",
      unit_divisor = 1,
      y_label = "Runtime (s)"
    ))
  }

  if (comparison == "eczemapred_vs_patient_specific") {
    return(list(
      comparison = comparison,
      model_a = "EczemaPred",
      model_b = "Patient-specific XemaPred",
      model_order = c("EczemaPred", "Patient-specific XemaPred"),
      out_stub = "eczemapred_vs_patient_specific",
      title = "Patient-specific XemaPred vs EczemaPred",
      unit = "min",
      unit_divisor = 60,
      y_label = "Runtime (min)"
    ))
  }

  if (comparison == "all_three_models") {
    return(list(
      comparison = comparison,
      model_a = "EczemaPred",
      model_b = "Patient-specific XemaPred",
      model_order = c("EczemaPred", "XemaPred", "Patient-specific XemaPred"),
      out_stub = "all_three_models",
      title = "EczemaPred vs XemaPred vs patient-specific XemaPred",
      unit = "min",
      unit_divisor = 60,
      y_label = "Runtime (min)"
    ))
  }

  stop(
    "Unknown scalability comparison: ", comparison,
    "\nUse 'xemapred_vs_eczemapred', 'xemapred_vs_patient_specific', 'eczemapred_vs_patient_specific', or 'all_three_models'."
  )
}

sc_dir_create <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

sc_comp_time_dir <- function(ds, result_root = "results") {
  here::here(result_root, ds, "ALL_MODELS", "comp_time")
}

sc_infer_step_days_from_filename <- function(path) {
  b <- basename(path)

  if (grepl("every4days", b, ignore.case = TRUE)) return(4L)
  if (grepl("everyday", b, ignore.case = TRUE)) return(1L)

  m <- regmatches(b, regexpr("every[0-9]+days", b, ignore.case = TRUE))

  if (length(m) && nzchar(m)) {
    as.integer(gsub("[^0-9]", "", m))
  } else {
    1L
  }
}

sc_add_iter_plot <- function(df, step_days = 1L) {
  use_iter_day <- "iter_day" %in% names(df) &&
    any(!is.na(suppressWarnings(as.integer(df$iter_day))))

  df$iter_plot <- if (use_iter_day) {
    suppressWarnings(as.integer(df$iter_day))
  } else {
    suppressWarnings(as.integer(df$iter)) * as.integer(step_days)
  }

  df
}

sc_max_plot_iter <- function(ds) {
  if (ds == "PFDC") return(82L)
  if (ds == "Derexyl") return(117L)

  stop("Unknown dataset: ", ds)
}

sc_pick_col <- function(df, candidates) {
  nm <- names(df)
  hit <- candidates[tolower(candidates) %in% tolower(nm)]

  if (!length(hit)) return(NULL)

  nm[match(tolower(hit[[1]]), tolower(nm))]
}

sc_load_rdata_dataframe <- function(path) {
  e <- new.env(parent = emptyenv())
  load(path, envir = e)

  if (exists("comp_time", envir = e)) {
    return(tibble::as_tibble(get("comp_time", envir = e)))
  }

  for (nm in ls(e)) {
    x <- get(nm, envir = e)
    if (inherits(x, c("data.frame", "tbl_df", "tbl"))) {
      return(tibble::as_tibble(x))
    }
  }

  stop("No data frame found in: ", path)
}

sc_find_comp_time_file <- function(
    ds,
    model_prefixes,
    result_root = "results",
    schedule = "everyday",
    strict_schedule = TRUE
) {
  dir_ct <- sc_comp_time_dir(ds, result_root = result_root)

  if (!dir.exists(dir_ct)) {
    return(NA_character_)
  }

  for (prefix in model_prefixes) {
    f_default <- file.path(
      dir_ct,
      paste0("comp_time_", prefix, "_", schedule, "_", ds, ".RData")
    )

    if (file.exists(f_default)) {
      return(f_default)
    }
  }

  if (isTRUE(strict_schedule)) {
    return(NA_character_)
  }

  all_files <- list.files(dir_ct, "^comp_time_.*\\.RData$", full.names = TRUE)

  if (!length(all_files)) {
    return(NA_character_)
  }

  base <- basename(all_files)

  for (prefix in model_prefixes) {
    pattern <- paste0("^comp_time_", prefix, "_")
    hit <- all_files[grepl(pattern, base)]

    if (length(hit)) {
      return(hit[[1]])
    }
  }

  NA_character_
}

sc_load_cohort_runtime <- function(
    ds,
    model_label,
    model_prefixes,
    result_root = "results",
    items = sc_default_items(),
    time_priority = c("compute_time", "run_time", "runtime", "elapsed")
) {
  f <- sc_find_comp_time_file(
    ds = ds,
    model_prefixes = model_prefixes,
    result_root = result_root
  )

  if (is.na(f)) {
    warning("[WARN] No comp_time file found for ", model_label, " / ", ds)
    return(tibble::tibble())
  }

  cat("[LOAD][", model_label, "][", ds, "] ", f, "\n", sep = "")

  df <- sc_load_rdata_dataframe(f)

  score_col <- sc_pick_col(df, c("score", "item", "variable", "question", "score_name"))
  iter_col <- sc_pick_col(df, c("iter", "iteration", "iter_idx", "it", "day", "t", "k"))
  iter_day_col <- sc_pick_col(df, c("iter_day", "day_idx"))

  time_col <- NULL
  for (candidate in time_priority) {
    candidate_hit <- sc_pick_col(df, candidate)
    if (!is.null(candidate_hit)) {
      time_col <- candidate_hit
      break
    }
  }

  if (is.null(score_col) || is.null(iter_col) || is.null(time_col)) {
    warning("[WARN] Missing score/iter/runtime column in: ", f)
    return(tibble::tibble())
  }

  step_days <- if ("timing_step_days" %in% names(df)) {
    x <- suppressWarnings(as.integer(df$timing_step_days[[1]]))
    if (is.na(x) || x < 1) 1L else x
  } else {
    sc_infer_step_days_from_filename(f)
  }

  out <- tibble::tibble(
    Dataset = ds,
    Model = model_label,
    score = tolower(as.character(df[[score_col]])),
    iter = suppressWarnings(as.integer(df[[iter_col]])),
    iter_day = if (!is.null(iter_day_col)) {
      suppressWarnings(as.integer(df[[iter_day_col]]))
    } else {
      NA_integer_
    },
    runtime_sec = suppressWarnings(as.numeric(df[[time_col]]))
  ) %>%
    dplyr::filter(
      !is.na(.data$score),
      !is.na(.data$iter),
      is.finite(.data$runtime_sec),
      .data$runtime_sec > 0
    )

  if (!is.null(items) && length(items)) {
    out <- out %>% dplyr::filter(.data$score %in% tolower(items))
  }

  out <- sc_add_iter_plot(out, step_days = step_days)

  out %>%
    dplyr::group_by(.data$Dataset, .data$Model, .data$score, .data$iter, .data$iter_plot) %>%
    dplyr::summarise(
      runtime_sec = dplyr::first(.data$runtime_sec[is.finite(.data$runtime_sec) & .data$runtime_sec > 0]),
      .groups = "drop"
    )
}

sc_load_patient_specific_runtime <- function(
    ds,
    result_root = "results",
    items = sc_default_items(),
    model_prefixes = c("SoloFast", "PatientSpecificXemaPred"),
    model_label = "Patient-specific XemaPred"
) {
  dir_ct <- sc_comp_time_dir(ds, result_root = result_root)

  if (!dir.exists(dir_ct)) {
    warning("[WARN] comp_time directory missing: ", dir_ct)
    return(tibble::tibble())
  }

  all_files <- list.files(dir_ct, "^comp_time_.*\\.RData$", full.names = TRUE)

  if (!length(all_files)) {
    return(tibble::tibble())
  }

  base <- basename(all_files)

  keep <- rep(FALSE, length(all_files))

  for (prefix in model_prefixes) {
    keep <- keep | grepl(
      paste0("^comp_time_", prefix, "_(everyday|every[0-9]+days)_"),
      base
    )
  }

  files <- all_files[keep]

  if (!length(files)) {
    warning("[WARN] No patient-specific comp_time file found for ", ds)
    return(tibble::tibble())
  }

  cat("[LOAD][", model_label, "][", ds, "]\n", sep = "")
  print(basename(files))

  purrr::map_dfr(files, function(f) {
    df <- sc_load_rdata_dataframe(f)

    ds_col <- sc_pick_col(df, c("Dataset"))
    score_col <- sc_pick_col(df, c("score", "Score", "item", "Item"))
    iter_col <- sc_pick_col(df, c("iter", "Iter", "iteration", "Iteration", "t", "T", "k", "K", "day", "Day"))
    iter_day_col <- sc_pick_col(df, c("iter_day", "Iter_day", "day_idx", "Day_idx"))
    patient_col <- sc_pick_col(df, c("Patient", "patient", "id", "ID"))

    time_col <- sc_pick_col(df, c("run_time", "run_raw", "compute_time", "run_time_total", "runtime", "elapsed"))

    if (is.null(ds_col) || is.null(score_col) || is.null(iter_col) || is.null(patient_col) || is.null(time_col)) {
      warning("[WARN] Missing required patient-specific runtime columns in: ", f)
      return(tibble::tibble())
    }

    step_days <- if ("timing_step_days" %in% names(df)) {
      x <- suppressWarnings(as.integer(df$timing_step_days[[1]]))
      if (is.na(x) || x < 1) 1L else x
    } else {
      sc_infer_step_days_from_filename(f)
    }

    out <- tibble::tibble(
      Dataset = as.character(df[[ds_col]]),
      Model = model_label,
      score = tolower(as.character(df[[score_col]])),
      iter = suppressWarnings(as.integer(df[[iter_col]])),
      iter_day = if (!is.null(iter_day_col)) {
        suppressWarnings(as.integer(df[[iter_day_col]]))
      } else {
        NA_integer_
      },
      patient = as.character(df[[patient_col]]),
      runtime_sec = suppressWarnings(as.numeric(df[[time_col]]))
    ) %>%
      dplyr::filter(
        .data$Dataset == ds,
        .data$score %in% tolower(items),
        !is.na(.data$iter),
        !is.na(.data$patient),
        is.finite(.data$runtime_sec),
        .data$runtime_sec > 0
      )

    sc_add_iter_plot(out, step_days = step_days)
  })
}

sc_load_runtime_for_comparison <- function(
    comparison,
    datasets,
    result_root = "results",
    items = sc_default_items()
) {
  settings <- sc_comparison_settings(comparison)

  out <- list()

  if ("EczemaPred" %in% settings$model_order) {
    out$eczemapred <- purrr::map_dfr(datasets, function(ds) {
      sc_load_cohort_runtime(
        ds = ds,
        model_label = "EczemaPred",
        model_prefixes = c("EczemaPred"),
        result_root = result_root,
        items = items,
        time_priority = c("run_time", "compute_time", "runtime", "elapsed")
      )
    })
  }

  if ("XemaPred" %in% settings$model_order) {
    out$xemapred <- purrr::map_dfr(datasets, function(ds) {
      sc_load_cohort_runtime(
        ds = ds,
        model_label = "XemaPred",
        model_prefixes = c("XemaPred", "Fast"),
        result_root = result_root,
        items = items,
        time_priority = c("compute_time", "run_time", "runtime", "elapsed")
      )
    })
  }

  if ("Patient-specific XemaPred" %in% settings$model_order) {
    out$patient_specific <- purrr::map_dfr(datasets, function(ds) {
      sc_load_patient_specific_runtime(
        ds = ds,
        result_root = result_root,
        items = items
      )
    })
  }

  dplyr::bind_rows(out) %>%
    dplyr::filter(.data$Model %in% settings$model_order)
}

sc_cap_runtime_days <- function(runtime_df) {
  runtime_df %>%
    dplyr::rowwise() %>%
    dplyr::mutate(iter_cap = sc_max_plot_iter(as.character(.data$Dataset))) %>%
    dplyr::ungroup() %>%
    dplyr::filter(.data$iter_plot <= .data$iter_cap) %>%
    dplyr::select(-iter_cap)
}

sc_make_scalability_df <- function(runtime_df, settings) {
  if (!"patient" %in% names(runtime_df)) {
    runtime_df$patient <- NA_character_
  }

  cohort_models <- setdiff(settings$model_order, "Patient-specific XemaPred")

  cohort_latency <- runtime_df %>%
    dplyr::filter(.data$Model %in% cohort_models) %>%
    dplyr::group_by(.data$Dataset, .data$Model, .data$iter_plot) %>%
    dplyr::summarise(
      update_latency_sec = max(.data$runtime_sec, na.rm = TRUE),
      n_scores = dplyr::n_distinct(.data$score),
      .groups = "drop"
    )

  patient_latency <- tibble::tibble()

  if ("Patient-specific XemaPred" %in% settings$model_order) {
    patient_latency <- runtime_df %>%
      dplyr::filter(.data$Model == "Patient-specific XemaPred") %>%
      dplyr::group_by(.data$Dataset, .data$iter_plot, .data$patient) %>%
      dplyr::summarise(
        patient_latency_sec = max(.data$runtime_sec, na.rm = TRUE),
        n_scores = dplyr::n_distinct(.data$score),
        .groups = "drop"
      ) %>%
      dplyr::filter(
        is.finite(.data$patient_latency_sec),
        .data$patient_latency_sec > 0
      ) %>%
      dplyr::group_by(.data$Dataset, .data$iter_plot) %>%
      dplyr::summarise(
        update_latency_sec = max(.data$patient_latency_sec, na.rm = TRUE),
        median_patient_latency_sec = median(.data$patient_latency_sec, na.rm = TRUE),
        q25_patient_latency_sec = stats::quantile(.data$patient_latency_sec, 0.25, na.rm = TRUE, names = FALSE),
        q75_patient_latency_sec = stats::quantile(.data$patient_latency_sec, 0.75, na.rm = TRUE, names = FALSE),
        q95_patient_latency_sec = stats::quantile(.data$patient_latency_sec, 0.95, na.rm = TRUE, names = FALSE),
        n_patients = dplyr::n_distinct(.data$patient),
        .groups = "drop"
      ) %>%
      dplyr::mutate(Model = "Patient-specific XemaPred")
  }

  dplyr::bind_rows(
    cohort_latency,
    patient_latency
  ) %>%
    dplyr::mutate(
      Dataset = factor(.data$Dataset, levels = c("PFDC", "Derexyl")),
      Model = factor(.data$Model, levels = settings$model_order)
    ) %>%
    dplyr::filter(.data$Model %in% settings$model_order)
}

sc_make_item_scalability_df <- function(runtime_df, settings, items = sc_default_items()) {
  if (!"patient" %in% names(runtime_df)) {
    runtime_df$patient <- NA_character_
  }

  cohort_models <- setdiff(settings$model_order, "Patient-specific XemaPred")

  cohort_item_latency <- runtime_df %>%
    dplyr::filter(.data$Model %in% cohort_models) %>%
    dplyr::group_by(.data$Dataset, .data$score, .data$Model, .data$iter_plot) %>%
    dplyr::summarise(
      update_latency_sec = max(.data$runtime_sec, na.rm = TRUE),
      n_runs = dplyr::n(),
      .groups = "drop"
    )

  patient_item_latency <- tibble::tibble()

  if ("Patient-specific XemaPred" %in% settings$model_order) {
    patient_item_latency <- runtime_df %>%
      dplyr::filter(.data$Model == "Patient-specific XemaPred") %>%
      dplyr::group_by(.data$Dataset, .data$score, .data$iter_plot) %>%
      dplyr::summarise(
        update_latency_sec = max(.data$runtime_sec, na.rm = TRUE),
        median_patient_item_sec = median(.data$runtime_sec, na.rm = TRUE),
        q25_patient_item_sec = stats::quantile(.data$runtime_sec, 0.25, na.rm = TRUE, names = FALSE),
        q75_patient_item_sec = stats::quantile(.data$runtime_sec, 0.75, na.rm = TRUE, names = FALSE),
        q95_patient_item_sec = stats::quantile(.data$runtime_sec, 0.95, na.rm = TRUE, names = FALSE),
        n_patients = dplyr::n_distinct(.data$patient),
        .groups = "drop"
      ) %>%
      dplyr::mutate(Model = "Patient-specific XemaPred")
  }

  item_labels <- sc_item_labels()

  dplyr::bind_rows(
    cohort_item_latency,
    patient_item_latency
  ) %>%
    dplyr::mutate(
      Dataset = factor(.data$Dataset, levels = c("PFDC", "Derexyl")),
      score = factor(.data$score, levels = items),
      score_label = factor(
        dplyr::recode(as.character(.data$score), !!!item_labels),
        levels = item_labels[items]
      ),
      Model = factor(.data$Model, levels = settings$model_order)
    ) %>%
    dplyr::filter(.data$Model %in% settings$model_order)
}

sc_speedup_summary <- function(scalability_df, settings) {
  model_a <- settings$model_a
  model_b <- settings$model_b

  wide <- scalability_df %>%
    dplyr::select(Dataset, iter_plot, Model, update_latency_sec) %>%
    tidyr::pivot_wider(
      names_from = Model,
      values_from = update_latency_sec
    )

  if (!all(c(model_a, model_b) %in% names(wide))) {
    stop("Cannot compute speedup. Missing columns: ", paste(setdiff(c(model_a, model_b), names(wide)), collapse = ", "))
  }

  wide %>%
    dplyr::filter(
      !is.na(.data[[model_a]]),
      !is.na(.data[[model_b]]),
      .data[[model_a]] > 0,
      .data[[model_b]] > 0
    ) %>%
    dplyr::mutate(
      speedup_fold = .data[[model_a]] / .data[[model_b]],
      runtime_reduction_pct = 100 * (1 - .data[[model_b]] / .data[[model_a]])
    ) %>%
    dplyr::group_by(.data$Dataset) %>%
    dplyr::summarise(
      n_updates = dplyr::n(),
      first_iter = min(.data$iter_plot),
      final_iter = max(.data$iter_plot),
      first_model_a = .data[[model_a]][which.min(.data$iter_plot)],
      final_model_a = .data[[model_a]][which.max(.data$iter_plot)],
      first_model_b = .data[[model_b]][which.min(.data$iter_plot)],
      final_model_b = .data[[model_b]][which.max(.data$iter_plot)],
      max_model_a = max(.data[[model_a]], na.rm = TRUE),
      max_model_b = max(.data[[model_b]], na.rm = TRUE),
      median_speedup = median(.data$speedup_fold, na.rm = TRUE),
      mean_speedup = mean(.data$speedup_fold, na.rm = TRUE),
      max_speedup = max(.data$speedup_fold, na.rm = TRUE),
      final_speedup = .data$speedup_fold[which.max(.data$iter_plot)],
      median_runtime_reduction_pct = median(.data$runtime_reduction_pct, na.rm = TRUE),
      final_runtime_reduction_pct = .data$runtime_reduction_pct[which.max(.data$iter_plot)],
      model_a_growth_fold = .data$final_model_a / .data$first_model_a,
      model_b_growth_fold = .data$final_model_b / .data$first_model_b,
      max_runtime_speedup = .data$max_model_a / .data$max_model_b,
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      Dataset_label = dplyr::recode(as.character(.data$Dataset), !!!sc_dataset_labels()),
      model_a = model_a,
      model_b = model_b,
      unit = settings$unit,
      final_model_a_label = paste0(round(.data$final_model_a / settings$unit_divisor, 3), " ", settings$unit),
      final_model_b_label = paste0(round(.data$final_model_b / settings$unit_divisor, 3), " ", settings$unit),
      max_model_a_label = paste0(round(.data$max_model_a / settings$unit_divisor, 3), " ", settings$unit),
      max_model_b_label = paste0(round(.data$max_model_b / settings$unit_divisor, 3), " ", settings$unit),
      median_speedup_label = paste0(round(.data$median_speedup, 1), "x"),
      final_speedup_label = paste0(round(.data$final_speedup, 1), "x"),
      max_runtime_speedup_label = paste0(round(.data$max_runtime_speedup, 1), "x"),
      final_reduction_label = paste0(round(.data$final_runtime_reduction_pct, 1), "% lower runtime")
    )
}

sc_recent_speedup_summary <- function(
    scalability_df,
    settings,
    n_recent = 10L
) {

  model_a <- settings$model_a
  model_b <- settings$model_b

  wide <- scalability_df %>%
    dplyr::select(
      Dataset,
      iter_plot,
      N_training_observations,
      Model,
      update_latency_sec
    ) %>%
    tidyr::pivot_wider(
      names_from = Model,
      values_from = update_latency_sec
    ) %>%
    dplyr::filter(
      !is.na(.data[[model_a]]),
      !is.na(.data[[model_b]]),
      .data[[model_a]] > 0,
      .data[[model_b]] > 0
    ) %>%
    dplyr::mutate(
      runtime_ratio = .data[[model_a]] / .data[[model_b]]
    )

  recent <- wide %>%
    dplyr::group_by(.data$Dataset) %>%
    dplyr::arrange(.data$iter_plot, .by_group = TRUE) %>%
    dplyr::slice_tail(n = n_recent) %>%
    dplyr::ungroup()

  summary <- recent %>%
    dplyr::group_by(.data$Dataset) %>%
    dplyr::summarise(
      n_updates = dplyr::n(),

      first_day = min(.data$iter_plot),
      last_day  = max(.data$iter_plot),

      median_model_a_sec = median(.data[[model_a]], na.rm = TRUE),
      q25_model_a_sec = stats::quantile(
        .data[[model_a]], 0.25, na.rm = TRUE, names = FALSE
      ),
      q75_model_a_sec = stats::quantile(
        .data[[model_a]], 0.75, na.rm = TRUE, names = FALSE
      ),

      median_model_b_sec = median(.data[[model_b]], na.rm = TRUE),
      q25_model_b_sec = stats::quantile(
        .data[[model_b]], 0.25, na.rm = TRUE, names = FALSE
      ),
      q75_model_b_sec = stats::quantile(
        .data[[model_b]], 0.75, na.rm = TRUE, names = FALSE
      ),

      median_runtime_ratio = median(.data$runtime_ratio, na.rm = TRUE),
      q25_runtime_ratio = stats::quantile(
        .data$runtime_ratio, 0.25, na.rm = TRUE, names = FALSE
      ),
      q75_runtime_ratio = stats::quantile(
        .data$runtime_ratio, 0.75, na.rm = TRUE, names = FALSE
      ),

      .groups = "drop"
    ) %>%
    dplyr::mutate(
      Dataset_label =
        dplyr::recode(as.character(.data$Dataset), !!!sc_dataset_labels())
    )

  list(
    raw = recent,
    summary = summary
  )
}

sc_fmt_runtime <- function(x_sec, settings, digits = 3) {
  out <- x_sec / settings$unit_divisor
  paste0(round(out, digits), " ", settings$unit)
}

sc_fmt_fold <- function(x, digits = 1) {
  paste0(round(x, digits), "x")
}

sc_fmt_pct_lower <- function(x, digits = 1) {
  paste0(round(x, digits), "% lower")
}

sc_speedup_summary_display <- function(speedup_summary, settings) {
  out <- speedup_summary %>%
    dplyr::arrange(factor(.data$Dataset_label, levels = c("Dataset 1", "Dataset 2"))) %>%
    dplyr::transmute(
      Dataset = .data$Dataset_label,
      `Updates` = .data$n_updates,
      `Training-day range` = paste0(.data$first_iter, "-", .data$final_iter),
      final_model_a = sc_fmt_runtime(.data$final_model_a, settings),
      final_model_b = sc_fmt_runtime(.data$final_model_b, settings),
      `Final speed-up` = sc_fmt_fold(.data$final_speedup),
      `Final reduction` = sc_fmt_pct_lower(.data$final_runtime_reduction_pct),
      `Median speed-up` = sc_fmt_fold(.data$median_speedup),
      max_model_a = sc_fmt_runtime(.data$max_model_a, settings),
      max_model_b = sc_fmt_runtime(.data$max_model_b, settings),
      `Peak runtime ratio` = sc_fmt_fold(.data$max_runtime_speedup),
      model_a_growth = sc_fmt_fold(.data$model_a_growth_fold),
      model_b_growth = sc_fmt_fold(.data$model_b_growth_fold)
    )

  names(out)[names(out) == "final_model_a"] <- paste0("Final ", settings$model_a)
  names(out)[names(out) == "final_model_b"] <- paste0("Final ", settings$model_b)
  names(out)[names(out) == "max_model_a"] <- paste0("Peak ", settings$model_a)
  names(out)[names(out) == "max_model_b"] <- paste0("Peak ", settings$model_b)
  names(out)[names(out) == "model_a_growth"] <- paste0("Runtime growth: ", settings$model_a)
  names(out)[names(out) == "model_b_growth"] <- paste0("Runtime growth: ", settings$model_b)

  out
}

sc_max_runtime_timing_display <- function(max_runtime_timing, settings) {
  max_runtime_timing %>%
    dplyr::arrange(
      factor(.data$Dataset_label, levels = c("Dataset 1", "Dataset 2")),
      factor(.data$Model, levels = settings$model_order)
    ) %>%
    dplyr::transmute(
      Dataset = .data$Dataset_label,
      Model = as.character(.data$Model),
      `Maximum runtime` = .data$max_runtime_label,
      `Timing day` = .data$max_runtime_day,
      `Timing rule` = .data$timing_note
    )
}

sc_max_runtime_timing <- function(
    scalability_df,
    settings,
    plateau_model = "XemaPred",
    plateau_tolerance_pct = 1
) {
  scalability_df %>%
    dplyr::group_by(.data$Dataset, .data$Model) %>%
    dplyr::group_modify(function(df, key) {
      max_runtime_sec <- max(df$update_latency_sec, na.rm = TRUE)

      df <- df %>%
        dplyr::mutate(
          max_runtime_sec = max_runtime_sec,
          pct_below_max = 100 * (.data$max_runtime_sec - .data$update_latency_sec) / .data$max_runtime_sec
        )

      if (as.character(key$Model[[1]]) %in% plateau_model) {
        chosen <- df %>%
          dplyr::filter(.data$pct_below_max <= plateau_tolerance_pct) %>%
          dplyr::slice_min(.data$iter_plot, n = 1, with_ties = FALSE)
      } else {
        chosen <- df %>%
          dplyr::slice_max(.data$update_latency_sec, n = 1, with_ties = FALSE)
      }

      chosen
    }) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      Dataset_label = dplyr::recode(as.character(.data$Dataset), !!!sc_dataset_labels()),
      max_runtime = .data$max_runtime_sec / settings$unit_divisor,
      max_runtime_label = paste0(round(.data$max_runtime, 3), " ", settings$unit),
      timing_note = dplyr::if_else(
        as.character(.data$Model) %in% plateau_model & .data$pct_below_max > 0,
        paste0("earliest day within ", plateau_tolerance_pct, "% of max"),
        "observed maximum"
      )
    ) %>%
    dplyr::select(
      Dataset,
      Dataset_label,
      Model,
      max_runtime,
      max_runtime_label,
      max_runtime_day = iter_plot,
      pct_below_max,
      timing_note
    )
}

sc_nice_runtime_label <- function(x_sec) {
  dplyr::case_when(
    is.na(x_sec) ~ NA_character_,
    x_sec >= 3600 ~ paste0(round(x_sec / 3600, 2), " h"),
    x_sec >= 60 ~ paste0(round(x_sec / 60, 2), " min"),
    x_sec >= 1 ~ paste0(round(x_sec, 2), " s"),
    TRUE ~ paste0(round(1000 * x_sec, 1), " ms")
  )
}

sc_runtime_plateau_timing <- function(
    scalability_df,
    model_filter = "XemaPred",
    tolerance_pct = 1
) {
  scalability_df %>%
    dplyr::filter(as.character(.data$Model) %in% model_filter) %>%
    dplyr::group_by(.data$Dataset, .data$Model) %>%
    dplyr::mutate(
      max_runtime_sec = max(.data$update_latency_sec, na.rm = TRUE),
      pct_below_max = 100 * (.data$max_runtime_sec - .data$update_latency_sec) / .data$max_runtime_sec,
      within_plateau = .data$pct_below_max <= tolerance_pct
    ) %>%
    dplyr::filter(.data$within_plateau) %>%
    dplyr::slice_min(.data$iter_plot, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      Dataset_label = dplyr::recode(as.character(.data$Dataset), !!!sc_dataset_labels()),
      plateau_runtime_label = sc_nice_runtime_label(.data$update_latency_sec),
      max_runtime_label = sc_nice_runtime_label(.data$max_runtime_sec),
      pct_below_max_label = paste0(round(.data$pct_below_max, 3), "% below max")
    ) %>%
    dplyr::select(
      Dataset = .data$Dataset_label,
      Model = .data$Model,
      `Plateau day` = .data$iter_plot,
      `Runtime at plateau` = .data$plateau_runtime_label,
      `Observed maximum` = .data$max_runtime_label,
      `Difference from maximum` = .data$pct_below_max_label
    )
}

sc_dual_axis_spec <- function(ds, it_min, it_max) {
  if (ds == "PFDC") {
    top_days <- c(1, 9, 17, 25, 37, 45, 53, 65, 73, 81)
    bottom_labs <- c(0, 300, 600, 900)
  } else if (ds == "Derexyl") {
    top_days <- c(1, 9, 21, 29, 41, 53, 65, 77, 89, 102, 115)
    bottom_labs <- c(0, 2500, 5000, 7500)
  } else {
    top_days <- numeric(0)
    bottom_labs <- numeric(0)
  }

  top_breaks <- top_days - 1L
  keep <- top_breaks >= it_min & top_breaks <= it_max

  bottom_breaks <- as.integer(round(seq(it_min, it_max, length.out = length(bottom_labs))))
  bottom_breaks <- unique(pmin(pmax(bottom_breaks, it_min), it_max))

  list(
    breaks_bottom = bottom_breaks,
    labels_bottom = bottom_labs[seq_along(bottom_breaks)],
    breaks_top = top_breaks[keep],
    labels_top = top_days[keep]
  )
}

sc_runtime_theme <- function(base_size = 14) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.placement = "outside",
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "top",
      legend.title = ggplot2::element_blank(),
      legend.box = "horizontal",
      axis.text.x.top = ggplot2::element_text(size = 10),
      axis.title.x.top = ggplot2::element_text(size = 11, face = "bold"),
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      plot.tag = ggplot2::element_text(face = "bold", size = 14),
      plot.tag.position = c(0.02, 0.90)
    )
}

sc_running_max <- function(df) {
  df %>%
    dplyr::arrange(.data$iter_plot) %>%
    dplyr::mutate(runtime_running_max = cummax(.data$runtime_plot))
}

sc_offset_running_max_line <- function(df) {
  df %>%
    dplyr::mutate(
      iter_plot_step = .data$N_training_observations,
      runtime_running_max_step = .data$runtime_running_max
    )
}

sc_get_y_limits <- function(df_plot, settings, logy = TRUE, y_axis_mode = c("fixed", "free")) {
  y_axis_mode <- match.arg(y_axis_mode)

  if (y_axis_mode == "free") {
    return(NULL)
  }

  vals <- df_plot$update_latency_sec / settings$unit_divisor
  vals <- vals[is.finite(vals) & vals > 0]

  if (!length(vals)) {
    return(NULL)
  }

  if (logy) {
    c(10 ^ floor(log10(min(vals))), 10 ^ ceiling(log10(max(vals))))
  } else {
    c(0, max(vals) * 1.05)
  }
}

sc_log_grid <- function(ds = NULL) {
  if (!is.null(ds) && ds == "PFDC") {
    return(c(0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100, 300, 1000, 3000))
  }

  c(0.001, 0.01, 0.1, 1, 10, 100, 1000)
}

sc_log10_limits_dynamic <- function(values, ds = NULL, min_break = NULL) {
  values <- values[is.finite(values) & values > 0]
  grid <- sc_log_grid(ds)

  if (!length(values)) {
    lower_fallback <- if (is.null(min_break) || !is.finite(min_break)) min(grid) else min_break
    return(c(lower_fallback, 100))
  }

  min_y <- min(values, na.rm = TRUE)
  max_y <- max(values, na.rm = TRUE)

  if (is.null(min_break) || !is.finite(min_break)) {
    lower_break <- max(grid[grid <= min_y])
    if (!is.finite(lower_break)) {
      lower_break <- 10 ^ floor(log10(min_y))
    }
  } else {
    lower_break <- min_break
  }

  upper_break <- min(grid[grid >= max_y])
  if (!is.finite(upper_break)) {
    upper_break <- 10 ^ ceiling(log10(max_y))
  }

  c(lower_break, upper_break)
}

sc_log10_breaks_dynamic <- function(limits, ds = NULL, min_break = NULL) {
  values <- limits[is.finite(limits) & limits > 0]
  grid <- sc_log_grid(ds)

  if (!length(values)) {
    lower_fallback <- if (is.null(min_break) || !is.finite(min_break)) min(grid) else min_break
    return(grid[grid >= lower_fallback & grid <= 100])
  }

  min_y <- min(values, na.rm = TRUE)
  max_y <- max(values, na.rm = TRUE)

  if (is.null(min_break) || !is.finite(min_break)) {
    lower_break <- max(grid[grid <= min_y])
    if (!is.finite(lower_break)) {
      lower_break <- 10 ^ floor(log10(min_y))
    }
  } else {
    lower_break <- min_break
  }

  upper_break <- min(grid[grid >= max_y])
  if (!is.finite(upper_break)) {
    upper_break <- 10 ^ ceiling(log10(max_y))
  }

  grid[grid >= lower_break & grid <= upper_break]
}

sc_log_axis_labels <- function(x) {
  vapply(x, function(z) {
    if (is.na(z) || !is.finite(z)) return(NA_character_)

    if (z < 1) {
      lab <- formatC(z, format = "f", digits = 3)
      lab <- sub("0+$", "", lab)
      lab <- sub("\\.$", "", lab)
      return(lab)
    }

    formatC(z, format = "f", digits = 1)
  }, character(1))
}

sc_use_custom_runtime_log_axis <- function(settings) {
  settings$comparison %in% c(
    "xemapred_vs_eczemapred",
    "xemapred_vs_patient_specific",
    "all_three_models"
  )
}


sc_max_across_items_log_limits <- function(settings) {
  if (settings$comparison == "xemapred_vs_eczemapred") {
    return(c(0.1, 1000))
  }

  if (settings$comparison == "xemapred_vs_patient_specific") {
    return(c(0.1, 1000))
  }

  if (settings$comparison == "all_three_models") {
    return(c(0.01, 1000))
  }

  NULL
}

sc_max_across_items_log_breaks <- function(settings) {
  if (settings$comparison == "xemapred_vs_eczemapred") {
    return(c(0.1, 1, 10, 100, 1000))
  }

  if (settings$comparison == "xemapred_vs_patient_specific") {
    return(c(0.1, 1, 10, 100, 1000))
  }

  if (settings$comparison == "all_three_models") {
    return(c(0.01, 0.1, 1, 10, 100, 1000))
  }

  ggplot2::waiver()
}

sc_item_level_log_breaks <- function(ds, limits, settings = NULL) {
  if (!is.null(settings) && settings$comparison == "all_three_models") {
    grid <- c(0.01, 0.1, 1, 10, 100, 1000)

  } else if (!is.null(settings) && settings$comparison == "xemapred_vs_patient_specific") {
    grid <- c(0.1, 1, 10, 100, 1000)

  } else if (!is.null(settings) && settings$comparison == "xemapred_vs_eczemapred") {
    grid <- if (ds == "Derexyl") {
      c(0.1, 1, 10, 100, 1000)
    } else if (ds == "PFDC") {
      c(0.1, 0.3, 1, 3, 10)
    } else {
      c(0.1, 1, 10, 100, 1000)
    }

  } else {
    grid <- c(0.1, 1, 10, 100, 1000)
  }

  limits <- limits[is.finite(limits) & limits > 0]

  if (!length(limits)) {
    return(grid)
  }

  upper <- max(limits, na.rm = TRUE)

  out <- grid[grid <= upper * 1.000001]

  if (!length(out)) {
    out <- min(grid)
  }

  out
}

sc_plot_scalability_dataset <- function(
    df_plot,
    ds,
    settings,
    logy = TRUE,
    y_axis_mode = c("fixed", "free"),
    show_y_label = TRUE,
    y_limits = NULL
) {
  y_axis_mode <- match.arg(y_axis_mode)

  dsd <- df_plot %>%
    dplyr::filter(
      as.character(.data$Dataset) == ds
    )

  if (!nrow(dsd)) {
    return(NULL)
  }



  dsd <- dsd %>%
    dplyr::mutate(
      runtime_value = .data$update_latency_sec / settings$unit_divisor,
      runtime_plot = if (logy) pmax(.data$runtime_value, 1e-8) else .data$runtime_value,
      Model = factor(as.character(.data$Model), levels = settings$model_order)
    )

  dsd_env <- dsd %>%
    dplyr::group_by(.data$Model) %>%
    dplyr::group_modify(~ sc_running_max(.x)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(.data$Model) %>%
    dplyr::group_modify(~ sc_offset_running_max_line(.x)) %>%
    dplyr::ungroup()

  # -----------------------------------------------------------------------
  # Runtime-ratio annotations for Dataset 1
  # -----------------------------------------------------------------------

  ratio_annotations <- NULL

  if (settings$comparison == "xemapred_vs_eczemapred") {

    reference_days <- if (ds == "Derexyl") {
      c(29L, 53L, 89L)
    } else if (ds == "PFDC") {
      c(25L, 53L, 73L)
    } else {
      integer(0)
    }

    ratio_annotations <- dsd %>%
      dplyr::filter(
        .data$iter_plot %in% reference_days,
        as.character(.data$Model) %in% c("EczemaPred", "XemaPred")
      ) %>%
      dplyr::select(
        iter_plot,
        N_training_observations,
        Model,
        runtime_plot
      ) %>%
      tidyr::pivot_wider(
        names_from = Model,
        values_from = runtime_plot
      ) %>%
      dplyr::filter(
        is.finite(EczemaPred),
        is.finite(XemaPred),
        EczemaPred > 0,
        XemaPred > 0
      ) %>%
      dplyr::mutate(
        ratio = EczemaPred / XemaPred,
        y_low = pmin(EczemaPred, XemaPred),
        y_high = pmax(EczemaPred, XemaPred),

        # Geometric midpoint for log-scale y-axis
        y_label = sqrt(y_low * y_high),

        label = paste0(round(ratio), "\u00d7")
      )
  }

  p <- ggplot2::ggplot(dsd, ggplot2::aes(x = .data$N_training_observations, colour = .data$Model)) +
    ggplot2::geom_line(
      ggplot2::aes(y = .data$runtime_plot),
      linewidth = 0.35,
      alpha = 0.25
    ) +
    ggplot2::geom_point(
      ggplot2::aes(y = .data$runtime_plot),
      size = 1.1,
      alpha = 0.45
    ) +
    ggplot2::geom_step(
      data = dsd_env,
      ggplot2::aes(
        x = .data$iter_plot_step,
        y = .data$runtime_running_max_step
      ),
      linewidth = 1.05,
      direction = "hv",
      alpha = 0.95
    ) +
    ggplot2::scale_colour_manual(
      values = sc_model_colours(),
      limits = settings$model_order,
      breaks = settings$model_order,
      drop = FALSE
    ) +
    sc_observation_axis(ds) +
    ggplot2::labs(
      tag = sc_dataset_labels()[[ds]],
      x = "Number of training observations",
      y = if (show_y_label) settings$y_label else NULL,
      colour = NULL
    ) +
    sc_runtime_theme()

  if (!is.null(ratio_annotations) && nrow(ratio_annotations) > 0) {

    x_offset <- 0.005 * diff(range(dsd$N_training_observations, na.rm = TRUE))

    p <- p +
      ggplot2::geom_segment(
        data = ratio_annotations,
        ggplot2::aes(
          x = .data$N_training_observations,
          xend = .data$N_training_observations,
          y = .data$y_low,
          yend = .data$y_high
        ),
        inherit.aes = FALSE,
        linewidth = 0.55,
        colour = "black",
        arrow = grid::arrow(
          ends = "both",
          type = "closed",
          length = grid::unit(0.08, "inches")
        )
      ) +
      ggplot2::geom_label(
        data = ratio_annotations,
        ggplot2::aes(
          x = .data$N_training_observations + x_offset,
          y = .data$y_label,
          label = .data$label
        ),
        inherit.aes = FALSE,
        size = 3.4,
        fontface = "bold",
        label.size = 0,
        fill = "white"
      )
  }

  if (logy) {
    if (sc_use_custom_runtime_log_axis(settings)) {
      ds_name <- as.character(ds)

      axis_values <- if (ds_name == "PFDC") {
        dsd$runtime_value
      } else if (!is.null(y_limits)) {
        y_limits
      } else {
        dsd$runtime_value
      }

      dynamic_limits <- sc_log10_limits_dynamic(
        axis_values,
        ds = ds_name,
        min_break = NULL
      )

      dynamic_breaks <- sc_log10_breaks_dynamic(
        dynamic_limits,
        ds = ds_name,
        min_break = NULL
      )

      if (!is.null(y_limits) && sc_use_custom_runtime_log_axis(settings)) {
        p <- p +
          ggplot2::scale_y_log10(
            breaks = sc_max_across_items_log_breaks(settings),
            labels = sc_log_axis_labels
          ) +
          ggplot2::coord_cartesian(
            ylim = y_limits
          )
      } else {
        p <- p +
          ggplot2::scale_y_log10(
            breaks = dynamic_breaks,
            labels = sc_log_axis_labels
          ) +
          ggplot2::coord_cartesian(
            ylim = dynamic_limits
          )
      }
    } else {
      p <- p +
        ggplot2::scale_y_log10(
          limits = y_limits,
          labels = scales::label_number(accuracy = 0.001)
        )
    }
  } else {
    p <- p +
      ggplot2::scale_y_continuous(
        limits = y_limits,
        labels = scales::label_number(accuracy = 0.01)
      )
  }

  p
}

sc_plot_scalability_combined <- function(
    df_plot,
    settings,
    logy = TRUE,
    y_axis_mode = c("fixed", "free")
) {
  y_axis_mode <- match.arg(y_axis_mode)

  y_limits <- sc_get_y_limits(
    df_plot = df_plot,
    settings = settings,
    logy = logy,
    y_axis_mode = y_axis_mode
  )

  if (
    isTRUE(logy) &&
    y_axis_mode == "fixed" &&
    sc_use_custom_runtime_log_axis(settings)
  ) {
    y_limits <- sc_max_across_items_log_limits(settings)
  }

  p1 <- sc_plot_scalability_dataset(
    df_plot,
    ds = "Derexyl",
    settings = settings,
    logy = logy,
    y_axis_mode = y_axis_mode,
    show_y_label = TRUE,
    y_limits = y_limits
  )

  p2 <- sc_plot_scalability_dataset(
    df_plot,
    ds = "PFDC",
    settings = settings,
    logy = logy,
    y_axis_mode = y_axis_mode,
    show_y_label = FALSE,
    y_limits = y_limits
  )

  p1 <- p1 +
    ggplot2::theme(
      plot.tag.position = c(0.15, 0.90)
    )

  p2 <- p2 +
    ggplot2::theme(
      plot.tag.position = c(0.15, 0.90)
    )

  if (is.null(p1) || is.null(p2)) {
    return(NULL)
  }

  (p1 | p2) +
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(
      legend.position = "top",
      legend.box = "horizontal"
    )
}

sc_observation_axis <- function(ds, t_horizon = 4) {

  training_axis <- make_training_axis_breaks(
    dataset = ds,
    t_horizon = t_horizon,
    n_breaks = 10
  )

  ggplot2::scale_x_continuous(
    sec.axis = ggplot2::dup_axis(
      breaks = training_axis$N,
      labels = training_axis$LastTime,
      name = "Training days"
    )
  )
}

sc_plot_item_scalability_dataset <- function(
    item_df,
    ds,
    settings,
    items = sc_default_items(),
    logy = TRUE
) {

  dsd <- item_df %>%
    dplyr::filter(
      as.character(.data$Dataset) == ds
    )

  if (!nrow(dsd)) {
    return(NULL)
  }

  dsd <- dsd %>%
    dplyr::mutate(
      runtime_value = .data$update_latency_sec / settings$unit_divisor,
      runtime_plot = if (logy) pmax(.data$runtime_value, 1e-8) else .data$runtime_value,
      score = factor(.data$score, levels = items),
      score_label = factor(.data$score_label, levels = sc_item_labels()[items]),
      Model = factor(as.character(.data$Model), levels = settings$model_order)
    )

  dsd_max <- dsd %>%
    dplyr::group_by(.data$score_label, .data$Model) %>%
    dplyr::arrange(.data$iter_plot, .by_group = TRUE) %>%
    dplyr::mutate(runtime_running_max = cummax(.data$runtime_plot)) %>%
    dplyr::group_modify(~ sc_offset_running_max_line(.x)) %>%
    dplyr::ungroup()

  p <- ggplot2::ggplot() +
    ggplot2::geom_line(
      data = dsd,
      ggplot2::aes(
        x = .data$N_training_observations,
        y = .data$runtime_plot,
        colour = .data$Model,
        group = .data$Model
      ),
      linewidth = 0.25,
      alpha = 0.18
    ) +
    ggplot2::geom_point(
      data = dsd,
      ggplot2::aes(
        x = .data$N_training_observations,
        y = .data$runtime_plot,
        colour = .data$Model,
        group = .data$Model
      ),
      size = 0.55,
      alpha = 0.35
    ) +
    ggplot2::geom_step(
      data = dsd_max,
      ggplot2::aes(
        x = .data$iter_plot_step,
        y = .data$runtime_running_max_step,
        colour = .data$Model,
        group = .data$Model
      ),
      linewidth = 0.85,
      alpha = 0.98,
      direction = "vh"
    ) +
    ggplot2::facet_wrap(
      ~ score_label,
      ncol = 3,
      scales = "free_y",
      drop = FALSE,
      strip.position = "top"
    ) +
    ggplot2::scale_colour_manual(
      values = sc_model_colours(),
      limits = settings$model_order,
      breaks = settings$model_order,
      drop = FALSE
    ) +
    sc_observation_axis(ds) +
    ggplot2::labs(
      tag = sc_dataset_labels()[[ds]],
      x = "Number of training observations",
      y = paste0("Per-update runtime (", settings$unit, ")"),
      colour = NULL
    ) +
    sc_runtime_theme() +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold", size = 11),
      legend.position = "top"
    )

    if (logy) {
      if (sc_use_custom_runtime_log_axis(settings)) {
        p <- p +
          ggplot2::scale_y_log10(
            breaks = function(x) sc_item_level_log_breaks(ds, x, settings),
            labels = sc_log_axis_labels
          )
      } else {
        p <- p +
          ggplot2::scale_y_log10(
            labels = scales::label_number(accuracy = 0.001)
          )
      }
    } else {
      p <- p +
        ggplot2::scale_y_continuous(
          labels = scales::label_number(accuracy = 0.01)
        )
    }

  p
}

sc_plot_item_scalability_combined <- function(
    item_df,
    settings,
    items = sc_default_items(),
    logy = TRUE
) {
  p1 <- sc_plot_item_scalability_dataset(
    item_df,
    ds = "Derexyl",
    settings = settings,
    items = items,
    logy = logy
  )

  p2 <- sc_plot_item_scalability_dataset(
    item_df,
    ds = "PFDC",
    settings = settings,
    items = items,
    logy = logy
  )

  if (is.null(p1) || is.null(p2)) {
    return(NULL)
  }

  p1 <- p1 +
    ggplot2::theme(
      plot.tag = ggplot2::element_text(face = "bold", size = 16),
      plot.tag.position = c(0.105, 0.95)
    )

  p2 <- p2 +
    ggplot2::theme(
      plot.tag = ggplot2::element_text(face = "bold", size = 16),
      plot.tag.position = c(0.105, 0.99)
    )

  (p1 / p2) +
    patchwork::plot_layout(ncol = 1, guides = "collect") &
    ggplot2::theme(
      legend.position = "top",
      legend.box = "horizontal"
    )
}

sc_save_plot <- function(plot, out_dir, filename_stem, width, height, formats = c("png", "pdf")) {
  if (is.null(plot)) {
    return(invisible(FALSE))
  }

  save_plot(
    plot = plot,
    out_dir = out_dir,
    filename_stem = filename_stem,
    width = width,
    height = height,
    dpi = 300,
    formats = formats,
    verbose = FALSE
  )

  cat("[SAVED] ", file.path(out_dir, filename_stem), ".*\n", sep = "")

  invisible(TRUE)
}

run_scalability_analysis <- function(
    comparison,
    datasets = c("PFDC", "Derexyl"),
    items = sc_default_items(),
    result_root = "results",
    plot_root = "plots",
    horizon = 1,
    logy = TRUE,
    y_axis_modes = c("fixed", "free"),
    save_item_level = TRUE,
    formats = c("png", "pdf")
) {
  settings <- sc_comparison_settings(comparison)

  out_dir <- get_plot_dir(
    "scalability",
    settings$out_stub,
    plot_root = plot_root
  )

  csv_dir <- file.path(out_dir, "csv")
  derexyl_dir <- file.path(out_dir, "Derexyl")
  pfdc_dir <- file.path(out_dir, "PFDC")
  combined_dir <- file.path(out_dir, "combined")

  purrr::walk(
    c(csv_dir, derexyl_dir, pfdc_dir, combined_dir),
    sc_dir_create
  )

  csv_dir <- file.path(out_dir, "csv")
  sc_dir_create(csv_dir)

  cat("============================================================\n")
  cat("[SCALABILITY] ", settings$title, "\n", sep = "")
  cat("============================================================\n")

    runtime_df <- sc_load_runtime_for_comparison(
    comparison = comparison,
    datasets = datasets,
    result_root = result_root,
    items = items
    ) %>%
    sc_cap_runtime_days()

    if (!"patient" %in% names(runtime_df)) {
    runtime_df$patient <- NA_character_
    }

  scalability_df <- sc_make_scalability_df(runtime_df, settings)
  item_scalability_df <- sc_make_item_scalability_df(runtime_df, settings, items = items)

  scalability_df <- rs_attach_training_observation_counts(
    speed_df = scalability_df,
    datasets = datasets,
    items = items,
    horizon_for_counts = horizon
  )

  item_scalability_df <- rs_attach_training_observation_counts(
    speed_df = item_scalability_df,
    datasets = datasets,
    items = items,
    horizon_for_counts = horizon
  )

  stopifnot(
    all(is.finite(scalability_df$N_training_observations)),
    all(is.finite(item_scalability_df$N_training_observations))
  )
  speedup_summary <- sc_speedup_summary(scalability_df, settings)
  max_runtime_timing <- sc_max_runtime_timing(scalability_df, settings)

  # -----------------------------------------------------------------------
  # Reference-day runtime ratios for figure annotations
  # -----------------------------------------------------------------------

  reference_days <- c(29L, 53L, 89L)

  reference_speedups <- scalability_df %>%
    dplyr::filter(
      as.character(.data$Dataset) == "Derexyl",
      .data$iter_plot %in% reference_days,
      as.character(.data$Model) %in% c("EczemaPred", "XemaPred")
    ) %>%
    dplyr::select(
      iter_plot,
      N_training_observations,
      Model,
      update_latency_sec
    ) %>%
    tidyr::pivot_wider(
      names_from = Model,
      values_from = update_latency_sec
    ) %>%
    dplyr::mutate(
      runtime_ratio = EczemaPred / XemaPred
    ) %>%
    dplyr::arrange(.data$iter_plot)

  cat("\n============================================================\n")
  cat("[REFERENCE-DAY RUNTIME RATIOS]\n")
  cat("============================================================\n")
  print(reference_speedups)
  cat("\n")

  recent_speedup <- sc_recent_speedup_summary(
    scalability_df = scalability_df,
    settings = settings,
    n_recent = 10L
  )

  recent_speedup_display <- recent_speedup$summary %>%
    dplyr::transmute(
      Dataset = .data$Dataset_label,
      `Late-stage days` = paste0(.data$first_day, "-", .data$last_day),

      `Median EczemaPred` =
        sc_fmt_runtime(.data$median_model_a_sec, settings),

      `Median XemaPred` =
        sc_fmt_runtime(.data$median_model_b_sec, settings),

      `Median runtime ratio` =
        paste0(round(.data$median_runtime_ratio, 1), "x"),

      `Runtime-ratio IQR` =
        paste0(
          round(.data$q25_runtime_ratio, 1),
          "-",
          round(.data$q75_runtime_ratio, 1),
          "x"
        )
    )

  cat("\n============================================================\n")
  cat("[LAST 10 UPDATES]\n")
  cat("Median values across the final 10 matched updates.\n")
  cat("Runtime ratio = EczemaPred / XemaPred.\n\n")

  sc_print_console_table(
    recent_speedup_display,
    width = 180
  )

  cat("\n")  

  readr::write_csv(
    runtime_df,
    file.path(csv_dir, paste0("runtime_raw_", settings$out_stub, "_H", horizon, ".csv"))
  )

  readr::write_csv(
    scalability_df,
    file.path(csv_dir, paste0("scalability_latency_", settings$out_stub, "_H", horizon, ".csv"))
  )

  readr::write_csv(
    speedup_summary,
    file.path(csv_dir, paste0("speedup_summary_", settings$out_stub, "_H", horizon, ".csv"))
  )

  readr::write_csv(
    max_runtime_timing,
    file.path(csv_dir, paste0("max_runtime_timing_", settings$out_stub, "_H", horizon, ".csv"))
  )

  readr::write_csv(
    item_scalability_df,
    file.path(csv_dir, paste0("item_level_scalability_", settings$out_stub, "_H", horizon, ".csv"))
  )

  speedup_summary_display <- sc_speedup_summary_display(speedup_summary, settings)
  max_runtime_timing_display <- sc_max_runtime_timing_display(max_runtime_timing, settings)

  readr::write_csv(
    speedup_summary_display,
    file.path(csv_dir, paste0("speedup_summary_readable_", settings$out_stub, "_H", horizon, ".csv"))
  )

  readr::write_csv(
    max_runtime_timing_display,
    file.path(csv_dir, paste0("max_runtime_timing_readable_", settings$out_stub, "_H", horizon, ".csv"))
  )

  cat("\n============================================================\n")
  cat("[SPEED-UP SUMMARY]\n")
  cat("============================================================\n")
  cat("Speed-up is calculated as ", settings$model_a, " runtime / ", settings$model_b, " runtime.\n\n", sep = "")

  speedup_final_display <- speedup_summary_display %>%
    dplyr::select(
      Dataset,
      Updates,
      `Training-day range`,
      dplyr::starts_with("Final "),
      `Final speed-up`,
      `Final reduction`
    )

  speedup_peak_display <- speedup_summary_display %>%
    dplyr::select(
      Dataset,
      `Median speed-up`,
      dplyr::starts_with("Peak "),
      `Peak runtime ratio`,
      dplyr::starts_with("Runtime growth:")
    )

  cat("[Final update]\n")
  sc_print_console_table(speedup_final_display, width = 180)
  cat("\n")

  cat("[Overall / peak behaviour]\n")
  sc_print_console_table(speedup_peak_display, width = 180)
  cat("\n")

  cat("\n============================================================\n")
  cat("[MAXIMUM RUNTIME TIMING]\n")
  cat("============================================================\n")
  cat("Maximum runtime is the highest per-update runtime reached by each model.\n\n")
  sc_print_console_table(max_runtime_timing_display, width = 180)
  cat("\n")

  scale_suffix <- ifelse(logy, "log", "linear")

  for (ds in datasets) {
    p_ds <- sc_plot_scalability_dataset(
      scalability_df,
      ds = ds,
      settings = settings,
      logy = logy,
      y_axis_mode = "free",
      show_y_label = TRUE
    ) +
    ggplot2::labs(tag = NULL)

    ds_out_dir <- if (ds == "Derexyl") derexyl_dir else pfdc_dir

    sc_save_plot(
      p_ds,
      out_dir = ds_out_dir,
      filename_stem = paste0("max_across_items_", settings$out_stub, "_", ds, "_H", horizon, "_", scale_suffix),
      width = 9,
      height = 5.5,
      formats = formats
    )
}

  for (ym in y_axis_modes) {
    p_combined <- sc_plot_scalability_combined(
      scalability_df,
      settings = settings,
      logy = logy,
      y_axis_mode = ym
    )

    sc_save_plot(
      p_combined,
      out_dir = combined_dir,
      filename_stem = paste0("max_across_items_", settings$out_stub, "_combined_H", horizon, "_", scale_suffix, "_", ym, "Y"),
      width = 14,
      height = 5.5,
      formats = formats
    )
  }

  if (isTRUE(save_item_level)) {

    for (ds in datasets) {
      p_item_ds <- sc_plot_item_scalability_dataset(
        item_scalability_df,
        ds = ds,
        settings = settings,
        items = items,
        logy = logy
      )  +
      ggplot2::labs(tag = NULL)

      ds_out_dir <- if (ds == "Derexyl") derexyl_dir else pfdc_dir

      sc_save_plot(
        p_item_ds,
        out_dir = ds_out_dir,
        filename_stem = paste0("item_level_", settings$out_stub, "_", ds, "_H", horizon, "_", scale_suffix),
        width = 12,
        height = 8,
        formats = formats
      )
    }

    p_item_combined <- sc_plot_item_scalability_combined(
      item_scalability_df,
      settings = settings,
      items = items,
      logy = logy
    )

    sc_save_plot(
      p_item_combined,
      out_dir = combined_dir,
      filename_stem = paste0("item_level_", settings$out_stub, "_combined_H", horizon, "_", scale_suffix),
      width = 12,
      height = 16,
      formats = formats
    )
  }

  cat("\n[DONE] Saved scalability outputs to: ", out_dir, "\n", sep = "")

  invisible(list(
    runtime_df = runtime_df,
    scalability_df = scalability_df,
    item_scalability_df = item_scalability_df,
    speedup_summary = speedup_summary,
    speedup_summary_display = speedup_summary_display,
    max_runtime_timing = max_runtime_timing,
    max_runtime_timing_display = max_runtime_timing_display,
    out_dir = out_dir
  ))
}