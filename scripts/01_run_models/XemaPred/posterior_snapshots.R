# -----------------------------------------------------------------------
# Population XemaPred posterior snapshots
# -----------------------------------------------------------------------

xemapred_posterior_file <- function(paths, it) {
  file.path(
    paths$posterior_dir,
    sprintf("posterior-%04d.rds", as.integer(it))
  )
}

get_xemapred_training_meta <- function(data_obj, it) {
  train <- lapply(data_obj$patients, function(pid) {
    split_xemapred_patient_iteration(data_obj, pid, it)$train
  }) |>
    dplyr::bind_rows()

  n_total <- nrow(data_obj$df_all2)

  list(
    max_training_day = if (nrow(train)) max(train$Time, na.rm = TRUE) else NA_real_,
    n_training_observations = nrow(train),
    n_total_observations = n_total,
    training_fraction = if (n_total > 0) nrow(train) / n_total else NA_real_,
    n_training_patients = dplyr::n_distinct(train$Patient)
  )
}

summarise_weighted_theta <- function(theta, w) {
  w <- normalize_weights(w)
  mu <- colSums(sweep(theta, 1, w, "*"))

  v <- colSums(
    sweep(sweep(theta, 2, mu, "-")^2, 1, w, "*")
  )

  data.frame(
    parameter = colnames(theta),
    mean = mu,
    sd = sqrt(v),
    stringsAsFactors = FALSE
  )
}

save_xemapred_posterior <- function(
    state,
    it,
    data_obj,
    run_info,
    paths
) {
  if (is.null(state$theta_u) || !NROW(state$theta_u)) {
    stop("[ERROR] XemaPred state has no theta_u.")
  }

  if (is.null(state$W) && is.null(state$w)) {
    stop("[ERROR] XemaPred state has no outer weights.")
  }

  theta <- as.matrix(state$theta_u)
  w_raw <- as.numeric(state$W %||% state$w)

  if (is.null(colnames(theta))) {
    stop("[ERROR] theta_u has no parameter names.")
  }

  if (length(w_raw) != nrow(theta)) {
    stop("[ERROR] theta_u and W dimensions do not match.")
  }

  w <- normalize_weights(w_raw)
  ess <- 1 / sum(w^2)

  meta <- get_xemapred_training_meta(data_obj, it)

  out <- list(
    dataset = run_info$dataset,
    score = run_info$score,
    model = run_info$mdl_name,
    model_family = run_info$model_family,
    iteration = as.integer(it),

    max_training_day = meta$max_training_day,
    n_training_observations = meta$n_training_observations,
    n_total_observations = meta$n_total_observations,
    training_fraction = meta$training_fraction,
    n_training_patients = meta$n_training_patients,

    n_theta = nrow(theta),
    n_x = state$Nx %||% NA_integer_,
    parameter_names = colnames(theta),

    theta_u = theta,
    weights = w,
    weights_raw = w_raw,

    ess = ess,
    ess_pct = 100 * ess / nrow(theta),
    unique_theta = nrow(unique(round(theta, 10))),

    raw_summary = summarise_weighted_theta(theta, w)
  )

  saveRDS(
    out,
    xemapred_posterior_file(paths, it)
  )

  invisible(out)
}