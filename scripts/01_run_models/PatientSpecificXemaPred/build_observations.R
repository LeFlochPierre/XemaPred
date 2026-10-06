# -----------------------------------------------------------------------
# Small observation/result builders
# -----------------------------------------------------------------------

empty_prediction_df <- function() {
  data.frame(
    Patient = integer(0),
    Time = integer(0),
    Horizon = integer(0),
    Iteration = integer(0),
    Score = numeric(0),
    lpd = numeric(0),
    RPS = numeric(0),
    y_pred = numeric(0),
    Samples = I(list()),
    stringsAsFactors = FALSE
  )
}

empty_diagnostic_df <- function() {
  data.frame(
    score = character(0),
    model = character(0),
    Patient = integer(0),
    iter = integer(0),
    run_time = numeric(0),
    compute_time = numeric(0),
    forecast_time = numeric(0),
    score_time = numeric(0),
    io_time = numeric(0),
    run_time_total = numeric(0),
    stringsAsFactors = FALSE
  )
}

error_prediction_df <- function(pid) {
  data.frame(
    Patient = pid,
    Time = NA_integer_,
    Horizon = NA_integer_,
    Iteration = NA_integer_,
    Score = NA_real_,
    lpd = NA_real_,
    RPS = NA_real_,
    y_pred = NA_real_,
    Samples = I(list(NULL)),
    stringsAsFactors = FALSE
  )
}

error_diagnostic_df <- function(pid, run_info) {
  data.frame(
    score = run_info$score,
    model = run_info$model,
    Patient = pid,
    iter = NA_integer_,
    run_time = NA_real_,
    compute_time = NA_real_,
    forecast_time = NA_real_,
    score_time = NA_real_,
    io_time = NA_real_,
    run_time_total = NA_real_,
    stringsAsFactors = FALSE
  )
}

get_existing_iterations <- function(iters_p) {
  if (is.null(iters_p) || !dir.exists(iters_p)) return(integer(0))

  existing_files <- list.files(iters_p, pattern = "^iter-[0-9]{4}\\.rds$", full.names = FALSE)
  as.integer(gsub("iter-([0-9]{4})\\.rds", "\\1", existing_files))
}

get_missing_iterations <- function(train_it, iters_p, save_mode, state) {
  if (identical(save_mode, "all")) {
    existing_its <- get_existing_iterations(iters_p)
    setdiff(train_it, existing_its)
  } else {
    train_it
  }
}
