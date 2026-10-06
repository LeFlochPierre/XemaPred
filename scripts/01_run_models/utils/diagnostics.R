# -----------------------------------------------------------------------
# Shared diagnostic utilities
#
# Purpose:
# - Check whether Stan returned valid samples.
# - Warn about poor Rhat or low effective sample size.
# - Store per-iteration runtime diagnostics in a consistent format.
# -----------------------------------------------------------------------

check_stan_fit_valid <- function(fit) {
  ok <- TRUE

  ok <- tryCatch({
    !is.null(fit@sim$samples) && length(fit@sim$samples) > 0
  }, error = function(e) {
    FALSE
  })

  if (!ok) {
    stop("[ERROR] Model failed to draw samples. Check Stan code and data quality.")
  }

  invisible(TRUE)
}

check_stan_diagnostics <- function(fit, rhat_threshold = 1.1, neff_threshold = 100) {
  summary_fit <- summary(fit)$summary

  if ("Rhat" %in% colnames(summary_fit)) {
    if (any(summary_fit[, "Rhat"] > rhat_threshold, na.rm = TRUE)) {
      warning(glue::glue("[WARNING] Rhat > {rhat_threshold} for some parameters."))
    }
  }

  if ("n_eff" %in% colnames(summary_fit)) {
    if (any(summary_fit[, "n_eff"] < neff_threshold, na.rm = TRUE)) {
      warning(glue::glue("[WARNING] n_eff < {neff_threshold} for some parameters."))
    }
  }

  invisible(summary_fit)
}

make_runtime_diag <- function(run_info, it, start, end) {
  data.frame(
    dataset = run_info$dataset,
    score = run_info$score,
    model = run_info$mdl_name,
    model_family = run_info$model_family,
    iter = it,
    run_time = as.numeric(end - start, units = "secs")
  )
}
