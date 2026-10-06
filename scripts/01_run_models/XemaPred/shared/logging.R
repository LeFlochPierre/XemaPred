# -----------------------------------------------------------------------
# XemaPred logging utilities
#
# Purpose:
# - Provide consistent, readable console logs.
# - Avoid glued-together messages caused by missing newlines.
# - Support optional debug logging controlled from YAML.
# -----------------------------------------------------------------------

xlog_time <- function() {
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
}

xlog_enabled <- function() {
  isTRUE(getOption("xemapred.verbose", TRUE))
}

xlog_debug_enabled <- function() {
  isTRUE(getOption("xemapred.debug", FALSE))
}

xlog <- function(level = "INFO", ..., sep = "") {
  if (!xlog_enabled() && level != "ERROR") {
    return(invisible(NULL))
  }

  msg <- paste(..., sep = sep)

  cat(
    "[",
    level,
    "] ",
    xlog_time(),
    " | ",
    msg,
    "\n",
    sep = ""
  )

  flush.console()

  invisible(NULL)
}

xlog_debug <- function(..., sep = "") {
  if (xlog_debug_enabled()) {
    xlog("DEBUG", ..., sep = sep)
  }

  invisible(NULL)
}

xlog_section <- function(title) {
  if (!xlog_enabled()) {
    return(invisible(NULL))
  }

  cat("\n")
  cat("=======================================================================\n")
  cat(title, "\n", sep = "")
  cat("=======================================================================\n")

  flush.console()

  invisible(NULL)
}

xlog_subsection <- function(title) {
  if (!xlog_enabled()) {
    return(invisible(NULL))
  }

  cat("\n")
  cat("-----------------------------------------------------------------------\n")
  cat(title, "\n", sep = "")
  cat("-----------------------------------------------------------------------\n")

  flush.console()

  invisible(NULL)
}

xlog_kv <- function(title, values) {
  if (!xlog_enabled()) {
    return(invisible(NULL))
  }

  xlog_subsection(title)

  names_values <- names(values)

  for (nm in names_values) {
    val <- values[[nm]]

    if (length(val) == 0 || is.null(val)) {
      val <- ""
    }

    if (is.logical(val)) {
      val <- ifelse(isTRUE(val), "true", "false")
    }

    if (is.infinite(val)) {
      val <- "Inf"
    }

    cat(sprintf("  %-22s %s\n", paste0(nm, ":"), as.character(val)))
  }

  flush.console()

  invisible(NULL)
}

xlog_paths <- function(paths) {
  xlog_kv(
    "Output directories",
    list(
      base = paths$base_dir,
      iters = paths$iters_dir,
      diag = paths$diag_dir,
      fits = paths$fits_dir,
      final = paths$final_dir,
      logs = paths$logs_dir
    )
  )
}

xlog_iteration_start <- function(it, total_it = NULL) {
  if (is.null(total_it)) {
    xlog_subsection(glue::glue("Iteration {it}"))
  } else {
    xlog_subsection(glue::glue("Iteration {it} / {total_it}"))
  }
}

xlog_iteration_summary <- function(
    it,
    obs_rows,
    obs_patients,
    update_time,
    pmmh_time,
    did_resample
) {
  xlog(
    "INFO",
    glue::glue(
      "iter={it} | obs={obs_rows} | patients={obs_patients} | ",
      "update={round(update_time, 2)}s | ",
      "pmmh={ifelse(is.finite(pmmh_time), round(pmmh_time, 2), NA)}s | ",
      "resample={did_resample}"
    )
  )
}