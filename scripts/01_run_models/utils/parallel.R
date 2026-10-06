# -----------------------------------------------------------------------
# Shared parallel execution utility
#
# Purpose:
# - Identify missing forward-chaining iterations.
# - Run only missing iterations in parallel.
# - Re-source required files inside each worker.
# -----------------------------------------------------------------------

run_missing_iterations_parallel <- function(
    data_obj,
    run_info,
    paths,
    fit_function_name,
    worker_sources
) {
  require_posterior <- identical(
    run_info$model_family,
    "eczemapred_item_models"
  )

  missing_its <- get_missing_iterations(
    data_obj,
    paths,
    require_posterior = require_posterior
  )

  cat(glue::glue(
    "[INFO] Resuming with {length(missing_its)} missing iterations...\n"
  ))

  if (length(missing_its) == 0) {
    cat("[INFO] No missing iterations. Nothing to run.\n")
    return(invisible(NULL))
  }

  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(purrr)
    library(glue)
    library(foreach)
    library(doParallel)
  })

  `%dopar%` <- foreach::`%dopar%`

  # Use a job-specific PSOCK port.
  psock_port <- suppressWarnings(
    as.integer(Sys.getenv("R_PARALLEL_PORT"))
  )

  if (is.na(psock_port)) {
    psock_port <- 20000 + (Sys.getpid() %% 30000)
  }

  cat(
    "[INFO] Creating PSOCK cluster with ",
    run_info$n_cluster,
    " workers on port ",
    psock_port,
    "\n",
    sep = ""
  )

  cl <- parallel::makePSOCKcluster(
    run_info$n_cluster,
    port = psock_port,
    outfile = ""
  )

  doParallel::registerDoParallel(cl)

  on.exit({
    try(parallel::stopCluster(cl), silent = TRUE)
  }, add = TRUE)

  foreach::foreach(
    it = rev(missing_its)
  ) %dopar% {
    source(here::here(
      "scripts",
      "00_setup",
      "00_init.R"
    ))

    suppressPackageStartupMessages({
      library(dplyr)
      library(tidyr)
      library(purrr)
      library(glue)
      library(foreach)
      library(doParallel)
    })

    for (src in worker_sources) {
      source(here::here(src))
    }

    fit_fun <- get(
      fit_function_name,
      mode = "function"
    )

    fit_fun(
      it,
      data_obj,
      run_info,
      paths
    )

    gc()
    NULL
  }

  cat("[INFO] Parallel processing complete.\n")

  invisible(NULL)
}