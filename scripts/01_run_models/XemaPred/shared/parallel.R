# -----------------------------------------------------------------------
# XemaPred parallel backend utilities
#
# Purpose:
# - Register a foreach backend for XemaPred SMC2 updates.
# - Use a job-specific PSOCK port when running on HPC.
# - Avoid port-collision errors between simultaneous PBS jobs.
# -----------------------------------------------------------------------

setup_xemapred_parallel <- function(run_info) {
  n_cluster <- as.integer(run_info$n_cluster %||% 1L)

  if (n_cluster <= 1L) {
    foreach::registerDoSEQ()
    cat("[INFO] Using sequential foreach backend.\n")
    return(NULL)
  }

  psock_port <- suppressWarnings(as.integer(Sys.getenv("R_PARALLEL_PORT")))

  if (is.na(psock_port)) {
    psock_port <- 20000 + (Sys.getpid() %% 30000)
  }

  cat(
    "[INFO] Creating XemaPred PSOCK cluster with ",
    n_cluster,
    " workers on port ",
    psock_port,
    "\n",
    sep = ""
  )

  cl <- parallel::makePSOCKcluster(
    n_cluster,
    port = psock_port,
    outfile = ""
  )

  doParallel::registerDoParallel(cl)

  cat("[INFO] Registered foreach workers: ", foreach::getDoParWorkers(), "\n", sep = "")

  cl
}

stop_xemapred_parallel <- function(cl) {
  if (!is.null(cl)) {
    try(parallel::stopCluster(cl), silent = TRUE)
  }

  foreach::registerDoSEQ()

  invisible(NULL)
}