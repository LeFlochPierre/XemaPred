# -----------------------------------------------------------------------
# Parallel setup helpers
# -----------------------------------------------------------------------

make_patient_cluster <- function(n_cluster) {
  if (n_cluster <= 1L) return(NULL)

  cl <- NULL
  try({ cl <- parallel::makeForkCluster(n_cluster, outfile = "") }, silent = TRUE)

  if (is.null(cl)) {
    cat("[WARN] makeForkCluster failed; falling back to PSOCK makeCluster.\n")
    cl <- parallel::makeCluster(n_cluster, outfile = "")
  }

  doParallel::registerDoParallel(cl)
  cl
}

stop_patient_cluster <- function(cl) {
  if (!is.null(cl)) parallel::stopCluster(cl)
}

get_parallel_exports <- function() {
  setdiff(ls(envir = .GlobalEnv), c("cl"))
}
