# -----------------------------------------------------------------------
# Lightweight logging helpers
# -----------------------------------------------------------------------

log_info <- function(...) {
  cat("[INFO] ", paste0(..., collapse = ""), "\n", sep = "")
}

log_warn <- function(...) {
  cat("[WARN] ", paste0(..., collapse = ""), "\n", sep = "")
}

log_error <- function(...) {
  cat("[ERROR] ", paste0(..., collapse = ""), "\n", sep = "")
}
