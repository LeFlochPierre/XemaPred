# -----------------------------------------------------------------------
# Shared plotting configuration
# -----------------------------------------------------------------------

DATASETS_DEFAULT <- c("Derexyl", "PFDC")

DATASET_LABELS <- c(
  Derexyl = "Dataset 1",
  PFDC    = "Dataset 2"
)

DATASET_COLS <- c(
  "Dataset 1" = "#7F6DFF",
  "Dataset 2" = "#35B873"
)

INTENSITY_ITEMS <- c(
  "dryness",
  "redness",
  "swelling",
  "oozing",
  "scratching",
  "thickening"
)

WIDE_ITEMS <- c(
  "extent",
  "itching",
  "sleep"
)

ITEM_ORDER_DISTRIBUTION <- c(
  INTENSITY_ITEMS,
  WIDE_ITEMS
)

ITEM_ORDER_MAIN <- c(
  "extent",
  "redness",
  "dryness",
  "swelling",
  "oozing",
  "scratching",
  "thickening",
  "itching",
  "sleep"
)

ITEM_ORDER_WITH_POSCORAD <- c(
  ITEM_ORDER_MAIN,
  "SCORAD"
)

ITEM_LABELS <- c(
  extent     = "Extent",
  redness    = "Redness",
  dryness    = "Dryness",
  swelling   = "Swelling",
  oozing     = "Oozing",
  scratching = "Scratching",
  thickening = "Thickening",
  itching    = "Itch",
  sleep      = "Sleep",
  SCORAD     = "PO-SCORAD"
)

dataset_label <- function(x) {
  dplyr::recode(x, !!!DATASET_LABELS)
}

item_label <- function(x) {
  dplyr::recode(x, !!!ITEM_LABELS)
}

parse_csv_arg <- function(x) {
  x <- gsub("\\s+", "", x)
  if (!nzchar(x)) character(0) else strsplit(x, ",")[[1]]
}

parse_logical_arg <- function(x) {
  x <- tolower(as.character(x))
  if (x %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (x %in% c("false", "f", "0", "no", "n")) return(FALSE)
  stop("Could not parse logical argument: ", x)
}
