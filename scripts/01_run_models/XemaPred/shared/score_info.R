# -----------------------------------------------------------------------
# XemaPred score metadata helpers
#
# Purpose:
# - Centralise PO-SCORAD score metadata for XemaPred.
# - Map scores to XemaPred model families:
#   intensity, subjective, extent.
# - Compute score scaling quantities used by the particle filters.
#
# Dependencies:
# - detail_POSCORAD() must be available from scripts/00_setup/00_init.R.
# -----------------------------------------------------------------------

get_xemapred_score_groups <- function() {
  list(
    extent = c("extent"),
    intensity = c(
      "dryness",
      "redness",
      "swelling",
      "oozing",
      "thickening",
      "scratching"
    ),
    subjective = c("itching", "sleep")
  )
}

get_xemapred_model_family <- function(score) {
  groups <- get_xemapred_score_groups()

  if (score %in% groups$intensity) {
    return("intensity")
  }

  if (score %in% groups$subjective) {
    return("subjective")
  }

  if (score %in% groups$extent) {
    return("extent")
  }

  stop("[ERROR] Unknown XemaPred score: ", score)
}

get_expected_xemapred_item_model <- function(score) {
  family <- get_xemapred_model_family(score)

  if (family == "intensity") {
    return("OrderedRW")
  }

  if (family == "subjective") {
    return("BinRW")
  }

  if (family == "extent") {
    return("BinMC")
  }

  stop("[ERROR] Unknown XemaPred model family: ", family)
}

get_xemapred_score_info <- function(score) {
  item_dict <- detail_POSCORAD()

  if (!score %in% item_dict[["Name"]]) {
    stop("[ERROR] Unknown PO-SCORAD score: ", score)
  }

  item_row <- item_dict |>
    dplyr::filter(.data$Name == score)

  item_label <- as.character(item_row[["Label"]])
  max_score <- as.numeric(item_row[["Maximum"]])

  reso <- if (score %in% c("SCORAD", "oSCORAD")) {
    1
  } else {
    as.numeric(item_row[["Resolution"]])
  }

  M_max <- as.integer(round(max_score / reso))
  K <- M_max + 1L

  list(
    score = score,
    item_label = item_label,
    model_family = get_xemapred_model_family(score),
    expected_model = get_expected_xemapred_item_model(score),
    max_score = max_score,
    reso = reso,
    M_max = M_max,
    K = K
  )
}

print_xemapred_score_info <- function(score_info) {
  cat(
    glue::glue(
      "[INFO] Score mapping:\n",
      "  score: {score_info$score}\n",
      "  family: {score_info$model_family}\n",
      "  item label: {score_info$item_label}\n",
      "  expected model: {score_info$expected_model}\n",
      "  max score: {score_info$max_score}\n",
      "  resolution: {score_info$reso}\n",
      "  M_max: {score_info$M_max}\n",
      "  K: {score_info$K}\n"
    )
  )

  invisible(score_info)
}