# -----------------------------------------------------------------------
# Score metadata and family dispatch
# -----------------------------------------------------------------------

intensity_scores <- function() {
  c("dryness", "redness", "swelling", "oozing", "thickening", "scratching")
}

subjective_scores <- function() {
  c("sleep", "itching")
}

get_score_family <- function(score) {
  if (score %in% intensity_scores()) {
    "intensity"
  } else if (score %in% subjective_scores()) {
    "subjective"
  } else if (identical(score, "extent")) {
    "extent"
  } else {
    stop(glue::glue("[ERROR] Unknown score '{score}' (family)."))
  }
}

get_score_mapping <- function(score) {
  item_dict <- detail_POSCORAD()
  if (!score %in% item_dict[["Name"]]) {
    stop(glue::glue("[ERROR] Unknown score '{score}'."))
  }

  item_row <- item_dict %>% dplyr::filter(.data$Name == score)
  item_lbl <- as.character(item_row[["Label"]])
  max_score <- as.numeric(item_row[["Maximum"]])
  reso <- if (score %in% c("SCORAD", "oSCORAD")) 1 else as.numeric(item_row[["Resolution"]])

  M_max <- as.integer(round(max_score / reso))
  K <- M_max + 1L

  list(
    item_label = item_lbl,
    max_score = max_score,
    reso = reso,
    M_max = M_max,
    K = K,
    model_family = get_score_family(score)
  )
}
