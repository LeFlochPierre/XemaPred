# -----------------------------------------------------------------------
# Shared PO-SCORAD score metadata
# -----------------------------------------------------------------------

fallback_score_metadata <- function() {
  tibble::tibble(
    Name = c(
      "dryness",
      "redness",
      "swelling",
      "oozing",
      "scratching",
      "thickening",
      "extent",
      "itching",
      "sleep",
      "SCORAD"
    ),
    Label = c(
      "Dryness",
      "Redness",
      "Swelling",
      "Scabs/Oozing",
      "Traces of scratching",
      "Thickening",
      "Extent",
      "Itching VAS",
      "Sleep disturbance VAS",
      "SCORAD"
    ),
    Maximum = c(3, 3, 3, 3, 3, 3, 100, 10, 10, 103),
    Resolution = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
  )
}

get_score_metadata <- function() {
  if (exists("detail_POSCORAD", mode = "function")) {
    out <- detail_POSCORAD() %>%
      mutate(
        Name = as.character(.data$Name),
        Label = as.character(.data$Label)
      )

    if (!"SCORAD" %in% out$Name) {
      out <- bind_rows(
        out,
        tibble::tibble(
          Name = "SCORAD",
          Label = "SCORAD",
          Maximum = 103,
          Resolution = 1
        )
      )
    }

    return(out)
  }

  fallback_score_metadata()
}

get_item_score_column <- function(item) {
  item_dict <- get_score_metadata()

  if (!item %in% item_dict$Name) {
    stop("Unknown item: ", item)
  }

  item_dict %>%
    filter(.data$Name == item) %>%
    slice(1) %>%
    pull(.data$Label) %>%
    as.character()
}

is_intensity_item <- function(item) {
  item %in% INTENSITY_ITEMS
}

get_item_maximum <- function(item) {
  if (item == "SCORAD") return(103)

  item_dict <- get_score_metadata()

  if (!item %in% item_dict$Name) {
    stop("Unknown item: ", item)
  }

  out <- item_dict %>%
    filter(.data$Name == item) %>%
    slice(1) %>%
    pull(.data$Maximum)

  if (length(out) == 0 || is.na(out)) {
    if (item %in% INTENSITY_ITEMS) return(3)
    if (item %in% c("itching", "sleep")) return(10)
    if (item == "extent") return(100)
  }

  as.numeric(out)
}

get_item_resolution <- function(item) {
  if (item == "SCORAD") return(1)

  item_dict <- get_score_metadata()

  out <- item_dict %>%
    filter(.data$Name == item) %>%
    slice(1) %>%
    pull(.data$Resolution)

  if (length(out) == 0 || is.na(out)) return(1)

  as.numeric(out)
}

get_item_ylim <- function(item, padded = FALSE) {
  if (item == "extent") {
    return(if (padded) c(-5, 105) else c(0, 100))
  }

  if (item %in% INTENSITY_ITEMS) {
    return(if (padded) c(-0.25, 3.25) else c(0, 3))
  }

  if (item %in% c("itching", "sleep")) {
    return(if (padded) c(-0.5, 10.5) else c(0, 10))
  }

  if (item == "SCORAD") {
    return(if (padded) c(-5, 105) else c(0, 103))
  }

  c(NA_real_, NA_real_)
}

make_item_scale_table <- function(items = ITEM_ORDER_DISTRIBUTION) {
  tibble::tibble(
    Item = items,
    xmin = 0,
    xmax = vapply(items, get_item_maximum, numeric(1)),
    step = vapply(items, get_item_resolution, numeric(1))
  ) %>%
    mutate(
      Item_display = factor(
        item_label(.data$Item),
        levels = unname(ITEM_LABELS[items])
      )
    )
}
