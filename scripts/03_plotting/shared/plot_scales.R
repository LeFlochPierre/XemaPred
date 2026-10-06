# -----------------------------------------------------------------------
# Shared plot scales and binning helpers
# -----------------------------------------------------------------------

round_to_step <- function(x, step = 1) {
  round(as.numeric(x) / step) * step
}

full_scale_breaks <- function(lims) {
  upper <- max(lims, na.rm = TRUE)

  if (upper <= 3.5) {
    return(seq(0, 3, by = 1))
  }

  if (upper <= 10.5) {
    return(seq(0, 10, by = 2))
  }

  if (upper <= 103.5) {
    return(c(seq(0, 100, by = 20), 103))
  }

  seq(0, upper, length.out = 6)
}

full_scale_labels <- function(x) {
  sprintf("%.1f", x)
}

get_x_break_interval <- function(dataset, plot_type = "patient_observed") {
  if (plot_type == "patient_observed") {
    if (dataset == "PFDC") return(10)
    if (dataset == "Derexyl") return(15)
  }

  if (dataset %in% c("PFDC", "Derexyl")) {
    return(20)
  }

  10
}

bin_score_column <- function(df, score_col = "Score") {
  df %>%
    mutate(
      max_bin_lower = floor((.data$xmax - 1e-9) / .data$bin_width) * .data$bin_width,
      bin_lower_raw = floor(.data[[score_col]] / .data$bin_width) * .data$bin_width,
      bin_lower = pmax(.data$xmin, pmin(.data$bin_lower_raw, .data$max_bin_lower)),
      bin_upper = pmin(.data$bin_lower + .data$bin_width, .data$xmax),
      bin_mid = (.data$bin_lower + .data$bin_upper) / 2
    ) %>%
    select(-max_bin_lower, -bin_lower_raw)
}

make_bin_grid <- function(scale_df, dataset_levels) {
  tidyr::expand_grid(
    Dataset_display = dataset_levels,
    scale_df
  ) %>%
    mutate(
      max_bin_lower = floor((.data$xmax - 1e-9) / .data$bin_width) * .data$bin_width
    ) %>%
    rowwise() %>%
    mutate(
      bin_lower = list(seq(.data$xmin, .data$max_bin_lower, by = .data$bin_width))
    ) %>%
    unnest(bin_lower) %>%
    ungroup() %>%
    mutate(
      bin_upper = pmin(.data$bin_lower + .data$bin_width, .data$xmax),
      bin_mid = (.data$bin_lower + .data$bin_upper) / 2
    ) %>%
    select(-max_bin_lower)
}
