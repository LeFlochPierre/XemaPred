#!/usr/bin/env Rscript
# -------------------------------------------------------------------------
# Observed PO-SCORAD data distribution
#
# Output:
#   plots/dataset_description/data_distribution.png
#   plots/dataset_description/data_distribution.pdf
#   plots/dataset_description/data_distribution.svg
#   plots/dataset_description/data_distribution_summary.csv
# -------------------------------------------------------------------------

# Shared plotting setup and helpers
source(here::here("scripts", "03_plotting", "shared", "00_plot_setup.R"))
source(here::here("scripts", "03_plotting", "shared", "plot_config.R"))
source(here::here("scripts", "03_plotting", "shared", "score_metadata.R"))
source(here::here("scripts", "03_plotting", "shared", "plot_paths.R"))
source(here::here("scripts", "03_plotting", "shared", "plot_scales.R"))
source(here::here("scripts", "03_plotting", "shared", "plot_themes.R"))
source(here::here("scripts", "03_plotting", "shared", "save_plots.R"))
source(here::here("scripts", "03_plotting", "shared", "observed_data.R"))

# -------------------------------------------------------------------------
# User options
# -------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

datasets <- DATASETS_DEFAULT
plot_root <- "plots"

# Binning is only used for variables with wider native scales.
bin_widths <- c(
  extent   = 5,
  itching  = 1,
  sleep    = 1,
  poscorad = 5
)

i <- 1
while (i <= length(args)) {
  key <- args[[i]]

  if (key %in% c("-h", "--help")) {
    cat(
"plot_score_distributions.R
=======================================

Usage:
  Rscript scripts/03_plotting/dataset_description/plot_score_distributions.R
    [--datasets Derexyl,PFDC]
    [--plot_root plots]
    [--extent_bin_width 5]
    [--itching_bin_width 1]
    [--sleep_bin_width 1]
    [--poscorad_bin_width 5]
"
    )
    quit(status = 0)

  } else if (key == "--datasets") {
    datasets <- parse_csv_arg(args[[i + 1]])
    i <- i + 2

  } else if (key == "--plot_root") {
    plot_root <- args[[i + 1]]
    i <- i + 2

  } else if (key == "--extent_bin_width") {
    bin_widths[["extent"]] <- as.numeric(args[[i + 1]])
    i <- i + 2

  } else if (key == "--itching_bin_width") {
    bin_widths[["itching"]] <- as.numeric(args[[i + 1]])
    i <- i + 2

  } else if (key == "--sleep_bin_width") {
    bin_widths[["sleep"]] <- as.numeric(args[[i + 1]])
    i <- i + 2

  } else if (key == "--poscorad_bin_width") {
    bin_widths[["poscorad"]] <- as.numeric(args[[i + 1]])
    i <- i + 2

  } else {
    stop("Unknown argument: ", key)
  }
}

cat("============================================================\n")
cat("[INIT] Observed score distributions\n")
cat("============================================================\n")
cat("[INFO] datasets  : ", paste(datasets, collapse = ","), "\n", sep = "")
cat("[INFO] plot_root : ", plot_root, "\n", sep = "")
cat("[INFO] bin widths:\n")
print(bin_widths)

out_dir <- get_plot_dir("dataset_description", plot_root = plot_root)

# -------------------------------------------------------------------------
# Load observed item-level and total PO-SCORAD data
# -------------------------------------------------------------------------

item_scales <- make_item_scale_table(ITEM_ORDER_DISTRIBUTION)

items_long <- load_observed_items_long(
  datasets = datasets,
  items = ITEM_ORDER_DISTRIBUTION
) %>%
  left_join(
    item_scales %>% select(Item, xmin, xmax, step),
    by = "Item"
  ) %>%
  mutate(
    Score_plot = round_to_step(Score, step),
    Score_plot = pmin(pmax(Score_plot, xmin), xmax)
  )

poscorad_long <- load_observed_poscorad_long(datasets = datasets) %>%
  mutate(
    Score_plot = round_to_step(Score, 1),
    Score_plot = pmin(pmax(Score_plot, 0), 103)
  )

# -------------------------------------------------------------------------
# Build native-scale distributions for the six intensity signs
# -------------------------------------------------------------------------

item_grid <- tidyr::expand_grid(
  Dataset_display = factor(
    unname(DATASET_LABELS[datasets]),
    levels = unname(DATASET_LABELS[datasets])
  ),
  item_scales
) %>%
  rowwise() %>%
  mutate(Score_plot = list(seq(xmin, xmax, by = step))) %>%
  unnest(Score_plot) %>%
  ungroup()

item_counts <- items_long %>%
  count(Dataset_display, Item, Item_display, Score_plot, name = "n")

item_totals <- items_long %>%
  count(Dataset_display, Item, Item_display, name = "n_total")

item_dist <- item_grid %>%
  left_join(
    item_counts,
    by = c("Dataset_display", "Item", "Item_display", "Score_plot")
  ) %>%
  left_join(
    item_totals,
    by = c("Dataset_display", "Item", "Item_display")
  ) %>%
  mutate(
    n = replace_na(n, 0L),
    prop = n / n_total
  )

# -------------------------------------------------------------------------
# Build binned distributions for extent, itch, and sleep
# -------------------------------------------------------------------------

wide_bin_scales <- item_scales %>%
  filter(Item %in% WIDE_ITEMS) %>%
  mutate(
    bin_width = case_when(
      Item == "extent"  ~ unname(bin_widths[["extent"]]),
      Item == "itching" ~ unname(bin_widths[["itching"]]),
      Item == "sleep"   ~ unname(bin_widths[["sleep"]]),
      TRUE ~ 1
    )
  )

wide_items_binned <- items_long %>%
  filter(Item %in% WIDE_ITEMS) %>%
  left_join(
    wide_bin_scales %>% select(Item, bin_width),
    by = "Item"
  ) %>%
  bin_score_column(score_col = "Score")

wide_item_grid <- make_bin_grid(
  scale_df = wide_bin_scales,
  dataset_levels = factor(
    unname(DATASET_LABELS[datasets]),
    levels = unname(DATASET_LABELS[datasets])
  )
)

wide_item_counts <- wide_items_binned %>%
  count(
    Dataset_display,
    Item,
    Item_display,
    bin_width,
    bin_lower,
    bin_upper,
    bin_mid,
    name = "n"
  )

wide_item_totals <- wide_items_binned %>%
  count(Dataset_display, Item, Item_display, name = "n_total")

wide_item_dist <- wide_item_grid %>%
  left_join(
    wide_item_counts,
    by = c(
      "Dataset_display",
      "Item",
      "Item_display",
      "bin_width",
      "bin_lower",
      "bin_upper",
      "bin_mid"
    )
  ) %>%
  left_join(
    wide_item_totals,
    by = c("Dataset_display", "Item", "Item_display")
  ) %>%
  mutate(
    n = replace_na(n, 0L),
    prop = n / n_total
  )

# -------------------------------------------------------------------------
# Build binned distribution for total PO-SCORAD
# -------------------------------------------------------------------------

poscorad_bin_width <- unname(bin_widths[["poscorad"]])

poscorad_binned <- poscorad_long %>%
  mutate(
    xmin = 0,
    xmax = 103,
    bin_width = poscorad_bin_width
  ) %>%
  bin_score_column(score_col = "Score")

poscorad_grid <- tibble(
  Dataset_display = factor(
    unname(DATASET_LABELS[datasets]),
    levels = unname(DATASET_LABELS[datasets])
  )
) %>%
  tidyr::expand_grid(
    bin_lower = seq(
      0,
      floor((103 - 1e-9) / poscorad_bin_width) * poscorad_bin_width,
      by = poscorad_bin_width
    )
  ) %>%
  mutate(
    xmin = 0,
    xmax = 103,
    bin_width = poscorad_bin_width,
    bin_upper = pmin(bin_lower + poscorad_bin_width, 103),
    bin_mid = (bin_lower + bin_upper) / 2
  )

poscorad_counts <- poscorad_binned %>%
  count(
    Dataset_display,
    bin_width,
    bin_lower,
    bin_upper,
    bin_mid,
    name = "n"
  )

poscorad_totals <- poscorad_binned %>%
  count(Dataset_display, name = "n_total")

poscorad_dist <- poscorad_grid %>%
  left_join(
    poscorad_counts,
    by = c(
      "Dataset_display",
      "bin_width",
      "bin_lower",
      "bin_upper",
      "bin_mid"
    )
  ) %>%
  left_join(poscorad_totals, by = "Dataset_display") %>%
  mutate(
    n = replace_na(n, 0L),
    prop = n / n_total
  )

# -------------------------------------------------------------------------
# Save a summary table
# -------------------------------------------------------------------------

item_summary <- items_long %>%
  group_by(Dataset_display, Item_display) %>%
  summarise(
    n_observations = n(),
    n_patients = n_distinct(Patient),
    mean = mean(Score, na.rm = TRUE),
    sd = sd(Score, na.rm = TRUE),
    median = median(Score, na.rm = TRUE),
    q1 = quantile(Score, 0.25, na.rm = TRUE),
    q3 = quantile(Score, 0.75, na.rm = TRUE),
    min = min(Score, na.rm = TRUE),
    max = max(Score, na.rm = TRUE),
    pct_zero = mean(Score == 0, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  mutate(Outcome = as.character(Item_display))

poscorad_summary <- poscorad_long %>%
  group_by(Dataset_display) %>%
  summarise(
    n_observations = n(),
    n_patients = n_distinct(Patient),
    mean = mean(Score, na.rm = TRUE),
    sd = sd(Score, na.rm = TRUE),
    median = median(Score, na.rm = TRUE),
    q1 = quantile(Score, 0.25, na.rm = TRUE),
    q3 = quantile(Score, 0.75, na.rm = TRUE),
    min = min(Score, na.rm = TRUE),
    max = max(Score, na.rm = TRUE),
    pct_zero = mean(Score == 0, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  mutate(Outcome = "PO-SCORAD total")

summary_table <- bind_rows(item_summary, poscorad_summary) %>%
  select(
    Dataset = Dataset_display,
    Outcome,
    n_observations,
    n_patients,
    mean,
    sd,
    median,
    q1,
    q3,
    min,
    max,
    pct_zero
  ) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

summary_path <- file.path(out_dir, "data_distribution_summary.csv")
readr::write_csv(summary_table, summary_path)

# -------------------------------------------------------------------------
# Plot A: six intensity signs on native 0-3 scale
# -------------------------------------------------------------------------

p_intensity <- ggplot(
  item_dist %>% filter(Item %in% INTENSITY_ITEMS),
  aes(x = Score_plot, y = prop, fill = Dataset_display)
) +
  geom_col(
    width = 0.92,
    colour = "white",
    linewidth = 0.08
  ) +
  facet_grid(
    Dataset_display ~ Item_display,
    scales = "free_x"
  ) +
  scale_fill_manual(values = DATASET_COLS, guide = "none") +
  scale_x_continuous(
    breaks = full_scale_breaks,
    labels = full_scale_labels,
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.06))
  ) +
  labs(
    title = "Intensity sign distributions",
    x = "Observed score",
    y = "Proportion of observations"
  ) +
  theme_distribution(base_size = 13)

# -------------------------------------------------------------------------
# Plot B: extent, itch, and sleep on binned native scales
# -------------------------------------------------------------------------

plot_binned_wide_item <- function(item_name) {
  plot_df <- wide_item_dist %>%
    filter(Item == item_name)

  item_xmax <- max(plot_df$xmax, na.rm = TRUE)

  ggplot(plot_df, aes(fill = Dataset_display)) +
    geom_rect(
      aes(
        xmin = bin_lower,
        xmax = bin_upper,
        ymin = 0,
        ymax = prop
      ),
      colour = "white",
      linewidth = 0.08
    ) +
    facet_grid(Dataset_display ~ ., scales = "free_y") +
    scale_fill_manual(values = DATASET_COLS, guide = "none") +
    scale_x_continuous(
      breaks = full_scale_breaks,
      labels = full_scale_labels,
      limits = c(0, item_xmax),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_y_continuous(
      labels = percent_format(accuracy = 1),
      expand = expansion(mult = c(0, 0.12))
    ) +
    labs(
      title = ITEM_LABELS[[item_name]],
      x = "Observed score",
      y = "Proportion"
    ) +
    theme_distribution(base_size = 13)
}

p_wide_items <- (
  plot_binned_wide_item("extent") |
    plot_binned_wide_item("itching") |
    plot_binned_wide_item("sleep")
) +
  plot_layout(widths = c(1.35, 1, 1))

# -------------------------------------------------------------------------
# Plot C: total PO-SCORAD distribution
# -------------------------------------------------------------------------

p_poscorad <- ggplot(poscorad_dist, aes(fill = Dataset_display)) +
  geom_rect(
    aes(
      xmin = bin_lower,
      xmax = bin_upper,
      ymin = 0,
      ymax = prop
    ),
    colour = "white",
    linewidth = 0.08
  ) +
  facet_wrap(~Dataset_display, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = DATASET_COLS, guide = "none") +
  scale_x_continuous(
    breaks = seq(0, 100, by = 20),
    labels = full_scale_labels,
    limits = c(0, 103),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.12))
  ) +
  labs(
    title = "Observed total PO-SCORAD distribution",
    x = "Observed PO-SCORAD score",
    y = "Proportion of observations"
  ) +
  theme_distribution(base_size = 13)

# -------------------------------------------------------------------------
# Combine and save only the final figure
# -------------------------------------------------------------------------

data_distribution_plot <- p_intensity / p_wide_items / p_poscorad +
  plot_layout(heights = c(1.05, 1.25, 1.1))

save_plot(
  data_distribution_plot,
  out_dir = out_dir,
  filename_stem = "data_distribution",
  width = 18,
  height = 13,
  formats = c("png", "pdf", "svg")
)

cat("\n[OK] Saved figure files to:\n")
cat("  ", file.path(out_dir, "data_distribution.png"), "\n", sep = "")
cat("  ", file.path(out_dir, "data_distribution.pdf"), "\n", sep = "")
cat("  ", file.path(out_dir, "data_distribution.svg"), "\n", sep = "")

cat("[OK] Saved summary table:\n")
cat("  ", summary_path, "\n", sep = "")
