# -----------------------------------------------------------------------
# Shared plotting themes
# -----------------------------------------------------------------------

theme_distribution <- function(base_size = 13) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0),
      strip.background = element_rect(fill = "grey95", colour = NA),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.spacing.x = unit(0.75, "lines"),
      panel.spacing.y = unit(0.75, "lines"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      legend.position = "none"
    )
}

theme_patient_base <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "white", colour = NA),
      axis.text = element_text(colour = "grey20"),
      axis.title = element_text(colour = "grey10")
    )
}

theme_patient_facets <- function(base_size = 11) {
  theme_patient_base(base_size = base_size) +
    theme(
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text.y.left = element_text(
        angle = 0,
        face = "bold",
        size = 15.5,
        colour = "grey10",
        hjust = 1,
        margin = margin(t = 0, r = 12, b = 0, l = 0)
      ),
      panel.background = element_rect(fill = "grey99", colour = NA),
      panel.border = element_rect(fill = NA, colour = "grey58", linewidth = 0.35),
      panel.grid.major.x = element_line(colour = "grey86", linewidth = 0.25),
      panel.grid.major.y = element_line(colour = "grey91", linewidth = 0.22)
    )
}

theme_learning_curve <- function(base_size = 14) {
  theme_bw(base_size = base_size) +
    theme(
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text = element_text(
        face = "bold",
        margin = margin(b = 2)
      ),
      axis.text.x.top = element_text(margin = margin(b = 1)),
      axis.title.x.top = element_text(margin = margin(b = 2)),
      axis.ticks.length.x.top = grid::unit(2, "pt"),
      panel.spacing.y = grid::unit(10, "pt"),
      axis.text.y = element_text(margin = margin(r = 10)),
      axis.title.y = element_text(margin = margin(r = 12)),
      legend.position = "top",
      legend.key.width = grid::unit(1.10, "cm"),
      legend.key.height = grid::unit(0.50, "cm"),
      legend.margin = margin(b = 2),
      legend.box.margin = margin(b = 2),
      plot.title = element_text(
        face = "bold",
        hjust = 0,
        margin = margin(b = 4)
      )
    )
}