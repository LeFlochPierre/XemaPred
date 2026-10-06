# -----------------------------------------------------------------------
# Save ggplot/patchwork objects using explicit graphics devices
# -----------------------------------------------------------------------

save_plot <- function(
    plot,
    out_dir,
    filename_stem,
    width,
    height,
    dpi = 320,
    formats = c("png"),
    bg = "white",
    verbose = TRUE,
    family = "sans"
) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # More stable on headless Linux / HPC than the default bitmap device.
  if (capabilities("cairo")) {
    options(bitmapType = "cairo")
  }

  # Force a simple font family. This helps avoid random grid/font crashes.
  plot <- tryCatch(
    {
      if (inherits(plot, "patchwork")) {
        plot & ggplot2::theme(text = ggplot2::element_text(family = family))
      } else if (inherits(plot, "ggplot")) {
        plot + ggplot2::theme(text = ggplot2::element_text(family = family))
      } else {
        plot
      }
    },
    error = function(e) plot
  )

  out_files <- character(0)

  for (ext in formats) {
    out_file <- file.path(out_dir, paste0(filename_stem, ".", ext))

    if (verbose) {
      cat("[SAVE] ", out_file, "\n", sep = "")
    }

    if (ext == "png") {
      if (!capabilities("cairo")) {
        stop("Cairo support is not available. PNG saving may be unstable in this environment.")
      }

      grDevices::png(
        filename = out_file,
        width = width,
        height = height,
        units = "in",
        res = dpi,
        type = "cairo-png",
        bg = bg
      )

      print(plot)
      grDevices::dev.off()

    } else if (ext == "pdf") {
      if (capabilities("cairo")) {
        grDevices::cairo_pdf(
          filename = out_file,
          width = width,
          height = height,
          family = family,
          bg = bg
        )
      } else {
        grDevices::pdf(
          file = out_file,
          width = width,
          height = height,
          family = family,
          bg = bg,
          useDingbats = FALSE
        )
      }

      print(plot)
      grDevices::dev.off()

    } else if (ext == "svg") {
      grDevices::svg(
        filename = out_file,
        width = width,
        height = height,
        family = family,
        bg = bg,
        onefile = FALSE
      )

      print(plot)
      grDevices::dev.off()

    } else {
      stop("Unknown plot format: ", ext)
    }

    out_files <- c(out_files, out_file)

    invisible(gc())
  }

  invisible(out_files)
}