suppressWarnings(
  suppressPackageStartupMessages({
    library(here)
    library(dplyr)
    library(tidyr)
    library(purrr)
    library(ggplot2)
    library(tibble)
    library(patchwork)
    library(scales)
    library(readr)
    library(glue)
    library(stringr)
  })
)

suppressWarnings(
  source(here::here("scripts", "00_setup", "00_init.R"))
)