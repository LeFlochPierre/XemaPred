#!/usr/bin/env Rscript
# =======================================================================
# smc_diagnostics_report.R
#
# Reads the instrumented per-iteration diagnostics for BOTH variants and
# produces everything needed for Supplementary Methods B.5 / B.6.
#
# Inputs (written by the runners):
#   results/<ds>/<item>/<Model>-XemaPred/H<h>/final/diagnostics.rds
#   results/<ds>/<item>/<Model>-PatientSpecificXemaPred/H<h>/final/diagnostics.rds
#
# Outputs:
#   smc_diag_by_item.csv     one row per dataset x item x variant
#   smc_diag_summary.csv     one row per variant x family  <- the paper table
#   smc_diag_settings.csv    settings actually used        <- B.5 table
#   smc_ess_trace.csv        ESS per iteration             <- for a figure
#
# Usage:
#   Rscript smc_diagnostics_report.R [horizon]
# =======================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(purrr); library(tidyr); library(readr)
})
options(width = 250)

a <- commandArgs(trailingOnly = TRUE)
getopt <- function(flag, default) {
  i <- match(flag, a); if (is.na(i) || i == length(a)) default else a[i + 1]
}
H      <- as.integer(getopt("--horizon", if (length(a) && !grepl("^--", a[1])) a[1] else "4"))
WHICH  <- tolower(getopt("--variant", "both"))   # pop | ps | both

DATASETS <- c("PFDC", "Derexyl")
MODEL <- c(extent = "BinMC", itching = "BinRW", sleep = "BinRW",
           dryness = "OrderedRW", redness = "OrderedRW", swelling = "OrderedRW",
           oozing = "OrderedRW", thickening = "OrderedRW", scratching = "OrderedRW")
FAMILY <- c(extent = "Extent", itching = "Subjective", sleep = "Subjective",
            dryness = "Intensity", redness = "Intensity", swelling = "Intensity",
            oozing = "Intensity", thickening = "Intensity", scratching = "Intensity")
VARIANT <- c(Population = "-XemaPred", `Patient-specific` = "-PatientSpecificXemaPred")
VARIANT <- switch(WHICH,
  pop  = VARIANT["Population"],
  ps   = VARIANT["Patient-specific"],
  both = VARIANT,
  stop("--variant must be pop, ps, or both"))
cat("Variant(s):", paste(names(VARIANT), collapse = ", "), " | horizon: H", H, "\n\n", sep = "")

`%||%` <- function(x, y) if (is.null(x)) y else x
col <- function(d, nm, default = NA) {
  if (nm %in% names(d)) return(d[[nm]])
  if (length(default) == nrow(d)) return(default)
  rep(default, length.out = nrow(d))
}
med    <- function(x) suppressWarnings(median(as.numeric(x), na.rm = TRUE))
mn     <- function(x) { x <- as.numeric(x); x <- x[is.finite(x)]; if (!length(x)) NA_real_ else min(x) }
mx     <- function(x) { x <- as.numeric(x); x <- x[is.finite(x)]; if (!length(x)) NA_real_ else max(x) }
rng    <- function(x) sprintf("%.0f (%.0f-%.0f)", med(x), mn(x), mx(x))

# -----------------------------------------------------------------------
# read every diagnostics file
# -----------------------------------------------------------------------
read_one <- function(ds, item, vlab) {
  fdir <- file.path("results", ds, item, paste0(MODEL[[item]], VARIANT[[vlab]]),
                    paste0("H", H), "final")
  cands <- c("diagnostics.rds", "diagnostics_ALL_PATIENTS.rds")
  f <- file.path(fdir, cands)
  f <- f[file.exists(f)]
  if (!length(f)) return(NULL)

  for (fi in f) {
    d <- try(readRDS(fi), silent = TRUE)
    if (inherits(d, "try-error") || !is.data.frame(d) || !nrow(d)) next
    return(
      as_tibble(d) |>
        mutate(Dataset = ds, Item = item, Family = FAMILY[[item]], Variant = vlab,
               .before = 1)
    )
  }
  NULL
}

raw <- list()
for (ds in DATASETS) for (item in names(MODEL)) for (vlab in names(VARIANT)) {
  d <- read_one(ds, item, vlab)
  if (is.null(d)) { cat(sprintf("[miss] %s / %-11s / %s\n", ds, item, vlab)); next }
  raw[[length(raw) + 1]] <- d
  cat(sprintf("[ok]   %s / %-11s / %-16s  %d rows\n", ds, item, vlab, nrow(d)))
}
if (!length(raw)) stop("No diagnostics.rds found. Check the horizon and that the runs finished.")

# -----------------------------------------------------------------------
# harmonise: the two variants record slightly different columns
# -----------------------------------------------------------------------
harm <- map_dfr(raw, function(d) {
  np <- as.numeric(col(d, "n_theta", col(d, "n_particles")))

  # pre-resampling ESS: ess_before (population) or ess_min (patient-specific)
  ess_pre <- as.numeric(col(d, "ess_before"))
  if (all(is.na(ess_pre))) ess_pre <- as.numeric(col(d, "ess_min"))

  tibble(
    Dataset = d$Dataset, Item = d$Item, Family = d$Family, Variant = d$Variant,
    Patient = as.character(col(d, "Patient")),
    iter    = as.integer(col(d, "iter")),

    n_particles = np,
    Nx          = as.numeric(col(d, "n_x", col(d, "pmmh_nx"))),
    np_sim      = as.numeric(col(d, "np_sim")),
    seed        = as.numeric(col(d, "seed")),

    ess_pre     = ess_pre,
    ess_pre_pct = 100 * ess_pre / np,
    ess_pct     = as.numeric(col(d, "ess_pct")),
    inner_ess   = as.numeric(col(d, "inner_ess_pct")),
    maxW        = as.numeric(col(d, "maxW")),

    unique_pct  = as.numeric(col(d, "unique_pct")),
    accept_cum  = as.numeric(col(d, "accept_cum")),
    attempt_cum = as.numeric(col(d, "attempt_cum")),

    resample    = as.numeric(col(d, "did_resample")),
    resamp_cum  = as.numeric(col(d, "n_resample_cum")),

    pmmh_window = as.numeric(col(d, "pmmh_window")),
    ess_thr     = as.numeric(col(d, "ess_theta_thr", col(d, "ess_frac_thr"))),
    proposal_sd = as.character(col(d, "proposal_sd")),

    update_time = as.numeric(col(d, "update_time", col(d, "compute_time")))
  )
})

write_csv(harm |> select(Dataset, Item, Variant, Patient, iter, ess_pre_pct, ess_pct),
          "smc_ess_trace.csv")

# -----------------------------------------------------------------------
# per dataset x item x variant
# -----------------------------------------------------------------------
by_item <- harm |>
  group_by(Dataset, Item, Family, Variant) |>
  summarise(
    n_iter      = n_distinct(iter),
    n_patients  = n_distinct(Patient[!is.na(Patient)]),
    n_particles = med(n_particles),
    Nx          = med(Nx),
    np_sim      = med(np_sim),
    seed        = med(seed),

    ess_pct_med = med(ess_pct),
    ess_pct_min = mn(ess_pct),
    ess_pre_med = med(ess_pre_pct),
    ess_pre_min = mn(ess_pre_pct),
    ess_pre_lo  = 100 * mean(ess_pre_pct < 25, na.rm = TRUE),
    ess_pre_late = med(ess_pre_pct[iter > 1]),
    inner_ess   = med(inner_ess),
    maxW_max    = mx(maxW),

    unique_med  = med(unique_pct),

    # accept/resample counters are cumulative -> take the final state
    accept_pct  = {
      key <- ifelse(is.na(Patient), "_pop_", Patient)
      idx <- unlist(lapply(split(seq_along(iter), key),
                           function(ii) ii[which.max(iter[ii])]))
      a <- sum(accept_cum[idx],  na.rm = TRUE)
      t <- sum(attempt_cum[idx], na.rm = TRUE)
      if (t > 0) 100 * a / t else NA_real_
    },
    n_resample  = {
      if (any(!is.na(resample))) sum(resample, na.rm = TRUE)
      else {
        key <- ifelse(is.na(Patient), "_pop_", Patient)
        idx <- unlist(lapply(split(seq_along(iter), key),
                             function(ii) ii[which.max(iter[ii])]))
        med(resamp_cum[idx])
      }
    },

    pmmh_window = med(pmmh_window),
    ess_thr     = med(ess_thr),
    proposal_sd = first(na.omit(proposal_sd)),
    update_s    = sum(update_time, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(by_item, "smc_diag_by_item.csv")

# -----------------------------------------------------------------------
# per variant x family  -> the supplementary table
# -----------------------------------------------------------------------
summ <- by_item |>
  group_by(Variant, Family) |>
  summarise(
    n_runs      = n(),
    particles   = paste(unique(na.omit(n_particles)), collapse = "/"),
    `ESS %`             = rng(ess_pre_med),
    `ESS % (iter>1)`    = rng(ess_pre_late),
    `ESS % min`         = sprintf("%.0f", mn(ess_pre_min)),
    `Iters ESS<25% (%)` = sprintf("%.0f", med(ess_pre_lo)),
    `Accept %`          = rng(accept_pct),
    `Distinct theta %`  = rng(unique_med),
    `Resampling events` = sprintf("%.0f", med(n_resample)),
    proposal_sd = first(na.omit(proposal_sd)),
    .groups = "drop"
  ) |>
  arrange(desc(Variant), Family)

write_csv(summ, "smc_diag_summary.csv")

# -----------------------------------------------------------------------
# same table, split by dataset
# -----------------------------------------------------------------------
summ_ds <- by_item |>
  group_by(Variant, Dataset, Family) |>
  summarise(
    n_runs      = n(),
    particles   = paste(unique(na.omit(n_particles)), collapse = "/"),
    n_iter      = med(n_iter),
    `ESS %`             = rng(ess_pre_med),
    `ESS % min`         = sprintf("%.0f", mn(ess_pre_min)),
    `Iters ESS<25% (%)` = sprintf("%.0f", med(ess_pre_lo)),
    `Accept %`          = rng(accept_pct),
    `Distinct theta %`  = rng(unique_med),
    `Resampling events` = sprintf("%.0f", med(n_resample)),
    .groups = "drop"
  ) |>
  arrange(desc(Variant), Dataset, Family)

write_csv(summ_ds, "smc_diag_summary_by_dataset.csv")

# -----------------------------------------------------------------------
# settings actually used -> B.5
# -----------------------------------------------------------------------
settings <- by_item |>
  distinct(Variant, Family, n_particles, Nx, np_sim, ess_thr, pmmh_window,
           seed, proposal_sd) |>
  arrange(desc(Variant), Family)

write_csv(settings, "smc_diag_settings.csv")

# -----------------------------------------------------------------------
cat("\n", strrep("=", 96), "\n",
    "SUPPLEMENTARY B.6 - diagnostics, median (range) across items and datasets\n",
    strrep("=", 96), "\n\n", sep = "")
print(as.data.frame(summ), right = FALSE)

cat("\n", strrep("=", 96), "\n",
    "SAME, SPLIT BY DATASET\n",
    strrep("=", 96), "\n\n", sep = "")
print(as.data.frame(summ_ds), right = FALSE)

cat("\n", strrep("=", 96), "\n",
    "SUPPLEMENTARY B.5 - settings actually used\n",
    strrep("=", 96), "\n\n", sep = "")
print(as.data.frame(settings), right = FALSE)

# -----------------------------------------------------------------------
# consistency checks
# -----------------------------------------------------------------------
cat("\n", strrep("-", 96), "\n", "CHECKS\n", strrep("-", 96), "\n", sep = "")

chk <- function(lbl, ok, detail = "")
  cat(sprintf("  %s  %-46s %s\n", if (isTRUE(ok)) "PASS" else "FAIL", lbl, detail))

s <- unique(na.omit(by_item$seed))
chk("single seed across all runs", length(s) == 1,
    paste(s, collapse = ", "))

ns <- unique(na.omit(by_item$np_sim))
chk("forecast draws identical across runs", length(ns) == 1,
    paste(ns, collapse = ", "))

chk("every run resampled at least once", all(by_item$n_resample > 0, na.rm = TRUE),
    sprintf("min = %s", mn(by_item$n_resample)))

chk("no missing ESS", !any(is.na(by_item$ess_pre_med)))
chk("no missing acceptance", !any(is.na(by_item$accept_pct)))

cat("\nWritten: smc_diag_by_item.csv, smc_diag_summary.csv,",
    "smc_diag_settings.csv, smc_ess_trace.csv\n\n")