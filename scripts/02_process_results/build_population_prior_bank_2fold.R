#!/usr/bin/env Rscript

# ======================================================================
# build_population_prior_bank_2fold.R
#
# Build fold-specific population-XemaPred prior banks for the held-out
# within-cohort experiment.
#
# IMPORTANT:
#   This script ONLY reads population XemaPred runs that were fitted on
#   one patient fold at a time:
#
#     PriorSourceFold1 -> prior bank for source fold 1
#     PriorSourceFold2 -> prior bank for source fold 2
#
#   Downstream use:
#     target patient fold 1 <- source prior fold 2
#     target patient fold 2 <- source prior fold 1
#
# Output:
#   priors/population_xemapred_2fold/<dataset>/fold_<k>/<item>.rds
#
# The selected posterior is the latest population posterior at or before
# the configured ~80% training day.
# ======================================================================

rm(list = ls())

# ======================================================================
# CONFIG
# ======================================================================

ROOT <- "results"
HORIZON <- "H4"

N_FOLDS <- 2L
SOURCE_FOLDS <- 1:2

# Keep consistent with the current H4 snapshot analysis.
TARGET_TRAINING_DAY <- c(
  Derexyl = 77,
  PFDC = 65
)

DATASETS <- c(
  "PFDC",
  "Derexyl"
)

ITEM_MODEL <- c(
  extent     = "BinMC",
  dryness    = "OrderedRW",
  redness    = "OrderedRW",
  swelling   = "OrderedRW",
  oozing     = "OrderedRW",
  thickening = "OrderedRW",
  scratching = "OrderedRW",
  itching    = "BinRW",
  sleep      = "BinRW"
)

ITEM_FAMILY <- c(
  extent     = "extent",
  dryness    = "intensity",
  redness    = "intensity",
  swelling   = "intensity",
  oozing     = "intensity",
  thickening = "intensity",
  scratching = "intensity",
  itching    = "subjective",
  sleep      = "subjective"
)

OUT_ROOT <- file.path(
  "priors",
  "population_xemapred_2fold"
)

`%||%` <- function(a, b) {
  if (!is.null(a)) a else b
}

# ======================================================================
# WEIGHT HELPERS
# ======================================================================

normalize_weights <- function(w) {
  w <- as.numeric(w)
  w[!is.finite(w) | w < 0] <- 0

  s <- sum(w)

  if (!is.finite(s) || s <= 0) {
    return(rep(1 / length(w), length(w)))
  }

  w / s
}

weighted_moments <- function(theta, w) {
  theta <- as.matrix(theta)
  w <- normalize_weights(w)

  if (nrow(theta) != length(w)) {
    stop("[ERROR] theta/weight dimension mismatch.")
  }

  mu <- colSums(
    sweep(theta, 1, w, "*")
  )

  z <- sweep(
    theta,
    2,
    mu,
    "-"
  )

  cov_mat <- crossprod(
    z,
    z * w
  )

  list(
    mean = mu,
    cov = cov_mat,
    sd = sqrt(
      pmax(
        diag(cov_mat),
        0
      )
    )
  )
}

# ======================================================================
# FIND THE CORRECT FOLD-SPECIFIC POPULATION POSTERIOR DIRECTORY
# ======================================================================

posterior_dir <- function(
    dataset,
    item,
    model,
    source_fold
) {

  # Accept both possible suffix renderings:
  #   OrderedRW-XemaPred-PriorSourceFold1
  #   OrderedRW-XemaPredPriorSourceFold1
  #
  # The recursive fallback below still requires the exact source-fold tag,
  # H4, item and model, so it cannot silently select a full-cohort run.
  expected_run_names <- c(
    paste0(
      model,
      "-XemaPred-PriorSourceFold",
      source_fold
    ),
    paste0(
      model,
      "-XemaPredPriorSourceFold",
      source_fold
    )
  )

  candidates <- unlist(
    lapply(
      expected_run_names,
      function(run_name) {
        c(
          file.path(
            ROOT,
            dataset,
            item,
            run_name,
            HORIZON,
            "posterior"
          ),
          file.path(
            ROOT,
            dataset,
            item,
            model,
            run_name,
            HORIZON,
            "posterior"
          )
        )
      }
    ),
    use.names = FALSE
  )

  found <- candidates[
    dir.exists(candidates)
  ]

  # Fallback discovery. Still require:
  #   - requested source-fold tag
  #   - requested model
  #   - requested horizon
  if (!length(found)) {
    item_root <- file.path(
      ROOT,
      dataset,
      item
    )

    if (dir.exists(item_root)) {
      all_dirs <- list.dirs(
        item_root,
        recursive = TRUE,
        full.names = TRUE
      )

      posterior_dirs <- all_dirs[
        basename(all_dirs) == "posterior"
      ]

      fold_tag <- paste0(
        "PriorSourceFold",
        source_fold
      )

      found <- posterior_dirs[
        grepl(fold_tag, posterior_dirs, fixed = TRUE) &
        grepl(model, posterior_dirs, fixed = TRUE) &
        grepl(
          paste0("/", HORIZON, "/"),
          posterior_dirs,
          fixed = TRUE
        )
      ]
    }
  }

  found <- unique(found)

  if (!length(found)) {
    stop(
      "[ERROR] Fold-specific posterior directory not found.\n",
      "  Dataset: ", dataset, "\n",
      "  Item: ", item, "\n",
      "  Source fold: ", source_fold, "\n",
      "  Required fold tag: PriorSourceFold", source_fold
    )
  }

  if (length(found) > 1L) {
    stop(
      "[ERROR] More than one posterior directory matched:\n  ",
      paste(found, collapse = "\n  ")
    )
  }

  # Leakage guard: path must explicitly identify the requested source fold.
  expected_tag <- paste0(
    "PriorSourceFold",
    source_fold
  )

  if (!grepl(expected_tag, found[[1]], fixed = TRUE)) {
    stop(
      "[ERROR] Fold mismatch in posterior path. Expected ",
      expected_tag,
      " but found:\n",
      found[[1]]
    )
  }

  found[[1]]
}

# ======================================================================
# SELECT ~80% POSTERIOR SNAPSHOT
# ======================================================================

select_target_day_posterior <- function(
    dir,
    target_day
) {

  files <- list.files(
    dir,
    pattern = "^posterior-[0-9]+\\.rds$",
    full.names = TRUE
  )

  if (!length(files)) {
    stop(
      "[ERROR] No posterior files in ",
      dir
    )
  }

  info <- do.call(
    rbind,
    lapply(
      files,
      function(f) {
        x <- readRDS(f)

        data.frame(
          file = f,
          iteration = x$iteration %||% NA_integer_,
          day = x$max_training_day %||% NA_real_,
          patients = x$n_training_patients %||% NA_real_,
          observations = x$n_training_observations %||% NA_real_,
          stringsAsFactors = FALSE
        )
      }
    )
  )

  info <- info[
    is.finite(info$day),
    ,
    drop = FALSE
  ]

  if (!nrow(info)) {
    stop(
      "[ERROR] No snapshots with valid max_training_day in ",
      dir
    )
  }

  # Use latest available training day at or before requested ~80% day.
  before <- info[
    info$day <= target_day,
    ,
    drop = FALSE
  ]

  if (nrow(before)) {
    selected_day <- max(
      before$day,
      na.rm = TRUE
    )

    candidates <- before[
      before$day == selected_day,
      ,
      drop = FALSE
    ]
  } else {
    # Fallback only if there is no earlier snapshot.
    nearest <- which.min(
      abs(info$day - target_day)
    )

    candidates <- info[
      nearest,
      ,
      drop = FALSE
    ]
  }

  # If multiple posterior files have the same training day,
  # use the largest iteration number.
  if (
    nrow(candidates) > 1L &&
    any(is.finite(candidates$iteration))
  ) {
    rank_iteration <- ifelse(
      is.finite(candidates$iteration),
      candidates$iteration,
      -Inf
    )

    selected <- candidates[
      which.max(rank_iteration),
      ,
      drop = FALSE
    ]
  } else {
    selected <- candidates[
      1,
      ,
      drop = FALSE
    ]
  }

  message(
    "[SELECT] target day ",
    target_day,
    " -> iteration ",
    selected$iteration,
    " / training day ",
    selected$day,
    " / n patients ",
    selected$patients
  )

  selected
}

# ======================================================================
# EXTRACT XEMAPRED POSTERIOR
# ======================================================================

extract_xemapred_posterior <- function(obj) {

  if (is.null(obj$theta_u)) {
    stop(
      "[ERROR] XemaPred posterior has no theta_u."
    )
  }

  theta <- as.matrix(
    obj$theta_u
  )

  raw_w <- obj$weights %||%
    obj$weights_raw

  if (is.null(raw_w)) {
    stop(
      "[ERROR] XemaPred posterior has no weights or weights_raw."
    )
  }

  w <- normalize_weights(
    raw_w
  )

  if (nrow(theta) != length(w)) {
    stop(
      "[ERROR] theta_u/weights dimension mismatch."
    )
  }

  if (any(!is.finite(theta))) {
    stop(
      "[ERROR] Non-finite theta_u values."
    )
  }

  list(
    theta = theta,
    weights = w
  )
}

# ======================================================================
# BUILD ONE FOLD-SPECIFIC PRIOR
# ======================================================================

build_prior <- function(
    dataset,
    item,
    source_fold
) {

  model <- unname(
    ITEM_MODEL[[item]]
  )

  family <- unname(
    ITEM_FAMILY[[item]]
  )

  dir <- posterior_dir(
    dataset = dataset,
    item = item,
    model = model,
    source_fold = source_fold
  )

  target_day <- unname(
    TARGET_TRAINING_DAY[[dataset]]
  )

  selected <- select_target_day_posterior(
    dir = dir,
    target_day = target_day
  )

  source_file <- selected$file[[1]]

  obj <- readRDS(
    source_file
  )

  post <- extract_xemapred_posterior(
    obj
  )

  theta <- post$theta
  w <- post$weights

  moments <- weighted_moments(
    theta,
    w
  )

  ess <- 1 / sum(
    w^2
  )

  # --------------------------------------------------------------------
  # Strict source-fold metadata
  # --------------------------------------------------------------------

  prior <- list(
    schema_version = 2L,

    prior_id = paste(
      "population_XemaPred_SMC2",
      dataset,
      paste0("fold", source_fold),
      item,
      paste0("day", obj$max_training_day %||% NA),
      sep = "__"
    ),

    prior_type =
      "heldout_within_cohort_population_posterior",

    prior_power = 1.0,

    source_method =
      "population_XemaPred_SMC2",

    source_dataset =
      dataset,

    n_patient_folds =
      N_FOLDS,

    source_fold =
      as.integer(source_fold),

    # Explicitly record which target fold this prior is allowed to inform.
    allowed_target_fold =
      as.integer(
        if (source_fold == 1L) 2L else 1L
      ),

    snapshot_mode =
      "target_day",

    target_training_day =
      target_day,

    score =
      item,

    model =
      model,

    model_family =
      family,

    horizon =
      HORIZON,

    source_iteration =
      obj$iteration %||%
      NA_integer_,

    source_training_day =
      obj$max_training_day %||%
      NA_real_,

    source_training_fraction =
      obj$training_fraction %||%
      NA_real_,

    source_n_training_observations =
      obj$n_training_observations %||%
      NA_real_,

    source_n_training_patients =
      obj$n_training_patients %||%
      NA_real_,

    source_file =
      normalizePath(
        source_file,
        winslash = "/",
        mustWork = FALSE
      ),

    # Authoritative prior representation.
    n_theta =
      nrow(theta),

    parameter_names =
      colnames(theta),

    theta_u =
      theta,

    weights =
      w,

    ess =
      ess,

    ess_pct =
      100 * ess / nrow(theta),

    unique_theta =
      nrow(
        unique(
          round(
            theta,
            10
          )
        )
      ),

    raw_weighted_mean =
      moments$mean,

    raw_weighted_cov =
      moments$cov,

    raw_weighted_sd =
      moments$sd,

    raw_summary =
      data.frame(
        parameter = colnames(theta),
        mean = as.numeric(moments$mean),
        sd = as.numeric(moments$sd),
        stringsAsFactors = FALSE
      )
  )

  # --------------------------------------------------------------------
  # Save:
  #
  # priors/population_xemapred_2fold/<dataset>/fold_<source>/<item>.rds
  # --------------------------------------------------------------------

  out_dir <- file.path(
    OUT_ROOT,
    dataset,
    paste0(
      "fold_",
      source_fold
    )
  )

  dir.create(
    out_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  out_file <- file.path(
    out_dir,
    paste0(
      item,
      ".rds"
    )
  )

  saveRDS(
    prior,
    out_file
  )

  message(
    "[OK] ",
    dataset,
    " / fold ",
    source_fold,
    " / ",
    item,
    " | day=",
    prior$source_training_day,
    " | patients=",
    prior$source_n_training_patients,
    " | theta=",
    nrow(theta),
    "x",
    ncol(theta),
    " | ESS=",
    round(
      ess,
      1
    ),
    " | allowed target fold=",
    prior$allowed_target_fold
  )

  data.frame(
    dataset = dataset,
    source_fold = source_fold,
    allowed_target_fold = prior$allowed_target_fold,
    item = item,
    model = model,
    family = family,
    target_training_day = target_day,
    iteration = prior$source_iteration,
    training_day = prior$source_training_day,
    n_training_patients = prior$source_n_training_patients,
    n_training_observations = prior$source_n_training_observations,
    n_theta = nrow(theta),
    n_parameters = ncol(theta),
    parameter_names = paste(
      colnames(theta),
      collapse = ";"
    ),
    ess = ess,
    source_file = prior$source_file,
    prior_file = normalizePath(
      out_file,
      winslash = "/",
      mustWork = FALSE
    ),
    status = "OK",
    stringsAsFactors = FALSE
  )
}

# ======================================================================
# RUN
# ======================================================================

dir.create(
  OUT_ROOT,
  recursive = TRUE,
  showWarnings = FALSE
)

manifest <- list()

for (dataset in DATASETS) {

  for (source_fold in SOURCE_FOLDS) {

    for (item in names(ITEM_MODEL)) {

      key <- paste(
        dataset,
        paste0("fold", source_fold),
        item,
        sep = "__"
      )

      manifest[[key]] <- tryCatch(
        build_prior(
          dataset = dataset,
          item = item,
          source_fold = source_fold
        ),
        error = function(e) {

          message(
            "[ERROR] ",
            dataset,
            " / fold ",
            source_fold,
            " / ",
            item,
            ": ",
            conditionMessage(e)
          )

          data.frame(
            dataset = dataset,
            source_fold = source_fold,
            allowed_target_fold =
              if (source_fold == 1L) 2L else 1L,
            item = item,
            model = unname(
              ITEM_MODEL[[item]]
            ),
            family = unname(
              ITEM_FAMILY[[item]]
            ),
            target_training_day =
              TARGET_TRAINING_DAY[[dataset]],
            iteration = NA_integer_,
            training_day = NA_real_,
            n_training_patients = NA_real_,
            n_training_observations = NA_real_,
            n_theta = NA_integer_,
            n_parameters = NA_integer_,
            parameter_names = NA_character_,
            ess = NA_real_,
            source_file = NA_character_,
            prior_file = NA_character_,
            status = conditionMessage(e),
            stringsAsFactors = FALSE
          )
        }
      )
    }
  }
}

manifest <- do.call(
  rbind,
  manifest
)

rownames(manifest) <- NULL

manifest_file <- file.path(
  OUT_ROOT,
  "population_prior_2fold_manifest.csv"
)

write.csv(
  manifest,
  manifest_file,
  row.names = FALSE
)

# ======================================================================
# SUMMARY
# ======================================================================

cat(
  "\n============================================================\n"
)

cat(
  "2-fold population XemaPred prior bank complete\n"
)

cat(
  "============================================================\n"
)

cat(
  "Source root: ",
  ROOT,
  "\n",
  sep = ""
)

cat(
  "Horizon:     ",
  HORIZON,
  "\n",
  sep = ""
)

cat(
  "Snapshot:    target-day (~80%)\n"
)

cat(
  "Output:      ",
  OUT_ROOT,
  "\n",
  sep = ""
)

cat(
  "Successful:  ",
  sum(
    manifest$status == "OK"
  ),
  " / ",
  nrow(manifest),
  "\n",
  sep = ""
)

cat(
  "Manifest:    ",
  manifest_file,
  "\n",
  sep = ""
)

cat(
  "============================================================\n"
)

if (
  any(
    manifest$status != "OK"
  )
) {

  cat(
    "\nFailures:\n"
  )

  print(
    manifest[
      manifest$status != "OK",
      c(
        "dataset",
        "source_fold",
        "item",
        "status"
      ),
      drop = FALSE
    ],
    row.names = FALSE
  )
}

# ======================================================================
# FINAL LEAKAGE-AUDIT SUMMARY
# ======================================================================

cat(
  "\n[INFO] Leakage audit mapping:\n"
)

cat(
  "  Source fold 1 priors -> allowed target fold 2 only\n"
)

cat(
  "  Source fold 2 priors -> allowed target fold 1 only\n"
)

cat(
  "\n[INFO] Prior files:\n"
)

ok_manifest <- manifest[
  manifest$status == "OK",
  ,
  drop = FALSE
]

if (nrow(ok_manifest)) {
  print(
    ok_manifest[
      ,
      c(
        "dataset",
        "source_fold",
        "allowed_target_fold",
        "item",
        "training_day",
        "n_training_patients",
        "prior_file"
      ),
      drop = FALSE
    ],
    row.names = FALSE
  )
}
