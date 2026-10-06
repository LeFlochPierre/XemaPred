# -----------------------------------------------------------------------
# Build population-level XemaPred observations
#
# Purpose:
# - Build the new observation stream for one forward-chaining iteration.
# - Include only observations that have not yet been assimilated by the
#   cached SMC2 state.
#
# Output:
# - intensity:  Patient, TimeRel, yc
# - subjective: Patient, TimeRel, y
# - extent:     Patient, TimeRel, y
# -----------------------------------------------------------------------

get_xemapred_patient_last_time <- function(state, pid, run_info = NULL) {
  pid_chr <- as.character(pid)

  t_last <- min(vapply(
    seq_len(state$J),
    function(j) {
      tl <- state$inner[[j]][[pid_chr]]$time_last

      if (is.null(tl) || is.na(tl)) {
        -1L
      } else {
        as.integer(tl)
      }
    },
    integer(1)
  ))

  t_last
}

build_xemapred_obs_for_patient <- function(pid, it, state, data_obj, run_info) {
  split <- split_xemapred_patient_iteration(
    data_obj = data_obj,
    pid = pid,
    it = it
  )

  train <- split$train

  if (nrow(train) == 0) {
    return(NULL)
  }

  t_last <- get_xemapred_patient_last_time(
    state = state,
    pid = pid,
    run_info = run_info
  )

  keep <- which(as.integer(train$TimeRel) > t_last)

  if (!length(keep)) {
    return(NULL)
  }

  family <- run_info$model_family

  if (family == "intensity") {
    train_vals <- as.integer(round(train$Score / run_info$reso))
    train_vals <- pmin(run_info$M_max, pmax(0L, train_vals))
    yc_train <- train_vals + 1L

    return(
      tibble::tibble(
        Patient = pid,
        TimeRel = as.integer(train$TimeRel[keep]),
        yc = as.integer(yc_train[keep])
      )
    )
  }

  if (family == "subjective") {
    M_binom <- as.integer(run_info$M_max)

    y_int <- as.integer(round(train$Score / run_info$reso))
    y_int <- pmin(M_binom, pmax(0L, y_int))

    return(
      tibble::tibble(
        Patient = pid,
        TimeRel = as.integer(train$TimeRel[keep]),
        y = as.integer(y_int[keep])
      )
    )
  }

  if (family == "extent") {
    y_int <- as.integer(round(train$Score))
    y_int <- pmin(100L, pmax(0L, y_int))

    return(
      tibble::tibble(
        Patient = pid,
        TimeRel = as.integer(train$TimeRel[keep]),
        y = as.integer(y_int[keep])
      )
    )
  }

  stop("[ERROR] Unknown XemaPred model family: ", family)
}

build_xemapred_cohort_obs <- function(it, state, data_obj, run_info) {
  obs_list <- vector("list", length(data_obj$patients))
  names(obs_list) <- as.character(data_obj$patients)

  for (pid in data_obj$patients) {
    obs_list[[as.character(pid)]] <- build_xemapred_obs_for_patient(
      pid = pid,
      it = it,
      state = state,
      data_obj = data_obj,
      run_info = run_info
    )
  }

  obs_df <- dplyr::bind_rows(obs_list)

  cat(
    glue::glue(
      "[INFO] Iter {it}: new obs rows = {nrow(obs_df)}",
      if (nrow(obs_df) > 0) {
        glue::glue(", patients = {dplyr::n_distinct(obs_df$Patient)}")
      } else {
        ""
      },
      "\n"
    )
  )

  obs_df
}