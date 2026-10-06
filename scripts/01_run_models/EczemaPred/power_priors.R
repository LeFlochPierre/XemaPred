# -----------------------------------------------------------------------
# Power-prior extraction for EczemaPred item-level models
#
# Purpose:
# - Extract posterior mean and standard deviation summaries from the final
#   EczemaPred fit.
# - Store model-specific summaries for BinMC, OrderedRW, and BinRW.
# -----------------------------------------------------------------------

mean_sd <- function(x) {
  c(mean(x), stats::sd(x))
}

extract_power_prior <- function(fit, mdl_name) {
  post <- rstan::extract(fit)

  cat("[DEBUG] Posterior variables extracted from Stan fit:\n")
  print(names(post))

  if (mdl_name == "BinMC") {
    return(extract_binmc_power_prior(post))
  }

  if (mdl_name == "OrderedRW") {
    return(extract_orderedrw_power_prior(post))
  }

  if (mdl_name == "BinRW") {
    return(extract_binrw_power_prior(post))
  }

  warning("[WARNING] No power-prior extraction rule for model: ", mdl_name)
  list()
}

extract_binmc_power_prior <- function(post) {
  power_prior <- list()

  if (!is.null(post$mu_logit_p10)) {
    power_prior$mu_logit_p10 <- mean_sd(post$mu_logit_p10)
  } else {
    warning("[WARNING] post$mu_logit_p10 is NULL.")
  }

  if (!is.null(post$sigma_logit_p10)) {
    power_prior$sigma_logit_p10 <- mean_sd(post$sigma_logit_p10)
  } else {
    warning("[WARNING] post$sigma_logit_p10 is NULL.")
  }

  if (!is.null(post$sigma)) {
    power_prior$sigma <- mean_sd(post$sigma)
  } else {
    warning("[WARNING] post$sigma is NULL.")
  }

  power_prior
}

extract_orderedrw_power_prior <- function(post) {
  power_prior <- list()

  if (!is.null(post$sigma_lat)) {
    power_prior$sigma_lat <- mean_sd(post$sigma_lat)
  } else {
    warning("[WARNING] post$sigma_lat is NULL.")
  }

  if (!is.null(post$sigma_meas)) {
    power_prior$sigma_meas <- mean_sd(post$sigma_meas)
  } else {
    warning("[WARNING] post$sigma_meas is NULL.")
  }

  if (!is.null(post$ct) && length(dim(post$ct)) == 2) {
    power_prior$ct <- apply(post$ct, 2, mean_sd)
  } else {
    warning("[WARNING] ct parameter missing or malformed in posterior.")
  }

  power_prior
}

extract_binrw_power_prior <- function(post) {
  power_prior <- list()

  if (!is.null(post$sigma)) {
    power_prior$sigma <- mean_sd(post$sigma)
  } else {
    warning("[WARNING] post$sigma is NULL.")
  }

  if (!is.null(post$mu_logit_y0)) {
    power_prior$mu_logit_y0 <- mean_sd(post$mu_logit_y0)
  } else {
    warning("[WARNING] post$mu_logit_y0 is NULL.")
  }

  if (!is.null(post$sigma_logit_y0)) {
    power_prior$sigma_logit_y0 <- mean_sd(post$sigma_logit_y0)
  } else {
    warning("[WARNING] post$sigma_logit_y0 is NULL.")
  }

  power_prior
}
