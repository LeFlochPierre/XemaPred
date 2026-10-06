# -----------------------------------------------------------------------
# EczemaPred item-model iteration wrapper
#
# Purpose:
# - Fit one forward-chaining iteration for BinMC, OrderedRW, or BinRW.
# - Use the shared Stan iteration fitter, because these models follow the
#   standard EczemaModel() / EczemaFit() validation workflow.
#
# Called by:
# - eczemapred_item_models/run_validation.R through shared/parallel.R
# -----------------------------------------------------------------------

fit_eczemapred_item_iteration <- function(it, data_obj, run_info, paths) {
  fit_stan_validation_iteration(it, data_obj, run_info, paths)
}