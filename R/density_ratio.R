
#' Estimate density ratios for a modified treatment policy (deprecated)
#'
#' @description
#' `density_ratio()` is a deprecated alias for [fit_g()]. All arguments are
#' forwarded unchanged. Use `fit_g()` in new code.
#'
#' @inheritParams fit_g
#' @param ... Passed to [fit_g()].
#' @return A `"g_fit"` object; see [fit_g()].
#' @seealso [fit_g()]
#' @export
density_ratio <- function(
    df, a_names, tmax,
    baseline = NULL, tv_names = character(0),
    no_lag_vars = character(0),
    policy_names = character(0),
    sl_g, k = 2, inner_v = 10L,
    cluster = NULL, v = 5L,
    seed = 1L,
    policy_spec_fun,
    id = "id", time = "trial_time",
    bounds = 1e-5,
    dr_sl = FALSE,
    drop_small_cluster_splits = TRUE,
    parallel = FALSE,
    t_workers = NULL,
    fold_workers = NULL,
    sl_workers = NULL,
    verbose = TRUE
) {
  .Deprecated("fit_g")
  fit_g(
    df = df, a_names = a_names, tmax = tmax,
    baseline = baseline, tv_names = tv_names,
    no_lag_vars = no_lag_vars, policy_names = policy_names,
    sl_g = sl_g, k = k, inner_v = inner_v,
    cluster = cluster, v = v, seed = seed,
    policy_spec_fun = policy_spec_fun,
    id = id, time = time,
    bounds = bounds,
    method = "classification",
    dr_sl  = dr_sl,
    drop_small_cluster_splits = drop_small_cluster_splits,
    parallel = parallel,
    t_workers = t_workers, fold_workers = fold_workers,
    sl_workers = sl_workers,
    verbose = verbose
  )
}
