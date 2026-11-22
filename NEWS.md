# CausalState 0.11.1

## Bug fixes

* `itmle()`: at each flexible targeting step, zero-weight at-risk training
  patients now receive the fitted fluctuation correction at their history,
  matching the rule already applied to validation patients. Previously these
  patients had `eps_tr = 0` at the current inner step and their propagated
  `Q` could enter the preceding regression stale. Point estimates may change,
  particularly under MTPs whose density ratios hit zero on at-risk patients;
  the fix is a no-op if no density ratio actually hits zero.

## Internal

* `itmle_target_fit()` now fits a single SuperLearner per flexible inner
  targeting step per outer fold (previously `v_target` separate SuperLearners
  whose outer-validation predictions were bagged). The candidate predictions
  feeding the training recursion are out of fold (from SuperLearner's `$Z`),
  combined via `fit$method$computePred()` with ensemble weights estimated
  across all eligible outer-training rows. Zero-weight at-risk training rows
  and outer-validation rows share the same full-refit ensemble prediction.
  This uses the same selection principle as Algorithm 1 of Luedtke et al.
  (2018): the candidate predictions are out of fold; the ensemble weights are
  pooled across eligible rows.
* When clustered CV is requested but fewer than `v_sl_inner_itmle` clusters
  survive filtering, the inner CV is now run with `V = number_of_clusters`
  rather than falling back to row-wise CV. With fewer than two clusters the
  targeting step contributes `eps = 0`.
* `reg_workers` and `sl_workers` are now combined into the single
  SuperLearner's `mc.cores` budget (there is no longer an outer v_target
  `par_lapply` loop to spend `reg_workers` on).

## Breaking changes

* `v_target_itmle` argument of `itmle()` is removed. Passing it now raises
  R's standard "unused argument" error. Existing scripts need only drop the
  argument; the default behaviour is unchanged by the removal.

# CausalState 0.11.0

## Bug fixes

* `itmle()`: the inner fluctuation design at math-time `t` was missing `A_t`,
  violating the Luedtke et al. (2017, Alg. 4) score equation
  `w * (Y - Q_t(H_t, A_t))`. `itmle_make_cols()` now includes `A_u` for
  `u = 1..t` (previously `u = 1..t-1`), matching the `tv_names` loop.
  Point estimates from `itmle()` under non-trivial MTPs will change.

## New features

* Categorical treatment is now declared by factor columns in the input data.
  Ordinal propensity estimation requires an ordered factor and uses its
  declared level order. Numeric binary 0/1 and continuous treatments retain
  their existing coding. No estimator has a separate `treatment_levels`
  argument; factor levels and ordered status are saved as `fit_g()` metadata
  and checked by SDR and iTMLE.
* Factor treatment and treatment history use fixed numeric reference
  indicators in treatment, transition, outcome, and iTMLE targeting models.
  Policies receive factor columns, and policy application preserves their
  labels, levels, and ordering. This permits separate GLM category effects.
* Ordinal density ratios sum natural probabilities over every policy
  preimage, including many-to-one dose caps and minimum-dose rules.
* Preserve exact zero density ratios when SDR and iTMLE assemble weights.
  Zero policy probability was previously treated as invalid and replaced
  with one.
* `sim_cat()` returns a panel with an ordered four-level treatment for use
  with `fit_g(method = "propensity")`.

# CausalState 0.10.2

## Bug fixes

* `density_ratio()`: crashed fold workers (returned as `try-error` by
  `mclapply`) now cause an immediate hard stop with a descriptive message,
  rather than being silently converted to `Rt_t = 1`. A second guard stops
  if any subject has `NA` density-ratio weight after the fold loop.

## Internal

* `R/diagnostics.R` merged into `R/helpers.R`; `diagnostics.R` removed.

# CausalState 0.10.1

## CRAN resubmission fixes

* Title shortened to comply with CRAN length requirement.
* DOI references added to `Description` field in DESCRIPTION.
* Added `\value` tag to `sl_tmle.Rd`.
* Removed dead `cat()`/`SL_WARN_FILE` logging block from `R/sl.R`.
* Removed `inst/sim_study/` (non-standard directory containing local
  simulation scripts with global environment side-effects).
* Replaced non-ASCII em dash characters in `R/helpers.R` with `—`
  escapes.
* Added documentation for exported `sim_bin()`, `sim_cont()`, and
  `sim_multi()` (missing man pages).
* Added `skip_on_cran()` to `tests/testthat/test-accuracy.R`.

# CausalState 0.10.0

## Breaking changes

* New `pool_time = c("spline", "linear", "factor")` argument on `sdr()`,
  `itmle()`, and `qreg()`, defaulting to `"spline"` (restricted cubic
  spline, `K = min(5, tmax)` equally-spaced knots, `df = K - 1`). Applies
  only when `pool_g_death = TRUE` or `pool_q_exit = TRUE`. Point estimates
  from pooled fits will change compared to 0.9.0; pass `pool_time = "linear"`
  to reproduce the previous behaviour.

* `Q_rem` at the terminal block (`tt == tmax`) now uses `sl_recursive`
  (routed through `sl_rec_early` when `tt <= rec_transition`), reverting
  the 0.9.0 behaviour where it used `sl_y`. The family switch (binomial
  when the outcome is binomial at `tmax`) is retained. Terminal-block
  estimates change.

* `sl_rec_simple` renamed to `sl_rec_early` in `sdr()`, `itmle()`, and
  internal helpers. Update any call that used the old argument name.

* `weight_diagnostics()` default `trim` changed from `0.99` to `1` (no
  trimming by default). Pass `trim = 0.99` explicitly for the previous
  output.

* `$diagnostics$recursion_diag` cleanup (`sdr()` and `itmle()`):
    * Renamed internal `pseudo_pre_ar_tr` to `Y_target_ar`; the diagnostic
      table now carries `Y_target_{mean,sd,min,max}` and
      `resid_{sd,q95_abs,max_abs}` columns.
    * Dropped `has_rem`, `has_dex`, `has_qexit`, `has_qrem` (superseded by
      the existing `used_const_*` flags).
    * Dropped `mean_pseudo_pre_ar`, `sd_pseudo_pre_ar`, `Q_pre_diff_*`,
      `delta_mean`, `corr_*`, `exit_contrib_*_mean`, `rem_contrib_*_mean`.
    * iTMLE additionally drops `Q_post_diff_*`, `delta_nat_target_mean`,
      `delta_shf_target_mean`, `delta_nat_vl_target_mean`, and
      `delta_shf_vl_target_mean`.

## New features

* Exported simulation helpers for vignettes, examples, and user
  experimentation:
    * `sim_bin()` -- single binary treatment.
    * `sim_cont()` -- single continuous treatment.
    * `sim_multi(n_binary, n_continuous)` -- arbitrary number of binary
      and/or continuous treatments (columns `A_b1..A_bK, A_c1..A_cM`).

* User-facing diagnostic helpers:
    * `weight_diagnostics()` -- per-time weight summary (instantaneous and
      cumulative-product means, quantiles, Kish ESS) from a density-ratio
      object.
    * `branch_cal_summary()` -- fold-averaged per-`(t, branch)`
      calibration table from an `sdr()` / `itmle()` fit.
    * `contrast()` -- risk difference, risk ratio, and odds ratio between
      two fits (or between a fit and the observed mean), with SEs derived
      from the influence curves.

## Documentation

* New `vignette("diagnostics")` walking through the recommended workflow:
  density-ratio weights and trim selection, natural-course Q calibration
  and SuperLearner tuning, final estimation, then post-estimation sanity
  checks.
* `?sdr` and `?itmle` gain a `Diagnostics` section listing all
  diagnostic slots (`branch_cal`, `recursion_diag`, `target_cal`,
  `target_sl`, `sl_summary`, `ic_df`) with one-sentence purposes and a
  pointer to the vignette.
* `?sdr` and `?itmle` `@examples` reworked to use the new exported
  simulation helpers -- one block each for binary, continuous, and
  mixed multi-treatment scenarios.
* `sl_workers` documentation simplified: no longer names specific
  learners as incompatible, just notes that fork-parallel evaluation is
  not compatible with all learners.

## Bug fixes

* `weight_diagnostics()` cumulative-product weights now computed over all
  `N` subjects (previously restricted, causing under-reporting of the
  mean cumulative weight).

# CausalState 0.9.0

* First public pre-release.
* Main estimators: `sdr()`, `itmle()`, `density_ratio()`, `qreg()`.
* Supports joint multi-treatment policies via `a_names = c("A1", "A2", ...)`.
* Gaussian outcome support with internal `[0, 1]` scaling via `y_bounds`.
* Pooled death model option (`pool_g_death`) and pooled exit-outcome model
  option (`pool_q_exit`) for sparse late time points.
* Custom SuperLearner wrappers for SDR (`SL.tgt.*`) and iTMLE (`SL.tmle_*`) targeting steps.
* Cross-fitted density ratio weights reusable across estimators and time horizons.
* Cluster-robust standard errors via `cluster` argument.
