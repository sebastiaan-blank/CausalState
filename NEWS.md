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
