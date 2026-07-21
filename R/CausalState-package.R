#' CausalState: SDR and iTMLE for State-Aware Longitudinal Modified Treatment Policies
#'
#' Two sequentially doubly robust estimators for longitudinal modified treatment
#' policies (MTPs) in settings with absorbing state transitions - such as ICU
#' stay, ward admission, or ED episodes - where treatment is structurally
#' inapplicable after transition to an absorbing state (e.g. discharge or
#' death).
#'
#' @section Sequential Doubly Robust estimator - [sdr()]:
#' Implements the LMTP-SDR estimator of Diaz et al. (2021), which extends the
#' SDR construction of Luedtke et al. (2017) to general modified treatment
#' policies via density-ratio weighting.  Starting from terminal pseudo-outcomes
#' and working backward in time, the estimator propagates an EIF-corrected
#' pseudo-outcome through a sequence of Q-regressions.  The final estimate is
#' the mean of these pseudo-outcomes at the first time point.  The estimator is
#' sequentially doubly robust (2^K-robust): consistent whenever, at each time
#' point, either the treatment model g_t or the outcome model Q_t is correctly
#' specified.  Under the natural-course policy the density ratios equal one and
#' the recursion collapses algebraically to `mean(Y)` regardless of model
#' quality - a natural-course SDR run therefore cannot serve as a model
#' calibration check.
#'
#' @section Infinite-dimensional TMLE - [itmle()]:
#' Implements the cross-validated infinite-dimensional TMLE of Luedtke et al.
#' (2017, Section 5 and Appendix 12), adapted to general MTPs via density-ratio
#' weighting following Diaz et al. (2021).  Like SDR, the estimator fits a
#' backward Q-regression; unlike SDR, it then applies an infinite-dimensional
#' fluctuation (targeting step) that solves the efficient influence equation
#' semiparametrically without restricting the fluctuation to a parametric
#' submodel.  The cross-validated variant is used throughout to prevent the
#' targeting step from overfitting the Q estimates.  iTMLE is also sequentially
#' doubly robust.  Unlike SDR, running under the natural-course policy does not
#' collapse to `mean(Y)` and therefore provides a genuine calibration check on
#' the Q models after fluctuation.
#'
#' @section Density ratio estimation - [density_ratio()]:
#' Estimates the density ratios needed by both estimators via a flexible
#' classification-based approach.  An alternative metalearner based on
#' Wu and Benkeser combines base learners by minimising a log density-ratio
#' loss rather than the default NNLS.  Density ratios can be computed once
#' and reused across both estimators and across multiple time horizons.
#'
#' @section State transition model and probabilistic mixing:
#' At each time point t, a subject in state (in_state = 1, alive = 1) faces
#' three mutually exclusive outcomes: remain in state, exit alive (discharge),
#' or die.  The backward Q-recursion assembles the expected outcome as a
#' probability-weighted mixture over these three branches:
#'
#' ```
#' Q_t = p_rem * Q_rem + (1 - p_rem) * [p_dex * Q_death + (1 - p_dex) * Q_dc]
#' ```
#'
#' Two regression components estimate the transition probabilities.
#' `g_remain` (fitted via `sl_remain`) estimates P(in_state = 1 at t+1 |
#' in_state = 1 at t, history) and is trained on all subjects currently in
#' state.  `g_death_exit` (fitted via `sl_death`) estimates P(death | exited,
#' history) and is trained only on the exit subset -- those subjects who leave
#' the state at each time point.
#'
#' The outcome components `Q_death` and `Q_dc` are typically overridden to
#' fixed absorbing-state constants via [absorb_rule()] (e.g. Y = 0 for death)
#' and require no separate regression.  `Q_rem` uses `sl_recursive` at
#' intermediate time points and `sl_y` at the terminal time point.
#'
#' **Event counts and the death regression.**  Because `g_death_exit` trains
#' on the exit subset, its effective sample size can be much smaller than the
#' full cohort -- and within that subset, deaths may be rare.  When per-time
#' death counts are low the SuperLearner fit can degenerate (near-empty
#' training sets, all-zero predictions, solver failures).
#'
#' Practical guidance:
#' \itemize{
#'   \item Inspect `diagnostics$branch_cal` from any estimator: the `g_death`
#'     rows show per-fold, per-time calibration and effective training sizes.
#'   \item As a rough floor, fewer than approximately 20 death events at a
#'     given time point makes per-time fitting unreliable; below 10, estimates
#'     at that step should be treated with caution.
#'   \item Set `pool_g_death = TRUE` to fit a single model across all time
#'     points (with time as a covariate), borrowing strength across follow-up.
#'     Use when per-time counts are sparse but the death hazard does not change
#'     markedly over time.
#'   \item When [absorb_rule()] fixes `Q_death` to a constant (the common
#'     case), a misspecified `g_death_exit` has limited impact: the death
#'     branch enters as P(death | exit) * constant, and any misspecification
#'     shifts the discharge branch by the complementary probability, which the
#'     DR correction via density ratios partially compensates.
#' }
#'
#' @section State machinery:
#' [absorb_rule()] defines outcome overrides at absorbing states (e.g. forcing
#' Y = 0 for subjects who die before the end of follow-up).
#'
#' @section Custom SuperLearner wrappers:
#' SDR-targeted learners ([SL.tgt.glm()], [SL.tgt.glmnet()],
#' [SL.tgt.xgboost()], etc.) for the Q-regression step, and iTMLE
#' fluctuation learners ([SL.tmle_glm()], [SL.tmle_glmnet_ridge()], etc.)
#' that carry the logit offset required for cross-validated TMLE as a data
#' column.
#'
#' @references
#' Diaz I, Williams N, Hoffman KL, Schenck EJ (2021). Nonparametric Causal
#' Effects Based on Longitudinal Modified Treatment Policies. *JASA*
#' 118(542):846-857.
#'
#' Luedtke AR, Sofrygin O, van der Laan MJ, Carone M (2017). Sequential
#' Double Robustness in Right-Censored Longitudinal Models. arXiv:1705.02459.
#'
#' @keywords package
"_PACKAGE"

utils::globalVariables(c(
  ".", "..a_names", "..cols", "..cols_hist", "..cols_t", "..L_now_cols",
  ".__cl", ".__id", ".__time", ".__y",
  "Rt_cum", "Rt_t", "S",
  "alive_last", "changed", "component", "delta", "fold",
  "global_fold", "ic", "id", "id_std", "in_state_last",
  "n", "natural", "p_raw", "r_raw", "row_type",
  "shifted", "time_std", "value", "weight", "y", "y1"
))
