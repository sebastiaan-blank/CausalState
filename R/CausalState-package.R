#' CausalState: Causal Inference in a Longitudinal Transitioning State Environment
#'
#' Two sequentially doubly robust estimators for longitudinal modified treatment
#' policies (MTPs) in settings with transitioning states - such as ICU
#' stay, ward admission, or ED episodes - where treatment is structurally
#' inapplicable after a state transition (e.g. discharge or death).
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
#' At each time point t, a subject in state (`in_state = 1`, `alive = 1`)
#' faces three mutually exclusive outcomes: remain in state, exit alive
#' (discharge), or die.  Four regression components handle this, corresponding
#' directly to the `sl_` arguments of [sdr()], [itmle()], and [qreg()]:
#'
#' \describe{
#'   \item{`g_remain` (`sl_remain`)}{Probability of remaining in state at
#'     the next time point, given current history.  Trained on all subjects
#'     currently in state.}
#'   \item{`g_death_exit` (`sl_death`)}{Probability of dying, conditional on
#'     having left the state (i.e. among exiters only).  Trained on the exit
#'     subset -- those subjects who leave the state at each time point.}
#'   \item{`Q_rem` (`sl_recursive` at intermediate time points; `sl_y` at
#'     the terminal time point)}{Expected outcome conditional on remaining in
#'     state.  This is the component propagated backward through the
#'     Q-recursion.}
#'   \item{`Q_exit` (`sl_y`)}{Expected outcome conditional on exiting.  In
#'     most applications this is not estimated but instead fixed to a constant
#'     for each exit type (death, discharge) via [absorb_rule()].}
#' }
#'
#' **Probabilistic mixture.**  The backward Q-recursion assembles the expected
#' outcome at time t as a probability-weighted mixture over the three
#' transition branches:
#'
#' ```
#' Q(t) = g_remain * Q_rem
#'       + (1 - g_remain) * [ g_death_exit * Q_death
#'                           + (1 - g_death_exit) * Q_discharge ]
#' ```
#'
#' where `g_remain` and `g_death_exit` are the model predictions at the
#' subject's current history, and `Q_death` / `Q_discharge` are the outcome
#' values assigned to each exit branch (constants when [absorb_rule()] is
#' used, or `sl_y` predictions otherwise).
#'
#' **Event counts and the death regression.**  Because `g_death_exit` is
#' trained on the exit subset only, its effective sample size is often much
#' smaller than the full cohort -- and within that subset, deaths may be rare.
#'
#' Before every SuperLearner fit, the code checks whether the training data
#' are sufficient to fit a model (internal function `can_fit_bin`): fitting
#' proceeds only when the training set has at least 30 observations and at
#' least 5 events of each class (deaths and discharges within the exit subset
#' for `g_death_exit`; remainers and exiters for `g_remain`).  When either
#' threshold is not met, the estimator substitutes the empirical proportion as
#' a time-constant prediction instead of fitting a SuperLearner.  The fallback
#' is recorded in `diagnostics$branch_cal` (the `p_*_const` columns are
#' non-`NA` when a constant was used).
#'
#' Even above the floor, with few events the SuperLearner ensemble will
#' typically collapse to `SL.mean` or a near-intercept logistic fit.
#' Setting `pool_g_death = TRUE` fits a single `g_death_exit` model across
#' all time points (with time included as a covariate), so the 30/5 thresholds
#' apply to the full follow-up pooled rather than each time point
#' individually.  This is the recommended remedy when per-time death counts
#' are sparse but the death hazard is roughly stable over time.
#'
#' To inspect event counts: `diagnostics$branch_cal` contains one row per
#' fold and time point; `n_tr_exit` is the exit-subset size, and the number
#' of deaths is approximately `n_tr_exit * mean_target` for the `g_death`
#' rows.
#'
#' When [absorb_rule()] fixes `Q_death` to a constant -- the common case --
#' the impact of a poorly fitted `g_death_exit` is limited: misspecification
#' shifts probability mass between the death and discharge branches of the
#' mixture, but the sequentially doubly robust correction via density-ratio
#' weighting partially compensates provided the treatment model is
#' well-specified.
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
#' 118(542):846-857. doi:10.1080/01621459.2021.1955691.
#'
#' Luedtke AR, Sofrygin O, van der Laan MJ, Carone M (2017). Sequential
#' Double Robustness in Right-Censored Longitudinal Models. arXiv:1705.02459.
#'
#' Rotnitzky A, Robins J, Babino L (2017). On the multiply robust estimation
#' of the mean of the g-functional. arXiv:1705.08582.
#'
#' Wu C, Benkeser D (2024). Nonparametric Efficient Estimation of Marginal
#' Structural Models using Targeted Machine Learning. arXiv:2408.10847.
#'
#' Williams NT, Diaz I (2023). lmtp: An R package for estimating the causal
#' effects of modified treatment policies. *Observational Studies*.
#'
#' Bang H, Robins JM (2005). Doubly robust estimation in missing data and
#' causal inference models. *Biometrics* 61(4):962-973.
#'
#' Haneuse S, Rotnitzky A (2013). Estimation of the effect of interventions
#' that modify the received treatment. *Statistics in Medicine*
#' 32(30):5260-5277.
#'
#' Diaz Munoz I, van der Laan MJ (2012). Population intervention causal
#' effects based on stochastic interventions. *Biometrics* 68(2):541-549.
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
