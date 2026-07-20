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
