#' CausalState: SDR and iTMLE for State-Aware Longitudinal Modified Treatment Policies
#'
#' Sequential Doubly Robust (SDR) and infinite-dimensional Targeted Minimum
#' Loss-based Estimation (iTMLE) estimators for longitudinal modified treatment
#' policies in settings with absorbing state transitions (e.g. ICU, ward, ED
#' care episodes).  Both estimators are sequentially doubly robust (SDR):
#' consistent whenever, at each time point, either the treatment model or the
#' outcome model is correctly specified.  Treatment is applicable only in active
#' states and becomes structurally inapplicable after transition to an absorbing
#' state (e.g. discharge or death).  Supports asymmetric g- and Q-model
#' regularisation, k-fold cross-fitting, and pluggable SuperLearner ensembles.
#'
#' @section Main estimators:
#' - [sdr()] — Sequential Doubly Robust estimator (EIF-based pseudo-outcome update)
#' - [itmle()] — Infinite-dimensional TMLE (infinite-dimensional targeting step; also SDR)
#'
#' @section State machinery:
#' - [absorb_rule()] — Define rules for entering absorbing states
#' - [sdr_competing()] — SDR estimator for competing-event outcomes under state transitions
#'
#' @section Custom SuperLearner wrappers:
#' SDR-targeted learners ([SL.tgt.glm()], [SL.tgt.glmnet()],
#' [SL.tgt.xgboost()], etc.) and TMLE learners ([SL.tmle_glm()],
#' [SL.tmle_glmnet_ridge()], etc.) for use in `SL.library` arguments.
#'
#' @keywords internal
"_PACKAGE"
