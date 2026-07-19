\name{itmle}
\alias{itmle}
\title{Infinite-dimensional TMLE (iTMLE) estimator for longitudinal MTPs}
\description{
Estimates the mean counterfactual outcome under a modified treatment policy
(MTP) using the infinite-dimensional targeted minimum loss-based estimator
(iTMLE) of Luedtke et al. (2017), Section 5 / Algorithm 4.  SDR and iTMLE
share the same backward Q-regression; they differ only in the update step:
iTMLE applies an infinite-dimensional TMLE fluctuation rather than the
EIF-based pseudo-outcome update used by \code{\link[=sdr]{sdr()}}.  Both are sequentially
doubly robust (SDR): consistent whenever, at each time point t, either the
treatment model \eqn{g_t} or the outcome model \eqn{Q_t} is consistently
estimated (2^\eqn{K}-robust, Luedtke et al. Definition 2).  The targeting
step uses a SuperLearner ensemble with custom learner wrappers (see
\link{sl_itmle}) that handle the logit-offset structure required for
cross-validated TMLE.  Requires pre-computed density ratio weights from
\code{\link[=density_ratio]{density_ratio()}}.
}
\usage{
itmle(
  df,
  weight_object,
  tmax,
  id,
  time,
  alive,
  in_state,
  y,
  baseline,
  tv_names = character(0),
  a_names = character(0),
  no_lag_vars = character(0),
  policy_names = character(0),
  sl_g = NULL,
  sl_recursive = NULL,
  sl_y = NULL,
  outcome_family = c("binomial", "gaussian"),
  y_bounds = NULL,
  bounds = 1e-05,
  trim = 0.99,
  absorb = list(),
  policy_spec_fun = function(D_block, t, a_names) NULL,
  k = 2,
  seed = 1,
  parallel = FALSE,
  fold_workers = NULL,
  reg_workers = NULL,
  sl_workers = NULL,
  inner_v = 5L,
  sl_tmle = NULL,
  v_target_itmle = 10L,
  v_sl_inner_itmle = 10L,
  cluster = NULL,
  cluster_se_only = FALSE,
  pool_g_death = FALSE
)
}
\arguments{
\item{df}{A \code{data.frame} in long format (one row per subject per time point).}

\item{weight_object}{Output from \code{\link[=density_ratio]{density_ratio()}}, or a \code{data.frame}
containing at least columns for subject id, time, \code{Rt_t}
(instantaneous density ratio), and \code{global_fold} (cross-fitting fold
assignment).}

\item{tmax}{Integer. Maximum follow-up time (number of time points).}

\item{id}{Character. Name of the subject identifier column.}

\item{time}{Character. Name of the integer time column.}

\item{alive}{Character. Name of the binary alive indicator column
(\code{1} = alive, \code{0} = dead).}

\item{in_state}{Character. Name of the binary active-state indicator column
(\code{1} = in active state e.g. ICU, \code{0} = transitioned out).}

\item{y}{Character. Name of the outcome column.}

\item{baseline}{Character vector of baseline (time-invariant) covariate
names.}

\item{tv_names}{Character vector of time-varying covariate names (excluding
treatment).}

\item{a_names}{Character vector of treatment variable names.}

\item{no_lag_vars}{Character vector of variables in \code{tv_names} or \code{a_names}
that should not be lagged (e.g. already-encoded temporal features).}

\item{policy_names}{Character vector of column names holding the shifted
treatment values under the MTP (one per treatment variable in \code{a_names}).}

\item{sl_g}{SuperLearner library for the treatment/censoring models
(\code{g_remain} and \code{g_death_exit}). Required.}

\item{sl_recursive}{SuperLearner library for the recursive Q-remain model.
Required.}

\item{sl_y}{SuperLearner library for the Q-exit (outcome-at-exit) model.
Required.}

\item{outcome_family}{\code{"binomial"} (default) or \code{"gaussian"}.}

\item{y_bounds}{Optional numeric vector of length 2 giving \code{c(min, max)}
for scaling a continuous outcome to \verb{[0, 1]}.  When \code{NULL} the
bounds are inferred from the observed range of \code{y} in \code{df}
(\code{min(y)}, \code{max(y)}).  Predictions that fall outside the scaling
range are clipped to \verb{[0, 1]} before back-transformation, so if the
outcome distribution under the MTP may extend beyond the training range,
supply explicit bounds wide enough to cover it.  Ignored for
\code{outcome_family = "binomial"}.}

\item{bounds}{Numeric. Probability clipping bound for g and Q predictions.
Default \code{1e-5}.}

\item{trim}{Quantile used to cap instantaneous density ratios before they
are assembled into the cumulative weight matrix.  Trimming is applied to
the \strong{full} \code{weights_dt} (all time points computed by \code{\link[=density_ratio]{density_ratio()}})
before any subsetting to the current \code{tmax}.  Consequently the trim
threshold is identical whether you call the estimator with \code{tmax = 2} or
\code{tmax = 14}, making results directly comparable across horizons that share
the same weight object.  The same \code{trim} value must be passed to
\code{\link[=sdr_competing]{sdr_competing()}} when reusing the same weight object, so
that the competing-event analysis operates on identically trimmed weights.
Default \code{0.99}.}

\item{absorb}{List of \code{\link[=absorb_rule]{absorb_rule()}} objects specifying outcome overrides
at absorbing states.}

\item{policy_spec_fun}{A function \verb{(D_block, t, a_names)} returning a
\code{data.table} with the shifted treatment values for time \code{t}. Used when
the policy cannot be pre-computed as static columns.}

\item{k}{Integer. Number of lags of time-varying covariates and treatment
to include in models. Default \code{2}.}

\item{seed}{Integer random seed. Default \code{1}.}

\item{parallel}{Logical. Enable parallel outer cross-fitting via
\code{\link[parallel:mclapply]{parallel::mclapply()}}. Default \code{FALSE}.}

\item{fold_workers}{Integer number of worker processes for outer folds.
\code{NULL} uses \code{\link[parallel:detectCores]{parallel::detectCores()}}.}

\item{reg_workers}{Integer number of worker processes for within-fold
regression parallelism.}

\item{sl_workers}{Integer. Workers for parallel learner evaluation within
each SuperLearner call via \code{\link[SuperLearner:mcSuperLearner]{SuperLearner::mcSuperLearner()}} (fork-based;
Linux/Mac only). Ignored when \code{parallel = FALSE}. Default \code{NULL}
(sequential). \strong{Note:} incompatible with \code{dbarts}/\code{SL.dbarts} learners —
dbarts' C-level RNG state does not survive forking. Remove dbarts from
\code{sl_g}/\code{sl_y}/\code{sl_recursive} when \code{sl_workers} is set.}

\item{inner_v}{Integer. Number of inner cross-validation folds for
SuperLearner. Default \code{5L}.}

\item{sl_tmle}{SuperLearner library for the iTMLE targeting step. Should
consist of learners from \link{sl_itmle} that accept an offset column
(\code{._sl_offset}). Defaults to the package-provided \code{sl_tmle} vector.}

\item{v_target_itmle}{Integer. Number of cross-validation folds for the
outer targeting loop. Default \code{10L}.}

\item{v_sl_inner_itmle}{Integer. Number of inner CV folds inside the
targeting SuperLearner. Default \code{10L}.}

\item{cluster}{Character. Name of a cluster variable for cluster-robust
standard errors (e.g. hospital). If \code{NULL}, subject id is used.}

\item{cluster_se_only}{Logical. If \code{TRUE}, skip fitting and return only
the fold structure for SE computation. Default \code{FALSE}.}

\item{pool_g_death}{Logical. If \code{TRUE}, fit a single pooled
\code{g_death_exit} model across all time points (with time as a covariate)
rather than a separate model per time point. Default \code{FALSE}.}
}
\value{
A named list with:
\describe{
\item{\code{psi}}{Point estimate of \code{E[Y(d)]} under the MTP (after targeting).}
\item{\code{se}}{Standard error from the efficient influence curve.}
\item{\code{ci}}{95\% Wald confidence interval.}
\item{\code{ic}}{Per-subject influence curve values (length n).}
\item{\code{psi_onestep}}{One-step pre-targeting estimate (diagnostic).}
\item{\code{targeting_gap}}{Mean EIF after targeting; near zero indicates convergence.}
\item{\code{sl_summary}}{\code{data.table} of SuperLearner weights per fold/time/component.}
\item{\code{fold_diag}}{Per-fold diagnostic summaries.}
\item{\code{target_diag}}{Per-iteration targeting diagnostics.}
\item{\code{fluct_diag}}{Fluctuation model diagnostics from the targeting step.}
}
}
\references{
Luedtke AR, Sofrygin O, van der Laan MJ, Carone M (2017). Sequential
Double Robustness in Right-Censored Longitudinal Models. arXiv:1705.02459.

Díaz I, Williams N, Hoffman KL, Schenck EJ (2021). Nonparametric Causal
Effects Based on Longitudinal Modified Treatment Policies. \emph{JASA}
118(542):846–857.
}
\seealso{
\code{\link[=density_ratio]{density_ratio()}}, \code{\link[=sdr]{sdr()}}, \code{\link[=absorb_rule]{absorb_rule()}}, \link{sl_itmle}
}
