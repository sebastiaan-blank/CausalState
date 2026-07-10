\name{qreg}
\alias{qreg}
\title{Pure Q-recursion estimator for longitudinal MTPs}
\description{
Fits sequential Q-models backwards through time and returns the
plug-in (substitution) estimator \eqn{\hat\Psi = n^{-1}\sum_i \hat Q_1(H_{i1})}
under the modified treatment policy (MTP).  No EIF update is applied, so
the estimate does not depend on density-ratio weights.
}
\usage{
qreg(
  df,
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
  absorb = list(),
  policy_spec_fun = function(D_block, t, a_names) NULL,
  k = 2,
  seed = 1,
  parallel = FALSE,
  fold_workers = NULL,
  reg_workers = NULL,
  sl_workers = NULL,
  inner_v = 5L,
  cluster = NULL,
  pool_g_death = FALSE,
  weight_object = NULL,
  v = 5L,
  trim = 0.99
)
}
\arguments{
\item{df}{A \code{data.frame} in long format (one row per subject per time point).}

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
for scaling a continuous outcome to \verb{[0, 1]}. Inferred from data if
\code{NULL}.}

\item{bounds}{Numeric. Probability clipping bound for g and Q predictions.
Default \code{1e-5}.}

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

\item{cluster}{Character. Name of a cluster variable for cluster-robust
standard errors (e.g. hospital). If \code{NULL}, subject id is used.}

\item{pool_g_death}{Logical. If \code{TRUE}, fit a single pooled
\code{g_death_exit} model across all time points (with time as a covariate)
rather than a separate model per time point. Default \code{FALSE}.}

\item{weight_object}{Optional. Output from \code{\link[=density_ratio]{density_ratio()}}.  If supplied,
an EIF-based SE is computed in addition to the naive SE.  The fold
structure is inherited from \code{weight_object} when provided; otherwise
independent folds are created using \code{v} and \code{seed}.}

\item{v}{Integer. Number of cross-fitting folds used when
\code{weight_object = NULL}. Default \code{5L}.}

\item{trim}{Quantile for capping density ratios before EIF-based SE
computation (only used when \code{weight_object} is supplied).
Default \code{0.99}.}
}
\value{
A list with:
\describe{
\item{\code{estimate}}{Point estimate \eqn{\hat\Psi}.}
\item{\code{se_naive}}{Naive SE (always present).}
\item{\code{se_eif}}{EIF-based SE (\code{NULL} if no \code{weight_object}).}
\item{\code{se}}{Best available SE: \code{se_eif} if available, else \code{se_naive}.}
\item{\code{ci}}{95\% Wald CI using \code{se}.}
\item{\code{psi_natural}}{Plug-in estimate under the natural course.}
\item{\code{psi_shifted}}{Plug-in estimate under the MTP (= \code{estimate}).}
\item{\code{sl_summary}}{SuperLearner weights per fold/time/component.}
\item{\code{fold_diag}}{Per-fold diagnostics.}
\item{\code{diagnostics}}{List with \code{branch_cal} and \code{diag_table}.}
}
}
\details{
Two standard-error options are available:
\itemize{
\item \strong{Naive SE} (always computed): \eqn{s(\hat Q_1) / \sqrt{n}},
treating the fitted Q as fixed.  Under-covers because it ignores
Q-model estimation uncertainty.
\item \strong{EIF-based SE} (requires \code{weight_object}): applies the
efficient influence function formula using the cross-fitted Q predictions
and pre-computed density-ratio weights.  \eqn{\hat\Psi} itself is
unchanged; only the SE improves.
}

The primary use case is as a sensitivity check alongside \code{\link[=sdr]{sdr()}} or
\code{\link[=itmle]{itmle()}}: when density-ratio weights are extreme (ESS collapse),
comparing \eqn{\hat\Psi_Q} to the doubly-robust estimates shows how much
the answer is driven by the weighting step.
}
\seealso{
\code{\link[=sdr]{sdr()}}, \code{\link[=itmle]{itmle()}}, \code{\link[=density_ratio]{density_ratio()}}
}
