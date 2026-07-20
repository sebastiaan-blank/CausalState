\name{method.WB_dr}
\alias{method.WB_dr}
\title{Wu-Benkeser log-density-ratio metalearner for SuperLearner}
\usage{
method.WB_dr(dr_floor = 1e-10)
}
\arguments{
\item{dr_floor}{Numeric scalar. Floor applied to density-ratio predictions
from each base learner before computing the log-DR loss.  Prevents
\code{log(0)} during optimisation.  Default \code{1e-10}.}
}
\value{
A SuperLearner method list (with \code{computeCoef} and
\code{computePred} slots) suitable for passing to the \code{method}
argument of \code{\link[SuperLearner]{SuperLearner}}.
}
\description{
A SuperLearner metalearner that combines base learners by minimising a
log-density-ratio loss rather than the default non-negative least squares
(NNLS).  Weights are found by BFGS optimisation over the softmax-reparameterised
simplex; if BFGS fails to converge, equal weights are used as a fallback.

The metalearner operates in density-ratio space: base learners must return
density ratios (not probabilities).  Standard SuperLearner wrappers such as
\code{SL.glm} or \code{SL.xgboost} return probabilities and are not compatible
with this metalearner.  Custom DR-returning wrappers are required; built-in
wrappers for the WB pathway are planned for version 1.0.

\strong{Clipping note:} the standard \code{bounds} clipping applied to all
probability-scale predictions elsewhere in the pipeline does not apply here.
DR-scale flooring is handled by \code{dr_floor} inside this function.
}
\references{
Wu and Benkeser (2021). Improved inference for vaccine-induced immune
responses via shape-constrained methods.
}
\seealso{
\code{\link{density_ratio}}
}
