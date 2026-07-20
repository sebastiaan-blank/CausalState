\name{sl_tmle}
\alias{sl_tmle}
\title{Default SuperLearner library for the iTMLE targeting step}
\usage{
sl_tmle
}
\format{
A character vector of SuperLearner learner names.
}
\description{
A character vector naming the default library of learners used in the
iTMLE targeting step.  Contains empty, intercept-only, GLM, and XGBoost
wrappers at three tree depths.  All wrappers carry the logit offset
required for cross-validated TMLE as a data column (the
\dQuote{offset-as-column} convention) and clip internal probability
predictions at \code{1e-5}.

Pass this vector to the \code{sl_target} argument of \code{\link{itmle}}
or construct your own library from the \code{SL.tmle_*} wrappers exported
by this package.
}
\seealso{
\code{\link{itmle}}, \code{\link{SL.tmle_glm}},
\code{\link{SL.tmle_glmnet_ridge}}, \code{\link{SL.tmle_xgb_d1}}
}
