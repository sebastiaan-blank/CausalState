\name{print.sdr_fit}
\alias{print.sdr_fit}
\alias{print.itmle_fit}
\alias{print.qreg_fit}
\alias{print.branch_cal_summary}
\title{Print methods for CausalState output objects}
\usage{
\method{print}{sdr_fit}(x, digits = 3L, ...)

\method{print}{itmle_fit}(x, digits = 3L, ...)

\method{print}{qreg_fit}(x, digits = 3L, ...)

\method{print}{branch_cal_summary}(x, digits = 3L, ...)
}
\arguments{
\item{x}{A fitted object returned by \code{sdr()}, \code{itmle()},
\code{qreg()}, or \code{branch_cal_summary()}.}

\item{digits}{Number of decimal places for numeric output. Default \code{3L}.}

\item{...}{Ignored.}
}
\description{
Compact console summaries for the main CausalState output objects.

\code{print.sdr_fit} and \code{print.itmle_fit} display the point estimate,
standard error, 95\% Wald CI, and the observed outcome mean.

\code{print.qreg_fit} displays the shifted estimate, SE, 95\% CI, and the
natural-course plug-in for comparison.

\code{print.branch_cal_summary} renders a wide table with super-column
headers per branch (\code{g_remain}, \code{g_death}, \code{q_rem},
\code{q_exit}), showing fold-weighted \code{n_vl}, target mean,
predicted mean, and calibration slope for each time step.  The underlying
\code{data.table} (returned invisibly) contains all additional metrics.
}
\value{
The input object, invisibly.
}
\seealso{
\code{\link{sdr}}, \code{\link{itmle}}, \code{\link{qreg}},
\code{\link{branch_cal_summary}}
}
