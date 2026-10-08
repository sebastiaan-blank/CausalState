
#' Covariate balance diagnostics via standardised differences
#'
#' Computes unweighted and weighted standardised mean differences (SMD) at
#' each time point to assess whether the density ratio weights from [fit_g()]
#' achieve covariate balance between the treated and control groups.
#'
#' For binary treatment the SMD at time \eqn{t} for covariate \eqn{X} is:
#' \deqn{SMD_t = \frac{\bar{X}_{t,1} - \bar{X}_{t,0}}{s_{\text{pool}}}}
#' where \eqn{s_{\text{pool}}} is the unweighted pooled standard deviation
#' (used as a common reference scale for both unweighted and weighted SMDs).
#' The weighted version uses cumulative IPW weights
#' \eqn{W_i(t) = \prod_{s=1}^{t} r_s(i)} to compute weighted group means.
#'
#' For continuous treatment the reported statistic is the Pearson correlation
#' between treatment and each covariate (unweighted and IPW-weighted).
#'
#' A conventional threshold for acceptable balance is \eqn{|SMD| < 0.10}.
#'
#' @param fit A `"g_fit"` object returned by [fit_g()].
#' @param df The original data frame passed to [fit_g()] (long format, one row
#'   per subject per time point).
#' @param covariates Character vector of covariate names to assess. Defaults
#'   to `NULL`, which uses all variables passed to [fit_g()] (`baseline` and
#'   `tv_names` from the stored settings). Pass an explicit vector to check
#'   any variable in `df`, including variables not in the model.
#' @param tmax Integer. Maximum time point to evaluate. Defaults to the value
#'   stored in `fit$settings$tmax`.
#' @param id Character. Subject id column. Defaults to `fit$settings$id`.
#' @param time Character. Time column. Defaults to `fit$settings$time`.
#'
#' @return A `data.table` (class `"balance_check"`) with columns:
#'   \describe{
#'     \item{`t`}{Time point.}
#'     \item{`covariate`}{Covariate name.}
#'     \item{`smd_unweighted`}{Unweighted SMD (or correlation for continuous
#'       treatment).}
#'     \item{`smd_weighted`}{IPW-weighted SMD (or correlation).}
#'     \item{`abs_improvement`}{Absolute reduction in |SMD| after weighting.}
#'     \item{`n`}{Number of subjects at risk at time \eqn{t}.}
#'   }
#'
#' @seealso [fit_g()], [weight_diagnostics()]
#'
#' @examples
#' \donttest{
#' # After fitting: wr <- fit_g(df, ..., method = "propensity")
#' # bc <- balance_check(wr, df, covariates = c("age", "SOFA", "lactate"))
#' # print(bc)
#' }
#'
#' @export
balance_check <- function(fit, df, covariates = NULL, tmax = NULL, id = NULL, time = NULL) {

  `%||%` <- function(x, y) if (!is.null(x)) x else y

  if (!inherits(fit, "g_fit"))
    stop("`fit` must be a \"g_fit\" object returned by fit_g().", call. = FALSE)

  s    <- fit$settings
  id   <- id   %||% s$id
  time <- time %||% s$time
  tmax <- as.integer(tmax %||% s$tmax)

  a_names   <- s$a_names
  is_binary <- isTRUE(s$is_binary)

  if (is.null(covariates)) {
    covariates <- unique(c(
      as.character(s$baseline %||% character(0)),
      as.character(s$tv_names %||% character(0))
    ))
    if (!length(covariates))
      stop("No covariates found in fit$settings. Pass covariates explicitly.", call. = FALSE)
  }
  covariates <- as.character(covariates)

  # Build cumulative weights from instantaneous Rt_t
  wdt <- data.table::copy(extract_weights_dt(fit))
  data.table::setnames(wdt, c(id, time), c(".__id", ".__time"), skip_absent = TRUE)
  if (!all(c(".__id", ".__time") %in% names(wdt))) {
    wdt[, .__id   := as.character(get(id))]
    wdt[, .__time := as.integer(get(time))]
  } else {
    wdt[, .__id   := as.character(.__id)]
    wdt[, .__time := as.integer(.__time)]
  }
  data.table::setorderv(wdt, c(".__id", ".__time"))
  wdt[, Rt_t := pmax(0, Rt_t, na.rm = TRUE)]
  wdt[, cum_w := cumprod(Rt_t), by = .__id]
  wdt <- wdt[.__time >= 1L & .__time <= tmax]

  DT <- data.table::as.data.table(df)
  DT[, .__id   := as.character(get(id))]
  DT[, .__time := as.integer(get(time))]

  DT <- merge(
    DT,
    wdt[, .(.__id, .__time, cum_w)],
    by  = c(".__id", ".__time"),
    all.x = TRUE
  )
  DT[is.na(cum_w), cum_w := 1]

  a_col <- a_names[1L]
  if (!a_col %in% names(DT))
    stop(sprintf("Treatment column '%s' not found in df.", a_col), call. = FALSE)

  missing_cov <- setdiff(covariates, names(DT))
  if (length(missing_cov))
    warning(sprintf("Covariates not found in df, skipped: %s",
                    paste(missing_cov, collapse = ", ")), call. = FALSE)
  covariates <- intersect(covariates, names(DT))
  if (!length(covariates))
    stop("No valid covariates found in df.", call. = FALSE)

  results <- vector("list", tmax * length(covariates))
  idx <- 1L

  for (tt in seq_len(tmax)) {
    DT_t <- DT[.__time == tt]
    n_t  <- nrow(DT_t)
    if (n_t < 4L) next

    A_t <- as.numeric(DT_t[[a_col]])
    W_t <- as.numeric(DT_t$cum_w)
    W_t[!is.finite(W_t) | W_t < 0] <- 1

    for (cov in covariates) {
      X_t <- as.numeric(DT_t[[cov]])
      ok  <- is.finite(A_t) & is.finite(X_t) & is.finite(W_t)
      if (sum(ok) < 4L) next

      A_ok <- A_t[ok]; X_ok <- X_t[ok]; W_ok <- W_t[ok]

      if (is_binary) {
        grp1 <- A_ok == 1; grp0 <- A_ok == 0
        if (sum(grp1) < 2L || sum(grp0) < 2L) next

        m1_uw <- mean(X_ok[grp1]); m0_uw <- mean(X_ok[grp0])
        pooled_sd <- sqrt((var(X_ok[grp1]) + var(X_ok[grp0])) / 2)
        if (!is.finite(pooled_sd) || pooled_sd <= 0) next

        smd_uw <- (m1_uw - m0_uw) / pooled_sd

        w1 <- W_ok[grp1]; w0 <- W_ok[grp0]
        sw1 <- sum(w1); sw0 <- sum(w0)
        if (sw1 <= 0 || sw0 <= 0) next
        m1_w <- sum(w1 * X_ok[grp1]) / sw1
        m0_w <- sum(w0 * X_ok[grp0]) / sw0
        smd_w <- (m1_w - m0_w) / pooled_sd

      } else {
        # Continuous treatment: weighted Pearson correlation
        smd_uw <- stats::cor(A_ok, X_ok)

        sw  <- sum(W_ok)
        mAw <- sum(W_ok * A_ok) / sw
        mXw <- sum(W_ok * X_ok) / sw
        cov_w  <- sum(W_ok * (A_ok - mAw) * (X_ok - mXw)) / sw
        sdA_w  <- sqrt(sum(W_ok * (A_ok - mAw)^2) / sw)
        sdX_w  <- sqrt(sum(W_ok * (X_ok - mXw)^2) / sw)
        smd_w  <- if (sdA_w > 0 && sdX_w > 0) cov_w / (sdA_w * sdX_w) else NA_real_
      }

      results[[idx]] <- data.table::data.table(
        t              = as.integer(tt),
        covariate      = cov,
        smd_unweighted = round(smd_uw, 4L),
        smd_weighted   = round(smd_w,  4L),
        n              = n_t
      )
      idx <- idx + 1L
    }
  }

  if (idx == 1L) {
    warning("No balance statistics could be computed.", call. = FALSE)
    return(invisible(NULL))
  }

  out <- data.table::rbindlist(results[seq_len(idx - 1L)], use.names = TRUE)
  out[, abs_improvement := round(abs(smd_unweighted) - abs(smd_weighted), 4L)]

  attr(out, "is_binary") <- is_binary
  attr(out, "a_names")   <- a_names
  class(out) <- c("balance_check", class(out))
  out
}


#' @export
print.balance_check <- function(x, ...) {
  is_binary <- isTRUE(attr(x, "is_binary"))
  n_t   <- length(unique(x$t))
  n_cov <- length(unique(x$covariate))

  cat(sprintf(
    "Balance check  (%s,  %d time points,  %d covariates)\n",
    if (is_binary) "binary treatment — SMD" else "continuous treatment — weighted correlation",
    n_t, n_cov
  ))
  cat(sprintf(
    "  Max |unweighted|: %.3f    Max |weighted|: %.3f    (threshold |SMD| < 0.10)\n\n",
    max(abs(x$smd_unweighted), na.rm = TRUE),
    max(abs(x$smd_weighted),   na.rm = TRUE)
  ))

  # Per-time summary: max |SMD| and worst covariate at each t
  per_t <- x[, .(
    max_smd_uw  = round(max(abs(smd_unweighted), na.rm = TRUE), 3),
    max_smd_w   = round(max(abs(smd_weighted),   na.rm = TRUE), 3),
    worst_cov   = covariate[which.max(abs(smd_unweighted))],
    n_cov_gt10  = sum(abs(smd_weighted) > 0.10, na.rm = TRUE),
    n           = n[1L]
  ), by = t][order(t)]

  print(per_t, row.names = FALSE)
  cat("\n(Full covariate x time detail in the returned data.table)\n")
  invisible(x)
}
