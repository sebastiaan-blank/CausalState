# helpers.R \u2014 user-facing print methods, inference helpers, and diagnostics
#
# print.sdr_fit             \u2014 print method for sdr() output
# print.itmle_fit           \u2014 print method for itmle() output
# print.qreg_fit            \u2014 print method for qreg() output
# contrast()                \u2014 RD / RR / OR between two fitted estimators
# print.CausalState_contrast
# weight_diagnostics()      \u2014 per-time weight summary for a density_ratio object
# branch_cal_summary()      \u2014 fold-weighted calibration table per t x branch
# print.branch_cal_summary


# ------------------------------------------------------------------
# print methods for estimator output objects
# ------------------------------------------------------------------

.fmt_est <- function(label, est, se, ci, digits) {
  cat(sprintf("  %-12s  %.*f   SE %.*f   95%% CI [%.*f, %.*f]\n",
              label,
              digits, est,
              digits, se,
              digits, ci[1],
              digits, ci[2]))
}

#' @export
print.sdr_fit <- function(x, digits = 3L, ...) {
  vi <- x$settings$variable_info
  cat(sprintf("SDR  (outcome: %s, family: %s, tmax: %d)\n",
              vi$outcome, x$settings$outcome_family, vi$tmax))
  cat(sprintf("n = %d subjects\n\n", nrow(x$ic_df)))
  .fmt_est("Intervention", x$psi, x$se, x$ci, digits)
  if (!is.null(x$Y_obs))
    cat(sprintf("  %-12s  %.*f\n", "Observed mean", digits, x$Y_obs))
  invisible(x)
}

#' @export
print.itmle_fit <- function(x, digits = 3L, ...) {
  vi <- x$settings$variable_info
  cat(sprintf("iTMLE  (outcome: %s, family: %s, tmax: %d)\n",
              vi$outcome, x$settings$outcome_family, vi$tmax))
  cat(sprintf("n = %d subjects\n\n", nrow(x$ic_df)))
  .fmt_est("Intervention", x$psi, x$se, x$ci, digits)
  if (!is.null(x$Y_obs))
    cat(sprintf("  %-12s  %.*f\n", "Observed mean", digits, x$Y_obs))
  invisible(x)
}

#' @export
print.qreg_fit <- function(x, digits = 3L, ...) {
  vi <- x$settings$variable_info
  cat(sprintf("Q-regression  (outcome: %s, family: %s, tmax: %d)\n",
              vi$outcome, x$settings$outcome_family, vi$tmax))
  cat(sprintf("n = %d subjects\n\n", x$n))
  .fmt_est("Shifted", x$estimate, x$se, x$ci, digits)
  if (!is.null(x$psi_natural))
    cat(sprintf("  %-12s  %.*f\n", "Natural", digits, x$psi_natural))
  invisible(x)
}


# ------------------------------------------------------------------
# contrast -- risk difference, risk ratio, and odds ratio
# ------------------------------------------------------------------

#' Compute contrasts between two fitted estimators
#'
#' @description
#' Takes an intervention-arm fit and a reference (from [sdr()] or [itmle()], or the
#' crude observed mean) and returns the risk difference (RD), risk ratio (RR), and
#' odds ratio (OR) with standard errors derived from the per-subject influence curves.
#'
#' @details
#' All standard errors use the delta method on the efficient influence curves:
#' \itemize{
#'   \item RD: \eqn{IC_{RD,i} = IC_{1,i} - IC_{0,i}}; Wald CI on the natural scale.
#'   \item RR: \eqn{IC_{\log RR,i} = IC_{1,i}/\psi_1 - IC_{0,i}/\psi_0};
#'     SE and 95\% CI computed on the log scale, then exponentiated.
#'   \item OR: \eqn{IC_{\log OR,i} = IC_{1,i}/(\psi_1(1-\psi_1)) - IC_{0,i}/(\psi_0(1-\psi_0))};
#'     SE and 95\% CI computed on the log scale, then exponentiated.
#' }
#' RR and OR are omitted for Gaussian outcomes.
#'
#' When \code{fit0 = NULL} and \code{y_col} is supplied, the reference is the
#' crude observed mean of \code{y_col} (last non-missing value per subject).
#' The influence curve for the observed mean is \eqn{IC_{0,i} = Y_i - \bar{Y}}.
#' This gives a simple descriptive contrast against the raw data rather than a
#' causal comparison between two counterfactual estimates.
#'
#' @param fit1 Fitted object (intervention arm): output of [sdr()] or [itmle()].
#'   Must contain \code{$psi} and \code{$ic_df} with columns \code{id} and \code{ic}.
#' @param fit0 Fitted object (reference arm), same type as \code{fit1}, or
#'   \code{NULL} to use the crude observed mean as the reference (requires
#'   \code{df} and \code{y_col}).
#' @param df Optional long-format data frame. Required when \code{cluster} is
#'   specified, or when \code{fit0 = NULL}.
#' @param id_col Name of the subject-id column in \code{df}. Defaults to the
#'   value stored in \code{fit1$settings}.
#' @param cluster \code{NULL} (default, IID standard errors), or a single
#'   character string naming a column in \code{df} for cluster-robust SEs.
#' @param y_col Name of the outcome column in \code{df}. Required when
#'   \code{fit0 = NULL}; ignored otherwise.
#'
#' @return A list of class \code{"CausalState_contrast"} with:
#'   \describe{
#'     \item{\code{psi1}, \code{psi0}}{Point estimates.}
#'     \item{\code{obs_ref}}{Logical; \code{TRUE} when the reference is the
#'       observed mean rather than a second estimator fit.}
#'     \item{\code{RD}, \code{se_RD}, \code{ci_RD}}{Risk difference and 95 pct Wald CI.}
#'     \item{\code{RR}, \code{se_log_RR}, \code{ci_RR}}{Risk ratio, SE on log scale,
#'       and 95 pct CI (exponentiated). \code{NULL} for Gaussian outcomes.}
#'     \item{\code{OR}, \code{se_log_OR}, \code{ci_OR}}{Odds ratio, SE on log scale,
#'       and 95 pct CI (exponentiated). \code{NULL} for Gaussian outcomes.}
#'     \item{\code{n}}{Number of matched subjects.}
#'     \item{\code{table}}{Summary data frame printed by default.}
#'   }
#'
#' @seealso [sdr()], [itmle()]
#'
#' @examples
#' \donttest{
#' library(SuperLearner)
#'
#' sim_ex <- function(n = 2000L, tmax = 5L) {
#'   set.seed(42L)
#'   rows <- vector("list", n)
#'   for (i in seq_len(n)) {
#'     age <- round(rnorm(1, 65, 10)); L1 <- rnorm(1)
#'     pat <- list()
#'     for (t in seq_len(tmax)) {
#'       A     <- rbinom(1, 1, plogis(0.3 * L1 - 0.4))
#'       u     <- runif(1)
#'       p_die <- plogis(-4.0 + 0.2 * L1 - 0.1 * age / 10)
#'       p_dc  <- plogis(-2.5 + 0.5 * A)
#'       if (u < p_die) {
#'         alive <- 0L; in_state <- 0L
#'       } else if (u < p_die + p_dc) {
#'         alive <- 1L; in_state <- 0L
#'       } else {
#'         alive <- 1L; in_state <- 1L
#'       }
#'       Y <- rbinom(1, 1, plogis(-0.5 + 0.4 * A - 0.2 * L1))
#'       pat[[length(pat) + 1L]] <- data.frame(
#'         id = i, time = t, age = age, L1 = L1,
#'         A = A, alive = alive, in_state = in_state, Y = Y
#'       )
#'       if (in_state == 0L) break
#'       if (t < tmax) L1 <- L1 + rnorm(1, -0.1 * A, 0.3)
#'     }
#'     rows[[i]] <- do.call(rbind, pat)
#'   }
#'   do.call(rbind, rows)
#' }
#' df <- sim_ex()
#' sl_lib <- c("SL.mean", "SL.glm")
#'
#' policy_nat <- function(D_block, t, a_names) D_block[, ..a_names, drop = FALSE]
#' policy_sft <- function(D_block, t, a_names) {
#'   out <- D_block[, ..a_names, drop = FALSE]
#'   out[[a_names[1]]] <- pmin(D_block[[a_names[1]]] + 0.3, 1)
#'   out
#' }
#'
#' dr_args <- list(
#'   df = df, a_names = "A", tmax = 5L, baseline = "age", tv_names = "L1",
#'   sl_g = sl_lib, k = 1L, inner_v = 3L, v = 3L, seed = 1L,
#'   id = "id", time = "time"
#' )
#' wr_nat <- do.call(density_ratio, c(dr_args, list(policy_spec_fun = policy_nat)))
#' wr_sft <- do.call(density_ratio, c(dr_args, list(policy_spec_fun = policy_sft)))
#'
#' sdr_args <- list(
#'   df = df, tmax = 5L, id = "id", time = "time",
#'   alive = "alive", in_state = "in_state", y = "Y",
#'   baseline = "age", tv_names = "L1", a_names = "A",
#'   sl_remain = sl_lib, sl_death = sl_lib,
#'   sl_recursive = sl_lib, sl_y = sl_lib,
#'   k = 1L, inner_v = 3L, parallel = FALSE, seed = 1L
#' )
#' res_nat <- do.call(sdr, c(sdr_args,
#'   list(weight_object = wr_nat, policy_spec_fun = policy_nat)))
#' res_sft <- do.call(sdr, c(sdr_args,
#'   list(weight_object = wr_sft, policy_spec_fun = policy_sft)))
#'
#' # Two-estimator contrast (causal RD/RR/OR)
#' ctr <- contrast(res_sft, res_nat)
#' ctr$RD; ctr$ci_RD; ctr$RR; ctr$OR
#'
#' # Observed-mean reference (descriptive)
#' ctr_obs <- contrast(res_sft, df = df, y_col = "Y")
#' ctr_obs$table
#' }
#'
#' @export
contrast <- function(fit1, fit0 = NULL, df = NULL, id_col = NULL, cluster = NULL,
                     y_col = NULL) {
  psi1 <- fit1$psi

  vi <- fit1$settings$variable_info
  if (is.null(id_col))  id_col  <- if (!is.null(vi$id))      vi$id      else "id"
  if (is.null(cluster)) cluster <- if (!is.null(vi$cluster))  vi$cluster else NULL

  ic1 <- data.table::as.data.table(fit1$ic_df)[, .(id, ic1 = ic)]
  if (anyDuplicated(ic1$id))
    stop("contrast: fit1$ic_df has duplicate subject IDs.", call. = FALSE)

  obs_ref <- is.null(fit0)
  if (obs_ref) {
    if (is.null(df) || is.null(y_col))
      stop("contrast: provide fit0, or df + y_col for an observed-mean reference.",
           call. = FALSE)
    if (!id_col %in% names(df))
      stop(sprintf("contrast: id_col '%s' not found in df.", id_col), call. = FALSE)
    if (!y_col %in% names(df))
      stop(sprintf("contrast: y_col '%s' not found in df.", y_col), call. = FALSE)
    dt_y <- data.table::as.data.table(df)[
      !is.na(get(y_col)),
      .(y = as.numeric(get(y_col))[.N]),
      by = id_col
    ]
    data.table::setnames(dt_y, id_col, "id")
    if (anyDuplicated(dt_y$id))
      stop("contrast: df has duplicate subject IDs after extracting y_col.", call. = FALSE)
    psi0 <- mean(dt_y$y)
    ic0  <- dt_y[, .(id, ic0 = y - psi0)]
  } else {
    psi0 <- fit0$psi
    ic0  <- data.table::as.data.table(fit0$ic_df)[, .(id, ic0 = ic)]
    if (anyDuplicated(ic0$id))
      stop("contrast: fit0$ic_df has duplicate subject IDs.", call. = FALSE)
  }

  only1 <- sum(!ic1$id %in% ic0$id)
  only0 <- sum(!ic0$id %in% ic1$id)
  if (only1 > 0L || only0 > 0L)
    stop(sprintf(
      "contrast: ID mismatch -- %d subject(s) in fit1 only, %d in reference only. Both objects must cover exactly the same subjects.",
      only1, only0), call. = FALSE)

  merged <- merge(ic1, ic0, by = "id", all = FALSE)
  n <- nrow(merged)
  if (n == 0L) stop("contrast: no overlapping subject IDs between fit1 and fit0.", call. = FALSE)

  cl <- NULL
  if (!is.null(cluster)) {
    if (is.null(df))
      stop("contrast: `df` must be provided when a cluster column is set.", call. = FALSE)
    if (!is.character(cluster) || length(cluster) != 1L)
      stop("contrast: `cluster` must be a single column name string.", call. = FALSE)
    if (!id_col %in% names(df))
      stop(sprintf("contrast: id_col '%s' not found in df.", id_col), call. = FALSE)
    if (!cluster %in% names(df))
      stop(sprintf("contrast: cluster column '%s' not found in df.", cluster), call. = FALSE)
    cl_tbl <- unique(df[, c(id_col, cluster), drop = FALSE])
    cl_map <- setNames(cl_tbl[[cluster]], as.character(cl_tbl[[id_col]]))
    cl     <- cl_map[as.character(merged$id)]
  }

  .se <- function(ic_vec) {
    ok <- is.finite(ic_vec)
    if (is.null(cl)) {
      n_ok <- sum(ok)
      if (n_ok < 2L) return(NA_real_)
      return(stats::sd(ic_vec[ok]) / sqrt(n_ok))
    }
    ic_c  <- ic_vec[ok] - mean(ic_vec[ok])
    cl_ok <- cl[ok]
    dt    <- data.table::data.table(cl = cl_ok, ic = ic_c)
    S     <- dt[, .(S = sum(ic)), by = cl]
    G     <- nrow(S)
    n_ok  <- sum(ok)
    if (G < 2L) return(stats::sd(ic_vec[ok]) / sqrt(n_ok))
    sqrt((G / (G - 1L)) * sum(S$S^2) / n_ok^2)
  }

  se1   <- .se(merged$ic1)
  se0   <- .se(merged$ic0)

  ic_rd <- merged$ic1 - merged$ic0
  rd    <- psi1 - psi0
  se_rd <- .se(ic_rd)

  gaussian <- identical(fit1$settings$outcome_family, "gaussian") ||
    (!obs_ref && identical(fit0$settings$outcome_family, "gaussian"))

  rr <- se_log_rr <- or <- se_log_or <- NULL
  if (!gaussian) {
    ic_log_rr <- merged$ic1 / psi1 - merged$ic0 / psi0
    rr         <- psi1 / psi0
    se_log_rr  <- .se(ic_log_rr)

    ic_log_or <- merged$ic1 / (psi1 * (1 - psi1)) - merged$ic0 / (psi0 * (1 - psi0))
    or         <- (psi1 / (1 - psi1)) / (psi0 / (1 - psi0))
    se_log_or  <- .se(ic_log_or)
  }

  ref_lbl <- if (obs_ref) "Observed" else "Control"

  if (gaussian) {
    tbl <- data.frame(
      Estimate   = c(psi1, psi0, rd),
      `CI Lower` = c(psi1 - 1.96 * se1, psi0 - 1.96 * se0, rd - 1.96 * se_rd),
      `CI Upper` = c(psi1 + 1.96 * se1, psi0 + 1.96 * se0, rd + 1.96 * se_rd),
      check.names = FALSE,
      row.names   = c("Intervention", ref_lbl, "Mean Difference")
    )
  } else {
    tbl <- data.frame(
      Estimate   = c(psi1, psi0, rd, rr, or),
      `CI Lower` = c(
        psi1 - 1.96 * se1,
        psi0 - 1.96 * se0,
        rd   - 1.96 * se_rd,
        exp(log(rr) - 1.96 * se_log_rr),
        exp(log(or) - 1.96 * se_log_or)
      ),
      `CI Upper` = c(
        psi1 + 1.96 * se1,
        psi0 + 1.96 * se0,
        rd   + 1.96 * se_rd,
        exp(log(rr) + 1.96 * se_log_rr),
        exp(log(or) + 1.96 * se_log_or)
      ),
      check.names = FALSE,
      row.names   = c("Intervention", ref_lbl, "Risk Difference", "Risk Ratio", "Odds Ratio")
    )
  }

  out <- list(
    psi1      = psi1,
    psi0      = psi0,
    obs_ref   = obs_ref,
    se1       = se1,
    se0       = se0,
    RD        = rd,
    se_RD     = se_rd,
    ci_RD     = c(rd - 1.96 * se_rd, rd + 1.96 * se_rd),
    RR        = rr,
    se_log_RR = se_log_rr,
    ci_RR     = if (is.null(rr)) NULL else exp(log(rr) + c(-1.96, 1.96) * se_log_rr),
    OR        = or,
    se_log_OR = se_log_or,
    ci_OR     = if (is.null(or)) NULL else exp(log(or) + c(-1.96, 1.96) * se_log_or),
    n         = n,
    gaussian  = gaussian,
    table     = tbl
  )
  class(out) <- "CausalState_contrast"
  out
}

#' @export
print.CausalState_contrast <- function(x, digits = 3L, ...) {
  tbl <- x$table
  tbl[] <- lapply(tbl, round, digits = digits)
  cat(sprintf("n = %d subjects\n\n", x$n))
  print(tbl)
  invisible(x)
}



# ------------------------------------------------------------------
# weight_diagnostics -- per-time weight summary for a density_ratio object
# ------------------------------------------------------------------

#' Per-time weight diagnostics for a density ratio object
#'
#' Summarises the instantaneous density ratios and their cumulative products at
#' each time step, restricted to subjects at risk at that time.  This matches
#' how the weights enter the SDR and iTMLE corrections: the EIF term at time t
#' uses the product of ratios from t=1 through t, and only at-risk subjects
#' contribute non-zero correction terms.
#'
#' Call with different \code{trim} values to assess trim sensitivity before
#' passing the \code{weight_object} to \code{sdr()} or \code{itmle()}.
#'
#' @param weight_object Output of \code{density_ratio()}.
#' @param trim Quantile used to cap \code{Rt_t} before computing summaries.
#'   Default \code{1} (no trimming) so the returned summaries reflect the raw
#'   density ratios; pass a value \code{< 1} (e.g. \code{0.99}) to see how a
#'   given trim reshapes the weight distribution before choosing what to hand
#'   to \code{sdr()}/\code{itmle()}.
#'
#' @return A \code{data.table} with one row per time step and columns:
#'   \code{t}, \code{n_risk};
#'   instantaneous ratio stats (\code{Rt_mean}, \code{Rt_median},
#'   \code{Rt_p05}, \code{Rt_p95}, \code{Rt_max}, \code{Rt_ess});
#'   cumulative product stats (\code{cum_mean}, \code{cum_median},
#'   \code{cum_p05}, \code{cum_p95}, \code{cum_max}, \code{cum_ess}).
#'   ESS uses the Kish approximation: \eqn{(\sum w)^2 / \sum w^2}.
#'
#' @seealso [density_ratio()], [sdr()], [itmle()]
#'
#' @export
weight_diagnostics <- function(weight_object, trim = 1) {
  weights_dt <- extract_weights_dt(weight_object)
  s    <- weight_object$settings
  id   <- s$id
  time <- s$time
  tmax <- as.integer(s$tmax)

  WW <- data.table::copy(weights_dt)

  if (!is.null(trim) && trim < 1) {
    trim_val <- stats::quantile(WW$Rt_t[is.finite(WW$Rt_t)], trim, na.rm = TRUE)
    WW[, Rt_t := pmin(Rt_t, trim_val)]
  }

  WW[, .__id   := as.character(get(id))]
  WW[, .__time := as.integer(get(time))]
  WW <- WW[.__time >= 1L & .__time <= tmax]
  WW <- WW[, .(Rt_t = data.table::first(Rt_t)), by = .(.__id, .__time)]

  all_ids <- unique(WW$.__id)
  N <- length(all_ids)

  wide <- data.table::dcast(WW, .__id ~ .__time, value.var = "Rt_t", fill = NA_real_)

  dens <- matrix(NA_real_, nrow = N, ncol = tmax)
  time_cols <- setdiff(names(wide), ".__id")
  for (tc in time_cols) {
    tt <- suppressWarnings(as.integer(tc))
    if (!is.na(tt) && tt >= 1L && tt <= tmax)
      dens[, tt] <- as.numeric(wide[[tc]])
  }

  at_risk <- !is.na(dens)

  dens_filled <- dens
  dens_filled[is.na(dens_filled)] <- 1
  cum_dens <- t(apply(dens_filled, 1L, cumprod))

  .ess <- function(w) {
    w <- w[is.finite(w) & w > 0]
    if (length(w) < 2L) return(NA_real_)
    sum(w)^2 / sum(w^2)
  }

  rows <- lapply(seq_len(tmax), function(tt) {
    ar <- at_risk[, tt]
    rt <- dens[ar, tt]
    cp <- cum_dens[, tt]
    data.table::data.table(
      t          = tt,
      n_risk     = sum(ar),
      Rt_mean    = mean_or_na(rt),
      Rt_median  = q_or_na(rt, 0.5),
      Rt_p05     = q_or_na(rt, 0.05),
      Rt_p95     = q_or_na(rt, 0.95),
      Rt_max     = max_or_na(rt),
      Rt_ess     = .ess(rt),
      cum_mean   = mean_or_na(cp),
      cum_median = q_or_na(cp, 0.5),
      cum_p05    = q_or_na(cp, 0.05),
      cum_p95    = q_or_na(cp, 0.95),
      cum_max    = max_or_na(cp),
      cum_ess    = .ess(cp)
    )
  })

  data.table::rbindlist(rows)
}


# ------------------------------------------------------------------
# branch_cal_summary -- fold-averaged calibration per t x branch
# ------------------------------------------------------------------

#' Per-time, per-branch calibration summary from an SDR or iTMLE fit
#'
#' Aggregates the raw per-fold calibration diagnostics stored in
#' \code{diagnostics$branch_cal} to one row per time step per branch,
#' weighted by the validation-set size at each fold.  Returns a tidy
#' \code{data.table} suitable for further filtering, plotting, or
#' \code{data.table::dcast()}.
#'
#' @param fit Output of \code{sdr()} or \code{itmle()}.
#'
#' @return A \code{data.table} with columns \code{t}, \code{branch}
#'   (one of \code{"g_remain"}, \code{"g_death"}, \code{"q_rem"},
#'   \code{"q_exit"}), \code{n_tr}, \code{tgt_tr}, \code{pred_tr},
#'   \code{n_vl}, \code{tgt_vl}, \code{pred_vl}, \code{cal_slope_vl};
#'   and where applicable \code{brier_vl}, \code{auc_vl} (binary
#'   branches) or \code{rmse_vl}, \code{cor_vl} (Gaussian branches).
#'
#' @seealso [sdr()], [itmle()]
#'
#' @export
branch_cal_summary <- function(fit) {
  bc <- fit$diagnostics$branch_cal
  if (is.null(bc) || nrow(bc) == 0L)
    stop("branch_cal_summary: no branch_cal diagnostics found in fit object.",
         call. = FALSE)

  bc <- data.table::as.data.table(bc)

  .wagg <- function(vals, wts) {
    ok <- is.finite(vals) & is.finite(wts) & wts > 0
    if (!any(ok)) return(NA_real_)
    sum(vals[ok] * wts[ok]) / sum(wts[ok])
  }

  .agg_branch <- function(bc, branch,
                          n_tr_col, n_vl_col,
                          tgt_tr_col, pred_tr_col,
                          tgt_vl_col, pred_vl_col,
                          slope_vl_col,
                          brier_vl_col  = NULL, auc_vl_col  = NULL,
                          rmse_vl_col   = NULL, cor_vl_col  = NULL) {
    bc[, {
      w_vl <- get(n_vl_col)
      data.table::data.table(
        branch       = branch,
        n_tr         = sum(get(n_tr_col),    na.rm = TRUE),
        tgt_tr       = .wagg(get(tgt_tr_col),  get(n_tr_col)),
        pred_tr      = .wagg(get(pred_tr_col), get(n_tr_col)),
        n_vl         = sum(w_vl,               na.rm = TRUE),
        tgt_vl       = .wagg(get(tgt_vl_col),  w_vl),
        pred_vl      = .wagg(get(pred_vl_col), w_vl),
        cal_slope_vl = .wagg(get(slope_vl_col), w_vl),
        brier_vl     = if (!is.null(brier_vl_col)) .wagg(get(brier_vl_col), w_vl) else NA_real_,
        auc_vl       = if (!is.null(auc_vl_col))   .wagg(get(auc_vl_col),   w_vl) else NA_real_,
        rmse_vl      = if (!is.null(rmse_vl_col))  .wagg(get(rmse_vl_col),  w_vl) else NA_real_,
        cor_vl       = if (!is.null(cor_vl_col))   .wagg(get(cor_vl_col),   w_vl) else NA_real_
      )
    }, by = t]
  }

  qexit_binary <- "qexit_type" %in% names(bc) &&
    any(bc$qexit_type == "binomial", na.rm = TRUE)

  parts <- list(
    .agg_branch(bc, "g_remain",
      n_tr_col    = "n_tr_ar",      n_vl_col    = "n_vl_ar",
      tgt_tr_col  = "target_grem_tr", pred_tr_col = "pred_grem_tr",
      tgt_vl_col  = "target_grem_vl", pred_vl_col = "pred_grem_vl",
      slope_vl_col = "grem_vl_cal_slope",
      brier_vl_col = "grem_vl_brier", auc_vl_col  = "grem_vl_auc"),

    .agg_branch(bc, "g_death",
      n_tr_col    = "n_tr_exit",      n_vl_col    = "n_vl_exit",
      tgt_tr_col  = "target_gdeath_tr", pred_tr_col = "pred_gdeath_tr",
      tgt_vl_col  = "target_gdeath_vl", pred_vl_col = "pred_gdeath_vl",
      slope_vl_col = "gdeath_vl_cal_slope",
      brier_vl_col = "gdeath_vl_brier", auc_vl_col  = "gdeath_vl_auc"),

    .agg_branch(bc, "q_rem",
      n_tr_col    = "n_tr_rem",      n_vl_col    = "n_vl_rem",
      tgt_tr_col  = "target_qrem_tr", pred_tr_col = "pred_qrem_tr",
      tgt_vl_col  = "target_qrem_vl", pred_vl_col = "pred_qrem_vl",
      slope_vl_col = "qrem_vl_cal_slope",
      rmse_vl_col  = "qrem_vl_rmse", cor_vl_col  = "qrem_vl_cor"),

    .agg_branch(bc, "q_exit",
      n_tr_col    = "n_tr_exit",      n_vl_col    = "n_vl_exit",
      tgt_tr_col  = "target_qexit_tr", pred_tr_col = "pred_qexit_tr",
      tgt_vl_col  = "target_qexit_vl", pred_vl_col = "pred_qexit_vl",
      slope_vl_col = "qexit_vl_cal_slope",
      brier_vl_col = if (qexit_binary)  "qexit_vl_brier" else NULL,
      auc_vl_col   = if (qexit_binary)  "qexit_vl_auc"   else NULL,
      rmse_vl_col  = if (!qexit_binary) "qexit_vl_rmse"  else NULL,
      cor_vl_col   = if (!qexit_binary) "qexit_vl_cor"   else NULL)
  )

  out <- data.table::rbindlist(parts, use.names = TRUE, fill = TRUE)
  out[, branch := factor(branch, levels = c("g_remain", "g_death", "q_rem", "q_exit"))]
  data.table::setorderv(out, c("t", "branch"))
  class(out) <- c("branch_cal_summary", class(out))
  out
}

#' @export
print.branch_cal_summary <- function(x, digits = 3L, ...) {
  branches  <- c("g_remain", "g_death", "q_rem", "q_exit")
  metrics   <- c("n_vl", "tgt_vl", "pred_vl", "cal_slope_vl")
  col_labels <- c("n_vl", "tgt", "pred", "slope")

  w_n  <- 6L
  w_f  <- digits + 4L
  col_widths <- c(w_n, w_f, w_f, w_f)
  branch_sep <- " | "
  sep_w      <- nchar(branch_sep)

  branch_w <- sum(col_widths) + (length(col_widths) - 1L)

  t_vals  <- sort(unique(x$t))
  t_label <- "t"
  t_w     <- max(nchar(as.character(t_vals)), nchar(t_label))

  fmt_cell <- function(v, w, is_int = FALSE) {
    if (is.na(v)) return(formatC("\u2014", width = w, flag = "-"))
    if (is_int)   return(formatC(as.integer(v), width = w, format = "d"))
    formatC(round(v, digits), width = w, format = "f", digits = digits)
  }

  branch_spans <- vapply(branches, function(b) {
    span <- branch_w + sep_w * (length(col_widths) - 1L)
    pad  <- branch_w - nchar(b)
    lpad <- floor(pad / 2)
    rpad <- pad - lpad
    paste0(strrep("-", lpad), b, strrep("-", rpad))
  }, character(1L))

  t_pad <- strrep(" ", t_w)
  cat("\nBranch calibration summary  (fold-weighted)\n\n")
  cat(t_pad, " ", paste(branch_spans, collapse = branch_sep), "\n", sep = "")

  sub_heads <- vapply(seq_along(col_labels), function(i)
    formatC(col_labels[i], width = col_widths[i], flag = "-"), character(1L))
  sub_block <- paste(sub_heads, collapse = " ")
  cat(formatC(t_label, width = t_w), " ",
      paste(rep(sub_block, length(branches)), collapse = branch_sep),
      "\n", sep = "")

  x_dt <- data.table::as.data.table(x)
  for (tt in t_vals) {
    row_parts <- vapply(branches, function(b) {
      r <- x_dt[t == tt & branch == b]
      if (nrow(r) == 0L) {
        return(paste(rep(formatC("\u2014", width = max(col_widths)), length(metrics)),
                     collapse = " "))
      }
      cells <- c(
        fmt_cell(r$n_vl,         col_widths[1L], is_int = TRUE),
        fmt_cell(r$tgt_vl,       col_widths[2L]),
        fmt_cell(r$pred_vl,      col_widths[3L]),
        fmt_cell(r$cal_slope_vl, col_widths[4L])
      )
      paste(cells, collapse = " ")
    }, character(1L))

    cat(formatC(tt, width = t_w, format = "d"), " ",
        paste(row_parts, collapse = branch_sep),
        "\n", sep = "")
  }

  invisible(x)
}
