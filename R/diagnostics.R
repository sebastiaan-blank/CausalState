# diagnostics.R
#
# Exported diagnostic helpers for density_ratio() / sdr() / itmle() fits.
#
# weight_diagnostics()       -- per-time weight summary for a density_ratio object
# branch_cal_summary()       -- fold-weighted calibration per t x branch
# print.branch_cal_summary   -- S3 print method

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
