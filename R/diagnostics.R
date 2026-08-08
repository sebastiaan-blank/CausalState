#################################################################
### Create a table with density ratio diagnostics 
#############################################################
ess_weight <- function(w) {
  w <- w[is.finite(w) & !is.na(w)]
  
  if (length(w) == 0 || sum(w^2) == 0) {
    return(NA_real_)
  }
  
  sum(w)^2 / sum(w^2)
}


diagnose_density_weights <- function(
    weight_object,
    trim = 0.99,
    digits = 3
) {
  if (!"weights_dt" %in% names(weight_object)) {
    stop("`weight_object` must contain a `weights_dt` element.")
  }
  
  weights <- weight_object$weights_dt
  
  if (!is.data.frame(weights)) {
    stop("`weight_object$weights_dt` must be a data.frame or data.table.")
  }
  
  if (ncol(weights) < 3) {
    stop("`weights_dt` must contain at least id, time, and Rt_t columns.")
  }
  
  id   <- names(weights)[1]
  time <- names(weights)[2]
  
  rt_col <- if ("Rt_t" %in% names(weights)) {
    "Rt_t"
  } else {
    names(weights)[3]
  }
  
  d <- dplyr::as_tibble(weights) |>
    dplyr::transmute(
      id          = .data[[id]],
      trial_block = .data[[time]],
      Rt_raw      = .data[[rt_col]]
    )
  
  if (any(!is.finite(d$Rt_raw) | d$Rt_raw < 0, na.rm = TRUE)) {
    stop("Incremental density ratios contain non-finite or negative values.")
  }
  
  if (is.null(trim) || is.na(trim)) {
    trim_threshold <- NA_real_
    
    d <- d |>
      dplyr::mutate(
        Rt_t = Rt_raw,
        was_trimmed = FALSE
      )
  } else {
    if (!is.numeric(trim) || length(trim) != 1 || trim <= 0 || trim >= 1) {
      stop("`trim` must be NULL/NA or a single number between 0 and 1.")
    }
    
    trim_threshold <- as.numeric(
      stats::quantile(d$Rt_raw, probs = trim, na.rm = TRUE)
    )
    
    d <- d |>
      dplyr::mutate(
        Rt_t = pmin(Rt_raw, trim_threshold),
        was_trimmed = Rt_raw > trim_threshold
      )
  }
  
  d <- d |>
    dplyr::arrange(id, trial_block) |>
    dplyr::group_by(id) |>
    dplyr::mutate(Rt_cum = cumprod(Rt_t)) |>
    dplyr::ungroup()
  
  d |>
    dplyr::group_by(trial_block) |>
    dplyr::summarise(
      at_risk = dplyr::n(),
      
      mean_inc    = mean(Rt_t, na.rm = TRUE),
      median_inc  = stats::median(Rt_t, na.rm = TRUE),
      ESS_inc     = ess_weight(Rt_t),
      ESS_pct_inc = ESS_inc / at_risk * 100,
      Q95_inc     = as.numeric(stats::quantile(Rt_t, 0.95, na.rm = TRUE)),
      Max_inc     = max(Rt_t, na.rm = TRUE),
      
      mean_cum    = mean(Rt_cum, na.rm = TRUE),
      median_cum  = stats::median(Rt_cum, na.rm = TRUE),
      ESS_cum     = ess_weight(Rt_cum),
      ESS_pct_cum = ESS_cum / at_risk * 100,
      Q95_cum     = as.numeric(stats::quantile(Rt_cum, 0.95, na.rm = TRUE)),
      Max_cum     = max(Rt_cum, na.rm = TRUE),
      
      n_trimmed   = sum(was_trimmed, na.rm = TRUE),
      pct_trimmed = n_trimmed / at_risk * 100,
      trim_threshold = trim_threshold,
      
      .groups = "drop"
    ) |>
    dplyr::mutate(
      dplyr::across(
        where(is.numeric) &
          !dplyr::all_of(c("trial_block", "at_risk", "n_trimmed")),
        ~ round(.x, digits)
      )
    )
}


############################################
# Branch diagnostics 
#######################################

branch_calibration <- function(fit, digits = 3) {
  bd <- fit$diagnostics$branch_calibration
  if (is.null(bd)) bd <- fit$diagnostics$branch_cal
  if (is.null(bd)) bd <- fit$diagnostics$branch_diag
  if (is.null(bd)) stop("No branch calibration diagnostics found.")
  
  d <- dplyr::as_tibble(bd)
  
  wmean <- function(x, w) {
    ok <- is.finite(x) & is.finite(w) & !is.na(x) & !is.na(w) & w > 0
    if (!any(ok)) return(NA_real_)
    stats::weighted.mean(x[ok], w[ok])
  }
  
  d |>
    dplyr::group_by(t) |>
    dplyr::summarise(
      grem_target_tr = wmean(target_grem_tr, n_tr_ar),
      grem_pred_tr   = wmean(pred_grem_tr,   n_tr_ar),
      grem_error_tr  = grem_target_tr - grem_pred_tr,
      grem_brier_tr  = wmean(grem_tr_brier, n_tr_ar),
      grem_auc_tr    = wmean(grem_tr_auc, n_tr_ar),
      
      grem_target_vl = wmean(target_grem_vl, n_vl_ar),
      grem_pred_vl   = wmean(pred_grem_vl,   n_vl_ar),
      grem_error_vl  = grem_target_vl - grem_pred_vl,
      grem_brier_vl  = wmean(grem_vl_brier, n_vl_ar),
      grem_auc_vl    = wmean(grem_vl_auc, n_vl_ar),
      
      gdeath_target_tr = wmean(target_gdeath_tr, n_tr_exit),
      gdeath_pred_tr   = wmean(pred_gdeath_tr,   n_tr_exit),
      gdeath_error_tr  = gdeath_target_tr - gdeath_pred_tr,
      gdeath_brier_tr  = wmean(gdeath_tr_brier, n_tr_exit),
      gdeath_auc_tr    = wmean(gdeath_tr_auc, n_tr_exit),
      
      gdeath_target_vl = wmean(target_gdeath_vl, n_vl_exit),
      gdeath_pred_vl   = wmean(pred_gdeath_vl,   n_vl_exit),
      gdeath_error_vl  = gdeath_target_vl - gdeath_pred_vl,
      gdeath_brier_vl  = wmean(gdeath_vl_brier, n_vl_exit),
      gdeath_auc_vl    = wmean(gdeath_vl_auc, n_vl_exit),
      
      qrem_target_tr = wmean(target_qrem_tr, n_tr_rem),
      qrem_pred_tr   = wmean(pred_qrem_tr,   n_tr_rem),
      qrem_error_tr  = qrem_target_tr - qrem_pred_tr,
      qrem_mse_tr    = wmean(qrem_tr_mse, n_tr_rem),
      qrem_rmse_tr   = sqrt(qrem_mse_tr),
      qrem_mae_tr    = wmean(qrem_tr_mae, n_tr_rem),
      qrem_cor_tr    = wmean(qrem_tr_cor, n_tr_rem),
      
      qrem_target_vl = wmean(target_qrem_vl, n_vl_rem),
      qrem_pred_vl   = wmean(pred_qrem_vl,   n_vl_rem),
      qrem_error_vl  = qrem_target_vl - qrem_pred_vl,
      qrem_mse_vl    = wmean(qrem_vl_mse, n_vl_rem),
      qrem_rmse_vl   = sqrt(qrem_mse_vl),
      qrem_mae_vl    = wmean(qrem_vl_mae, n_vl_rem),
      qrem_cor_vl    = wmean(qrem_vl_cor, n_vl_rem),
      
      qexit_target_tr = wmean(target_qexit_tr, n_tr_exit),
      qexit_pred_tr   = wmean(pred_qexit_tr,   n_tr_exit),
      qexit_error_tr  = qexit_target_tr - qexit_pred_tr,
      qexit_brier_tr  = wmean(qexit_tr_brier, n_tr_exit),
      qexit_auc_tr    = wmean(qexit_tr_auc, n_tr_exit),
      
      qexit_target_vl = wmean(target_qexit_vl, n_vl_exit),
      qexit_pred_vl   = wmean(pred_qexit_vl,   n_vl_exit),
      qexit_error_vl  = qexit_target_vl - qexit_pred_vl,
      qexit_brier_vl  = wmean(qexit_vl_brier, n_vl_exit),
      qexit_auc_vl    = wmean(qexit_vl_auc, n_vl_exit),
      
      .groups = "drop"
    ) |>
    dplyr::mutate(
      grem_brier_ratio   = grem_brier_vl / grem_brier_tr,
      gdeath_brier_ratio = gdeath_brier_vl / gdeath_brier_tr,
      qrem_rmse_ratio    = qrem_rmse_vl / qrem_rmse_tr,
      qexit_brier_ratio  = qexit_brier_vl / qexit_brier_tr,
      
      grem_error_gap   = abs(grem_error_vl)   - abs(grem_error_tr),
      gdeath_error_gap = abs(gdeath_error_vl) - abs(gdeath_error_tr),
      qrem_error_gap   = abs(qrem_error_vl)   - abs(qrem_error_tr),
      qexit_error_gap  = abs(qexit_error_vl)  - abs(qexit_error_tr),
      
      grem_auc_drop   = grem_auc_tr - grem_auc_vl,
      gdeath_auc_drop = gdeath_auc_tr - gdeath_auc_vl,
      qexit_auc_drop  = qexit_auc_tr - qexit_auc_vl,
      qrem_cor_drop   = qrem_cor_tr - qrem_cor_vl
    ) |>
    dplyr::mutate(
      dplyr::across(where(is.numeric), ~ round(.x, digits))
    ) |>
    dplyr::arrange(t)
}

################################################
## recursion diagnostics 
#############################################


recursion_diagnostics <- function(fit, digits = 3) {
  d <- fit$diagnostics$recursion_diag
  if (is.null(d)) d <- fit$diagnostics$fold_diag
  if (is.null(d)) stop("No recursion diagnostics found.")
  
  d <- dplyr::as_tibble(d)
  
  wmean <- function(x, w) {
    ok <- is.finite(x) & is.finite(w) & !is.na(x) & !is.na(w) & w > 0
    if (!any(ok)) return(NA_real_)
    stats::weighted.mean(x[ok], w[ok])
  }
  
  d |>
    dplyr::group_by(t) |>
    dplyr::summarise(
      n_train_total = sum(n_train_at_risk, na.rm = TRUE),
      n_valid_total = sum(n_valid_at_risk, na.rm = TRUE),
      
      Q_nat_pre_mean = wmean(Q_nat_pre_mean, n_train_at_risk),
      Q_shf_pre_mean = wmean(Q_shf_pre_mean, n_train_at_risk),
      Q_pre_diff_mean = wmean(Q_pre_diff_mean, n_train_at_risk),
      Q_pre_diff_q95_abs = wmean(Q_pre_diff_q95_abs, n_train_at_risk),
      
      pseudo_pre_mean  = wmean(pseudo_pre_mean, n_train_at_risk),
      pseudo_post_mean = wmean(pseudo_post_mean, n_train_at_risk),
      
      pseudo_pre_min  = min(pseudo_pre_min, na.rm = TRUE),
      pseudo_pre_max  = max(pseudo_pre_max, na.rm = TRUE),
      pseudo_post_min = min(pseudo_post_min, na.rm = TRUE),
      pseudo_post_max = max(pseudo_post_max, na.rm = TRUE),
      
      delta_mean    = wmean(delta_mean, n_train_at_risk),
      delta_sd      = wmean(delta_sd, n_train_at_risk),
      delta_q95_abs = wmean(delta_q95_abs, n_train_at_risk),
      delta_max_abs = max(delta_max_abs, na.rm = TRUE),
      
      corr_mean    = wmean(corr_mean, n_train_at_risk),
      corr_sd      = wmean(corr_sd, n_train_at_risk),
      corr_q95_abs = wmean(corr_q95_abs, n_train_at_risk),
      corr_max_abs = max(corr_max_abs, na.rm = TRUE),
      
      n_post_below_0 = sum(n_post_below_0, na.rm = TRUE),
      n_post_above_1 = sum(n_post_above_1, na.rm = TRUE),
      
      .groups = "drop"
    ) |>
    dplyr::mutate(
      prop_post_below_0 = n_post_below_0 / n_train_total,
      prop_post_above_1 = n_post_above_1 / n_train_total
    ) |>
    dplyr::mutate(
      dplyr::across(where(is.numeric), ~ round(.x, digits))
    ) |>
    dplyr::arrange(t)
}


# ------------------------------------------------------------------
# Targeting SL diagnostics (future development)
# ------------------------------------------------------------------

# overview of targeting SL 
plot_itmle_targeting_sl <- function(fit,
                                    top_n = 8L,
                                    min_mean_weight = 0.01) {
  stopifnot(requireNamespace("dplyr", quietly = TRUE))
  stopifnot(requireNamespace("ggplot2", quietly = TRUE))
  stopifnot(requireNamespace("stringr", quietly = TRUE))
  stopifnot(requireNamespace("scales", quietly = TRUE))
  
  target_sl <- fit$diagnostics$target_sl
  
  df <- target_sl %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(
      tt = dplyr::coalesce(.data$tt, .data$t),
      learner = as.character(.data$learner),
      weight = as.numeric(.data$weight)
    )
  
  w <- df %>%
    dplyr::group_by(tt, learner) %>%
    dplyr::summarise(
      mean_weight = mean(weight, na.rm = TRUE),
      .groups = "drop"
    )
  
  keep_learners <- w %>%
    dplyr::group_by(learner) %>%
    dplyr::summarise(overall = mean(mean_weight, na.rm = TRUE), .groups = "drop") %>%
    dplyr::filter(overall >= min_mean_weight) %>%
    dplyr::slice_max(overall, n = top_n, with_ties = FALSE) %>%
    dplyr::pull(learner)
  
  w_plot <- w %>%
    dplyr::mutate(
      learner_plot = dplyr::if_else(learner %in% keep_learners, learner, "Other")
    ) %>%
    dplyr::group_by(tt, learner_plot) %>%
    dplyr::summarise(mean_weight = sum(mean_weight, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(
      learner_plot = stringr::str_replace_all(learner_plot, "^SL\\.", ""),
      learner_plot = stringr::str_replace_all(learner_plot, "_All$", "")
    )
  
  ggplot2::ggplot(
    w_plot,
    ggplot2::aes(
      x = tt,
      y = mean_weight,
      colour = learner_plot,
      linetype = learner_plot,
      shape = learner_plot,
      group = learner_plot
    )
  ) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 2.4) +
    ggplot2::scale_x_reverse(breaks = sort(unique(w_plot$tt))) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      labels = scales::percent_format(accuracy = 1)
    ) +
    ggplot2::labs(
      title = "iTMLE targeting SuperLearner weights",
      subtitle = "Mean metalearner weights across folds and targeting steps",
      x = "Time point",
      y = "Mean SL weight",
      colour = "Learner",
      linetype = "Learner",
      shape = "Learner"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank()
    )
}


## plots per s 

plot_itmle_targeting_sl_by_s <- function(fit,
                                         top_n = 6L,
                                         min_mean_weight = 0.01,
                                         ncol = 4) {
  stopifnot(requireNamespace("dplyr", quietly = TRUE))
  stopifnot(requireNamespace("ggplot2", quietly = TRUE))
  stopifnot(requireNamespace("stringr", quietly = TRUE))
  stopifnot(requireNamespace("scales", quietly = TRUE))
  
  df <- fit$diagnostics$target_sl %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(
      s = as.integer(.data$s),
      t = as.integer(.data$t),
      learner = as.character(.data$learner),
      weight = as.numeric(.data$weight)
    )
  
  w <- df %>%
    dplyr::group_by(s, t, learner) %>%
    dplyr::summarise(
      mean_weight = mean(weight, na.rm = TRUE),
      .groups = "drop"
    )
  
  keep_learners <- w %>%
    dplyr::group_by(learner) %>%
    dplyr::summarise(
      overall = mean(mean_weight, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(overall >= min_mean_weight) %>%
    dplyr::slice_max(overall, n = top_n, with_ties = FALSE) %>%
    dplyr::pull(learner)
  
  w_plot <- w %>%
    dplyr::mutate(
      learner_plot = dplyr::if_else(learner %in% keep_learners, learner, "Other"),
      learner_plot = stringr::str_replace_all(learner_plot, "^SL\\.", ""),
      learner_plot = stringr::str_replace_all(learner_plot, "_All$", "")
    ) %>%
    dplyr::group_by(s, t, learner_plot) %>%
    dplyr::summarise(
      mean_weight = sum(mean_weight, na.rm = TRUE),
      .groups = "drop"
    )
  
  ggplot2::ggplot(
    w_plot,
    ggplot2::aes(
      x = t,
      y = mean_weight,
      colour = learner_plot,
      linetype = learner_plot,
      shape = learner_plot,
      group = learner_plot
    )
  ) +
    ggplot2::geom_line(linewidth = 0.75) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::facet_wrap(~ s, ncol = ncol) +
    ggplot2::scale_x_reverse(breaks = sort(unique(w_plot$t))) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      labels = scales::percent_format(accuracy = 1)
    ) +
    ggplot2::labs(
      title = "iTMLE targeting SuperLearner weights",
      subtitle = "Each panel is outer s; x-axis is targeting t",
      x = "Targeting t",
      y = "Mean SL weight",
      colour = "Learner",
      linetype = "Learner",
      shape = "Learner"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank()
    )
}



#############################
## Branch diagnostics 
##########################################
branch_check <- function(fit, digits = 3) {
  bd <- fit$diagnostics$branch_cal
  if (is.null(bd)) stop("No `fit$diagnostics$branch_cal` found.")
  
  bd |>
    dplyr::as_tibble() |>
    dplyr::group_by(t) |>
    dplyr::summarise(
      n_ar_tr = sum(n_tr_ar, na.rm = TRUE),
      n_ar_vl = sum(n_vl_ar, na.rm = TRUE),
      
      Q_rem_tr_on_rem   = weighted.mean(Q_rem_tr_on_rem, n_rem_tr, na.rm = TRUE),
      Q_rem_vl_on_rem   = weighted.mean(Q_rem_vl_on_rem, n_rem_vl, na.rm = TRUE),
      Q_rem_tr_all      = weighted.mean(Q_rem_tr_all, n_tr_ar, na.rm = TRUE),
      Q_rem_vl_all      = weighted.mean(Q_rem_vl_all, n_vl_ar, na.rm = TRUE),
      
      Q_exit_tr_on_exit = weighted.mean(Q_exit_tr_on_exit, n_exit_tr, na.rm = TRUE),
      Q_exit_vl_on_exit = weighted.mean(Q_exit_vl_on_exit, n_exit_vl, na.rm = TRUE),
      Q_exit_tr_all     = weighted.mean(Q_exit_tr_all, n_tr_ar, na.rm = TRUE),
      Q_exit_vl_all     = weighted.mean(Q_exit_vl_all, n_vl_ar, na.rm = TRUE),
      
      prem_Qrem_tr      = weighted.mean(prem_Qrem_tr, n_tr_ar, na.rm = TRUE),
      prem_Qrem_vl      = weighted.mean(prem_Qrem_vl, n_vl_ar, na.rm = TRUE),
      pexit_Qexit_tr    = weighted.mean(pexit_Qexit_tr, n_tr_ar, na.rm = TRUE),
      pexit_Qexit_vl    = weighted.mean(pexit_Qexit_vl, n_vl_ar, na.rm = TRUE),
      
      Q_mix_tr          = weighted.mean(Q_mix_tr, n_tr_ar, na.rm = TRUE),
      Q_mix_vl          = weighted.mean(Q_mix_vl, n_vl_ar, na.rm = TRUE),
      
      Y_rem_tr = weighted.mean(Y_rem_tr, n_rem_tr, na.rm = TRUE),
      Y_exit_tr = weighted.mean(Y_exit_tr, n_exit_tr, na.rm = TRUE),
      Y_train_ar = weighted.mean(Y_train_ar, n_tr_ar, na.rm = TRUE),
      Y_valid_next = weighted.mean(Y_valid_next, n_vl_ar, na.rm = TRUE),
      
      Q_rem_vl_vs_Y_rem_tr = Q_rem_vl_on_rem - Y_rem_tr,
      Q_exit_vl_vs_Y_exit_tr = Q_exit_vl_on_exit - Y_exit_tr,
      Q_mix_vl_vs_Y_valid_next = Q_mix_vl - Y_valid_next,
      
      .groups = "drop"
    ) |>
    dplyr::mutate(
      Q_rem_on_rem_diff   = Q_rem_vl_on_rem - Q_rem_tr_on_rem,
      Q_rem_all_diff      = Q_rem_vl_all - Q_rem_tr_all,
      Q_exit_on_exit_diff = Q_exit_vl_on_exit - Q_exit_tr_on_exit,
      Q_exit_all_diff     = Q_exit_vl_all - Q_exit_tr_all,
      prem_Qrem_diff      = prem_Qrem_vl - prem_Qrem_tr,
      pexit_Qexit_diff    = pexit_Qexit_vl - pexit_Qexit_tr,
      mix_diff            = Q_mix_vl - Q_mix_tr
    ) |>
    dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, digits))) |>
    dplyr::arrange(t)
}

#########################################################
# recusrsion diagnostic 
##########################################################

recursion_check <- function(fit, digits = 3) {
  fit$diagnostics$recursion_diag |>
    dplyr::as_tibble() |>
    dplyr::group_by(t) |>
    dplyr::summarise(
      n_train_at_risk_rep = sum(n_train_at_risk, na.rm = TRUE),
      
      Q_shf_pre_mean    = weighted.mean(Q_shf_pre_mean, n_train_at_risk, na.rm = TRUE),
      pseudo_pre_mean   = weighted.mean(pseudo_pre_mean, n_train_at_risk, na.rm = TRUE),
      pseudo_post_mean  = weighted.mean(pseudo_post_mean, n_train_at_risk, na.rm = TRUE),
      delta_mean        = weighted.mean(delta_mean, n_train_at_risk, na.rm = TRUE),
      corr_mean         = weighted.mean(corr_mean, n_train_at_risk, na.rm = TRUE),
      
      n_below_0 = sum(n_post_below_0, na.rm = TRUE),
      n_above_1 = sum(n_post_above_1, na.rm = TRUE),
      
      pct_below_0 = n_below_0 / n_train_at_risk_rep,
      pct_above_1 = n_above_1 / n_train_at_risk_rep,
      
      .groups = "drop"
    ) |>
    dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, digits))) |>
    dplyr::arrange(t)
}





