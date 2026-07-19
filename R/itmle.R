itmle_make_cols <- function(t, baseline, tv_names, a_names, all_names) {
  cols <- intersect(baseline, all_names)
  
  if (t <= 0L) {
    return(unique(cols))
  }
  
  for (u in seq_len(t)) {
    cols <- c(cols, intersect(paste0(tv_names, "_t", u), all_names))
  }
  
  if (t > 1L) {
    for (u in seq_len(t - 1L)) {
      cols <- c(cols, intersect(paste0(a_names, "_t", u), all_names))
    }
  }
  
  unique(cols)
}


itmle_target_diag_one <- function(Y, Q0, Q1, w, s, t, bounds = 1e-5) {
  ok <- is.finite(Y) & is.finite(Q0) & is.finite(Q1) &
    is.finite(w) & w >= 0
  
  Y  <- Y[ok]
  Q0 <- Q0[ok]
  Q1 <- Q1[ok]
  w  <- w[ok]
  
  if (!length(Y) || sum(w) <= 0 || sum(w^2) <= 0) {
    return(data.table::data.table(
      s = s, t = t, n = length(Y),
      ess_r = NA_real_, w_q99 = NA_real_, w_max = NA_real_,
      sc0 = NA_real_, sc1 = NA_real_, sc_r = NA_real_,
      d_q95 = NA_real_, d_max = NA_real_,
      b_lo = NA_real_, b_hi = NA_real_
    ))
  }
  
  sw  <- sum(w)
  ess <- sw^2 / sum(w^2)
  
  sc0 <- sum(w * (Y - Q0)) / sw
  sc1 <- sum(w * (Y - Q1)) / sw
  d   <- Q1 - Q0
  
  data.table::data.table(
    s = s,
    t = t,
    n = length(Y),
    
    q0 = mean(Q0),
    q1 = mean(Q1),
    d_mean = mean(d),
    
    ess_r = ess / length(Y),
    w_q99 = as.numeric(stats::quantile(w, 0.99, na.rm = TRUE)),
    w_max = max(w, na.rm = TRUE),
    
    sc0  = sc0,
    sc1  = sc1,
    sc_r = abs(sc1) / pmax(abs(sc0), 1e-12),
    
    d_q95 = as.numeric(stats::quantile(abs(d), 0.95, na.rm = TRUE)),
    d_max = max(abs(d), na.rm = TRUE),
    
    b_lo = mean(Q1 <= bounds * 10),
    b_hi = mean(Q1 >= 1 - bounds * 10)
  )
}


itmle_weight <- function(density_ratios, t, s) {
  if (t == s) {
    return(rep(1, nrow(density_ratios)))
  }
  u_idx     <- seq(t + 1L, s)
  ratio_mat <- density_ratios[, u_idx, drop = FALSE]
  apply(ratio_mat, 1L, prod)
}

itmle_update <- function(Q, eps, bounds = 1e-5, clip_fun = NULL) {
  if (is.null(clip_fun)) clip_fun <- function(x) clip(x, bounds)
  out <- plogis(qlogis(clip_fun(Q)) + eps)
  clip_fun(out)
}

itmle_target_fit <- function(
    Y_train, X_train,
    X_pred_tr, X_pred_vl,
    offset_train, offset_pred_tr, offset_pred_vl,
    weights,
    sl_lib,
    v_target            = 10L,
    v_sl_inner          = 10L,
    parallel            = FALSE,
    reg_workers         = NULL,
    sl_workers          = NULL,
    cluster_inner       = NULL,
    seed                = 1L,
    min_n_train_fold    = 20L,
    min_obs_with_weight = 30L,
    clip_fun            = function(x) pmin(pmax(x, 1e-5), 1 - 1e-5)
) {
  N_tr <- length(Y_train)
  
  ok <- weights > 0 & is.finite(weights) &
    is.finite(Y_train) & is.finite(offset_train)
  
  if (sum(ok) < min_obs_with_weight) {
    Q_tr <- clip_fun(plogis(offset_pred_tr))
    Q_vl <- clip_fun(plogis(offset_pred_vl))
    
    return(list(
      eps_tr  = rep(0, nrow(X_pred_tr)),
      eps_vl  = rep(0, nrow(X_pred_vl)),
      Q_tr    = Q_tr,
      Q_vl    = Q_vl,
      sl_meta = data.table::data.table(
        learner = "empty_fallback_too_few_weighted",
        weight  = 1,
        fold    = NA_integer_
      ),
      
      mean_Q_tr_before = mean(Q_tr, na.rm = TRUE),
      mean_Q_tr_after  = mean(Q_tr, na.rm = TRUE),
      mean_eps_tr      = 0,
      sd_eps_tr        = 0,
      
      mean_Q_vl_before = mean(Q_vl, na.rm = TRUE),
      mean_Q_vl_after  = mean(Q_vl, na.rm = TRUE),
      mean_eps_vl      = 0,
      sd_eps_vl        = 0
    ))
  }
  
  X_train_aug   <- cbind(as.data.frame(X_train),   "._sl_offset" = offset_train)
  X_pred_tr_aug <- cbind(as.data.frame(X_pred_tr), "._sl_offset" = offset_pred_tr)
  X_pred_vl_aug <- cbind(as.data.frame(X_pred_vl), "._sl_offset" = offset_pred_vl)
  
  use_cluster <- !is.null(cluster_inner) &&
    length(unique(cluster_inner)) >= v_target
  
  outer_folds <- with_seed(seed, {
    if (use_cluster) {
      origami::make_folds(
        n           = N_tr,
        V           = v_target,
        cluster_ids = cluster_inner,
        fold_fun    = origami::folds_vfold
      )
    } else {
      origami::make_folds(
        n        = N_tr,
        V        = v_target,
        fold_fun = origami::folds_vfold
      )
    }
  })
  
  eps_tr_cv <- numeric(N_tr)
  Q_tr_cv   <- clip_fun(plogis(offset_pred_tr))
  
  eps_vl_mat <- matrix(
    NA_real_,
    nrow = nrow(X_pred_vl),
    ncol = length(outer_folds)
  )
  
  ok_idx <- which(ok)

  use_mc_tgt <- !is.null(sl_workers)
  if (use_mc_tgt) {
    old_cores_tgt <- getOption("mc.cores")
    options(mc.cores = as.integer(sl_workers))
    on.exit(options(mc.cores = old_cores_tgt), add = TRUE)
  }
  sl_fn_tgt <- if (use_mc_tgt) SuperLearner::mcSuperLearner else SuperLearner::SuperLearner

  fold_res <- par_lapply(
    seq_along(outer_folds),
    function(v) {
      
      tr_idx <- outer_folds[[v]]$training_set
      vl_idx <- outer_folds[[v]]$validation_set
      
      tr_v <- intersect(tr_idx, ok_idx)
      vl_v <- intersect(vl_idx, ok_idx)
      
      if (length(tr_v) < min_n_train_fold || length(vl_v) < 1L) {
        return(list(
          vl_idx = vl_idx,
          eps_vl = rep(0, length(vl_idx)),
          Q_vl   = clip_fun(plogis(offset_pred_tr[vl_idx])),
          eps_vl_full = rep(NA_real_, nrow(X_pred_vl)),
          sl_meta = data.table::data.table(
            learner = "fold_skipped_too_small",
            weight  = 1,
            fold    = v
          )
        ))
      }
      
      inner_use_cluster <- !is.null(cluster_inner) &&
        length(unique(cluster_inner[tr_v])) >= v_sl_inner
      
      inner_folds <- with_seed(seed + 1000L * v, {
        if (inner_use_cluster) {
          origami::make_folds(
            n           = length(tr_v),
            V           = v_sl_inner,
            cluster_ids = cluster_inner[tr_v],
            fold_fun    = origami::folds_vfold
          )
        } else {
          origami::make_folds(
            n        = length(tr_v),
            V        = v_sl_inner,
            fold_fun = origami::folds_vfold
          )
        }
      })
      
      sl_validRows <- lapply(inner_folds, function(f) f$validation_set)
      
      fit_v <- tryCatch(
        with_seed(
          seed + 2000L * v,
          suppressWarnings(
            sl_fn_tgt(
              Y          = clip_fun(Y_train[tr_v]),
              X          = X_train_aug[tr_v, , drop = FALSE],
              SL.library = sl_lib,
              family     = stats::binomial(),
              obsWeights = weights[tr_v],
              id         = if (!is.null(cluster_inner)) {
                cluster_inner[tr_v]
              } else {
                seq_along(tr_v)
              },
              cvControl  = list(
                V = v_sl_inner,
                validRows = sl_validRows
              )
            )
          )
        ),
        error = function(e) {
          message(sprintf(
            "[itmle_target_fit fold %d] SL error: %s",
            v, conditionMessage(e)
          ))
          NULL
        }
      )
      
      if (is.null(fit_v)) {
        return(list(
          vl_idx = vl_idx,
          eps_vl = rep(0, length(vl_idx)),
          Q_vl   = clip_fun(plogis(offset_pred_tr[vl_idx])),
          eps_vl_full = rep(NA_real_, nrow(X_pred_vl)),
          sl_meta = data.table::data.table(
            learner = "fold_sl_error",
            weight  = 1,
            fold    = v
          )
        ))
      }
      
      p_tr_v <- as.numeric(
        stats::predict(
          fit_v,
          X_pred_tr_aug[vl_v, , drop = FALSE]
        )$pred
      )
      
      p_tr_v <- clip_fun(p_tr_v)
      
      eps_v <- qlogis(p_tr_v) - offset_pred_tr[vl_v]
      eps_v[!is.finite(eps_v)] <- 0
      
      p_vl_v <- as.numeric(
        stats::predict(
          fit_v,
          X_pred_vl_aug
        )$pred
      )
      
      p_vl_v <- clip_fun(p_vl_v)
      
      eps_vl_full <- qlogis(p_vl_v) - offset_pred_vl
      eps_vl_full[!is.finite(eps_vl_full)] <- 0
      
      list(
        vl_idx = vl_v,
        eps_vl = eps_v,
        Q_vl   = p_tr_v,
        eps_vl_full = eps_vl_full,
        sl_meta = data.table::data.table(
          learner = fit_v$libraryNames,
          weight  = fit_v$coef,
          fold    = v
        )
      )
    },
    workers  = reg_workers,
    parallel = parallel,
    seed     = TRUE
  )
  
  sl_meta_all <- vector("list", length(fold_res))
  
  for (v in seq_along(fold_res)) {
    rr <- fold_res[[v]]
    
    eps_tr_cv[rr$vl_idx] <- rr$eps_vl
    Q_tr_cv[rr$vl_idx]   <- rr$Q_vl
    
    eps_vl_mat[, v] <- rr$eps_vl_full
    sl_meta_all[[v]] <- rr$sl_meta
  }
  
  eps_vl_cv <- rowMeans(eps_vl_mat, na.rm = TRUE)
  
  bad_vl <- !is.finite(eps_vl_cv)
  if (any(bad_vl)) {
    eps_vl_cv[bad_vl] <- 0
  }
  
  Q_vl_cv <- clip_fun(plogis(offset_pred_vl + eps_vl_cv))
  
  list(
    eps_tr  = eps_tr_cv,
    eps_vl  = eps_vl_cv,
    Q_tr    = Q_tr_cv,
    Q_vl    = Q_vl_cv,
    sl_meta = data.table::rbindlist(sl_meta_all, use.names = TRUE, fill = TRUE), 
    
    mean_Q_tr_before = mean(clip_fun(plogis(offset_pred_tr)), na.rm = TRUE),
    mean_Q_tr_after  = mean(Q_tr_cv, na.rm = TRUE),
    mean_eps_tr      = mean(eps_tr_cv, na.rm = TRUE),
    sd_eps_tr        = sd(eps_tr_cv, na.rm = TRUE),
    
    mean_Q_vl_before = mean(clip_fun(plogis(offset_pred_vl)), na.rm = TRUE),
    mean_Q_vl_after  = mean(Q_vl_cv, na.rm = TRUE),
    mean_eps_vl      = mean(eps_vl_cv, na.rm = TRUE),
    sd_eps_vl        = sd(eps_vl_cv, na.rm = TRUE)
  )
}

itmle_inner_target <- function(
    Q_mix_nat_tr, Q_mix_shf_tr,
    Q_mix_nat_vl, Q_mix_shf_vl,
    Y_pseudo_tr,
    D_wide,
    id_tr_ar, id_vl_ar,
    density_ratios_tr, density_ratios_vl,
    s, t0, K,
    ids,
    sl_lib_target,
    cluster_tr   = NULL,
    v_target     = 10L,
    v_sl_inner   = 10L,
    parallel     = FALSE,
    reg_workers  = NULL,
    sl_workers   = NULL,
    seed         = 1L,
    baseline, tv_names, a_names,
    bounds       = 1e-5,
    clip_fun     = NULL
) {
  if (is.null(clip_fun)) {
    clip_fun <- function(x) clip(x, bounds)
  }
  
  Q_run_nat_tr <- Q_mix_nat_tr
  Q_run_shf_tr <- Q_mix_shf_tr
  Q_run_nat_vl <- Q_mix_nat_vl
  Q_run_shf_vl <- Q_mix_shf_vl
  
  n_steps <- s - t0 + 1L
  
  eps_log     <- vector("list", n_steps)
  target_diag <- vector("list", n_steps)
  sl_diag     <- vector("list", n_steps)
  
  row_pos_tr <- match(ids[id_tr_ar], D_wide$.__id)
  row_pos_vl <- match(ids[id_vl_ar], D_wide$.__id)
  
  if (anyNA(row_pos_tr) || anyNA(row_pos_vl)) {
    stop("D_wide row match failed: ids[id_tr_ar] or ids[id_vl_ar] not found in D_wide$.__id.")
  }
  
  

  for (t in seq(s, t0, by = -1L)) {

    idx <- s - t + 1L

    if (t == 1L) {
      weight_tr <- itmle_weight(density_ratios_tr, t = 0L, s = s)

      off_tr <- qlogis(clip_fun(Q_run_nat_tr))
      eps_4a <- 0
      ok_w   <- is.finite(weight_tr) & weight_tr > 0
      if (sum(ok_w) >= 2L) {
        fit_4a <- tryCatch(
          stats::glm.fit(
            x       = matrix(1L, length(Y_pseudo_tr), 1L),
            y       = clip_fun(Y_pseudo_tr),
            offset  = off_tr,
            weights = weight_tr,
            family  = stats::quasibinomial()
          ),
          error = function(e) NULL
        )
        if (!is.null(fit_4a) && is.finite(fit_4a$coefficients[1L])) {
          eps_4a <- fit_4a$coefficients[1L]
        }
      }

      Q0    <- Q_run_nat_tr
      Q_vl0 <- Q_run_nat_vl
      Q_run_nat_tr <- itmle_update(Q_run_nat_tr, eps_4a, bounds, clip_fun)
      Q_run_shf_tr <- itmle_update(Q_run_shf_tr, eps_4a, bounds, clip_fun)
      Q_run_nat_vl <- itmle_update(Q_run_nat_vl, eps_4a, bounds, clip_fun)
      Q_run_shf_vl <- itmle_update(Q_run_shf_vl, eps_4a, bounds, clip_fun)

      target_diag[[idx]] <- itmle_target_diag_one(
        Y = Y_pseudo_tr, Q0 = Q0, Q1 = Q_run_nat_tr,
        w = weight_tr, s = s, t = 1L, bounds = bounds
      )
      target_diag[[idx]][, `:=`(
        mean_Q_tr_before = mean(clip_fun(plogis(off_tr))),
        mean_Q_tr_after  = mean(Q_run_nat_tr),
        mean_eps_tr      = eps_4a,
        sd_eps_tr        = 0,
        mean_Q_vl_before = mean(clip_fun(Q_vl0)),
        mean_Q_vl_after  = mean(Q_run_nat_vl),
        mean_eps_vl      = eps_4a,
        sd_eps_vl        = 0
      )]

      eps_log[[idx]] <- list(
        s      = s,
        t      = 1L,
        eps_tr = rep(eps_4a, length(Q_run_nat_tr)),
        eps_vl = rep(eps_4a, length(Q_run_nat_vl)),
        Q_tr   = Q_run_nat_tr,
        Q_vl   = Q_run_nat_vl
      )

    } else {
      cols_t <- itmle_make_cols(
        t         = t - 1L,
        baseline  = baseline,
        tv_names  = tv_names,
        a_names   = a_names,
        all_names = names(D_wide)
      )

      H_t_tr <- D_wide[row_pos_tr, ..cols_t]
      H_t_vl <- D_wide[row_pos_vl, ..cols_t]

      weight_tr <- itmle_weight(density_ratios_tr, t = t - 1L, s = s)

      fit_pred <- itmle_target_fit(
        Y_train        = Y_pseudo_tr,
        X_train        = H_t_tr,
        X_pred_tr      = H_t_tr,
        X_pred_vl      = H_t_vl,

        offset_train   = qlogis(clip_fun(Q_run_nat_tr)),
        offset_pred_tr = qlogis(clip_fun(Q_run_nat_tr)),
        offset_pred_vl = qlogis(clip_fun(Q_run_nat_vl)),

        weights        = weight_tr,
        sl_lib         = sl_lib_target,
        v_target       = v_target,
        v_sl_inner     = v_sl_inner,
        reg_workers    = reg_workers,
        parallel       = parallel,
        sl_workers     = sl_workers,
        cluster_inner  = cluster_tr,
        seed           = seed + t,
        clip_fun       = clip_fun
      )

      Q0 <- Q_run_nat_tr
      Q_run_nat_tr <- itmle_update(Q_run_nat_tr, fit_pred$eps_tr, bounds, clip_fun)
      Q_run_shf_tr <- itmle_update(Q_run_shf_tr, fit_pred$eps_tr, bounds, clip_fun)
      Q_run_nat_vl <- itmle_update(Q_run_nat_vl, fit_pred$eps_vl, bounds, clip_fun)
      Q_run_shf_vl <- itmle_update(Q_run_shf_vl, fit_pred$eps_vl, bounds, clip_fun)

      target_diag[[idx]] <- itmle_target_diag_one(
        Y      = Y_pseudo_tr,
        Q0     = Q0,
        Q1     = Q_run_nat_tr,
        w      = weight_tr,
        s      = s,
        t      = t,
        bounds = bounds
      )

      target_diag[[idx]][
        , `:=`(
          mean_Q_tr_before = fit_pred$mean_Q_tr_before,
          mean_Q_tr_after  = fit_pred$mean_Q_tr_after,
          mean_eps_tr      = fit_pred$mean_eps_tr,
          sd_eps_tr        = fit_pred$sd_eps_tr,

          mean_Q_vl_before = fit_pred$mean_Q_vl_before,
          mean_Q_vl_after  = fit_pred$mean_Q_vl_after,
          mean_eps_vl      = fit_pred$mean_eps_vl,
          sd_eps_vl        = fit_pred$sd_eps_vl
        )
      ]

      if (!is.null(fit_pred$sl_meta)) {
        sl_diag[[idx]] <- data.table::as.data.table(fit_pred$sl_meta)[
          , `:=`(s = s, t = t)
        ]
      }

      eps_log[[idx]] <- list(
        s      = s,
        t      = t,
        eps_tr = fit_pred$eps_tr,
        eps_vl = fit_pred$eps_vl,
        Q_tr   = fit_pred$Q_tr,
        Q_vl   = fit_pred$Q_vl
      )
    }
  }
  
  list(
    Q_star_nat_tr = Q_run_nat_tr,
    Q_star_shf_tr = Q_run_shf_tr,
    Q_star_nat_vl = Q_run_nat_vl,
    Q_star_shf_vl = Q_run_shf_vl,
    eps_log       = eps_log,
    target_diag   = data.table::rbindlist(target_diag, use.names = TRUE, fill = TRUE),
    sl_diag       = data.table::rbindlist(Filter(Negate(is.null), sl_diag), use.names = TRUE, fill = TRUE)
  )
}

itmle_eic_se <- function(ic, cluster_by_id = NULL, use_second_moment = TRUE) {
  ic <- as.numeric(ic)
  ok_ic <- is.finite(ic)
  
  if (is.null(cluster_by_id)) {
    n_eff <- sum(ok_ic)
    if (n_eff < 2L) {
      stop("SE: <2 finite IC values.", call. = FALSE)
    }
    
    ic_ok <- ic[ok_ic]
    
    se <- if (isTRUE(use_second_moment)) {
      sqrt(mean(ic_ok^2) / n_eff)
    } else {
      stats::sd(ic_ok) / sqrt(n_eff)
    }
    
    return(list(
      se = se,
      type = if (isTRUE(use_second_moment)) {
        "iid_eic_second_moment"
      } else {
        "iid_centered"
      },
      n_eff = n_eff,
      n_clusters = NA_integer_,
      df = n_eff - 1L,
      eic_mean = mean(ic_ok),
      eic_second_moment = mean(ic_ok^2),
      se_centered = stats::sd(ic_ok) / sqrt(n_eff),
      se_second_moment = sqrt(mean(ic_ok^2) / n_eff)
    ))
  }
  
  cl <- cluster_by_id
  if (length(cl) != length(ic)) {
    stop("SE: cluster_by_id must have same length as ic.", call. = FALSE)
  }
  
  ok <- ok_ic & !is.na(cl)
  
  dt <- data.table::data.table(
    cl = cl[ok],
    ic = ic[ok]
  )
  
  n_eff <- nrow(dt)
  if (n_eff < 2L) {
    stop("SE: <2 finite IC values after clustering filter.", call. = FALSE)
  }
  
  S <- dt[, .(S = sum(ic), n = .N), by = cl]
  G <- nrow(S)
  
  if (G < 2L) {
    stop("SE: need >=2 clusters for cluster-robust SE.", call. = FALSE)
  }
  
  Sbar <- mean(S$S)
  
  se_second_moment <- sqrt(
    (G / (G - 1)) * sum(S$S^2) / (n_eff^2)
  )
  
  se_centered <- sqrt(
    (G / (G - 1)) * sum((S$S - Sbar)^2) / (n_eff^2)
  )
  
  se <- if (isTRUE(use_second_moment)) {
    se_second_moment
  } else {
    se_centered
  }
  
  list(
    se = se,
    type = if (isTRUE(use_second_moment)) {
      "cluster_eic_second_moment"
    } else {
      "cluster_centered"
    },
    n_eff = n_eff,
    n_clusters = G,
    df = G - 1L,
    eic_mean = mean(dt$ic),
    eic_second_moment = mean(dt$ic^2),
    mean_cluster_score = Sbar,
    max_abs_cluster_score = max(abs(S$S)),
    se_centered = se_centered,
    se_second_moment = se_second_moment#,
  )
}



#' Infinite-dimensional TMLE (iTMLE) estimator for longitudinal MTPs
#'
#' Estimates the mean counterfactual outcome under a modified treatment policy
#' (MTP) using the infinite-dimensional targeted minimum loss-based estimator
#' (iTMLE) of Luedtke et al. (2017), Section 5 / Algorithm 4.  The companion
#' [sdr()] function implements the LMTP-SDR estimator of Díaz et al. (2021),
#' which applies the Luedtke et al. SDR construction to the density-ratio / MTP
#' setting.  The two estimators share the same backward Q-regression; they
#' differ only in the update step: iTMLE applies an infinite-dimensional TMLE
#' fluctuation rather than the EIF-based pseudo-outcome update used by [sdr()].
#' Both are sequentially doubly robust (2^\eqn{K}-robust, Luedtke et al.
#' Definition 2): consistent whenever, at each time point t, either the
#' treatment model \eqn{g_t} or the outcome model \eqn{Q_t} is consistently
#' estimated.  The targeting step uses a SuperLearner ensemble with custom
#' learner wrappers (see [sl_itmle]) that handle the logit-offset structure
#' required for cross-validated TMLE.  Requires pre-computed density ratio
#' weights from [density_ratio()].
#'
#' @inheritParams sdr
#' @param sl_tmle SuperLearner library for the iTMLE targeting step. Should
#'   consist of learners from [sl_itmle] that accept an offset column
#'   (`._sl_offset`). Defaults to the package-provided `sl_tmle` vector.
#' @param v_target_itmle Integer. Number of cross-validation folds for the
#'   outer targeting loop. Default `10L`.
#' @param v_sl_inner_itmle Integer. Number of inner CV folds inside the
#'   targeting SuperLearner. Default `10L`.
#'
#' @return A named list with:
#'   \describe{
#'     \item{`psi`}{Point estimate of `E[Y(d)]` under the MTP (after targeting).}
#'     \item{`se`}{Standard error from the efficient influence curve.}
#'     \item{`ci`}{95\% Wald confidence interval.}
#'     \item{`ic`}{Per-subject influence curve values (length n).}
#'     \item{`psi_onestep`}{One-step pre-targeting estimate (diagnostic).}
#'     \item{`targeting_gap`}{Mean EIF after targeting; near zero indicates convergence.}
#'     \item{`sl_summary`}{`data.table` of SuperLearner weights per fold/time/component.}
#'     \item{`fold_diag`}{Per-fold diagnostic summaries.}
#'     \item{`target_diag`}{Per-iteration targeting diagnostics.}
#'     \item{`fluct_diag`}{Fluctuation model diagnostics from the targeting step.}
#'   }
#'
#' @references
#' Luedtke AR, Sofrygin O, van der Laan MJ, Carone M (2017). Sequential
#' Double Robustness in Right-Censored Longitudinal Models. arXiv:1705.02459.
#'
#' Díaz I, Williams N, Hoffman KL, Schenck EJ (2021). Nonparametric Causal
#' Effects Based on Longitudinal Modified Treatment Policies. *JASA*
#' 118(542):846–857.
#'
#' @seealso [density_ratio()], [sdr()], [absorb_rule()], [sl_itmle]
#'
#' @export
itmle <- function(
    df, weight_object, tmax,
    id, time, alive, in_state, y,
    baseline,
    tv_names = character(0),
    a_names = character(0),
    no_lag_vars = character(0),
    policy_names = character(0),
    sl_remain = NULL,
    sl_death  = NULL,
    sl_recursive = NULL,
    sl_rec_simple = NULL,
    rec_transition = NULL,
    sl_y = NULL,
    outcome_family = c("binomial", "gaussian"),
    y_bounds = NULL,
    bounds = 1e-5,
    trim = 0.99,
    absorb = list(),
    policy_spec_fun = function(D_block, t, a_names) NULL,
    k = 2,
    seed = 1,
    parallel     = FALSE,
    fold_workers = NULL,
    reg_workers  = NULL,
    sl_workers   = NULL,
    inner_v = 5L,
    sl_tmle = NULL,
    v_target_itmle = 10L,
    v_sl_inner_itmle = 10L,
    cluster = NULL,
    cluster_se_only = FALSE,
    pool_g_death = FALSE
) {
  stopifnot(requireNamespace("data.table", quietly = TRUE))
  stopifnot(requireNamespace("SuperLearner", quietly = TRUE))
  if (is.null(sl_tmle)) sl_tmle <- get("sl_tmle", envir = asNamespace("CausalState"))

  if (!isTRUE(parallel)) sl_workers <- NULL

  weights_dt <- extract_weights_dt(weight_object)
  
  keep_cols <- unique(c(
    id, time, a_names, y, alive, in_state,
    baseline, tv_names, 
    no_lag_vars, policy_names, cluster
  ))
  keep_cols <- intersect(keep_cols, names(df))
  df <- as.data.frame(df)[, keep_cols, drop = FALSE]
  
  
  data_check(data = df, id = id, time = time,
             a_names = a_names, tv_names = tv_names,
             bs = baseline, y_col = y, alive_col = alive,
             in_state_col = in_state, tmin = 1, tmax = tmax)
  
  set.seed(seed)
  
  DT <- data.table::as.data.table(df)
  data.table::setorderv(DT, c(id, time))
  rm(df); gc()
  t_min <- DT[, .SD[1L], by = id][[time]] |> max()
  
  if (is.null(sl_remain)) stop("`sl_remain` must be provided (SL library for g_remain).", call. = FALSE)
  if (is.null(sl_death))  stop("`sl_death` must be provided (SL library for g_death_exit).", call. = FALSE)
  if (is.null(sl_y)) stop("`sl_y` must be provided (SL library for Q_exit).", call. = FALSE)
  if (is.null(sl_recursive))  stop("`sl_recursive` must be provided (SL library for Q_rem).", call. = FALSE)
  if (!is.null(sl_rec_simple) && is.null(rec_transition)) {
    stop("`sl_rec_simple` requires `rec_transition` to be specified.", call. = FALSE)
  }

  outcome_family <- match.arg(outcome_family, c("binomial", "gaussian"))

  prep <- prepare_long(DT, id, time, y, cluster = cluster)
  D      <- prep$D

  cols_to_pivot <- intersect(c(tv_names, a_names), names(D))
  
  D_wide <- data.table::dcast(
    D[.__time >= 1L],
    formula   = .__id ~ .__time,
    value.var = cols_to_pivot,
    sep       = "_t"
  )
  
  D_baseline <- unique(D[, .SD, .SDcols = c(".__id", baseline)], by = ".__id")
  D_wide <- D_baseline[D_wide, on = ".__id"]
  
  data.table::setorderv(D_wide, ".__id")
  
  hist_out <- prepare_history_lags(
    DT            = D,
    a_names       = a_names,
    tv_names       = tv_names,
    k             = k,
    id        = ".__id",
    no_lag_vars   = no_lag_vars
  )
  D <- hist_out$DT
  
  scale_info <- scale_y(D, y, outcome_family, y_bounds, bounds)
  is_binom   <- scale_info$is_binom
  
  Y_info <- build_Y_by_id(D, scale_info)
  ids    <- Y_info$ids
  Y_init <- Y_info$Y
  N      <- length(ids)
  
  row_index <- build_id_time_index(D, ids, tmax)
  
  if (!is.null(cluster)) {
    cl_by_id <- D[, .(cl = data.table::first(na.omit(.__cl))), by = .__id]
    cluster_by_id <- cl_by_id$cl[match(ids, cl_by_id$.__id)]
  } else {
    cluster_by_id <- ids
  }
  
  w_tmax <- max(as.integer(weights_dt[[time]]), na.rm = TRUE)
  if (w_tmax < tmax) {
    stop(sprintf(
      "weight_object covers time points up to %d but tmax = %d. Re-run density_ratio() with tmax >= %d.",
      w_tmax, tmax, tmax
    ), call. = FALSE)
  }
  density_ratios <- build_density_ratios(
    weights_dt, id, time, ids, tmax,
    row_index = row_index,
    trim = trim
  )
  
  
  fold_info <- get_folds_from_weights(
    weights_dt      = weights_dt,
    id              = id,
    ids             = ids,
    cluster_by_id   = cluster_by_id,
    cluster_se_only = cluster_se_only
  )
  folds <- fold_info$folds
  
  all_names <- names(D)
  cols_by_t <- lapply(seq_len(tmax), function(tt) {
    make_cols(
      tt        = tt,
      baseline  = baseline,
      tv_names  = tv_names,
      a_names   = a_names,
      k         = k,
      all_names = all_names, 
      t_min = t_min
    )
  })
  cols_by_t_Q <- cols_by_t
  
  D_shifted <- vector("list", length = tmax)
  if (length(a_names)) {
    for (tt in seq_len(tmax)) {
      D_shifted[[tt]] <- overwrite_policy_history_for_Q(
        D,
        t               = tt,
        a_names         = a_names,
        policy_spec_fun = policy_spec_fun,
        id              = id,
        time            = time,
        k               = k, 
        t_min           = t_min 
      )
    }
  } else {
    for (tt in seq_len(tmax)) {
      D_shifted[[tt]] <- D
    }
  }
  
  pred_nat_all <- matrix(NA_real_, nrow = N, ncol = tmax + 1L)
  pred_shf_all <- matrix(NA_real_, nrow = N, ncol = tmax + 1L)
  
  pred_nat_all[, tmax + 1L] <- Y_init
  pred_shf_all[, tmax + 1L] <- Y_init
  
  sl_chunks <- list()
  
  fold_worker <- function(f_idx) {
    
    is_binom <- scale_info$is_binom
    
    fold   <- folds[[f_idx]]
    tr_ids <- fold$training_set
    vl_ids <- fold$validation_set
    
    Y_train <- Y_init[tr_ids]
    
    n_tr <- length(tr_ids)
    n_vl <- length(vl_ids)
    
    nat_train <- matrix(NA_real_, nrow = n_tr, ncol = tmax + 1L)
    shf_train <- matrix(NA_real_, nrow = n_tr, ncol = tmax + 1L)
    nat_valid <- matrix(NA_real_, nrow = n_vl, ncol = tmax + 1L)
    shf_valid <- matrix(NA_real_, nrow = n_vl, ncol = tmax + 1L)
    
    nat_train[, tmax + 1L] <- Y_train
    shf_train[, tmax + 1L] <- Y_train
    nat_valid[, tmax + 1L] <- Y_init[vl_ids]
    shf_valid[, tmax + 1L] <- Y_init[vl_ids]
    
    
    sl_chunks_here  <- list()
    fold_diag_here  <- list()
    target_diag <- list()
    target_sl <- list()
    branch_diag_here <- list()
    
    eic_tplus_valid <- numeric(length(vl_ids))
    
    pooled_dex <- NULL
    
    if (isTRUE(pool_g_death)) {
      pooled_dex <- fit_pooled_g_death_exit(
        D              = D,
        tr_ids         = tr_ids,
        row_index      = row_index,
        tmax           = tmax,
        baseline       = baseline,
        tv_names       = tv_names,
        a_names        = a_names,
        all_names      = all_names,
        k              = k,
        t_min          = t_min,
        alive          = alive,
        in_state       = in_state,
        cluster_by_id  = cluster_by_id,
        sl_death       = sl_death,
        inner_v        = inner_v,
        seed           = seed,
        f_idx          = f_idx,
        sl_workers     = sl_workers
      )
      
      if (is.null(pooled_dex$fit)) {
        message(sprintf(
          "[fold %d] pool_g_death=TRUE but pooled g_death_exit fit failed (%s); falling back to per-time fits",
          f_idx, pooled_dex$reason
        ))
        fit_dex_pooled <- NULL
        sl_dex_pooled  <- NULL
        cn_pool        <- NULL
        cols_pool      <- NULL
        all_Y_pool     <- NULL
        pool_g_death   <- FALSE
      } else {
        fit_dex_pooled <- pooled_dex$fit
        sl_dex_pooled  <- pooled_dex$sl
        cn_pool        <- pooled_dex$cols
        cols_pool      <- pooled_dex$cols_pool
        all_Y_pool     <- pooled_dex$Y
      }
    } else {
      fit_dex_pooled <- NULL
      sl_dex_pooled  <- NULL
      cn_pool        <- NULL
      cols_pool      <- NULL
      all_Y_pool     <- NULL
    }
    
    
    for (tt in rev(seq_len(tmax))) {
      
      message(sprintf("[fold %d][t=%d][pid=%d] %s", f_idx, tt, Sys.getpid(), format(Sys.time(), "%H:%M:%S")))
      
      
      sl_dex <- sl_rem <- sl_qexit <- sl_qrem <- NULL
      fit_dex <- fit_rem <- fit_qexit <- fit_qrem <- NULL
      cn_qexit <- cn_qrem <- NULL
      
      at_risk_tr <- !is.na(row_index[tr_ids, tt])
      at_risk_vl <- !is.na(row_index[vl_ids, tt])
      
      if (!any(at_risk_tr)) {
        nat_train[, tt] <- nat_train[, tt + 1L]
        shf_train[, tt] <- shf_train[, tt + 1L]
        nat_valid[, tt] <- nat_valid[, tt + 1L]
        shf_valid[, tt] <- shf_valid[, tt + 1L]
        next
      }
      
      idx_tr_ar <- which(at_risk_tr)
      idx_vl_ar <- which(at_risk_vl)
      
      id_tr_ar <- tr_ids[idx_tr_ar]
      id_vl_ar <- vl_ids[idx_vl_ar]
      
      rows_tr   <- row_index[id_tr_ar, tt]
      cols_base <- cols_by_t_Q[[tt]]
      
      X_nat_tr_base <- make_design(D, rows_tr, cols_base)
      X_shf_tr_base <- patch_shifted_design(
        X_nat        = X_nat_tr_base,
        D_shifted_tt = D_shifted[[tt]],
        rows         = rows_tr,
        a_names      = a_names,
        tt           = tt,
        k            = k, 
        t_min         = t_min
      )
      
      X_nat_tr_haz <- X_nat_tr_base
      X_shf_tr_haz <- X_shf_tr_base
      
      alive_end <- as.integer(D[[alive]][rows_tr])
      in_state_end   <- as.integer(D[[in_state]][rows_tr])
      
      D_tr <- as.integer(alive_end == 0L)
      R_tr <- as.integer(alive_end == 1L & in_state_end == 1L)
      C_tr <- as.integer(R_tr == 0L & D_tr == 0L)
      
      exit_idx <- which(R_tr == 0L)
      rem_idx  <- which(R_tr == 1L)
      
      cl_tr_ar   <- cluster_by_id[id_tr_ar]
      cl_exit_ar <- cluster_by_id[id_tr_ar[exit_idx]]
      cl_rem_ar  <- cluster_by_id[id_tr_ar[rem_idx]]
      
      Y_pseudo_ar <- Y_train[at_risk_tr]
      n_ar        <- length(idx_tr_ar)
      
      
      reg_res <- fit_transition_regressions(
        X_nat_tr_base = X_nat_tr_base,
        X_shf_tr_base = X_shf_tr_base,
        X_nat_tr_haz = X_nat_tr_haz,
        X_shf_tr_haz = X_shf_tr_haz,
        Y_pseudo_ar = Y_pseudo_ar,
        Y_init_ar = Y_init[id_tr_ar],
        R_tr = R_tr,
        D_tr = D_tr,
        exit_idx = exit_idx,
        rem_idx = rem_idx,
        n_ar = n_ar,
        cl_tr_ar = cl_tr_ar,
        cl_exit_ar = cl_exit_ar,
        cl_rem_ar = cl_rem_ar,
        pool_g_death = pool_g_death,
        fit_dex_pooled = fit_dex_pooled,
        cn_pool = cn_pool,
        cols_pool = cols_pool,
        D = D,
        rows_tr = rows_tr,
        D_shifted_tt = D_shifted[[tt]],
        a_names = a_names,
        tt = tt,
        k = k,
        t_min = t_min,
        is_binom = is_binom,
        tmax = tmax,
        scale_info = scale_info,
        sl_remain = sl_remain,
        sl_death = sl_death,
        sl_y = sl_y,
        sl_recursive = sl_recursive,
        sl_rec_simple = sl_rec_simple,
        rec_transition = rec_transition,
        inner_v = inner_v,
        seed = seed,
        f_idx = f_idx,
        parallel = parallel,
        reg_workers = reg_workers,
        sl_workers = sl_workers
      )

      r_rem  <- reg_res$g_remain
      r_dex  <- reg_res$g_death_exit
      r_exit <- reg_res$Q_exit
      r_qrem <- reg_res$Q_rem

      sl_rem      <- r_rem$sl_rem
      fit_rem     <- r_rem$fit_rem
      p_rem_const <- r_rem$p_rem_const
      p_rem_nat   <- r_rem$p_rem_nat
      p_rem_shf   <- r_rem$p_rem_shf

      sl_dex      <- r_dex$sl_dex
      fit_dex     <- r_dex$fit_dex
      p_dex_const <- r_dex$p_dex_const
      p_dex_nat   <- r_dex$p_dex_nat
      p_dex_shf   <- r_dex$p_dex_shf
      
      sl_qexit    <- r_exit$sl_qexit
      fit_qexit   <- r_exit$fit_qexit
      cn_qexit    <- r_exit$cn_qexit
      muY         <- r_exit$muY
      q_death_nat <- r_exit$q_death_nat
      q_dc_nat    <- r_exit$q_dc_nat
      q_death_shf <- r_exit$q_death_shf
      q_dc_shf    <- r_exit$q_dc_shf
      
      sl_qrem     <- r_qrem$sl_qrem
      fit_qrem    <- r_qrem$fit_qrem
      cn_qrem     <- r_qrem$cn_qrem
      muP         <- r_qrem$muP
      q_rem_nat   <- r_qrem$q_rem_nat
      q_rem_shf   <- r_qrem$q_rem_shf
      
      rm(reg_res, r_rem, r_dex, r_exit, r_qrem)
      gc()
      
      q_death_nat_raw <- q_death_nat
      q_dc_nat_raw    <- q_dc_nat
      q_death_shf_raw <- q_death_shf
      q_dc_shf_raw    <- q_dc_shf
      
      
      if (length(absorb)) {
        D_block_tr <- D[rows_tr]
        
        q_death_nat <- apply_absorb_branch(q_death_nat, D_block_tr, tt, "death", scale_info, absorb)
        q_death_shf <- apply_absorb_branch(q_death_shf, D_block_tr, tt, "death", scale_info, absorb)
        q_dc_nat    <- apply_absorb_branch(q_dc_nat,    D_block_tr, tt, "dc",    scale_info, absorb)
        q_dc_shf    <- apply_absorb_branch(q_dc_shf,    D_block_tr, tt, "dc",    scale_info, absorb)
      }
      
      Q_nat_ar <- p_rem_nat * q_rem_nat +
        (1 - p_rem_nat) * (p_dex_nat * q_death_nat + (1 - p_dex_nat) * q_dc_nat)
      
      Q_shf_ar <- p_rem_shf * q_rem_shf +
        (1 - p_rem_shf) * (p_dex_shf * q_death_shf + (1 - p_dex_shf) * q_dc_shf)
      
      nat_train[at_risk_tr, tt] <- Q_nat_ar
      shf_train[at_risk_tr, tt] <- Q_shf_ar
      
      n_val_tt <- sum(at_risk_vl)
      
      meta <- sl_meta(sl_rem,   f_idx, tt, "g_remain",     n_ar,             n_val_tt)
      if (!is.null(meta)) sl_chunks_here[[length(sl_chunks_here) + 1L]] <- meta
      
      if (!is.null(sl_dex)) {
        meta <- sl_meta(sl_dex, f_idx, tt, "g_death_exit", length(exit_idx), n_val_tt)
      } else if (isTRUE(pool_g_death) && !is.null(fit_dex_pooled) && tt == tmax) {
        meta <- sl_meta(sl_dex_pooled, f_idx, 0L, "g_death_exit_pooled",
                            length(all_Y_pool), n_val_tt)
      } else {
        meta <- NULL
      }
      if (!is.null(meta)) sl_chunks_here[[length(sl_chunks_here) + 1L]] <- meta  
      
      meta <- sl_meta(sl_qexit, f_idx, tt, "Q_exit",       length(exit_idx), n_val_tt)
      if (!is.null(meta)) sl_chunks_here[[length(sl_chunks_here) + 1L]] <- meta
      
      meta <- sl_meta(sl_qrem,  f_idx, tt, "Q_rem",        length(rem_idx),  n_val_tt)
      if (!is.null(meta)) sl_chunks_here[[length(sl_chunks_here) + 1L]] <- meta
      
      if (any(at_risk_vl)) {
        
        rows_vl <- row_index[id_vl_ar, tt]
        
        X_nat_vl_base <- make_design(D, rows_vl, cols_base)
        X_shf_vl_base <- patch_shifted_design(
          X_nat        = X_nat_vl_base,
          D_shifted_tt = D_shifted[[tt]],
          rows         = rows_vl,
          a_names      = a_names,
          tt           = tt,
          k            = k, 
          t_min         = t_min 
        )
        
        cn_base       <- colnames(X_nat_tr_base)
        X_nat_vl_base <- align_cols(X_nat_vl_base, cn_base)
        X_shf_vl_base <- align_cols(X_shf_vl_base, cn_base)
        
        X_nat_vl_haz <- X_nat_vl_base
        X_shf_vl_haz <- X_shf_vl_base
        
        if (!is.null(fit_rem)) {
          p_rem_nat_vl <- scale_info$clip(sl_predict(fit_rem, X_nat_vl_haz))
          p_rem_shf_vl <- scale_info$clip(sl_predict(fit_rem, X_shf_vl_haz))
        } else {
          p_rem_nat_vl <- rep(p_rem_const, nrow(X_nat_vl_haz))
          p_rem_shf_vl <- rep(p_rem_const, nrow(X_nat_vl_haz))
        }
        
        if (!is.null(fit_dex)) {
          if (isTRUE(pool_g_death)) {
            X_nat_dex_vl      <- make_design(D, rows_vl, cols_pool)
            X_nat_dex_vl$.__t <- tt
            X_shf_dex_vl <- patch_shifted_design(
              X_nat        = X_nat_dex_vl,
              D_shifted_tt = D_shifted[[tt]],
              rows         = rows_vl,
              a_names      = a_names,
              tt           = tt,
              k            = k,
              t_min        = t_min
            )
            
            X_shf_dex_vl$.__t <- tt
            X_nat_dex_vl <- align_cols(X_nat_dex_vl, cn_pool)
            X_shf_dex_vl <- align_cols(X_shf_dex_vl, cn_pool)
            p_dex_nat_vl <- scale_info$clip(sl_predict(fit_dex, X_nat_dex_vl))
            p_dex_shf_vl <- scale_info$clip(sl_predict(fit_dex, X_shf_dex_vl))
          } else {
            p_dex_nat_vl <- scale_info$clip(sl_predict(fit_dex, X_nat_vl_haz))
            p_dex_shf_vl <- scale_info$clip(sl_predict(fit_dex, X_shf_vl_haz))
          }
        } else {
          p_dex_nat_vl <- rep(p_dex_const, nrow(X_nat_vl_haz))
          p_dex_shf_vl <- rep(p_dex_const, nrow(X_nat_vl_haz))
        }
        
        if (!is.null(fit_qexit)) {
          X_nat_d_vl <- X_nat_vl_base; X_nat_d_vl$exit_status <- 1
          X_nat_c_vl <- X_nat_vl_base; X_nat_c_vl$exit_status <- 0
          X_shf_d_vl <- X_shf_vl_base; X_shf_d_vl$exit_status <- 1
          X_shf_c_vl <- X_shf_vl_base; X_shf_c_vl$exit_status <- 0
          
          X_nat_d_vl <- align_cols(X_nat_d_vl, cn_qexit)
          X_nat_c_vl <- align_cols(X_nat_c_vl, cn_qexit)
          X_shf_d_vl <- align_cols(X_shf_d_vl, cn_qexit)
          X_shf_c_vl <- align_cols(X_shf_c_vl, cn_qexit)
          
          q_death_nat_vl <- scale_info$clip(sl_predict(fit_qexit, X_nat_d_vl))
          q_dc_nat_vl    <- scale_info$clip(sl_predict(fit_qexit, X_nat_c_vl))
          q_death_shf_vl <- scale_info$clip(sl_predict(fit_qexit, X_shf_d_vl))
          q_dc_shf_vl    <- scale_info$clip(sl_predict(fit_qexit, X_shf_c_vl))
        } else {
          q_death_nat_vl <- rep(muY, nrow(X_nat_vl_base))
          q_dc_nat_vl    <- rep(muY, nrow(X_nat_vl_base))
          q_death_shf_vl <- rep(muY, nrow(X_nat_vl_base))
          q_dc_shf_vl    <- rep(muY, nrow(X_nat_vl_base))
        }
        
        if (!is.null(fit_qrem)) {
          X_nat_r_vl <- align_cols(X_nat_vl_base, cn_qrem)
          X_shf_r_vl <- align_cols(X_shf_vl_base, cn_qrem)
          
          q_rem_nat_vl <- scale_info$clip(sl_predict(fit_qrem, X_nat_r_vl))
          q_rem_shf_vl <- scale_info$clip(sl_predict(fit_qrem, X_shf_r_vl))
        } else {
          q_rem_nat_vl <- rep(muP, nrow(X_nat_vl_base))
          q_rem_shf_vl <- rep(muP, nrow(X_nat_vl_base))
        }
        
        
        q_death_nat_vl_raw <- q_death_nat_vl
        q_dc_nat_vl_raw    <- q_dc_nat_vl
        q_death_shf_vl_raw <- q_death_shf_vl
        q_dc_shf_vl_raw    <- q_dc_shf_vl
        
        if (length(absorb)) {
          D_block_vl <- D[rows_vl]
          
          q_death_nat_vl <- apply_absorb_branch(q_death_nat_vl, D_block_vl, tt, "death", scale_info, absorb)
          q_death_shf_vl <- apply_absorb_branch(q_death_shf_vl, D_block_vl, tt, "death", scale_info, absorb)
          q_dc_nat_vl    <- apply_absorb_branch(q_dc_nat_vl,    D_block_vl, tt, "dc",    scale_info, absorb)
          q_dc_shf_vl    <- apply_absorb_branch(q_dc_shf_vl,    D_block_vl, tt, "dc",    scale_info, absorb)
        }
        
        
        Q_nat_vl <- p_rem_nat_vl * q_rem_nat_vl +
          (1 - p_rem_nat_vl) * (p_dex_nat_vl * q_death_nat_vl + (1 - p_dex_nat_vl) * q_dc_nat_vl)
        
        Q_shf_vl <- p_rem_shf_vl * q_rem_shf_vl +
          (1 - p_rem_shf_vl) * (p_dex_shf_vl * q_death_shf_vl + (1 - p_dex_shf_vl) * q_dc_shf_vl)
        
        nat_valid[idx_vl_ar, tt] <- Q_nat_vl
        shf_valid[idx_vl_ar, tt] <- Q_shf_vl
      }
      
      if (any(!at_risk_tr)) {
        nat_train[!at_risk_tr, tt] <- nat_train[!at_risk_tr, tt + 1L]
        shf_train[!at_risk_tr, tt] <- shf_train[!at_risk_tr, tt + 1L]
      }
      if (any(!at_risk_vl)) {
        nat_valid[!at_risk_vl, tt] <- nat_valid[!at_risk_vl, tt + 1L]
        shf_valid[!at_risk_vl, tt] <- shf_valid[!at_risk_vl, tt + 1L]
      }
      
      
      D_vl <- R_vl <- integer(0)
      rem_idx_vl <- exit_idx_vl <- integer(0)
      Y_next_vl <- numeric(0)
      q_exit_obs_tr_raw <- ifelse(D_tr == 1L, q_death_nat_raw, q_dc_nat_raw)
      q_exit_obs_vl_raw <- numeric(0)
      
      if (any(at_risk_vl)) {
        alive_end_vl <- as.integer(D[[alive]][rows_vl])
        icu_end_vl   <- as.integer(D[[in_state]][rows_vl])
        
        D_vl <- as.integer(alive_end_vl == 0L)
        R_vl <- as.integer(alive_end_vl == 1L & icu_end_vl == 1L)
        
        rem_idx_vl  <- which(R_vl == 1L)
        exit_idx_vl <- which(R_vl == 0L)
        
        q_exit_obs_vl_raw <- ifelse(D_vl == 1L, q_death_nat_vl_raw, q_dc_nat_vl_raw)
      }

      m_grem_tr <- binom_metric(R_tr, p_rem_nat)
      m_grem_vl <- binom_metric(R_vl, p_rem_nat_vl)

      m_gdeath_tr <- binom_metric(D_tr[exit_idx], p_dex_nat[exit_idx])
      m_gdeath_vl <- binom_metric(D_vl[exit_idx_vl], p_dex_nat_vl[exit_idx_vl])

      m_qrem_tr <- gauss_metric(Y_init[id_tr_ar[rem_idx]], q_rem_nat[rem_idx])
      m_qrem_vl <- gauss_metric(Y_init[id_vl_ar[rem_idx_vl]], q_rem_nat_vl[rem_idx_vl])
      
      if (isTRUE(is_binom)) {
        m_qexit_tr <- binom_metric(
          Y_init[id_tr_ar[exit_idx]],
          q_exit_obs_tr_raw[exit_idx]
        )
        
        m_qexit_vl <- binom_metric(
          Y_init[id_vl_ar[exit_idx_vl]],
          q_exit_obs_vl_raw[exit_idx_vl]
        )
      } else {
        m_qexit_tr <- gauss_metric(
          Y_init[id_tr_ar[exit_idx]],
          q_exit_obs_tr_raw[exit_idx]
        )
        
        m_qexit_vl <- gauss_metric(
          Y_init[id_vl_ar[exit_idx_vl]],
          q_exit_obs_vl_raw[exit_idx_vl]
        )
      }
      
      branch_diag_row <- data.table::data.table(
        fold = f_idx,
        t    = tt,
        
        n_tr_ar   = length(id_tr_ar),
        n_vl_ar   = length(id_vl_ar),
        n_tr_rem  = length(rem_idx),
        n_vl_rem  = length(rem_idx_vl),
        n_tr_exit = length(exit_idx),
        n_vl_exit = length(exit_idx_vl),
        
        target_grem_tr = m_grem_tr$target,
        pred_grem_tr   = m_grem_tr$pred,
        target_grem_vl = m_grem_vl$target,
        pred_grem_vl   = m_grem_vl$pred,
        
        target_gdeath_tr = m_gdeath_tr$target,
        pred_gdeath_tr   = m_gdeath_tr$pred,
        target_gdeath_vl = m_gdeath_vl$target,
        pred_gdeath_vl   = m_gdeath_vl$pred,
        
        target_qrem_tr = m_qrem_tr$target,
        pred_qrem_tr   = m_qrem_tr$pred,
        target_qrem_vl = m_qrem_vl$target,
        pred_qrem_vl   = m_qrem_vl$pred,
        
        target_qexit_tr = m_qexit_tr$target,
        pred_qexit_tr   = m_qexit_tr$pred,
        target_qexit_vl = m_qexit_vl$target,
        pred_qexit_vl   = m_qexit_vl$pred,
        
        qexit_type = if (isTRUE(is_binom)) "binomial" else "gaussian",
        
        grem_tr_brier = m_grem_tr$brier,
        grem_vl_brier = m_grem_vl$brier,
        grem_tr_auc   = m_grem_tr$auc,
        grem_vl_auc   = m_grem_vl$auc,
        grem_tr_cal_int   = m_grem_tr$cal_int,
        grem_vl_cal_int   = m_grem_vl$cal_int,
        grem_tr_cal_slope = m_grem_tr$cal_slope,
        grem_vl_cal_slope = m_grem_vl$cal_slope,
        
        gdeath_tr_brier = m_gdeath_tr$brier,
        gdeath_vl_brier = m_gdeath_vl$brier,
        gdeath_tr_auc   = m_gdeath_tr$auc,
        gdeath_vl_auc   = m_gdeath_vl$auc,
        gdeath_tr_cal_int   = m_gdeath_tr$cal_int,
        gdeath_vl_cal_int   = m_gdeath_vl$cal_int,
        gdeath_tr_cal_slope = m_gdeath_tr$cal_slope,
        gdeath_vl_cal_slope = m_gdeath_vl$cal_slope,
        
        qrem_tr_mse  = m_qrem_tr$mse,
        qrem_vl_mse  = m_qrem_vl$mse,
        qrem_tr_rmse = m_qrem_tr$rmse,
        qrem_vl_rmse = m_qrem_vl$rmse,
        qrem_tr_mae  = m_qrem_tr$mae,
        qrem_vl_mae  = m_qrem_vl$mae,
        qrem_tr_cor  = m_qrem_tr$cor,
        qrem_vl_cor  = m_qrem_vl$cor,
        qrem_tr_cal_int   = m_qrem_tr$cal_int,
        qrem_vl_cal_int   = m_qrem_vl$cal_int,
        qrem_tr_cal_slope = m_qrem_tr$cal_slope,
        qrem_vl_cal_slope = m_qrem_vl$cal_slope,
        
        qexit_tr_brier = if (isTRUE(is_binom)) m_qexit_tr$brier else NA_real_,
        qexit_vl_brier = if (isTRUE(is_binom)) m_qexit_vl$brier else NA_real_,
        qexit_tr_auc   = if (isTRUE(is_binom)) m_qexit_tr$auc else NA_real_,
        qexit_vl_auc   = if (isTRUE(is_binom)) m_qexit_vl$auc else NA_real_,
        
        qexit_tr_mse  = if (!isTRUE(is_binom)) m_qexit_tr$mse else NA_real_,
        qexit_vl_mse  = if (!isTRUE(is_binom)) m_qexit_vl$mse else NA_real_,
        qexit_tr_rmse = if (!isTRUE(is_binom)) m_qexit_tr$rmse else NA_real_,
        qexit_vl_rmse = if (!isTRUE(is_binom)) m_qexit_vl$rmse else NA_real_,
        qexit_tr_mae  = if (!isTRUE(is_binom)) m_qexit_tr$mae else NA_real_,
        qexit_vl_mae  = if (!isTRUE(is_binom)) m_qexit_vl$mae else NA_real_,
        qexit_tr_cor  = if (!isTRUE(is_binom)) m_qexit_tr$cor else NA_real_,
        qexit_vl_cor  = if (!isTRUE(is_binom)) m_qexit_vl$cor else NA_real_,
        
        qexit_tr_cal_int   = m_qexit_tr$cal_int,
        qexit_vl_cal_int   = m_qexit_vl$cal_int,
        qexit_tr_cal_slope = m_qexit_tr$cal_slope,
        qexit_vl_cal_slope = m_qexit_vl$cal_slope
      )
      
      branch_diag_here[[length(branch_diag_here) + 1L]] <- branch_diag_row

      pseudo_pre_ar_tr <- shf_train[at_risk_tr, tt + 1L]

      targ <- itmle_inner_target(
        Q_mix_nat_tr      = Q_nat_ar,
        Q_mix_shf_tr      = Q_shf_ar,
        Q_mix_nat_vl      = Q_nat_vl,
        Q_mix_shf_vl      = Q_shf_vl,
        Y_pseudo_tr       = shf_train[at_risk_tr, tt + 1L],

        D_wide            = D_wide,
        id_tr_ar          = id_tr_ar,
        id_vl_ar          = id_vl_ar,

        density_ratios_tr = density_ratios[id_tr_ar, , drop = FALSE],
        density_ratios_vl = density_ratios[id_vl_ar, , drop = FALSE],

        s                 = tt,
        t0                = 1L,
        K                 = tmax,
        ids               = ids,

        sl_lib_target     = sl_tmle,
        cluster_tr        = cluster_by_id[id_tr_ar],

        v_target          = v_target_itmle,
        v_sl_inner        = v_sl_inner_itmle,
        parallel          = parallel,
        reg_workers       = reg_workers,
        sl_workers        = sl_workers,

        seed              = est_seed_for(seed, fold = f_idx, t = tt,
                                         component = "target", phase = "fit"),

        baseline          = baseline,
        tv_names          = tv_names,
        a_names           = a_names,

        bounds            = bounds,
        clip_fun          = scale_info$clip
      )

      target_diag[[length(target_diag) + 1L]] <- data.table::copy(targ$target_diag)[
        , `:=`(
          fold = f_idx,
          tt   = tt
        )
      ]

      target_sl[[length(target_sl) + 1L]] <- data.table::copy(targ$sl_diag)[
        , `:=`(
          fold = f_idx,
          tt   = tt
        )
      ]

      Q_nat_post  <- targ$Q_star_nat_tr
      Q_shf_post  <- targ$Q_star_shf_tr
      targ_nat_vl <- targ$Q_star_nat_vl
      targ_shf_vl <- targ$Q_star_shf_vl

      nat_train[at_risk_tr, tt] <- targ$Q_star_nat_tr
      shf_train[at_risk_tr, tt] <- targ$Q_star_shf_tr

      if (length(idx_vl_ar) > 0L) {
        nat_valid[idx_vl_ar, tt] <- targ$Q_star_nat_vl
        shf_valid[idx_vl_ar, tt] <- targ$Q_star_shf_vl

        dr_4a_vl <- density_ratios[id_vl_ar, seq_len(tt), drop = FALSE]
        w4a_vl   <- apply(dr_4a_vl, 1L, prod)
        Qnext_vl <- shf_valid[idx_vl_ar, tt + 1L]
        eic_tplus_valid[idx_vl_ar] <- eic_tplus_valid[idx_vl_ar] +
          w4a_vl * (Qnext_vl - nat_valid[idx_vl_ar, tt])
      }

      Y_train <- shf_train[, tt]
      
      exit_branch_nat_tr <- p_dex_nat * q_death_nat +
        (1 - p_dex_nat) * q_dc_nat
      
      exit_branch_shf_tr <- p_dex_shf * q_death_shf +
        (1 - p_dex_shf) * q_dc_shf
      
      diag_row <- data.table::data.table(
        fold = f_idx,
        t    = tt,
        
        n_train_at_risk = sum(at_risk_tr),
        n_valid_at_risk = sum(at_risk_vl),
        
        n_train_death  = sum(D_tr == 1L),
        n_train_dc     = sum(C_tr == 1L),
        n_train_remain = sum(R_tr == 1L),
        
        has_rem   = !is.null(sl_rem),
        has_dex   = !is.null(sl_dex) ||
          (isTRUE(pool_g_death) && !is.null(fit_dex_pooled)),
        has_qexit = !is.null(sl_qexit),
        has_qrem  = !is.null(sl_qrem),
        
        used_const_rem = if (!is.null(sl_rem)) isTRUE(sl_rem$used_const) else NA,
        used_const_dex = if (!is.null(sl_dex)) {
          isTRUE(sl_dex$used_const)
        } else if (isTRUE(pool_g_death) && !is.null(fit_dex_pooled)) {
          FALSE
        } else if (!is.na(p_dex_const)) {
          TRUE
        } else NA,
        used_const_qexit = if (!is.null(sl_qexit)) isTRUE(sl_qexit$used_const) else NA,
        used_const_qrem  = if (!is.null(sl_qrem))  isTRUE(sl_qrem$used_const)  else NA,
        
        p_rem_nat_mean = mean_or_na(p_rem_nat),
        p_rem_nat_sd   = sd_or_na(p_rem_nat),
        p_rem_shf_mean = mean_or_na(p_rem_shf),
        p_rem_shf_sd   = sd_or_na(p_rem_shf),
        
        p_dex_nat_mean = mean_or_na(p_dex_nat),
        p_dex_nat_sd   = sd_or_na(p_dex_nat),
        p_dex_shf_mean = mean_or_na(p_dex_shf),
        p_dex_shf_sd   = sd_or_na(p_dex_shf),
        
        q_rem_nat_mean = mean_or_na(q_rem_nat),
        q_rem_nat_sd   = sd_or_na(q_rem_nat),
        q_rem_shf_mean = mean_or_na(q_rem_shf),
        q_rem_shf_sd   = sd_or_na(q_rem_shf),
        
        q_exit_nat_mean = mean_or_na(exit_branch_nat_tr),
        q_exit_nat_sd   = sd_or_na(exit_branch_nat_tr),
        q_exit_shf_mean = mean_or_na(exit_branch_shf_tr),
        q_exit_shf_sd   = sd_or_na(exit_branch_shf_tr),
        
        rem_contrib_nat_mean = mean_or_na(p_rem_nat * q_rem_nat),
        rem_contrib_shf_mean = mean_or_na(p_rem_shf * q_rem_shf),
        
        exit_contrib_nat_mean = mean_or_na((1 - p_rem_nat) * exit_branch_nat_tr),
        exit_contrib_shf_mean = mean_or_na((1 - p_rem_shf) * exit_branch_shf_tr),
        
        target_pre_mean = mean_or_na(pseudo_pre_ar_tr),
        target_pre_sd   = sd_or_na(pseudo_pre_ar_tr),
        
        Q_nat_pre_mean = mean_or_na(Q_nat_ar),
        Q_nat_pre_sd   = sd_or_na(Q_nat_ar),
        Q_nat_pre_min  = min_or_na(Q_nat_ar),
        Q_nat_pre_max  = max_or_na(Q_nat_ar),
        
        Q_shf_pre_mean = mean_or_na(Q_shf_ar),
        Q_shf_pre_sd   = sd_or_na(Q_shf_ar),
        Q_shf_pre_min  = min_or_na(Q_shf_ar),
        Q_shf_pre_max  = max_or_na(Q_shf_ar),
        
        Q_pre_diff_mean    = mean_or_na(Q_shf_ar - Q_nat_ar),
        Q_pre_diff_sd      = sd_or_na(Q_shf_ar - Q_nat_ar),
        Q_pre_diff_min     = min_or_na(Q_shf_ar - Q_nat_ar),
        Q_pre_diff_max     = max_or_na(Q_shf_ar - Q_nat_ar),
        Q_pre_diff_q95_abs = q_or_na(abs(Q_shf_ar - Q_nat_ar), 0.95),
        
        Q_nat_post_mean = mean_or_na(Q_nat_post),
        Q_nat_post_sd   = sd_or_na(Q_nat_post),
        Q_nat_post_min  = min_or_na(Q_nat_post),
        Q_nat_post_max  = max_or_na(Q_nat_post),
        
        Q_shf_post_mean = mean_or_na(Q_shf_post),
        Q_shf_post_sd   = sd_or_na(Q_shf_post),
        Q_shf_post_min  = min_or_na(Q_shf_post),
        Q_shf_post_max  = max_or_na(Q_shf_post),
        
        Q_post_diff_mean    = mean_or_na(Q_shf_post - Q_nat_post),
        Q_post_diff_sd      = sd_or_na(Q_shf_post - Q_nat_post),
        Q_post_diff_min     = min_or_na(Q_shf_post - Q_nat_post),
        Q_post_diff_max     = max_or_na(Q_shf_post - Q_nat_post),
        Q_post_diff_q95_abs = q_or_na(abs(Q_shf_post - Q_nat_post), 0.95),
        
        delta_nat_target_mean = mean_or_na(Q_nat_post - Q_nat_ar),
        delta_nat_target_sd   = sd_or_na(Q_nat_post - Q_nat_ar),
        delta_nat_target_q95_abs = q_or_na(abs(Q_nat_post - Q_nat_ar), 0.95),
        delta_nat_target_max_abs = max_or_na(abs(Q_nat_post - Q_nat_ar)),
        
        delta_shf_target_mean = mean_or_na(Q_shf_post - Q_shf_ar),
        delta_shf_target_sd   = sd_or_na(Q_shf_post - Q_shf_ar),
        delta_shf_target_q95_abs = q_or_na(abs(Q_shf_post - Q_shf_ar), 0.95),
        delta_shf_target_max_abs = max_or_na(abs(Q_shf_post - Q_shf_ar)),
        
        n_nat_post_below_0 = sum(Q_nat_post < 0, na.rm = TRUE),
        n_nat_post_above_1 = sum(Q_nat_post > 1, na.rm = TRUE),
        n_shf_post_below_0 = sum(Q_shf_post < 0, na.rm = TRUE),
        n_shf_post_above_1 = sum(Q_shf_post > 1, na.rm = TRUE),
        
        Q_nat_vl_pre_mean  = mean_or_na(Q_nat_vl),
        Q_shf_vl_pre_mean  = mean_or_na(Q_shf_vl),
        Q_nat_vl_post_mean = mean_or_na(targ_nat_vl),
        Q_shf_vl_post_mean = mean_or_na(targ_shf_vl),

        delta_nat_vl_target_mean = mean_or_na(targ_nat_vl - Q_nat_vl),
        delta_shf_vl_target_mean = mean_or_na(targ_shf_vl - Q_shf_vl),

        delta_nat_vl_target_q95_abs =
          q_or_na(abs(targ_nat_vl - Q_nat_vl), 0.95),

        delta_shf_vl_target_q95_abs =
          q_or_na(abs(targ_shf_vl - Q_shf_vl), 0.95)
      )
      
      fold_diag_here[[length(fold_diag_here) + 1L]] <- diag_row
      
      rm(sl_rem, sl_dex, sl_qexit, sl_qrem,
         fit_rem, fit_dex, fit_qexit, fit_qrem)
      gc()
      
    } # end tt loop
    
    
    fold_diag_dt <- if (length(fold_diag_here)) {
      data.table::rbindlist(fold_diag_here, use.names = TRUE, fill = TRUE)
    } else NULL
    
    branch_diag_dt <- if (length(branch_diag_here)) {
      data.table::rbindlist(branch_diag_here, use.names = TRUE, fill = TRUE)
    } else NULL
    
    target_diag <- if (length(target_diag)) {
      data.table::rbindlist(target_diag, use.names = TRUE, fill = TRUE)
    } else NULL
    
    target_sl <- if (length(target_sl)) {
      data.table::rbindlist(target_sl, use.names = TRUE, fill = TRUE)
    } else NULL
    
    list(
      fold        = f_idx,
      valid_ids   = vl_ids,
      nat_valid   = nat_valid,
      shf_valid   = shf_valid,
      eic_tplus   = eic_tplus_valid,
      sl_meta     = if (length(sl_chunks_here)) {
        data.table::rbindlist(sl_chunks_here, use.names = TRUE, fill = TRUE)
      } else NULL,
      fold_diag   = fold_diag_dt,
      branch_diag = branch_diag_dt,
      target_diag = target_diag,
      target_sl   = target_sl
    )
  } # end fold_worker
  
  res_by_fold <- par_lapply(
    X        = seq_along(folds),
    FUN      = fold_worker,
    workers  = fold_workers,
    parallel = parallel,
    seed     = TRUE
  )
  
  fold_diag_chunks <- list()
  branch_diag_chunks <- list()
  target_diag_chunks <- list()
  target_sl_chunks <- list()

  eic_tplus_all <- numeric(N)

  failed <- sapply(res_by_fold, function(x) inherits(x, "error") || is.character(x))
  if (any(failed)) {
    for (i in which(failed)) message(sprintf("fold %d: %s", i, as.character(res_by_fold[[i]])))
    stop("SDR: fold worker(s) failed - see messages above")
  }

  for (fr in res_by_fold) {
    if (inherits(fr, "error") || is.character(fr)) {
      message(sprintf("fold failed: %s", as.character(fr)))
      next
    }

    vl_ids <- fr$valid_ids
    pred_nat_all[vl_ids, ] <- fr$nat_valid
    pred_shf_all[vl_ids, ] <- fr$shf_valid
    eic_tplus_all[vl_ids]  <- eic_tplus_all[vl_ids] + fr$eic_tplus
    
    if (!is.null(fr$sl_meta)) {
      sl_chunks[[length(sl_chunks) + 1L]] <- fr$sl_meta
    }
    
    if (!is.null(fr$fold_diag)) {
      fold_diag_chunks[[length(fold_diag_chunks) + 1L]] <- fr$fold_diag
    }
    
    if (!is.null(fr$branch_diag)) {
      branch_diag_chunks[[length(branch_diag_chunks) + 1L]] <- fr$branch_diag
    }
    
    if (!is.null(fr$target_diag)) {
      target_diag_chunks[[length(target_diag_chunks) + 1L]] <- fr$target_diag
    }
    
    if (!is.null(fr$target_sl)) {
      target_sl_chunks[[length(target_sl_chunks) + 1L]] <- fr$target_sl
    }
  }
  
  fold_diag <- if (length(fold_diag_chunks)) {
    data.table::rbindlist(fold_diag_chunks, use.names = TRUE, fill = TRUE)
  } else NULL
  
  branch_diag <- if (length(branch_diag_chunks)) {
    data.table::rbindlist(branch_diag_chunks, use.names = TRUE, fill = TRUE)
  } else NULL
  
  target_diag <- if (length(target_diag_chunks)) {
    data.table::rbindlist(target_diag_chunks, use.names = TRUE, fill = TRUE)
  } else NULL
  
  target_sl <- if (length(target_sl_chunks)) {
    data.table::rbindlist(target_sl_chunks, use.names = TRUE, fill = TRUE)
  } else NULL
  
  t_vec      <- seq_len(tmax)
  n_at_risk  <- integer(tmax)
  mean_Q_nat <- numeric(tmax)
  mean_Q_int <- numeric(tmax)
  delta      <- numeric(tmax)
  mean_Y_obs <- numeric(tmax)
  
  for (tt in t_vec) {
    at_risk <- !is.na(row_index[, tt])
    n_risk  <- sum(at_risk)
    
    if (!n_risk) {
      n_at_risk[tt]  <- 0L
      mean_Q_nat[tt] <- NA_real_
      mean_Q_int[tt] <- NA_real_
      delta[tt]      <- NA_real_
      mean_Y_obs[tt] <- NA_real_
    } else {
      n_at_risk[tt]  <- n_risk
      mean_Q_nat[tt] <- mean(pred_nat_all[at_risk, tt], na.rm = TRUE)
      mean_Q_int[tt] <- mean(pred_shf_all[at_risk, tt], na.rm = TRUE)
      delta[tt]      <- mean_Q_int[tt] - mean_Q_nat[tt]
      mean_Y_obs[tt] <- mean(Y_init[at_risk], na.rm = TRUE)
    }
  }
  
  diag_table <- data.frame(
    t          = t_vec,
    n_at_risk  = n_at_risk,
    mean_Q_nat = mean_Q_nat,
    mean_Q_int = mean_Q_int,
    delta      = delta,
    mean_Y_obs = mean_Y_obs
  )
  
  diag_table <- diag_table[diag_table$n_at_risk > 0L, , drop = FALSE]
  

  stopifnot(all(!is.na(row_index[, 1L])))

  start_t_global <- 1L

  psi_shifted_scaled <- mean(pred_shf_all[, start_t_global], na.rm = TRUE)
  psi_scaled         <- psi_shifted_scaled

  eic_t0_scaled <- pred_shf_all[, start_t_global] - psi_shifted_scaled

  ic_shifted_scaled <- as.numeric(eic_t0_scaled + eic_tplus_all)
  ic_scaled         <- ic_shifted_scaled
  
  targeting_gap_scaled <- mean(ic_scaled[is.finite(ic_scaled)], na.rm = TRUE)
  psi_shifted_onestep_scaled <- psi_shifted_scaled + targeting_gap_scaled
  
  se_fit <- itmle_eic_se(
    ic = ic_scaled,
    cluster_by_id = cluster_by_id,
    use_second_moment = TRUE
  )
  
  se_fit$targeting_gap_scaled <- targeting_gap_scaled
  se_fit$targeting_gap_z      <- if (is.finite(se_fit$se) && se_fit$se > 0) {
    targeting_gap_scaled / se_fit$se
  } else {
    NA_real_
  }
  
  se_scaled <- se_fit$se
  ci_scaled <- psi_scaled + c(-1, 1) * 1.96 * se_scaled
  
  
  if (scale_info$bounded) {
    psi <- scale_info$from_unit(psi_scaled)
    
    ic <- scale_info$y_rng * ic_scaled
    se <- scale_info$y_rng * se_scaled
    
    ci <- psi + c(-1, 1) * 1.96 * se
    
    psi_shifted_onestep <- scale_info$from_unit(psi_shifted_onestep_scaled)
    targeting_gap       <- scale_info$y_rng * targeting_gap_scaled
    
    if (!is.null(diag_table) && nrow(diag_table)) {
      diag_table$mean_Q_nat <- scale_info$from_unit(diag_table$mean_Q_nat)
      diag_table$mean_Q_int <- scale_info$from_unit(diag_table$mean_Q_int)
      diag_table$delta      <- scale_info$y_rng * diag_table$delta
      diag_table$mean_Y_obs <- scale_info$from_unit(diag_table$mean_Y_obs)
    }
    
  } else {
    psi <- psi_scaled
    
    ic <- ic_scaled
    se <- se_scaled
    ci <- ci_scaled
    
    psi_shifted_onestep <- psi_shifted_onestep_scaled
    targeting_gap       <- targeting_gap_scaled
  }
  
  
  sl_summary <- NULL
  if (length(sl_chunks)) {
    sl_summary <- data.table::rbindlist(sl_chunks, fill = TRUE)
    sl_summary <- sl_summary[order(component, t, fold, -weight)]
  }
  
  out <- list(
    psi = psi,
    se  = se,
    ci  = ci,
    ic_df = data.table::data.table(id = ids, ic = ic),

    psi_shifted_onestep = psi_shifted_onestep,
    targeting_gap       = targeting_gap,

    start_t = start_t_global,

    predictions = list(
      natural = pred_nat_all,
      shifted = pred_shf_all
    ),
    
    diagnostics = list(
      diag_table     = diag_table,
      recursion_diag = fold_diag,
      branch_diag    = branch_diag,
      target_cal     = target_diag,
      target_sl      = target_sl,
      sl_summary     = sl_summary,
      se_info        = se_fit
    ),
    
    weights = list(
      weights_dt = {
        x <- data.table::copy(data.table::as.data.table(weights_dt))
        if ("Rt_cum" %in% names(x)) x[, Rt_cum := NULL]
        x
      },
      trim = trim
    ),
    
    settings = list(
      pool_g_death = pool_g_death,
      trim = trim,
      variable_info = list(
        id           = id,
        time         = time,
        alive        = alive,
        in_state     = in_state,
        cluster      = cluster,
        baseline     = baseline,
        time_varying = tv_names,
        treatment    = a_names,
        outcome      = y,
        tmax         = tmax,
        k            = k
      )
    ),
    
    call = match.call()
  )
  
  class(out) <- c("itmle_fit", "list")
  out
}

itmle_competing <- function(
    DT,
    weights_dt,
    policy_spec_fun,
    id       = "id",
    time     = "trial_time",
    alive    = "alive",
    in_state = "in_state",
    cluster  = NULL,
    baseline,
    tv_names,
    a_names,
    k = 2L,
    tmax,
    n_folds  = 5L,
    seed     = 1L,
    sl_q,
    sl_tmle = NULL,
    inner_v          = 5L,
    v_target_itmle   = 10L,
    v_sl_inner_itmle = 10L,
    bounds   = 1e-5,
    trim     = 0.99,
    parallel     = FALSE,
    fold_workers = NULL,
    reg_workers  = NULL,
    sl_workers   = NULL,
    cluster_se_only = FALSE
) {
  stopifnot(requireNamespace("data.table", quietly = TRUE))
  if (is.null(sl_tmle)) sl_tmle <- get("sl_tmle", envir = asNamespace("CausalState"))

  if (!isTRUE(parallel)) sl_workers <- NULL

  set.seed(seed)

  prep <- prepare_long(DT, id = id, time = time, y = alive, cluster = cluster)
  D    <- prep$D
  D    <- D[.__time >= 1L & .__time <= tmax]
  data.table::setorderv(D, c(".__id", ".__time"))

  D <- expand_to_horizon(
    D, id = ".__id", time = ".__time",
    alive = alive, in_state = in_state, tmax = tmax
  )

  lag_vars <- unique(c(tv_names, a_names))
  if (k > 0L && length(lag_vars)) {
    D <- make_lags(D, vars = lag_vars, k = k, id = ".__id", time = ".__time", overwrite = FALSE)
  }

  ids       <- unique(D$.__id)
  N         <- length(ids)
  row_index <- build_id_time_index(D, ids = ids, tmax = tmax)

  cols_to_pivot <- intersect(c(tv_names, a_names), names(D))
  D_wide <- data.table::dcast(
    D[.__time >= 1L],
    formula   = .__id ~ .__time,
    value.var = cols_to_pivot,
    sep       = "_t"
  )
  D_baseline <- unique(D[, .SD, .SDcols = c(".__id", baseline)], by = ".__id")
  D_wide     <- D_baseline[D_wide, on = ".__id"]
  data.table::setorderv(D_wide, ".__id")

  dens <- build_density_ratios(
    weights_dt = weights_dt,
    id         = id,
    time       = time,
    ids        = ids,
    tmax       = tmax,
    row_index  = row_index,
    trim       = trim
  )

  cluster_by_id <- NULL
  if (!is.null(cluster)) {
    cl_dt         <- D[!is.na(.__cl), .(cl = data.table::first(na.omit(.__cl))), by = .__id]
    cluster_by_id <- cl_dt$cl[match(ids, cl_dt$.__id)]
  }

  fold_info <- get_folds_from_weights(
    weights_dt      = weights_dt,
    id              = id,
    ids             = ids,
    cluster_by_id   = cluster_by_id,
    cluster_se_only = cluster_se_only
  )
  folds <- fold_info$folds

  D_shifted_by_t <- vector("list", tmax)
  if (length(a_names)) {
    for (tt in seq_len(tmax)) {
      D_shifted_by_t[[tt]] <- overwrite_policy_history_for_Q(
        D               = D,
        t               = tt,
        a_names         = a_names,
        policy_spec_fun = policy_spec_fun,
        id              = id,
        time            = time,
        k               = k
      )
    }
  } else {
    for (tt in seq_len(tmax)) D_shifted_by_t[[tt]] <- D
  }

  all_names <- names(D)
  cols_by_t <- lapply(seq_len(tmax), function(tt) {
    make_cols(
      tt        = tt,
      baseline  = baseline,
      tv_names  = tv_names,
      a_names   = a_names,
      k         = k,
      all_names = all_names
    )
  })

  alive_nm   <- alive
  died_by_id <- D[, .(Y = as.integer(any(get(alive_nm) == 0L, na.rm = TRUE))), by = .__id]
  Y_by_id    <- died_by_id$Y[match(ids, died_by_id$.__id)]
  stopifnot(all(Y_by_id %in% c(0L, 1L)))

  scale_info <- scale_y(
    D              = data.table::data.table(tmp_y = Y_by_id),
    y              = "tmp_y",
    outcome_family = "binomial",
    bounds         = bounds
  )

  pred_nat_all <- matrix(NA_real_, nrow = N, ncol = tmax + 1L)
  pred_shf_all <- matrix(NA_real_, nrow = N, ncol = tmax + 1L)
  pred_nat_all[, tmax + 1L] <- Y_by_id
  pred_shf_all[, tmax + 1L] <- Y_by_id

  sl_chunks          <- list()
  fold_diag_chunks   <- list()
  target_diag_chunks <- list()
  target_sl_chunks   <- list()

  fold_worker <- function(f_idx) {
    set.seed(seed + f_idx)

    tr   <- folds[[f_idx]]$training_set
    vl   <- folds[[f_idx]]$validation_set
    n_tr <- length(tr)
    n_vl <- length(vl)

    nat_train <- matrix(NA_real_, nrow = n_tr, ncol = tmax + 1L)
    shf_train <- matrix(NA_real_, nrow = n_tr, ncol = tmax + 1L)
    nat_valid <- matrix(NA_real_, nrow = n_vl, ncol = tmax + 1L)
    shf_valid <- matrix(NA_real_, nrow = n_vl, ncol = tmax + 1L)

    nat_train[, tmax + 1L] <- Y_by_id[tr]
    shf_train[, tmax + 1L] <- Y_by_id[tr]
    nat_valid[, tmax + 1L] <- Y_by_id[vl]
    shf_valid[, tmax + 1L] <- Y_by_id[vl]

    Y_train <- as.numeric(Y_by_id[tr])

    sl_chunks_here   <- list()
    fold_diag_here   <- list()
    target_diag_here <- list()
    target_sl_here   <- list()

    .apply_gating <- function(mat, y1, d0, tt_col) {
      mat[!y1, tt_col]      <- 1
      mat[y1 & !d0, tt_col] <- 0
      still_na <- is.na(mat[, tt_col])
      if (any(still_na)) {
        mat[still_na, tt_col] <- mat[still_na, tt_col + 1L]
      }
      mat
    }

    for (tt in rev(seq_len(tmax))) {
      c1_tr <- !is.na(row_index[tr, tt])
      c1_vl <- !is.na(row_index[vl, tt])

      if (tt == 1L) {
        y1_tr <- rep(TRUE, n_tr); d0_tr <- rep(TRUE, n_tr)
        y1_vl <- rep(TRUE, n_vl); d0_vl <- rep(TRUE, n_vl)
      } else {
        prev_tr          <- row_index[tr, tt - 1L]
        okp_tr           <- !is.na(prev_tr)
        alive_prev_tr    <- integer(n_tr)
        in_state_prev_tr <- integer(n_tr)
        alive_prev_tr[okp_tr]    <- as.integer(D[[alive]][prev_tr[okp_tr]])
        in_state_prev_tr[okp_tr] <- as.integer(D[[in_state]][prev_tr[okp_tr]])
        y1_tr <- okp_tr & (alive_prev_tr == 1L)
        d0_tr <- okp_tr & (in_state_prev_tr == 1L)

        prev_vl          <- row_index[vl, tt - 1L]
        okp_vl           <- !is.na(prev_vl)
        alive_prev_vl    <- integer(n_vl)
        in_state_prev_vl <- integer(n_vl)
        alive_prev_vl[okp_vl]    <- as.integer(D[[alive]][prev_vl[okp_vl]])
        in_state_prev_vl[okp_vl] <- as.integer(D[[in_state]][prev_vl[okp_vl]])
        y1_vl <- okp_vl & (alive_prev_vl == 1L)
        d0_vl <- okp_vl & (in_state_prev_vl == 1L)
      }

      i_tr <- c1_tr & y1_tr & d0_tr
      i_vl <- c1_vl & y1_vl & d0_vl

      if (!any(i_tr)) {
        nat_train <- .apply_gating(nat_train, y1_tr, d0_tr, tt)
        shf_train <- .apply_gating(shf_train, y1_tr, d0_tr, tt)
        nat_valid <- .apply_gating(nat_valid, y1_vl, d0_vl, tt)
        shf_valid <- .apply_gating(shf_valid, y1_vl, d0_vl, tt)
        Y_train   <- shf_train[, tt]
        next
      }

      rows_tr_ar <- row_index[tr[i_tr], tt]
      cols_tt    <- cols_by_t[[tt]]

      X_nat_tr <- make_design(D, rows_tr_ar, cols_tt)
      X_shf_tr <- patch_shifted_design(
        X_nat        = X_nat_tr,
        D_shifted_tt = D_shifted_by_t[[tt]],
        rows         = rows_tr_ar,
        a_names      = a_names,
        tt           = tt,
        k            = k
      )

      cl_tr_ar <- if (!is.null(cluster_by_id)) cluster_by_id[tr[i_tr]] else NULL
      id_tr_ar <- tr[which(i_tr)]
      id_vl_ar <- vl[which(i_vl)]

      fam <- if (tt == tmax) stats::binomial() else stats::gaussian()

      sl_fit <- sl_block_fit(
        Y             = Y_train[i_tr],
        X             = X_nat_tr,
        family        = fam,
        sl_lib        = sl_q,
        v_inner       = max(2L, min(as.integer(inner_v), sum(i_tr))),
        id_inner      = if (is.null(cl_tr_ar)) ids[id_tr_ar] else NULL,
        cluster_inner = cl_tr_ar,
        seed_inner    = est_seed_for(seed, f_idx, tt, "Q_remain", "fit"),
        sl_workers    = sl_workers
      )
      fit <- sl_fit$fit

      n_val_tt <- sum(i_vl)
      meta     <- sl_meta(sl_fit, f_idx, tt, "Q", sum(i_tr), n_val_tt)
      if (!is.null(meta)) sl_chunks_here[[length(sl_chunks_here) + 1L]] <- meta

      cn_fit   <- colnames(X_nat_tr)
      q_nat_tr <- scale_info$clip(as.numeric(sl_predict(fit, X_nat_tr)))
      q_shf_tr <- scale_info$clip(as.numeric(sl_predict(fit, X_shf_tr)))

      Q_nat_vl <- numeric(0)
      Q_shf_vl <- numeric(0)
      if (any(i_vl)) {
        rows_vl_ar <- row_index[vl[i_vl], tt]
        X_nat_vl   <- make_design(D, rows_vl_ar, cols_tt)
        X_shf_vl   <- patch_shifted_design(
          X_nat        = X_nat_vl,
          D_shifted_tt = D_shifted_by_t[[tt]],
          rows         = rows_vl_ar,
          a_names      = a_names,
          tt           = tt,
          k            = k
        )
        X_nat_vl <- align_cols(X_nat_vl, cn_fit)
        X_shf_vl <- align_cols(X_shf_vl, cn_fit)
        Q_nat_vl <- scale_info$clip(as.numeric(sl_predict(fit, X_nat_vl)))
        Q_shf_vl <- scale_info$clip(as.numeric(sl_predict(fit, X_shf_vl)))
      }

      targ_nat_post <- q_nat_tr
      targ_shf_post <- q_shf_tr

      if (tt > 1L) {
        targ <- itmle_inner_target(
          Q_mix_nat_tr      = q_nat_tr,
          Q_mix_shf_tr      = q_shf_tr,
          Q_mix_nat_vl      = Q_nat_vl,
          Q_mix_shf_vl      = Q_shf_vl,
          Y_pseudo_tr       = shf_train[i_tr, tt + 1L],
          D_wide            = D_wide,
          id_tr_ar          = id_tr_ar,
          id_vl_ar          = id_vl_ar,
          density_ratios_tr = dens[id_tr_ar, , drop = FALSE],
          density_ratios_vl = dens[id_vl_ar, , drop = FALSE],
          s                 = tt,
          t0                = 0L,
          K                 = tmax,
          ids               = ids,
          sl_lib_target     = sl_tmle,
          cluster_tr        = if (!is.null(cluster_by_id)) cluster_by_id[id_tr_ar] else NULL,
          v_target          = v_target_itmle,
          v_sl_inner        = v_sl_inner_itmle,
          parallel          = parallel,
          reg_workers       = reg_workers,
          sl_workers        = sl_workers,
          seed              = est_seed_for(seed, fold = f_idx, t = tt,
                                           component = "target", phase = "fit"),
          baseline          = baseline,
          tv_names          = tv_names,
          a_names           = a_names,
          bounds            = bounds,
          clip_fun          = scale_info$clip
        )

        targ_nat_post <- targ$Q_star_nat_tr
        targ_shf_post <- targ$Q_star_shf_tr

        nat_train[i_tr, tt] <- targ$Q_star_nat_tr
        shf_train[i_tr, tt] <- targ$Q_star_shf_tr

        if (any(i_vl)) {
          nat_valid[i_vl, tt] <- targ$Q_star_nat_vl
          shf_valid[i_vl, tt] <- targ$Q_star_shf_vl
        }

        target_diag_here[[length(target_diag_here) + 1L]] <-
          data.table::copy(targ$target_diag)[, `:=`(fold = f_idx, tt = tt)]

        if (!is.null(targ$sl_diag) && nrow(targ$sl_diag) > 0L) {
          target_sl_here[[length(target_sl_here) + 1L]] <-
            data.table::copy(targ$sl_diag)[, `:=`(fold = f_idx, tt = tt)]
        }
      } else {
        if (any(i_vl)) {
          nat_valid[i_vl, tt] <- Q_nat_vl
          shf_valid[i_vl, tt] <- Q_shf_vl
        }
      }

      nat_train <- .apply_gating(nat_train, y1_tr, d0_tr, tt)
      shf_train <- .apply_gating(shf_train, y1_tr, d0_tr, tt)
      nat_valid <- .apply_gating(nat_valid, y1_vl, d0_vl, tt)
      shf_valid <- .apply_gating(shf_valid, y1_vl, d0_vl, tt)

      Y_train <- shf_train[, tt]

      fold_diag_here[[length(fold_diag_here) + 1L]] <- data.table::data.table(
        fold            = f_idx,
        t               = tt,
        n_train_at_risk = sum(i_tr),
        n_valid_at_risk = sum(i_vl),
        mean_Y_pseudo   = mean(Y_train[i_tr], na.rm = TRUE),
        mean_Q_nat_post = mean(targ_nat_post, na.rm = TRUE),
        mean_Q_shf_post = mean(targ_shf_post, na.rm = TRUE),
        used_const      = isTRUE(sl_fit$used_const)
      )
    } # end tt loop

    fold_diag_dt   <- if (length(fold_diag_here)) {
      data.table::rbindlist(fold_diag_here, use.names = TRUE, fill = TRUE)
    } else NULL
    target_diag_dt <- if (length(target_diag_here)) {
      data.table::rbindlist(target_diag_here, use.names = TRUE, fill = TRUE)
    } else NULL
    target_sl_dt   <- if (length(target_sl_here)) {
      data.table::rbindlist(target_sl_here, use.names = TRUE, fill = TRUE)
    } else NULL

    list(
      fold      = f_idx,
      vl        = vl,
      nat_valid = nat_valid,
      shf_valid = shf_valid,
      sl_meta   = if (length(sl_chunks_here)) {
        data.table::rbindlist(sl_chunks_here, use.names = TRUE, fill = TRUE)
      } else NULL,
      fold_diag   = fold_diag_dt,
      target_diag = target_diag_dt,
      target_sl   = target_sl_dt
    )
  } # end fold_worker

  res_by_fold <- par_lapply(
    X        = seq_along(folds),
    FUN      = fold_worker,
    workers  = fold_workers,
    parallel = parallel,
    seed     = TRUE
  )

  failed <- sapply(res_by_fold, function(x) inherits(x, "error") || is.character(x))
  if (any(failed)) {
    for (i in which(failed)) message(sprintf("fold %d: %s", i, as.character(res_by_fold[[i]])))
    stop("itmle_competing: fold worker(s) failed — see messages above", call. = FALSE)
  }

  for (res in res_by_fold) {
    vl <- res$vl
    pred_nat_all[vl, ] <- res$nat_valid
    pred_shf_all[vl, ] <- res$shf_valid

    if (!is.null(res$sl_meta))
      sl_chunks[[length(sl_chunks) + 1L]] <- res$sl_meta
    if (!is.null(res$fold_diag))
      fold_diag_chunks[[length(fold_diag_chunks) + 1L]] <- res$fold_diag
    if (!is.null(res$target_diag))
      target_diag_chunks[[length(target_diag_chunks) + 1L]] <- res$target_diag
    if (!is.null(res$target_sl))
      target_sl_chunks[[length(target_sl_chunks) + 1L]] <- res$target_sl
  }

  fold_diag <- if (length(fold_diag_chunks)) {
    data.table::rbindlist(fold_diag_chunks, use.names = TRUE, fill = TRUE)
  } else NULL
  target_diag <- if (length(target_diag_chunks)) {
    data.table::rbindlist(target_diag_chunks, use.names = TRUE, fill = TRUE)
  } else NULL
  target_sl <- if (length(target_sl_chunks)) {
    data.table::rbindlist(target_sl_chunks, use.names = TRUE, fill = TRUE)
  } else NULL

  backcumdr_all <- matrix(1, nrow = N, ncol = tmax)
  if (tmax > 1L) {
    for (t_bc in rev(seq_len(tmax - 1L))) {
      backcumdr_all[, t_bc] <- backcumdr_all[, t_bc + 1L] * dens[, t_bc + 1L]
    }
  }

  eic_tplus_all <- numeric(N)
  fluct_diag    <- vector("list", tmax)

  for (t in rev(seq_len(tmax))) {
    ar <- !is.na(row_index[, t])
    if (!any(ar)) next

    Qnext    <- pred_shf_all[ar, t + 1L]
    Q_nat_t  <- pred_nat_all[ar, t]
    Q_shf_t  <- pred_shf_all[ar, t]
    w_t      <- backcumdr_all[ar, t]

    Q_nat_cl <- scale_info$clip(Q_nat_t)
    Q_shf_cl <- scale_info$clip(Q_shf_t)
    Qnext_cl <- scale_info$clip(Qnext)

    eps_t <- 0
    ok_w  <- is.finite(w_t) & w_t > 0
    if (sum(ok_w) >= 2L) {
      fit_t <- tryCatch(
        stats::glm.fit(
          x       = matrix(1L, sum(ar), 1L),
          y       = Qnext_cl,
          offset  = qlogis(Q_nat_cl),
          weights = w_t,
          family  = stats::quasibinomial()
        ),
        error = function(e) NULL
      )
      if (!is.null(fit_t) && is.finite(fit_t$coefficients[1L])) {
        eps_t <- fit_t$coefficients[1L]
      }
    }

    pred_nat_all[ar, t] <- scale_info$clip(plogis(qlogis(Q_nat_cl) + eps_t))
    pred_shf_all[ar, t] <- scale_info$clip(plogis(qlogis(Q_shf_cl) + eps_t))

    eic_tplus_all[ar] <- eic_tplus_all[ar] +
      w_t * (Qnext - pred_nat_all[ar, t])

    sw  <- sum(w_t)
    ess <- if (sw > 0) sw^2 / sum(w_t^2) else NA_real_
    fluct_diag[[t]] <- data.table::data.table(
      t               = t,
      n_ar            = sum(ar),
      eps             = eps_t,
      ess_r           = ess / sum(ar),
      mean_Q_nat_pre  = mean(Q_nat_t),
      mean_Q_shf_pre  = mean(Q_shf_t),
      mean_Q_nat_post = mean(pred_nat_all[ar, t]),
      mean_Q_shf_post = mean(pred_shf_all[ar, t]),
      delta_nat       = mean(pred_nat_all[ar, t] - Q_nat_t),
      delta_shf       = mean(pred_shf_all[ar, t] - Q_shf_t)
    )
  }

  fluct_diag <- data.table::rbindlist(fluct_diag, use.names = TRUE, fill = TRUE)

  psi    <- mean(pred_shf_all[, 1L], na.rm = TRUE)
  eic_t0 <- pred_shf_all[, 1L] - psi
  ic     <- as.numeric(eic_t0 + eic_tplus_all)

  targeting_gap <- mean(ic[is.finite(ic)], na.rm = TRUE)
  psi_onestep   <- psi + targeting_gap

  se_fit <- itmle_eic_se(
    ic                = ic,
    cluster_by_id     = cluster_by_id,
    use_second_moment = TRUE
  )
  se_fit$targeting_gap <- targeting_gap
  se_fit$targeting_gap_z <- if (is.finite(se_fit$se) && se_fit$se > 0) {
    targeting_gap / se_fit$se
  } else NA_real_

  se <- se_fit$se
  ci <- psi + c(-1, 1) * 1.96 * se

  sl_summary <- if (length(sl_chunks)) {
    dt <- data.table::rbindlist(sl_chunks, fill = TRUE)
    dt[order(component, t, fold, -weight)]
  } else NULL

  list(
    psi           = psi,
    se            = se,
    ci            = ci,
    ic_df         = data.table::data.table(id = ids, ic = ic),
    psi_onestep   = psi_onestep,
    targeting_gap = targeting_gap,
    predictions   = list(
      natural = pred_nat_all,
      shifted = pred_shf_all
    ),
    sl_summary    = sl_summary,
    fold_diag     = fold_diag,
    target_diag   = target_diag,
    target_sl     = target_sl,
    fluct_diag    = fluct_diag
  )
}







