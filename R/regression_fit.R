fit_pooled_g_death_exit <- function(
    D,
    tr_ids,
    row_index,
    tmax,
    baseline,
    tv_names,
    a_names,
    all_names,
    k,
    t_min,
    alive,
    in_state,
    cluster_by_id,
    sl_death,
    inner_v,
    seed,
    f_idx,
    pool_time_basis = NULL,
    min_n = 30L,
    min_events = 5L,
    sl_workers = NULL
) {
  if (is.null(pool_time_basis))
    pool_time_basis <- make_pool_time_basis(tmax, mode = "linear")
  cols_pool <- make_cols(
    tt        = 1L,
    baseline  = baseline,
    tv_names  = tv_names,
    a_names   = a_names,
    k         = k,
    all_names = all_names,
    t_min     = t_min
  )

  pool_X  <- vector("list", tmax)
  pool_Y  <- vector("list", tmax)
  pool_cl <- vector("list", tmax)

  for (tt in seq_len(tmax)) {
    at_risk_p <- !is.na(row_index[tr_ids, tt])
    if (!any(at_risk_p)) next

    id_tr_ar_p <- tr_ids[at_risk_p]
    rows_tr_p  <- row_index[id_tr_ar_p, tt]

    alive_p <- as.integer(D[[alive]][rows_tr_p])
    in_state_p   <- as.integer(D[[in_state]][rows_tr_p])

    R_p <- as.integer(alive_p == 1L & in_state_p == 1L)
    D_p <- as.integer(alive_p == 0L)

    exit_idx_p <- which(R_p == 0L)
    if (!length(exit_idx_p)) next

    X_p <- make_design(D, rows_tr_p[exit_idx_p], cols_pool)
    X_p <- apply_pool_time_basis(X_p, tt, pool_time_basis)

    pool_X[[tt]]  <- X_p
    pool_Y[[tt]]  <- D_p[exit_idx_p]
    pool_cl[[tt]] <- cluster_by_id[id_tr_ar_p[exit_idx_p]]
  }

  pool_X  <- Filter(Negate(is.null), pool_X)
  pool_Y  <- Filter(Negate(is.null), pool_Y)
  pool_cl <- Filter(Negate(is.null), pool_cl)

  if (!length(pool_X)) {
    return(list(
      fit = NULL,
      sl = NULL,
      cols = NULL,
      cols_pool = cols_pool,
      X = NULL,
      Y = integer(0),
      cluster = NULL,
      reason = "no_exit_rows"
    ))
  }

  all_X_pool  <- data.table::rbindlist(pool_X, use.names = TRUE, fill = TRUE)
  all_Y_pool  <- unlist(pool_Y, use.names = FALSE)
  all_cl_pool <- unlist(pool_cl, use.names = FALSE)

  chk <- can_fit_bin(
    all_Y_pool,
    min_n      = min_n,
    min_events = min_events
  )

  if (!chk$ok || nrow(all_X_pool) < 2L) {
    return(list(
      fit = NULL,
      sl = NULL,
      cols = colnames(all_X_pool),
      cols_pool = cols_pool,
      X = all_X_pool,
      Y = all_Y_pool,
      cluster = all_cl_pool,
      reason = "sparse_or_one_class"
    ))
  }

  sl <- sl_block_fit(
    Y             = all_Y_pool,
    X             = all_X_pool,
    family        = stats::binomial(),
    sl_lib        = sl_death,
    v_inner       = max(2L, min(as.integer(inner_v), nrow(all_X_pool))),
    cluster_inner = all_cl_pool,
    seed_inner    = est_seed_for(seed, fold = f_idx, t = 0L,
                                 component = "g_death_exit_pooled", phase = "fit"),
    seed_cv_rows  = est_seed_for(seed, fold = f_idx, t = 0L,
                                 component = "g_death_exit_pooled", phase = "cv_rows"),
    seed_sl_fit   = est_seed_for(seed, fold = f_idx, t = 0L,
                                 component = "g_death_exit_pooled", phase = "sl_fit"),
    sl_workers    = sl_workers
  )

  list(
    fit = sl$fit,
    sl = sl,
    cols = colnames(all_X_pool),
    cols_pool = cols_pool,
    X = all_X_pool,
    Y = all_Y_pool,
    cluster = all_cl_pool,
    reason = "fit"
  )
}

fit_pooled_q_exit <- function(
    D,
    tr_ids,
    row_index,
    tmax,
    baseline,
    tv_names,
    a_names,
    all_names,
    k,
    t_min,
    alive,
    in_state,
    Y_init,
    cluster_by_id,
    sl_y,
    inner_v,
    seed,
    f_idx,
    is_binom,
    pool_time_basis = NULL,
    min_n = 30L,
    min_events = 5L,
    sl_workers = NULL
) {
  if (is.null(pool_time_basis))
    pool_time_basis <- make_pool_time_basis(tmax, mode = "linear")
  cols_pool <- make_cols(
    tt        = 1L,
    baseline  = baseline,
    tv_names  = tv_names,
    a_names   = a_names,
    k         = k,
    all_names = all_names,
    t_min     = t_min
  )

  pool_X  <- vector("list", tmax)
  pool_Y  <- vector("list", tmax)
  pool_cl <- vector("list", tmax)

  for (tt in seq_len(tmax)) {
    at_risk_p <- !is.na(row_index[tr_ids, tt])
    if (!any(at_risk_p)) next

    id_tr_ar_p <- tr_ids[at_risk_p]
    rows_tr_p  <- row_index[id_tr_ar_p, tt]

    alive_p    <- as.integer(D[[alive]][rows_tr_p])
    in_state_p <- as.integer(D[[in_state]][rows_tr_p])
    R_p <- as.integer(alive_p == 1L & in_state_p == 1L)
    D_p <- as.integer(alive_p == 0L)

    exit_idx_p <- which(R_p == 0L)
    if (!length(exit_idx_p)) next

    X_p <- make_design(D, rows_tr_p[exit_idx_p], cols_pool)
    X_p$exit_status <- D_p[exit_idx_p]
    X_p <- apply_pool_time_basis(X_p, tt, pool_time_basis)

    pool_X[[tt]]  <- X_p
    pool_Y[[tt]]  <- Y_init[id_tr_ar_p[exit_idx_p]]
    pool_cl[[tt]] <- cluster_by_id[id_tr_ar_p[exit_idx_p]]
  }

  pool_X  <- Filter(Negate(is.null), pool_X)
  pool_Y  <- Filter(Negate(is.null), pool_Y)
  pool_cl <- Filter(Negate(is.null), pool_cl)

  if (!length(pool_X)) {
    return(list(
      fit = NULL, sl = NULL, cols = NULL, cols_pool = cols_pool,
      X = NULL, Y = numeric(0), cluster = NULL, reason = "no_exit_rows"
    ))
  }

  all_X_pool  <- data.table::rbindlist(pool_X, use.names = TRUE, fill = TRUE)
  all_Y_pool  <- unlist(pool_Y, use.names = FALSE)
  all_cl_pool <- unlist(pool_cl, use.names = FALSE)

  if (isTRUE(is_binom)) {
    chk <- can_fit_bin(as.integer(round(all_Y_pool)),
                       min_n = min_n, min_events = min_events)
  } else {
    chk <- list(ok = nrow(all_X_pool) >= min_n)
  }

  if (!chk$ok || nrow(all_X_pool) < 2L) {
    return(list(
      fit = NULL, sl = NULL, cols = colnames(all_X_pool), cols_pool = cols_pool,
      X = all_X_pool, Y = all_Y_pool, cluster = all_cl_pool,
      reason = "sparse_or_one_class"
    ))
  }

  sl <- sl_block_fit(
    Y             = all_Y_pool,
    X             = all_X_pool,
    family        = if (isTRUE(is_binom)) stats::binomial() else stats::gaussian(),
    sl_lib        = sl_y,
    v_inner       = max(2L, min(as.integer(inner_v), nrow(all_X_pool))),
    cluster_inner = all_cl_pool,
    seed_inner    = est_seed_for(seed, fold = f_idx, t = 0L,
                                 component = "Q_exit_pooled", phase = "fit"),
    seed_cv_rows  = est_seed_for(seed, fold = f_idx, t = 0L,
                                 component = "Q_exit_pooled", phase = "cv_rows"),
    seed_sl_fit   = est_seed_for(seed, fold = f_idx, t = 0L,
                                 component = "Q_exit_pooled", phase = "sl_fit"),
    sl_workers    = sl_workers
  )

  list(
    fit       = sl$fit,
    sl        = sl,
    cols      = colnames(all_X_pool),
    cols_pool = cols_pool,
    X         = all_X_pool,
    Y         = all_Y_pool,
    cluster   = all_cl_pool,
    reason    = "fit"
  )
}

slim_sl <- function(sl_res) {
  if (is.null(sl_res)) return(NULL)
  list(
    fit        = list(coef = sl_res$fit$coef, cvRisk = sl_res$fit$cvRisk),
    v_eff      = sl_res$v_eff,
    used_const = sl_res$used_const
  )
}

build_pooled_pred_cache <- function(
  tmax, tr_ids, vl_ids, row_index,
  D, D_shifted,
  pool_g_death, fit_dex_pooled, cols_pool, cn_pool,
  pool_q_exit, fit_qexit_pooled, cols_pool_qexit, cn_qexit_pool,
  pool_time_basis, a_names, k, t_min, scale_info
) {
  tr_cache <- vector("list", tmax)
  vl_cache <- vector("list", tmax)

  for (tt in seq_len(tmax)) {
    at_risk_tr <- !is.na(row_index[tr_ids, tt])
    at_risk_vl <- !is.na(row_index[vl_ids, tt])
    if (!any(at_risk_tr)) next

    id_tr_ar <- tr_ids[which(at_risk_tr)]
    id_vl_ar <- vl_ids[which(at_risk_vl)]
    rows_tr  <- row_index[id_tr_ar, tt]
    rows_vl  <- row_index[id_vl_ar, tt]
    has_vl   <- length(rows_vl) > 0L

    tc <- list()
    vc <- list()

    if (isTRUE(pool_g_death) && !is.null(fit_dex_pooled)) {
      make_dex_X <- function(rows) {
        X_nat <- make_design(D, rows, cols_pool)
        X_nat <- apply_pool_time_basis(X_nat, tt, pool_time_basis$g_death)
        X_shf <- patch_shifted_design(
          X_nat = X_nat, D_shifted_tt = D_shifted[[tt]],
          rows = rows, a_names = a_names, tt = tt, k = k, t_min = t_min)
        X_shf <- apply_pool_time_basis(X_shf, tt, pool_time_basis$g_death)
        list(nat = align_cols(X_nat, cn_pool), shf = align_cols(X_shf, cn_pool))
      }
      Xtr <- make_dex_X(rows_tr)
      tc$p_dex_nat <- scale_info$clip(sl_predict(fit_dex_pooled, Xtr$nat))
      tc$p_dex_shf <- scale_info$clip(sl_predict(fit_dex_pooled, Xtr$shf))
      if (has_vl) {
        Xvl <- make_dex_X(rows_vl)
        vc$p_dex_nat <- scale_info$clip(sl_predict(fit_dex_pooled, Xvl$nat))
        vc$p_dex_shf <- scale_info$clip(sl_predict(fit_dex_pooled, Xvl$shf))
      }
    }

    if (isTRUE(pool_q_exit) && !is.null(fit_qexit_pooled)) {
      make_qex_X <- function(rows) {
        X_nat <- make_design(D, rows, cols_pool_qexit)
        X_shf <- patch_shifted_design(
          X_nat = X_nat, D_shifted_tt = D_shifted[[tt]],
          rows = rows, a_names = a_names, tt = tt, k = k, t_min = t_min)
        X_nat_d <- apply_pool_time_basis(X_nat, tt, pool_time_basis$q_exit); X_nat_d$exit_status <- 1
        X_nat_c <- apply_pool_time_basis(X_nat, tt, pool_time_basis$q_exit); X_nat_c$exit_status <- 0
        X_shf_d <- apply_pool_time_basis(X_shf, tt, pool_time_basis$q_exit); X_shf_d$exit_status <- 1
        X_shf_c <- apply_pool_time_basis(X_shf, tt, pool_time_basis$q_exit); X_shf_c$exit_status <- 0
        list(
          nd = align_cols(X_nat_d, cn_qexit_pool), nc = align_cols(X_nat_c, cn_qexit_pool),
          sd = align_cols(X_shf_d, cn_qexit_pool), sc = align_cols(X_shf_c, cn_qexit_pool)
        )
      }
      Xtr <- make_qex_X(rows_tr)
      tc$q_death_nat <- scale_info$clip(sl_predict(fit_qexit_pooled, Xtr$nd))
      tc$q_dc_nat    <- scale_info$clip(sl_predict(fit_qexit_pooled, Xtr$nc))
      tc$q_death_shf <- scale_info$clip(sl_predict(fit_qexit_pooled, Xtr$sd))
      tc$q_dc_shf    <- scale_info$clip(sl_predict(fit_qexit_pooled, Xtr$sc))
      if (has_vl) {
        Xvl <- make_qex_X(rows_vl)
        vc$q_death_nat <- scale_info$clip(sl_predict(fit_qexit_pooled, Xvl$nd))
        vc$q_dc_nat    <- scale_info$clip(sl_predict(fit_qexit_pooled, Xvl$nc))
        vc$q_death_shf <- scale_info$clip(sl_predict(fit_qexit_pooled, Xvl$sd))
        vc$q_dc_shf    <- scale_info$clip(sl_predict(fit_qexit_pooled, Xvl$sc))
      }
    }

    tr_cache[[tt]] <- if (length(tc)) tc else NULL
    vl_cache[[tt]] <- if (length(vc)) vc else NULL
  }

  list(tr = tr_cache, vl = vl_cache)
}

fit_transition_regressions <- function(
  X_nat_tr_base,
  X_shf_tr_base,
  X_nat_tr_haz,
  X_shf_tr_haz,
  Y_pseudo_ar,
  Y_init_ar,
  R_tr,
  D_tr,
  exit_idx,
  rem_idx,
  n_ar,

  cl_tr_ar,
  cl_exit_ar,
  cl_rem_ar,

  pool_g_death = FALSE,
  fit_dex_pooled = NULL,
  cn_pool = NULL,
  cols_pool = NULL,
  pool_q_exit = FALSE,
  fit_qexit_pooled = NULL,
  cn_qexit_pool = NULL,
  cols_pool_qexit = NULL,
  pool_time_basis = NULL,
  D = NULL,
  rows_tr = NULL,
  D_shifted_tt = NULL,
  a_names = NULL,
  tt = NULL,
  k = NULL,
  t_min = NULL,

  pooled_pred_tr = NULL,
  pooled_pred_vl = NULL,
  rows_vl = NULL,
  X_nat_vl_base = NULL,
  X_shf_vl_base = NULL,

  is_binom,
  tmax,
  scale_info,
  sl_remain,
  sl_death,
  sl_y,
  sl_recursive,
  sl_rec_early = NULL,
  rec_transition = NULL,
  inner_v,
  seed,
  f_idx,

  parallel = FALSE,
  reg_workers = NULL,
  sl_workers = NULL
) {
  fit_g_remain_one <- function() {
    p_rem_const <- NA_real_
    fit_rem <- NULL
    sl_rem <- NULL

    chkR <- can_fit_bin(R_tr, min_n = 30L, min_events = 5L)

    if (!chkR$ok) {
      p_rem_const <- scale_info$clip(chkR$p)
      p_rem_nat <- rep_const(p_rem_const, n_ar)
      p_rem_shf <- rep_const(p_rem_const, n_ar)
    } else {
      sl_rem <- sl_block_fit(
        Y = R_tr,
        X = X_nat_tr_haz,
        family = stats::binomial(),
        sl_lib = sl_remain,
        v_inner = max(2L, min(as.integer(inner_v), n_ar)),
        cluster_inner = cl_tr_ar,
        seed_inner = est_seed_for(seed, f_idx, tt, "g_remain", "fit"),
        seed_cv_rows = est_seed_for(seed, f_idx, tt, "g_remain", "cv_rows"),
        seed_sl_fit = est_seed_for(seed, f_idx, tt, "g_remain", "sl_fit"),
        sl_workers = sl_workers
      )
      fit_rem <- sl_rem$fit
      p_rem_nat <- scale_info$clip(sl_predict(fit_rem, X_nat_tr_haz))
      p_rem_shf <- scale_info$clip(sl_predict(fit_rem, X_shf_tr_haz))
    }

    list(
      sl_rem = sl_rem,
      fit_rem = fit_rem,
      p_rem_const = p_rem_const,
      p_rem_nat = p_rem_nat,
      p_rem_shf = p_rem_shf
    )
  }

  fit_g_death_exit_one <- function() {
    if (!is.null(pooled_pred_tr$p_dex_nat)) {
      return(list(
        sl_dex      = NULL,
        fit_dex     = NULL,
        p_dex_const = NA_real_,
        p_dex_nat   = pooled_pred_tr$p_dex_nat,
        p_dex_shf   = pooled_pred_tr$p_dex_shf
      ))
    }

    p_dex_const <- NA_real_
    p_dex_nat <- rep_const(0, n_ar)
    p_dex_shf <- rep_const(0, n_ar)
    fit_dex <- NULL
    sl_dex <- NULL

    if (isTRUE(pool_g_death)) {
      if (is.null(fit_dex_pooled)) {
        stop(sprintf("[fold %d][t=%d] pooled g_death fit is NULL", f_idx, tt), call. = FALSE)
      }

      X_nat_dex_t <- make_design(D, rows_tr, cols_pool)
      X_nat_dex_t <- apply_pool_time_basis(X_nat_dex_t, tt, pool_time_basis$g_death)

      X_shf_dex_t <- patch_shifted_design(
        X_nat = X_nat_dex_t,
        D_shifted_tt = D_shifted_tt,
        rows = rows_tr,
        a_names = a_names,
        tt = tt,
        k = k,
        t_min = t_min
      )
      X_shf_dex_t <- apply_pool_time_basis(X_shf_dex_t, tt, pool_time_basis$g_death)

      X_nat_dex_t <- align_cols(X_nat_dex_t, cn_pool)
      X_shf_dex_t <- align_cols(X_shf_dex_t, cn_pool)

      set.seed(est_seed_for(seed, f_idx, 0L, "g_death_exit_pooled", "sl_fit"))
      p_dex_nat <- scale_info$clip(sl_predict(fit_dex_pooled, X_nat_dex_t))
      p_dex_shf <- scale_info$clip(sl_predict(fit_dex_pooled, X_shf_dex_t))
      fit_dex <- fit_dex_pooled

    } else if (length(exit_idx) >= 2L) {
      chkDex <- can_fit_bin(D_tr[exit_idx], min_n = 30L, min_events = 5L)

      if (!chkDex$ok) {
        p_dex_const <- scale_info$clip(chkDex$p)
        p_dex_nat <- rep_const(p_dex_const, n_ar)
        p_dex_shf <- rep_const(p_dex_const, n_ar)
      } else {
        sl_dex <- sl_block_fit(
          Y = D_tr[exit_idx],
          X = X_nat_tr_haz[exit_idx, , drop = FALSE],
          family = stats::binomial(),
          sl_lib = sl_death,
          v_inner = max(2L, min(as.integer(inner_v), length(exit_idx))),
          cluster_inner = cl_exit_ar,
          seed_inner = est_seed_for(seed, f_idx, tt, "g_death_exit", "fit"),
          seed_cv_rows = est_seed_for(seed, f_idx, tt, "g_death_exit", "cv_rows"),
          seed_sl_fit = est_seed_for(seed, f_idx, tt, "g_death_exit", "sl_fit"),
          sl_workers = sl_workers
        )
        fit_dex <- sl_dex$fit
        p_dex_nat <- scale_info$clip(sl_predict(fit_dex, X_nat_tr_haz))
        p_dex_shf <- scale_info$clip(sl_predict(fit_dex, X_shf_tr_haz))
      }

    } else {
      if (length(exit_idx) >= 1L) {
        p_dex_const <- scale_info$clip(mean(D_tr[exit_idx], na.rm = TRUE))
      } else {
        p_dex_const <- 0
      }
      p_dex_nat <- rep_const(p_dex_const, n_ar)
      p_dex_shf <- rep_const(p_dex_const, n_ar)
    }

    list(
      sl_dex = sl_dex,
      fit_dex = fit_dex,
      p_dex_const = p_dex_const,
      p_dex_nat = p_dex_nat,
      p_dex_shf = p_dex_shf
    )
  }

  fit_Q_exit_one <- function() {
    muY <- if (length(exit_idx)) {
      mean(Y_init_ar[exit_idx], na.rm = TRUE)
    } else {
      mean(Y_init_ar, na.rm = TRUE)
    }

    if (!is.null(pooled_pred_tr$q_death_nat)) {
      return(list(
        sl_qexit    = NULL,
        fit_qexit   = NULL,
        cn_qexit    = NULL,
        muY         = muY,
        q_death_nat = pooled_pred_tr$q_death_nat,
        q_dc_nat    = pooled_pred_tr$q_dc_nat,
        q_death_shf = pooled_pred_tr$q_death_shf,
        q_dc_shf    = pooled_pred_tr$q_dc_shf
      ))
    }

    q_death_nat <- rep(muY, n_ar)
    q_death_shf <- rep(muY, n_ar)
    q_dc_nat <- rep(muY, n_ar)
    q_dc_shf <- rep(muY, n_ar)
    fit_qexit <- NULL
    cn_qexit <- NULL
    sl_qexit <- NULL

    if (isTRUE(pool_q_exit) && !is.null(fit_qexit_pooled)) {
      X_nat_pool <- make_design(D, rows_tr, cols_pool_qexit)
      X_shf_pool <- patch_shifted_design(
        X_nat        = X_nat_pool,
        D_shifted_tt = D_shifted_tt,
        rows         = rows_tr,
        a_names      = a_names,
        tt           = tt,
        k            = k,
        t_min        = t_min
      )

      X_nat_d <- apply_pool_time_basis(X_nat_pool, tt, pool_time_basis$q_exit); X_nat_d$exit_status <- 1
      X_nat_c <- apply_pool_time_basis(X_nat_pool, tt, pool_time_basis$q_exit); X_nat_c$exit_status <- 0
      X_shf_d <- apply_pool_time_basis(X_shf_pool, tt, pool_time_basis$q_exit); X_shf_d$exit_status <- 1
      X_shf_c <- apply_pool_time_basis(X_shf_pool, tt, pool_time_basis$q_exit); X_shf_c$exit_status <- 0

      X_nat_d <- align_cols(X_nat_d, cn_qexit_pool)
      X_nat_c <- align_cols(X_nat_c, cn_qexit_pool)
      X_shf_d <- align_cols(X_shf_d, cn_qexit_pool)
      X_shf_c <- align_cols(X_shf_c, cn_qexit_pool)

      set.seed(est_seed_for(seed, f_idx, 0L, "Q_exit_pooled", "sl_fit"))
      q_death_nat <- scale_info$clip(sl_predict(fit_qexit_pooled, X_nat_d))
      q_dc_nat    <- scale_info$clip(sl_predict(fit_qexit_pooled, X_nat_c))
      q_death_shf <- scale_info$clip(sl_predict(fit_qexit_pooled, X_shf_d))
      q_dc_shf    <- scale_info$clip(sl_predict(fit_qexit_pooled, X_shf_c))

      fit_qexit <- fit_qexit_pooled
      cn_qexit  <- cn_qexit_pool

    } else if (length(exit_idx) >= 2L) {
      X_exit <- X_nat_tr_base[exit_idx, , drop = FALSE]
      X_exit$exit_status <- as.numeric(D_tr[exit_idx])
      Y_exit <- Y_init_ar[exit_idx]

      sl_qexit <- sl_block_fit(
        Y = Y_exit,
        X = X_exit,
        family = if (is_binom) stats::binomial() else stats::gaussian(),
        sl_lib = sl_y,
        v_inner = max(2L, min(as.integer(inner_v), length(exit_idx))),
        cluster_inner = cl_exit_ar,
        seed_inner = est_seed_for(seed, f_idx, tt, "Q_exit", "fit"),
        seed_cv_rows = est_seed_for(seed, f_idx, tt, "Q_exit", "cv_rows"),
        seed_sl_fit = est_seed_for(seed, f_idx, tt, "Q_exit", "sl_fit"),
        sl_workers = sl_workers
      )
      fit_qexit <- sl_qexit$fit
      cn_qexit <- colnames(X_exit)

      X_nat_d <- X_nat_tr_base; X_nat_d$exit_status <- 1
      X_nat_c <- X_nat_tr_base; X_nat_c$exit_status <- 0
      X_shf_d <- X_shf_tr_base; X_shf_d$exit_status <- 1
      X_shf_c <- X_shf_tr_base; X_shf_c$exit_status <- 0

      X_nat_d <- align_cols(X_nat_d, cn_qexit)
      X_nat_c <- align_cols(X_nat_c, cn_qexit)
      X_shf_d <- align_cols(X_shf_d, cn_qexit)
      X_shf_c <- align_cols(X_shf_c, cn_qexit)

      q_death_nat <- scale_info$clip(sl_predict(fit_qexit, X_nat_d))
      q_dc_nat    <- scale_info$clip(sl_predict(fit_qexit, X_nat_c))
      q_death_shf <- scale_info$clip(sl_predict(fit_qexit, X_shf_d))
      q_dc_shf    <- scale_info$clip(sl_predict(fit_qexit, X_shf_c))
    }

    list(
      sl_qexit = sl_qexit,
      fit_qexit = fit_qexit,
      cn_qexit = cn_qexit,
      muY = muY,
      q_death_nat = q_death_nat,
      q_dc_nat = q_dc_nat,
      q_death_shf = q_death_shf,
      q_dc_shf = q_dc_shf
    )
  }

  fit_Q_rem_one <- function() {
    muP <- if (length(rem_idx)) {
      mean(Y_pseudo_ar[rem_idx], na.rm = TRUE)
    } else {
      mean(Y_pseudo_ar, na.rm = TRUE)
    }

    q_rem_nat <- rep(muP, n_ar)
    q_rem_shf <- rep(muP, n_ar)
    fit_qrem <- NULL
    cn_qrem <- NULL
    sl_qrem <- NULL

    if (length(rem_idx) >= 2L) {
      X_rem <- X_nat_tr_base[rem_idx, , drop = FALSE]
      Y_rem <- Y_pseudo_ar[rem_idx]

      ok_to_fit <- if (is_binom && tt == tmax) can_fit_bin(Y_rem)$ok else TRUE

      if (ok_to_fit) {
        sl_qrem <- sl_block_fit(
          Y = Y_rem,
          X = X_rem,
          family = if (is_binom && tt == tmax) stats::binomial() else stats::gaussian(),
          sl_lib = if (!is.null(sl_rec_early) && tt <= rec_transition) {
            sl_rec_early
          } else {
            sl_recursive
          },
          v_inner = max(2L, min(as.integer(inner_v), length(rem_idx))),
          cluster_inner = cl_rem_ar,
          seed_inner = est_seed_for(seed, f_idx, tt, "Q_remain", "fit"),
          seed_cv_rows = est_seed_for(seed, f_idx, tt, "Q_remain", "cv_rows"),
          seed_sl_fit = est_seed_for(seed, f_idx, tt, "Q_remain", "sl_fit"),
          sl_workers = sl_workers
        )
        fit_qrem <- sl_qrem$fit
        cn_qrem <- colnames(X_rem)

        X_nat_r <- align_cols(X_nat_tr_base, cn_qrem)
        X_shf_r <- align_cols(X_shf_tr_base, cn_qrem)

        q_rem_nat <- scale_info$clip(sl_predict(fit_qrem, X_nat_r))
        q_rem_shf <- scale_info$clip(sl_predict(fit_qrem, X_shf_r))
      }
    }

    list(
      sl_qrem = sl_qrem,
      fit_qrem = fit_qrem,
      cn_qrem = cn_qrem,
      muP = muP,
      q_rem_nat = q_rem_nat,
      q_rem_shf = q_rem_shf
    )
  }

  reg_jobs <- list(
    g_remain     = fit_g_remain_one,
    g_death_exit = fit_g_death_exit_one,
    Q_exit       = fit_Q_exit_one,
    Q_rem        = fit_Q_rem_one
  )

  reg_res <- par_lapply(
    X        = reg_jobs,
    FUN      = function(f) f(),
    workers  = reg_workers,
    parallel = parallel,
    seed     = TRUE
  )

  r_rem  <- reg_res$g_remain
  r_dex  <- reg_res$g_death_exit
  r_exit <- reg_res$Q_exit
  r_qrem <- reg_res$Q_rem

  valid_res <- NULL
  if (!is.null(rows_vl) && length(rows_vl) > 0L &&
      !is.null(X_nat_vl_base) && !is.null(X_shf_vl_base)) {
    n_vl <- length(rows_vl)

    if (!is.null(r_rem$fit_rem)) {
      p_rem_nat_vl <- scale_info$clip(sl_predict(r_rem$fit_rem, X_nat_vl_base))
      p_rem_shf_vl <- scale_info$clip(sl_predict(r_rem$fit_rem, X_shf_vl_base))
    } else {
      p_rem_nat_vl <- rep(r_rem$p_rem_const, n_vl)
      p_rem_shf_vl <- rep(r_rem$p_rem_const, n_vl)
    }

    if (!is.null(pooled_pred_vl$p_dex_nat)) {
      p_dex_nat_vl <- pooled_pred_vl$p_dex_nat
      p_dex_shf_vl <- pooled_pred_vl$p_dex_shf
    } else if (!is.null(r_dex$fit_dex)) {
      if (isTRUE(pool_g_death)) {
        X_nat_dex_vl <- make_design(D, rows_vl, cols_pool)
        X_nat_dex_vl <- apply_pool_time_basis(X_nat_dex_vl, tt, pool_time_basis$g_death)
        X_shf_dex_vl <- patch_shifted_design(
          X_nat = X_nat_dex_vl, D_shifted_tt = D_shifted_tt,
          rows = rows_vl, a_names = a_names, tt = tt, k = k, t_min = t_min)
        X_shf_dex_vl <- apply_pool_time_basis(X_shf_dex_vl, tt, pool_time_basis$g_death)
        p_dex_nat_vl <- scale_info$clip(sl_predict(r_dex$fit_dex, align_cols(X_nat_dex_vl, cn_pool)))
        p_dex_shf_vl <- scale_info$clip(sl_predict(r_dex$fit_dex, align_cols(X_shf_dex_vl, cn_pool)))
      } else {
        p_dex_nat_vl <- scale_info$clip(sl_predict(r_dex$fit_dex, X_nat_vl_base))
        p_dex_shf_vl <- scale_info$clip(sl_predict(r_dex$fit_dex, X_shf_vl_base))
      }
    } else {
      p_dex_nat_vl <- rep(r_dex$p_dex_const, n_vl)
      p_dex_shf_vl <- rep(r_dex$p_dex_const, n_vl)
    }

    if (!is.null(pooled_pred_vl$q_death_nat)) {
      q_death_nat_vl <- pooled_pred_vl$q_death_nat
      q_dc_nat_vl    <- pooled_pred_vl$q_dc_nat
      q_death_shf_vl <- pooled_pred_vl$q_death_shf
      q_dc_shf_vl    <- pooled_pred_vl$q_dc_shf
    } else if (!is.null(r_exit$fit_qexit)) {
      if (isTRUE(pool_q_exit)) {
        X_nat_qvl <- make_design(D, rows_vl, cols_pool_qexit)
        X_shf_qvl <- patch_shifted_design(
          X_nat = X_nat_qvl, D_shifted_tt = D_shifted_tt,
          rows = rows_vl, a_names = a_names, tt = tt, k = k, t_min = t_min)
        X_nat_d_vl <- apply_pool_time_basis(X_nat_qvl, tt, pool_time_basis$q_exit); X_nat_d_vl$exit_status <- 1
        X_nat_c_vl <- apply_pool_time_basis(X_nat_qvl, tt, pool_time_basis$q_exit); X_nat_c_vl$exit_status <- 0
        X_shf_d_vl <- apply_pool_time_basis(X_shf_qvl, tt, pool_time_basis$q_exit); X_shf_d_vl$exit_status <- 1
        X_shf_c_vl <- apply_pool_time_basis(X_shf_qvl, tt, pool_time_basis$q_exit); X_shf_c_vl$exit_status <- 0
      } else {
        X_nat_d_vl <- X_nat_vl_base; X_nat_d_vl$exit_status <- 1
        X_nat_c_vl <- X_nat_vl_base; X_nat_c_vl$exit_status <- 0
        X_shf_d_vl <- X_shf_vl_base; X_shf_d_vl$exit_status <- 1
        X_shf_c_vl <- X_shf_vl_base; X_shf_c_vl$exit_status <- 0
      }
      cn_qex <- r_exit$cn_qexit
      q_death_nat_vl <- scale_info$clip(sl_predict(r_exit$fit_qexit, align_cols(X_nat_d_vl, cn_qex)))
      q_dc_nat_vl    <- scale_info$clip(sl_predict(r_exit$fit_qexit, align_cols(X_nat_c_vl, cn_qex)))
      q_death_shf_vl <- scale_info$clip(sl_predict(r_exit$fit_qexit, align_cols(X_shf_d_vl, cn_qex)))
      q_dc_shf_vl    <- scale_info$clip(sl_predict(r_exit$fit_qexit, align_cols(X_shf_c_vl, cn_qex)))
    } else {
      q_death_nat_vl <- rep(r_exit$muY, n_vl)
      q_dc_nat_vl    <- rep(r_exit$muY, n_vl)
      q_death_shf_vl <- rep(r_exit$muY, n_vl)
      q_dc_shf_vl    <- rep(r_exit$muY, n_vl)
    }

    if (!is.null(r_qrem$fit_qrem)) {
      q_rem_nat_vl <- scale_info$clip(sl_predict(r_qrem$fit_qrem, align_cols(X_nat_vl_base, r_qrem$cn_qrem)))
      q_rem_shf_vl <- scale_info$clip(sl_predict(r_qrem$fit_qrem, align_cols(X_shf_vl_base, r_qrem$cn_qrem)))
    } else {
      q_rem_nat_vl <- rep(r_qrem$muP, n_vl)
      q_rem_shf_vl <- rep(r_qrem$muP, n_vl)
    }

    valid_res <- list(
      p_rem_nat   = p_rem_nat_vl,   p_rem_shf   = p_rem_shf_vl,
      p_dex_nat   = p_dex_nat_vl,   p_dex_shf   = p_dex_shf_vl,
      p_dex_const = r_dex$p_dex_const,
      q_death_nat = q_death_nat_vl, q_dc_nat    = q_dc_nat_vl,
      q_death_shf = q_death_shf_vl, q_dc_shf    = q_dc_shf_vl,
      q_rem_nat   = q_rem_nat_vl,   q_rem_shf   = q_rem_shf_vl
    )
  }

  list(
    g_remain     = r_rem,
    g_death_exit = r_dex,
    Q_exit       = r_exit,
    Q_rem        = r_qrem,
    valid        = valid_res
  )
}
