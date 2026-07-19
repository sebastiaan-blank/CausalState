#' Sequential Doubly Robust (SDR) estimator for longitudinal MTPs
#'
#' Estimates the mean counterfactual outcome under a modified treatment policy
#' (MTP) using the sequentially doubly robust estimator of Díaz et al. (2021),
#' which generalises the SDR construction of Luedtke et al. (2017) to the
#' density-ratio / MTP setting. Requires pre-computed density ratio weights
#' from [density_ratio()].
#'
#' @param df A `data.frame` in long format (one row per subject per time point).
#' @param weight_object Output from [density_ratio()], or a `data.frame`
#'   containing at least columns for subject id, time, `Rt_t`
#'   (instantaneous density ratio), and `global_fold` (cross-fitting fold
#'   assignment).
#' @param tmax Integer. Maximum follow-up time (number of time points).
#' @param id Character. Name of the subject identifier column.
#' @param time Character. Name of the integer time column.
#' @param alive Character. Name of the binary alive indicator column
#'   (`1` = alive, `0` = dead).
#' @param in_state Character. Name of the binary active-state indicator column
#'   (`1` = in active state e.g. ICU, `0` = transitioned out).
#' @param y Character. Name of the outcome column.
#' @param baseline Character vector of baseline (time-invariant) covariate
#'   names.
#' @param tv_names Character vector of time-varying covariate names (excluding
#'   treatment).
#' @param a_names Character vector of treatment variable names.
#' @param no_lag_vars Character vector of variables in `tv_names` or `a_names`
#'   that should not be lagged (e.g. already-encoded temporal features).
#' @param policy_names Character vector of column names holding the shifted
#'   treatment values under the MTP (one per treatment variable in `a_names`).
#' @param sl_remain SuperLearner library for the `g_remain` model
#'   (probability of remaining in state at each time point). Required.
#' @param sl_death SuperLearner library for the `g_death_exit` model
#'   (probability of death among exiters). Required.
#' @param sl_recursive SuperLearner library for the recursive Q-remain model.
#'   Required.
#' @param sl_y SuperLearner library for the Q-exit (outcome-at-exit) model.
#'   Required.
#' @param outcome_family `"binomial"` (default) or `"gaussian"`.
#' @param y_bounds Optional numeric vector of length 2 giving `c(min, max)`
#'   for scaling a continuous outcome to `[0, 1]`. When `NULL` the bounds are
#'   inferred from the observed range of `y` in `df` (`min(y)`, `max(y)`).
#'   Predictions outside the scaling range are clipped to `[0, 1]` before
#'   back-transformation, so supply explicit bounds if the outcome distribution
#'   under the MTP may extend beyond the training range. Ignored for
#'   `outcome_family = "binomial"`.
#' @param bounds Numeric. Probability clipping bound for g and Q predictions.
#'   Default `1e-5`.
#' @param trim Quantile used to cap instantaneous density ratios before they
#'   are assembled into the cumulative weight matrix.  Trimming is applied to
#'   the **full** `weights_dt` (all time points computed by [density_ratio()])
#'   before any subsetting to the current `tmax`.  Consequently the trim
#'   threshold is identical whether you call the estimator with `tmax = 2` or
#'   `tmax = 14`, making results directly comparable across horizons that share
#'   the same weight object.  The same `trim` value must be passed to
#'   [sdr_competing()] when reusing the same weight object, so
#'   that the competing-event analysis operates on identically trimmed weights.
#'   Default `0.99`.
#' @param absorb List of [absorb_rule()] objects specifying outcome overrides
#'   at absorbing states.
#' @param policy_spec_fun A function `(D_block, t, a_names)` returning a
#'   `data.table` with the shifted treatment values for time `t`. Used when
#'   the policy cannot be pre-computed as static columns.
#' @param k Integer. Number of lags of time-varying covariates and treatment
#'   to include in models. Default `2`.
#' @param seed Integer random seed. Default `1`.
#' @param parallel Logical. Enable parallel outer cross-fitting via
#'   [parallel::mclapply()]. Default `FALSE`.
#' @param fold_workers Integer number of worker processes for outer folds.
#'   `NULL` uses [parallel::detectCores()].
#' @param reg_workers Integer number of worker processes for within-fold
#'   regression parallelism.
#' @param sl_workers Integer. Workers for parallel learner evaluation within
#'   each SuperLearner call via [SuperLearner::mcSuperLearner()] (fork-based;
#'   Linux/Mac only). Ignored when `parallel = FALSE`. Default `NULL`
#'   (sequential). **Note:** incompatible with `dbarts`/`SL.dbarts` learners —
#'   dbarts' C-level RNG state does not survive forking. Remove dbarts from
#'   `sl_remain`/`sl_death`/`sl_y`/`sl_recursive` when `sl_workers` is set.
#'
#' @section Parallelism:
#'   `fold_workers`, `reg_workers`, and `sl_workers` can be used independently
#'   or together. Using a single level is robust with all learners. Combining
#'   two or more creates nested `mclapply` calls; in that case multi-threaded
#'   or GPU-based learners (e.g. xgboost with CUDA, OpenMP-based methods) may
#'   crash in the child processes and should be avoided or limited to one thread.
#'
#' @param inner_v Integer. Number of inner cross-validation folds for
#'   SuperLearner. Default `5L`.
#' @param cluster Character. Name of a cluster variable for cluster-robust
#'   standard errors (e.g. hospital). If `NULL`, subject id is used.
#' @param cluster_se_only Logical. If `TRUE`, skip fitting and return only
#'   the fold structure for SE computation. Default `FALSE`.
#' @param pool_g_death Logical. If `TRUE`, fit a single pooled
#'   `g_death_exit` model across all time points (with time as a covariate)
#'   rather than a separate model per time point. Default `FALSE`.
#'
#' @return A named list with:
#'   \describe{
#'     \item{`psi`}{Point estimate of `E[Y(d)]` under the MTP.}
#'     \item{`psi_nat`}{Plug-in estimate under the natural course.}
#'     \item{`psi_shf`}{Plug-in estimate under the MTP (equals `psi`).}
#'     \item{`rd`}{Risk difference `psi_shf - psi_nat`.}
#'     \item{`se`}{Standard error from the efficient influence curve.}
#'     \item{`ic`}{Per-subject influence curve values (length n).}
#'     \item{`sl_summary`}{`data.table` of SuperLearner weights per fold/time/component.}
#'     \item{`fold_diag`}{Per-fold diagnostic summaries.}
#'   }
#'
#' @references
#' Díaz I, Williams N, Hoffman KL, Schenck EJ (2021). Nonparametric Causal
#' Effects Based on Longitudinal Modified Treatment Policies. *JASA*
#' 118(542):846–857.
#'
#' Luedtke AR, Sofrygin O, van der Laan MJ, Carone M (2017). Sequential
#' Double Robustness in Right-Censored Longitudinal Models. arXiv:1705.02459.
#'
#' @seealso [density_ratio()], [itmle()], [absorb_rule()]
#'
#' @export
sdr <- function(
    df, weight_object, tmax,
    id, time, alive, in_state, y,
    baseline,
    tv_names = character(0),
    a_names = character(0),
    no_lag_vars = character(0),
    policy_names = character(0),
    sl_remain = NULL,   # SL library for g_remain
    sl_death  = NULL,   # SL library for g_death_exit
    sl_recursive = NULL, # SL fro REM 
    sl_y = NULL, # SL library for Q_exit
    outcome_family = c("binomial", "gaussian"),
    y_bounds = NULL,
    bounds = 1e-5,
    trim = 0.99,
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
    cluster_se_only = FALSE,
    pool_g_death = FALSE              # if TRUE, fit g_death once pooled over time with t as covariate
) {
  stopifnot(requireNamespace("data.table", quietly = TRUE))
  stopifnot(requireNamespace("SuperLearner", quietly = TRUE))

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
  if (is.null(sl_y)) stop("`sl_exit` must be provided (SL library for Q_exit).", call. = FALSE)
  if (is.null(sl_recursive))  stop("`sl_rem` must be provided (SL library for Q_rem).", call. = FALSE)
  
  outcome_family <- match.arg(outcome_family, c("binomial", "gaussian"))
  
  prep <- prepare_long(DT, id, time, y, cluster = cluster)
  D      <- prep$D
  
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
    branch_diag_here <- list()
    
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

      # Q_rem target is observed Y (terminal outcome) for remainers.
      # Calibration is only directly interpretable under a natural course run;
      # under a shifted policy the target is still natural-course Y, not the
      # counterfactual, so treat this metric as a model-quality check only.
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
      
      pseudo_pre_ar_tr <- Y_train[at_risk_tr]
      mean_pseudo_pre_ar <- mean_or_na(pseudo_pre_ar_tr)
      sd_pseudo_pre_ar   <- sd_or_na(pseudo_pre_ar_tr)
      
      dens_tr     <- density_ratios[tr_ids, , drop = FALSE]
      Y_train_new <- eif(
        density_ratios = dens_tr,
        shifted        = shf_train,
        natural        = nat_train,
        time           = tt,
        time_horizon   = tmax
      )
      
      pseudo_post_ar_tr   <- Y_train_new[at_risk_tr]
      mean_pseudo_post_ar <- mean_or_na(pseudo_post_ar_tr)
      sd_pseudo_post_ar   <- sd_or_na(pseudo_post_ar_tr)
      
      delta_vec         <- pseudo_post_ar_tr - pseudo_pre_ar_tr
      mean_delta_eif    <- mean_or_na(delta_vec)
      sd_delta_eif      <- sd_or_na(delta_vec)
      max_abs_delta_eif <- if (any(is.finite(delta_vec))) max(abs(delta_vec), na.rm = TRUE) else NA_real_
      q95_abs_delta_eif <- if (any(is.finite(delta_vec))) stats::quantile(abs(delta_vec), 0.95, na.rm = TRUE, names = FALSE) else NA_real_
      
      qshf_vec      <- shf_train[at_risk_tr, tt]
      corr_vec      <- pseudo_post_ar_tr - qshf_vec
      mean_corr_eif <- mean_or_na(corr_vec)
      sd_corr_eif   <- sd_or_na(corr_vec)
      
      Y_train <- Y_train_new
      
      
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
        
        Q_pre_diff_mean = mean_or_na(Q_shf_ar - Q_nat_ar),
        Q_pre_diff_sd   = sd_or_na(Q_shf_ar - Q_nat_ar),
        Q_pre_diff_min  = min_or_na(Q_shf_ar - Q_nat_ar),
        Q_pre_diff_max  = max_or_na(Q_shf_ar - Q_nat_ar),
        Q_pre_diff_q95_abs = q_or_na(abs(Q_shf_ar - Q_nat_ar), 0.95),
        
        Q_nat_pre_mean = mean_or_na(Q_nat_ar),
        Q_nat_pre_sd   = sd_or_na(Q_nat_ar),
        Q_nat_pre_min  = min_or_na(Q_nat_ar),
        Q_nat_pre_max  = max_or_na(Q_nat_ar),
        
        Q_shf_pre_mean = mean_or_na(Q_shf_ar),
        Q_shf_pre_sd   = sd_or_na(Q_shf_ar),
        Q_shf_pre_min  = min_or_na(Q_shf_ar),
        Q_shf_pre_max  = max_or_na(Q_shf_ar),
        
        pseudo_pre_mean = mean_pseudo_pre_ar,
        pseudo_pre_sd   = sd_pseudo_pre_ar,
        pseudo_pre_min  = min_or_na(pseudo_pre_ar_tr),
        pseudo_pre_max  = max_or_na(pseudo_pre_ar_tr),
        
        pseudo_post_mean = mean_pseudo_post_ar,
        pseudo_post_sd   = sd_pseudo_post_ar,
        pseudo_post_min  = min_or_na(pseudo_post_ar_tr),
        pseudo_post_max  = max_or_na(pseudo_post_ar_tr),
        
        delta_mean = mean_delta_eif,
        delta_sd   = sd_delta_eif,
        delta_q95_abs = q95_abs_delta_eif,
        delta_max_abs = max_abs_delta_eif,
        
        n_post_below_0 = sum(pseudo_post_ar_tr < 0, na.rm = TRUE),
        n_post_above_1 = sum(pseudo_post_ar_tr > 1, na.rm = TRUE),
        
        corr_mean = mean_corr_eif,
        corr_sd   = sd_corr_eif,
        corr_q95_abs = q_or_na(abs(corr_vec), 0.95),
        corr_max_abs = max_or_na(abs(corr_vec))
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
    } else {
      NULL
    }
    
    list(
      fold      = f_idx,
      valid_ids = vl_ids,
      nat_valid = nat_valid,
      shf_valid = shf_valid,
      sl_meta   = if (length(sl_chunks_here)) {
        data.table::rbindlist(sl_chunks_here, use.names = TRUE, fill = TRUE)
      } else NULL,
      fold_diag = fold_diag_dt, 
      branch_diag = branch_diag_dt
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
    if (!is.null(fr$sl_meta)) {
      sl_chunks[[length(sl_chunks) + 1]] <- fr$sl_meta
    }
    if (!is.null(fr$fold_diag)) {
      fold_diag_chunks[[length(fold_diag_chunks) + 1]] <- fr$fold_diag
    }
    
    if (!is.null(fr$branch_diag)) {
      branch_diag_chunks[[length(branch_diag_chunks) + 1L]] <- fr$branch_diag
    }
  }
  
  fold_diag <- if (length(fold_diag_chunks)) {
    data.table::rbindlist(fold_diag_chunks, use.names = TRUE, fill = TRUE)
  } else NULL
  
  branch_diag <- if (length(branch_diag_chunks)) {
    data.table::rbindlist(branch_diag_chunks, use.names = TRUE, fill = TRUE)
  } else {
    NULL
  }
  
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
  
  eif_vec <- eif(
    density_ratios = density_ratios,
    shifted        = pred_shf_all,
    natural        = pred_nat_all,
    time           = 1L,
    time_horizon   = tmax
  )
  
  psi_scaled <- mean(eif_vec, na.rm = TRUE)
  ic_scaled  <- eif_vec - psi_scaled
  
  
  psi_plugin_nat_scaled  <- mean(pred_nat_all[, start_t_global], na.rm = TRUE)
  psi_plugin_shf_scaled  <- mean(pred_shf_all[, start_t_global], na.rm = TRUE)
  psi_plugin_diff_scaled <- psi_plugin_shf_scaled - psi_plugin_nat_scaled
  
  psi_eif_gap_scaled <- psi_scaled - psi_plugin_shf_scaled
  
  psi_natural_scaled <- psi_plugin_nat_scaled
  psi_shifted_scaled <- psi_plugin_shf_scaled
  psi_diff_scaled    <- psi_plugin_diff_scaled
  
  ic_scaled <- as.numeric(ic_scaled)
  ok_ic     <- is.finite(ic_scaled)
  
  if (is.null(cluster_by_id)) {
    n_eff <- sum(ok_ic)
    if (n_eff < 2L) stop("SE: <2 finite IC values.", call. = FALSE)
    se_scaled <- stats::sd(ic_scaled[ok_ic]) / sqrt(n_eff)
  } else {
    cl  <- cluster_by_id
    ok  <- ok_ic & !is.na(cl)
    dt  <- data.table::data.table(cl = cl[ok], ic = ic_scaled[ok])
    n_eff <- nrow(dt)
    if (n_eff < 2L) stop("SE: <2 finite IC values after clustering filter.", call. = FALSE)
    S   <- dt[, .(S = sum(ic)), by = cl]
    G   <- nrow(S)
    if (G < 2L) stop("SE: need >=2 clusters for cluster-robust SE.", call. = FALSE)
    Sbar      <- mean(S$S)
    se_scaled <- sqrt((G / (G - 1)) * sum((S$S - Sbar)^2) / (n_eff^2))
    df        <- G - 1
  }
  
  if (scale_info$bounded) {
    psi         <- scale_info$from_unit(psi_scaled)
    psi_natural <- scale_info$from_unit(psi_natural_scaled)
    psi_shifted <- scale_info$from_unit(psi_shifted_scaled)
    psi_diff    <- psi_shifted - psi_natural
    
    psi_plugin_nat  <- scale_info$from_unit(psi_plugin_nat_scaled)
    psi_plugin_shf  <- scale_info$from_unit(psi_plugin_shf_scaled)
    psi_plugin_diff <- psi_plugin_shf - psi_plugin_nat
    psi_eif_gap     <- scale_info$y_rng * psi_eif_gap_scaled
    
    ic <- scale_info$y_rng * ic_scaled
    se <- scale_info$y_rng * se_scaled
    
    if (!is.null(diag_table) && nrow(diag_table)) {
      diag_table$mean_Q_nat <- scale_info$from_unit(diag_table$mean_Q_nat)
      diag_table$mean_Q_int <- scale_info$from_unit(diag_table$mean_Q_int)
      diag_table$delta      <- scale_info$y_rng * diag_table$delta
      diag_table$mean_Y_obs <- scale_info$from_unit(diag_table$mean_Y_obs)
    }
  } else {
    psi         <- psi_scaled
    psi_natural <- psi_natural_scaled
    psi_shifted <- psi_shifted_scaled
    psi_diff    <- psi_diff_scaled
    
    psi_plugin_nat  <- psi_plugin_nat_scaled
    psi_plugin_shf  <- psi_plugin_shf_scaled
    psi_plugin_diff <- psi_plugin_diff_scaled
    psi_eif_gap     <- psi_eif_gap_scaled
    
    ic <- ic_scaled
    se <- se_scaled
  }
  
  sl_summary <- NULL
  if (length(sl_chunks)) {
    sl_summary <- data.table::rbindlist(sl_chunks, fill = TRUE)
    sl_summary <- sl_summary[order(component, t, fold, -weight)]
  }
  
  out <- list(
    psi = psi,
    se = se,
    ic_df = data.table::data.table(id = ids, ic = ic),
    start_t = start_t_global,

    predictions = list(
      natural = pred_nat_all,
      shifted = pred_shf_all
    ),
    
    diagnostics = list(
      diag_table      = diag_table,      # cross-fitted trajectory
      recursion_diag  = fold_diag,       # fold/time SDR mechanics
      branch_cal      = branch_diag,     # branch calibration
      sl_summary = sl_summary
    ),
    
    decomposition = list(
      psi_eif         = psi,
      psi_plugin_nat  = psi_plugin_nat,
      psi_plugin_shf  = psi_plugin_shf,
      psi_plugin_diff = psi_plugin_diff,
      psi_eif_gap     = psi_eif_gap
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
  
  class(out) <- c("sdr_fit", "list")
  out
}


#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#







sdr_competing <- function(
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
    inner_v  = 5L,
    bounds   = 1e-5,
    trim     = 0.99,
    parallel     = FALSE,
    fold_workers = NULL,
    sl_workers   = NULL,
    cluster_se_only = FALSE
) {
  stopifnot(requireNamespace("data.table", quietly = TRUE))

  if (!isTRUE(parallel)) sl_workers <- NULL

  set.seed(seed)

  seed_for <- function(seed, f_idx, tt, component_code) {
    as.integer(seed + 1000L * f_idx + 10L * tt + component_code)
  }
  
  prep <- prepare_long(DT, id = id, time = time, y = alive, cluster = cluster)
  D    <- prep$D
  D <- D[.__time >= 1L & .__time <= tmax]
  data.table::setorderv(D, c(".__id", ".__time"))
  
  D <- expand_to_horizon(
    D, id = ".__id", time = ".__time",
    alive = alive, in_state = in_state, tmax = tmax
  )
  
  lag_vars <- unique(c(tv_names, a_names))
  if (k > 0L && length(lag_vars)) {
    D <- make_lags(D, vars = lag_vars, k = k, id = ".__id", time = ".__time", overwrite = FALSE)
  }
  
  ids <- unique(D$.__id)
  N   <- length(ids)
  row_index <- build_id_time_index(D, ids = ids, tmax = tmax)
  
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
    cl_dt <- D[!is.na(.__cl), .(cl = data.table::first(na.omit(.__cl))), by = .__id]
    cluster_by_id <- cl_dt$cl[match(ids, cl_dt$.__id)]
  }
  
  fold_info <- get_folds_from_weights(
    weights_dt    = weights_dt,
    id            = id,
    ids           = ids,
    cluster_by_id = cluster_by_id,
    cluster_se_only = cluster_se_only
  )
  folds <- fold_info$folds
  
  D_shifted_by_t <- vector("list", tmax)
  if (length(a_names)) {
    for (tt in seq_len(tmax)) {
      D_shifted_by_t[[tt]] <- overwrite_policy_history_for_Q(
        D       = D,
        t       = tt,
        a_names = a_names,
        policy_spec_fun = policy_spec_fun,
        id      = id,
        time    = time,
        k       = k
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
  Y_by_id <- died_by_id$Y[match(ids, died_by_id$.__id)]
  stopifnot(all(Y_by_id %in% c(0L, 1L)))
  
  scale_info <- scale_y(
    D = data.table::data.table(tmp_y = Y_by_id),
    y = "tmp_y",
    outcome_family = "binomial",
    bounds = bounds
  )
  
  pred_nat_all <- matrix(NA_real_, nrow = N, ncol = tmax + 1L)
  pred_shf_all <- matrix(NA_real_, nrow = N, ncol = tmax + 1L)
  
  pred_nat_all[, tmax + 1L] <- Y_by_id
  pred_shf_all[, tmax + 1L] <- Y_by_id
  
  sl_chunks <- list()
  fold_diag_chunks <- list()
  
  fold_worker <- function(f_idx) {
    set.seed(seed + f_idx)
    
    tr <- folds[[f_idx]]$training_set
    vl <- folds[[f_idx]]$validation_set
    
    n_tr <- length(tr)
    n_vl <- length(vl)
    
    pred_nat_tr <- matrix(NA_real_, nrow = n_tr, ncol = tmax + 1L)
    pred_shf_tr <- matrix(NA_real_, nrow = n_tr, ncol = tmax + 1L)
    pred_nat_vl <- matrix(NA_real_, nrow = n_vl, ncol = tmax + 1L)
    pred_shf_vl <- matrix(NA_real_, nrow = n_vl, ncol = tmax + 1L)
    
    pred_nat_tr[, tmax + 1L] <- Y_by_id[tr]
    pred_shf_tr[, tmax + 1L] <- Y_by_id[tr]
    pred_nat_vl[, tmax + 1L] <- Y_by_id[vl]
    pred_shf_vl[, tmax + 1L] <- Y_by_id[vl]
    
    dens_tr <- dens[tr, , drop = FALSE]
    
    Y_tr <- as.numeric(Y_by_id[tr])
    
    sl_chunks_here <- list()
    fold_diag_here <- list()
    
    for (tt in rev(seq_len(tmax))) {
      
      c1_tr <- !is.na(row_index[tr, tt])
      c1_vl <- !is.na(row_index[vl, tt])
      
      if (tt == 1L) {
        y1_tr <- rep(TRUE, n_tr)
        d0_tr <- rep(TRUE, n_tr)
        y1_vl <- rep(TRUE, n_vl)
        d0_vl <- rep(TRUE, n_vl)
      } else {
        prev_tr <- row_index[tr, tt - 1L]
        okp_tr  <- !is.na(prev_tr)
        alive_prev_tr <- integer(n_tr)
        in_state_prev_tr   <- integer(n_tr)
        alive_prev_tr[okp_tr] <- as.integer(D[[alive]][prev_tr[okp_tr]])
        in_state_prev_tr[okp_tr]   <- as.integer(D[[in_state]][prev_tr[okp_tr]])
        y1_tr <- okp_tr & (alive_prev_tr == 1L)
        d0_tr <- okp_tr & (in_state_prev_tr == 1L)
        
        prev_vl <- row_index[vl, tt - 1L]
        okp_vl  <- !is.na(prev_vl)
        alive_prev_vl <- integer(n_vl)
        in_state_prev_vl   <- integer(n_vl)
        alive_prev_vl[okp_vl] <- as.integer(D[[alive]][prev_vl[okp_vl]])
        in_state_prev_vl[okp_vl]   <- as.integer(D[[in_state]][prev_vl[okp_vl]])
        y1_vl <- okp_vl & (alive_prev_vl == 1L)
        d0_vl <- okp_vl & (in_state_prev_vl == 1L)
      }
      
      i_tr <- c1_tr & y1_tr & d0_tr
      i_vl <- c1_vl & y1_vl & d0_vl
      
      .apply_gating <- function(pred_mat, y1, d0, tt_col) {
        pred_mat[!y1, tt_col]          <- 1
        pred_mat[y1 & !d0, tt_col]     <- 0
        still_na <- is.na(pred_mat[, tt_col])
        if (any(still_na)) {
          pred_mat[still_na, tt_col] <- pred_mat[still_na, tt_col + 1L]
        }
        pred_mat
      }
      
      if (!any(i_tr)) {
        pred_nat_tr <- .apply_gating(pred_nat_tr, y1_tr, d0_tr, tt)
        pred_shf_tr <- .apply_gating(pred_shf_tr, y1_tr, d0_tr, tt)
        pred_nat_vl <- .apply_gating(pred_nat_vl, y1_vl, d0_vl, tt)
        pred_shf_vl <- .apply_gating(pred_shf_vl, y1_vl, d0_vl, tt)
        
        Y_tr <- eif(
          density_ratios = dens_tr,
          shifted        = pred_shf_tr,
          natural        = pred_nat_tr,
          time           = tt,
          time_horizon   = tmax
        )
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
      
      cl_tr_ar <- if (!is.null(cluster_by_id)) {
        cluster_by_id[tr[i_tr]]
      } else NULL
      
      id_tr_ar <- ids[tr[i_tr]]
      
      fam <- if (tt == tmax) stats::binomial() else stats::gaussian()
      
      sl_fit <- sl_block_fit(
        Y             = Y_tr[i_tr],
        X             = X_nat_tr,
        family        = fam,
        sl_lib        = sl_q,
        v_inner       = max(2L, min(as.integer(inner_v), sum(i_tr))),
        id_inner      = if (is.null(cl_tr_ar)) id_tr_ar else NULL,
        cluster_inner = cl_tr_ar,
        seed_inner    = seed_for(seed, f_idx, tt, 1L),
        sl_workers    = sl_workers
      )
      fit <- sl_fit$fit

      n_val_tt <- sum(i_vl)
      meta <- sl_meta(sl_fit, f_idx, tt, "Q", sum(i_tr), n_val_tt)
      if (!is.null(meta)) sl_chunks_here[[length(sl_chunks_here) + 1L]] <- meta
      
      cn_fit <- colnames(X_nat_tr)
      
      q_nat_tr <- scale_info$clip(as.numeric(sl_predict(fit, X_nat_tr)))
      q_shf_tr <- scale_info$clip(as.numeric(sl_predict(fit, X_shf_tr)))
      
      pred_nat_tr[i_tr, tt] <- q_nat_tr
      pred_shf_tr[i_tr, tt] <- q_shf_tr
      
      pred_nat_tr <- .apply_gating(pred_nat_tr, y1_tr, d0_tr, tt)
      pred_shf_tr <- .apply_gating(pred_shf_tr, y1_tr, d0_tr, tt)
      
      if (any(i_vl)) {
        rows_vl_ar <- row_index[vl[i_vl], tt]
        
        X_nat_vl <- make_design(D, rows_vl_ar, cols_tt)
        X_shf_vl <- patch_shifted_design(
          X_nat        = X_nat_vl,
          D_shifted_tt = D_shifted_by_t[[tt]],
          rows         = rows_vl_ar,
          a_names      = a_names,
          tt           = tt,
          k            = k
        )
        
        X_nat_vl <- align_cols(X_nat_vl, cn_fit)
        X_shf_vl <- align_cols(X_shf_vl, cn_fit)
        
        q_nat_vl <- scale_info$clip(as.numeric(sl_predict(fit, X_nat_vl)))
        q_shf_vl <- scale_info$clip(as.numeric(sl_predict(fit, X_shf_vl)))
        
        pred_nat_vl[i_vl, tt] <- q_nat_vl
        pred_shf_vl[i_vl, tt] <- q_shf_vl
      }
      
      pred_nat_vl <- .apply_gating(pred_nat_vl, y1_vl, d0_vl, tt)
      pred_shf_vl <- .apply_gating(pred_shf_vl, y1_vl, d0_vl, tt)
      
      Y_tr <- eif(
        density_ratios = dens_tr,
        shifted        = pred_shf_tr,
        natural        = pred_nat_tr,
        time           = tt,
        time_horizon   = tmax
      )
      
      diag_row <- data.table::data.table(
        fold            = f_idx,
        t               = tt,
        n_train_at_risk = sum(i_tr),
        n_valid_at_risk = sum(i_vl),
        mean_Y_pseudo   = mean(Y_tr[i_tr], na.rm = TRUE),
        mean_Q_nat      = mean(q_nat_tr, na.rm = TRUE),
        mean_Q_shf      = mean(q_shf_tr, na.rm = TRUE),
        used_const       = isTRUE(sl_fit$used_const)
      )
      fold_diag_here[[length(fold_diag_here) + 1L]] <- diag_row
      
    } # end tt loop
    
    fold_diag_dt <- if (length(fold_diag_here)) {
      data.table::rbindlist(fold_diag_here, use.names = TRUE, fill = TRUE)
    } else NULL
    
    list(
      f_idx       = f_idx,
      vl          = vl,
      pred_nat_vl = pred_nat_vl,
      pred_shf_vl = pred_shf_vl,
      sl_meta     = if (length(sl_chunks_here)) {
        data.table::rbindlist(sl_chunks_here, use.names = TRUE, fill = TRUE)
      } else NULL,
      fold_diag   = fold_diag_dt
    )
  } # end fold_worker
  
  res_by_fold <- par_lapply(
    X        = seq_along(folds),
    FUN      = fold_worker,
    workers  = fold_workers,
    parallel = parallel,
    seed     = TRUE
  )
  
  for (res in res_by_fold) {
    vl <- res$vl
    pred_nat_all[vl, ] <- res$pred_nat_vl
    pred_shf_all[vl, ] <- res$pred_shf_vl
    
    if (!is.null(res$sl_meta)) {
      sl_chunks[[length(sl_chunks) + 1L]] <- res$sl_meta
    }
    if (!is.null(res$fold_diag)) {
      fold_diag_chunks[[length(fold_diag_chunks) + 1L]] <- res$fold_diag
    }
  }
  
  fold_diag <- if (length(fold_diag_chunks)) {
    data.table::rbindlist(fold_diag_chunks, use.names = TRUE, fill = TRUE)
  } else NULL
  
  sl_summary <- if (length(sl_chunks)) {
    data.table::rbindlist(sl_chunks, fill = TRUE)
  } else NULL
  
  eif_vec <- eif(
    density_ratios = dens,
    shifted        = pred_shf_all,
    natural        = pred_nat_all,
    time           = 1L,
    time_horizon   = tmax
  )
  
  psi     <- mean(eif_vec, na.rm = TRUE)
  ic      <- eif_vec - psi
  ok_ic   <- is.finite(ic)
  
  if (!is.null(cluster_by_id)) {
    dt_ic <- data.table::data.table(
      .__id = ids[ok_ic],
      ic    = ic[ok_ic],
      cl    = cluster_by_id[ok_ic]
    )
    se <- cluster_se(dt_ic, id = ".__id", cluster_col = "cl")
  } else {
    n_eff <- sum(ok_ic)
    se <- stats::sd(ic[ok_ic]) / sqrt(n_eff)
  }
  
  psi_nat <- mean(pred_nat_all[, 1L], na.rm = TRUE)
  psi_shf <- mean(pred_shf_all[, 1L], na.rm = TRUE)
  
  list(
    psi         = psi,
    psi_nat     = psi_nat,
    psi_shf     = psi_shf,
    rd          = psi_shf - psi_nat,
    se          = se,
    ic_df       = data.table::data.table(id = ids, ic = ic),
    predictions = list(
      natural = pred_nat_all,
      shifted = pred_shf_all
    ),
    sl_summary  = sl_summary,
    fold_diag   = fold_diag
  )
}




