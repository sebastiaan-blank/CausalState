
make_subject_folds_qreg <- function(ids, v, cluster_by_id, seed) {
  cl  <- as.character(cluster_by_id)
  fvec <- make_group_folds(cl, v = v, seed = seed)
  fold_ids <- sort(unique(fvec))
  lapply(fold_ids, function(k) {
    list(training_set   = which(fvec != k),
         validation_set = which(fvec == k))
  })
}


#' Pure Q-recursion estimator -- diagnostic use only
#'
#' Fits sequential Q-models backwards through time and returns the
#' plug-in (substitution) estimator \eqn{\hat\Psi = n^{-1}\sum_i \hat Q_1(H_{i1})}.
#' No EIF update is applied, so the estimate carries first-order bias and
#' \strong{should not be used as a real causal estimate}.
#'
#' The intended use is to pass the \strong{natural-course} (identity) policy --
#' that is, no shift -- and compare the resulting estimate to the observed mean
#' of \code{y}.  Close agreement indicates that the Q-models are
#' well-calibrated under the natural course, which is a prerequisite for the
#' doubly-robust estimates from [sdr()] or [itmle()] to be trustworthy.  Poor
#' agreement signals Q-model misspecification.
#'
#' \strong{Do not pass a shifted policy to \code{qreg} and interpret the
#' result as a causal estimate.}  For policy evaluation use [sdr()] or
#' [itmle()].
#'
#' The reported \code{se_naive} is \eqn{s(\hat Q_1) / \sqrt{n}}, which treats
#' the fitted Q as fixed and should be read only as a rough guide to Monte
#' Carlo variability of the plug-in.
#'
#' @inheritParams sdr
#' @param v Integer. Number of cross-fitting folds. Default \code{5L}.
#' @param n_boot Integer. Number of bootstrap replicates for SE estimation.
#'   \code{0L} (default) returns only the naive SE (\eqn{s(\hat Q_1)/\sqrt{n}}).
#'   When \code{> 0}, subjects (or clusters when \code{cluster} is set) are
#'   resampled with replacement, the full Q-recursion is refit for each
#'   replicate, and \code{se = sd(psi_boot)}.  The bootstrap distribution is
#'   stored in \code{decomposition$boot_psi_dist}.
#'
#' @return A list with:
#'   \describe{
#'     \item{\code{psi}}{Point estimate \eqn{\hat\Psi}.}
#'     \item{\code{se}}{Best available SE: bootstrap SE if \code{n_boot > 0},
#'       otherwise naive SE.}
#'     \item{\code{ci}}{95% Wald CI using \code{se}.}
#'     \item{\code{Y_obs}}{Observed mean of \code{y} (natural-course reference).}
#'     \item{\code{decomposition$se_naive}}{Naive SE \eqn{s(\hat Q_1)/\sqrt{n}}.}
#'     \item{\code{decomposition$se_boot}}{Bootstrap SE (\code{NULL} if
#'       \code{n_boot = 0}).}
#'     \item{\code{decomposition$boot_psi_dist}}{Vector of per-replicate psi
#'       values (\code{NULL} if \code{n_boot = 0}).}
#'     \item{\code{sl_summary}}{`data.table` of SuperLearner learner weights per
#'       (fold, time-point, model component). Same structure as in [sdr()].}
#'     \item{\code{diagnostics$branch_cal}}{Per-fold, per-time branch calibration
#'       table: empirical mean targets vs. predictions for `g_remain`, `g_death`,
#'       and `Q_remain`. Use this to check whether the Q-models are
#'       well-calibrated under the natural course -- the primary diagnostic
#'       purpose of [qreg()].}
#'     \item{\code{diagnostics$diag_table}}{Additional per-fold, per-time
#'       summary statistics (Q prediction ranges, pseudo-outcome statistics).}
#'   }
#'
#' @section Sample size requirements:
#'   This package targets large longitudinal datasets -- thousands of subjects
#'   -- typical of ICU, ward, or emergency department cohorts.  The Q-models
#'   require adequate observations within each branch (alive, in-state, exited)
#'   at every time point to be stable.  `pool_g_death` and `pool_q_exit` borrow
#'   strength across time steps when exit events are sparse, but there must
#'   still be sufficient events across the pooled structure.  All built-in
#'   examples use a minimum of 2,000 subjects.
#'
#' @seealso [sdr()], [itmle()], [density_ratio()]
#'
#' @examples
#' \donttest{
#' library(SuperLearner)
#'
#' # ICU-like DGP: patients die, discharge, or remain in state each time point
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
#'
#' # Natural-course policy (no shift) -- the correct use of qreg
#' policy_nat <- function(D_block, t, a_names) D_block[, ..a_names, drop = FALSE]
#'
#' sl_lib <- c("SL.mean", "SL.glm")
#'
#' # Q-model calibration check: estimate should be close to mean(df$Y)
#' res <- qreg(
#'   df              = df,
#'   tmax            = 5L,
#'   id              = "id",
#'   time            = "time",
#'   alive           = "alive",
#'   in_state        = "in_state",
#'   y               = "Y",
#'   baseline        = "age",
#'   tv_names        = "L1",
#'   a_names         = "A",
#'   sl_remain       = sl_lib,
#'   sl_death        = sl_lib,
#'   sl_recursive    = sl_lib,
#'   sl_y            = sl_lib,
#'   k               = 1L,
#'   inner_v         = 3L,
#'   parallel        = FALSE,
#'   seed            = 1L,
#'   policy_spec_fun = policy_nat
#' )
#'
#' # Compare plug-in to observed mean: close agreement = well-calibrated Q
#' res$estimate
#' mean(df$Y)
#' }
#'
#' # Multi-treatment: plug-in estimate under a joint policy over two treatments
#' \donttest{
#' library(SuperLearner)
#'
#' sim_multi <- function(n = 2000L, tmax = 5L) {
#'   set.seed(42L)
#'   rows <- vector("list", n)
#'   for (i in seq_len(n)) {
#'     age <- round(rnorm(1, 65, 10)); L1 <- rnorm(1)
#'     pat <- list()
#'     for (t in seq_len(tmax)) {
#'       A1    <- rbinom(1, 1, plogis(0.3 * L1 - 0.4))
#'       A2    <- rbinom(1, 1, plogis(0.2 * L1 + 0.3 * A1 - 0.3))
#'       u     <- runif(1)
#'       p_die <- plogis(-4.0 + 0.2 * L1 - 0.1 * age / 10)
#'       p_dc  <- plogis(-2.5 + 0.4 * A1 + 0.3 * A2)
#'       if (u < p_die) {
#'         alive <- 0L; in_state <- 0L
#'       } else if (u < p_die + p_dc) {
#'         alive <- 1L; in_state <- 0L
#'       } else {
#'         alive <- 1L; in_state <- 1L
#'       }
#'       Y <- rbinom(1, 1, plogis(-0.5 + 0.3 * A1 + 0.3 * A2 - 0.2 * L1))
#'       pat[[length(pat) + 1L]] <- data.frame(
#'         id = i, time = t, age = age, L1 = L1,
#'         A1 = A1, A2 = A2, alive = alive, in_state = in_state, Y = Y
#'       )
#'       if (in_state == 0L) break
#'       if (t < tmax) L1 <- L1 + rnorm(1, -0.1 * A1, 0.3)
#'     }
#'     rows[[i]] <- do.call(rbind, pat)
#'   }
#'   do.call(rbind, rows)
#' }
#' df2 <- sim_multi()
#'
#' policy_fn2 <- function(D_block, t, a_names) {
#'   out <- D_block[, ..a_names, drop = FALSE]
#'   out[[a_names[1]]] <- pmin(D_block[[a_names[1]]] + 0.3, 1)
#'   out[[a_names[2]]] <- pmin(D_block[[a_names[2]]] + 0.3, 1)
#'   out
#' }
#'
#' sl_lib <- c("SL.mean", "SL.glm")
#'
#' res2 <- qreg(
#'   df              = df2,
#'   tmax            = 5L,
#'   id              = "id",
#'   time            = "time",
#'   alive           = "alive",
#'   in_state        = "in_state",
#'   y               = "Y",
#'   baseline        = "age",
#'   tv_names        = "L1",
#'   a_names         = c("A1", "A2"),
#'   sl_remain       = sl_lib,
#'   sl_death        = sl_lib,
#'   sl_recursive    = sl_lib,
#'   sl_y            = sl_lib,
#'   k               = 1L,
#'   inner_v         = 3L,
#'   parallel        = FALSE,
#'   seed            = 1L,
#'   policy_spec_fun = policy_fn2
#' )
#' res2$estimate
#' res2$se
#' }
#'
#' @export
qreg <- function(
    df, tmax,
    id, time, alive, in_state, y,
    baseline,
    tv_names        = character(0),
    a_names         = character(0),
    no_lag_vars     = character(0),
    policy_names    = character(0),
    sl_remain       = NULL,
    sl_death        = NULL,
    sl_recursive    = NULL,
    sl_rec_early    = NULL,
    rec_transition  = NULL,
    sl_y            = NULL,
    outcome_family  = c("binomial", "gaussian"),
    y_bounds        = NULL,
    bounds          = 1e-5,
    absorb          = list(),
    policy_spec_fun = function(D_block, t, a_names) NULL,
    k               = 2,
    seed            = 1,
    parallel        = FALSE,
    fold_workers    = NULL,
    reg_workers     = NULL,
    sl_workers      = NULL,
    inner_v         = 5L,
    cluster         = NULL,
    pool_g_death    = FALSE,
    pool_q_exit     = FALSE,
    pool_time       = "spline",
    v               = 5L,
    n_boot          = 0L,
    verbose         = TRUE
) {
  stopifnot(requireNamespace("data.table", quietly = TRUE))
  stopifnot(requireNamespace("SuperLearner", quietly = TRUE))

  pool_time_basis <- make_pool_time_basis_list(pool_time, tmax, pool_g_death, pool_q_exit)

  if (!isTRUE(parallel)) {
    fold_workers <- NULL
    reg_workers  <- NULL
    sl_workers   <- NULL
  }

  keep_cols <- unique(c(
    id, time, a_names, y, alive, in_state,
    baseline, tv_names, no_lag_vars, policy_names, cluster
  ))
  keep_cols <- intersect(keep_cols, names(df))
  df <- as.data.frame(df)[, keep_cols, drop = FALSE]

  data_check(data = df, id = id, time = time,
             a_names = a_names, tv_names = tv_names,
             bs = baseline, y_col = y, alive_col = alive,
             in_state_col = in_state, tmin = 1, tmax = tmax,
             verbose = verbose)

  if (is.null(sl_remain))    stop("`sl_remain` must be provided (SL library for g_remain).", call. = FALSE)
  if (is.null(sl_death))     stop("`sl_death` must be provided (SL library for g_death_exit).", call. = FALSE)
  if (is.null(sl_y))         stop("`sl_y` must be provided.", call. = FALSE)
  if (is.null(sl_recursive)) stop("`sl_recursive` must be provided.", call. = FALSE)
  if (!is.null(sl_rec_early) && is.null(rec_transition))
    stop("`sl_rec_early` requires `rec_transition` to be specified.", call. = FALSE)

  set.seed(seed)

  DT <- data.table::as.data.table(df)
  data.table::setorderv(DT, c(id, time))
  rm(df); gc()
  t_min <- DT[, .SD[1L], by = id][[time]] |> max()

  outcome_family <- match.arg(outcome_family, c("binomial", "gaussian"))

  prep    <- prepare_long(DT, id, time, y, cluster = cluster)
  D       <- prep$D

  hist_out <- prepare_history_lags(
    DT          = D,
    a_names     = a_names,
    tv_names    = tv_names,
    k           = k,
    id          = ".__id",
    no_lag_vars = no_lag_vars
  )
  D <- hist_out$DT

  scale_info <- scale_y(D, y, outcome_family, y_bounds, bounds)
  is_binom   <- scale_info$is_binom

  Y_info  <- build_Y_by_id(D, scale_info)
  ids     <- Y_info$ids
  Y_init  <- Y_info$Y
  N       <- length(ids)

  row_index <- build_id_time_index(D, ids, tmax)

  if (!is.null(cluster)) {
    cl_by_id      <- D[, .(cl = data.table::first(na.omit(.__cl))), by = .__id]
    cluster_by_id <- cl_by_id$cl[match(ids, cl_by_id$.__id)]
  } else {
    cluster_by_id <- ids
  }

  folds <- make_subject_folds_qreg(ids, v = v,
                                    cluster_by_id = cluster_by_id,
                                    seed = seed)

  all_names  <- names(D)
  cols_by_t  <- lapply(seq_len(tmax), function(tt) {
    make_cols(tt = tt, baseline = baseline, tv_names = tv_names,
              a_names = a_names, k = k, all_names = all_names, t_min = t_min)
  })

  D_shifted <- vector("list", length = tmax)
  if (length(a_names)) {
    for (tt in seq_len(tmax)) {
      D_shifted[[tt]] <- overwrite_policy_history_for_Q(
        D, t = tt, a_names = a_names, policy_spec_fun = policy_spec_fun,
        id = id, time = time, k = k, t_min = t_min
      )
    }
  } else {
    for (tt in seq_len(tmax)) D_shifted[[tt]] <- D
  }

  pred_nat_all <- matrix(NA_real_, nrow = N, ncol = tmax + 1L)
  pred_shf_all <- matrix(NA_real_, nrow = N, ncol = tmax + 1L)
  pred_nat_all[, tmax + 1L] <- Y_init
  pred_shf_all[, tmax + 1L] <- Y_init

  sl_chunks <- list()

  fold_worker <- function(f_idx, folds, Y_init, row_index, cluster_by_id,
                          verbose_fold = verbose) {

    is_binom <- scale_info$is_binom

    fold   <- folds[[f_idx]]
    tr_ids <- fold$training_set
    vl_ids <- fold$validation_set

    n_tr <- length(tr_ids)
    n_vl <- length(vl_ids)

    nat_train <- matrix(NA_real_, nrow = n_tr, ncol = tmax + 1L)
    shf_train <- matrix(NA_real_, nrow = n_tr, ncol = tmax + 1L)
    nat_valid <- matrix(NA_real_, nrow = n_vl, ncol = tmax + 1L)
    shf_valid <- matrix(NA_real_, nrow = n_vl, ncol = tmax + 1L)

    nat_train[, tmax + 1L] <- Y_init[tr_ids]
    shf_train[, tmax + 1L] <- Y_init[tr_ids]
    nat_valid[, tmax + 1L] <- Y_init[vl_ids]
    shf_valid[, tmax + 1L] <- Y_init[vl_ids]

    sl_chunks_here   <- list()
    fold_diag_here   <- list()
    branch_diag_here <- list()

    .fit_pool_g <- function() fit_pooled_g_death_exit(
      D = D, tr_ids = tr_ids, row_index = row_index, tmax = tmax,
      baseline = baseline, tv_names = tv_names, a_names = a_names,
      all_names = all_names, k = k, t_min = t_min,
      alive = alive, in_state = in_state,
      cluster_by_id = cluster_by_id,
      sl_death = sl_death, inner_v = inner_v, seed = seed, f_idx = f_idx,
      pool_time_basis = pool_time_basis$g_death,
      sl_workers = sl_workers
    )
    .fit_pool_q <- function() fit_pooled_q_exit(
      D = D, tr_ids = tr_ids, row_index = row_index, tmax = tmax,
      baseline = baseline, tv_names = tv_names, a_names = a_names,
      all_names = all_names, k = k, t_min = t_min,
      alive = alive, in_state = in_state, Y_init = Y_init,
      cluster_by_id = cluster_by_id,
      sl_y = sl_y, inner_v = inner_v, seed = seed, f_idx = f_idx,
      is_binom = is_binom,
      pool_time_basis = pool_time_basis$q_exit,
      sl_workers = sl_workers
    )

    run_both_pooled <- isTRUE(pool_g_death) && isTRUE(pool_q_exit) &&
                       isTRUE(parallel) && isTRUE(reg_workers > 1L)

    if (run_both_pooled) {
      pool_res   <- par_lapply(
        X        = list(g_death = .fit_pool_g, q_exit = .fit_pool_q),
        FUN      = function(f) f(),
        workers  = reg_workers,
        parallel = TRUE,
        seed     = TRUE
      )
      pooled_dex <- pool_res$g_death
      pooled_qex <- pool_res$q_exit
    } else {
      pooled_dex <- if (isTRUE(pool_g_death)) .fit_pool_g() else NULL
      pooled_qex <- if (isTRUE(pool_q_exit))  .fit_pool_q() else NULL
    }

    if (isTRUE(pool_g_death)) {
      if (is.null(pooled_dex$fit)) {
        message(sprintf(
          "[fold %d] pool_g_death=TRUE but pooled g_death_exit failed (%s); per-time fallback",
          f_idx, pooled_dex$reason))
        fit_dex_pooled <- sl_dex_pooled <- cn_pool <- cols_pool <- all_Y_pool <- NULL
        pool_g_death <- FALSE
      } else {
        fit_dex_pooled <- pooled_dex$fit
        sl_dex_pooled  <- pooled_dex$sl
        cn_pool        <- pooled_dex$cols
        cols_pool      <- pooled_dex$cols_pool
        all_Y_pool     <- pooled_dex$Y
      }
    } else {
      fit_dex_pooled <- sl_dex_pooled <- cn_pool <- cols_pool <- all_Y_pool <- NULL
    }

    if (isTRUE(pool_q_exit)) {
      if (is.null(pooled_qex$fit)) {
        message(sprintf(
          "[fold %d] pool_q_exit=TRUE but pooled Q_exit failed (%s); per-time fallback",
          f_idx, pooled_qex$reason))
        fit_qexit_pooled <- sl_qexit_pooled <- cn_qexit_pool <- cols_pool_qexit <- all_Y_qexit_pool <- NULL
        pool_q_exit <- FALSE
      } else {
        fit_qexit_pooled  <- pooled_qex$fit
        sl_qexit_pooled   <- pooled_qex$sl
        cn_qexit_pool     <- pooled_qex$cols
        cols_pool_qexit   <- pooled_qex$cols_pool
        all_Y_qexit_pool  <- pooled_qex$Y
      }
    } else {
      fit_qexit_pooled <- sl_qexit_pooled <- cn_qexit_pool <- cols_pool_qexit <- all_Y_qexit_pool <- NULL
    }

    pooled_dex_used   <- isTRUE(pool_g_death) && !is.null(fit_dex_pooled)
    pooled_qexit_used <- isTRUE(pool_q_exit)  && !is.null(fit_qexit_pooled)

    if (pooled_dex_used || pooled_qexit_used) {
      pooled_pred_cache <- build_pooled_pred_cache(
        tmax             = tmax,
        tr_ids           = tr_ids,
        vl_ids           = vl_ids,
        row_index        = row_index,
        D                = D,
        D_shifted        = D_shifted,
        pool_g_death     = pool_g_death,
        fit_dex_pooled   = fit_dex_pooled,
        cols_pool        = cols_pool,
        cn_pool          = cn_pool,
        pool_q_exit      = pool_q_exit,
        fit_qexit_pooled = fit_qexit_pooled,
        cols_pool_qexit  = cols_pool_qexit,
        cn_qexit_pool    = cn_qexit_pool,
        pool_time_basis  = pool_time_basis,
        a_names          = a_names,
        k                = k,
        t_min            = t_min,
        scale_info       = scale_info
      )
      sl_dex_pooled   <- slim_sl(sl_dex_pooled)
      sl_qexit_pooled <- slim_sl(sl_qexit_pooled)
      rm(fit_dex_pooled, fit_qexit_pooled)
      gc(FALSE)
      fit_dex_pooled   <- NULL
      fit_qexit_pooled <- NULL
    } else {
      pooled_pred_cache <- NULL
    }

    for (tt in rev(seq_len(tmax))) {

      if (verbose_fold) message(sprintf("[qreg][fold %d][t=%d][pid=%d] %s",
                                   f_idx, tt, Sys.getpid(), format(Sys.time(), "%H:%M:%S")))

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
      id_tr_ar  <- tr_ids[idx_tr_ar]
      id_vl_ar  <- vl_ids[idx_vl_ar]
      rows_tr   <- row_index[id_tr_ar, tt]
      cols_base <- cols_by_t[[tt]]

      X_nat_tr_base <- make_design(D, rows_tr, cols_base)
      X_shf_tr_base <- patch_shifted_design(
        X_nat = X_nat_tr_base, D_shifted_tt = D_shifted[[tt]],
        rows = rows_tr, a_names = a_names, tt = tt, k = k, t_min = t_min
      )
      X_nat_tr_haz <- X_nat_tr_base
      X_shf_tr_haz <- X_shf_tr_base

      alive_end  <- as.integer(D[[alive]][rows_tr])
      in_state_end <- as.integer(D[[in_state]][rows_tr])
      D_tr <- as.integer(alive_end == 0L)
      R_tr <- as.integer(alive_end == 1L & in_state_end == 1L)
      C_tr <- as.integer(R_tr == 0L & D_tr == 0L)

      exit_idx <- which(R_tr == 0L)
      rem_idx  <- which(R_tr == 1L)

      cl_tr_ar   <- cluster_by_id[id_tr_ar]
      cl_exit_ar <- cluster_by_id[id_tr_ar[exit_idx]]
      cl_rem_ar  <- cluster_by_id[id_tr_ar[rem_idx]]

      Y_pseudo_ar <- shf_train[at_risk_tr, tt + 1L]
      n_ar        <- length(idx_tr_ar)

      rows_vl       <- if (any(at_risk_vl)) row_index[id_vl_ar, tt] else integer(0)
      cn_base_tr    <- colnames(X_nat_tr_base)
      X_nat_vl_base <- if (length(rows_vl)) {
        align_cols(make_design(D, rows_vl, cols_base), cn_base_tr)
      } else NULL
      X_shf_vl_base <- if (length(rows_vl)) {
        align_cols(patch_shifted_design(
          X_nat = make_design(D, rows_vl, cols_base), D_shifted_tt = D_shifted[[tt]],
          rows = rows_vl, a_names = a_names, tt = tt, k = k, t_min = t_min
        ), cn_base_tr)
      } else NULL

      reg_res <- fit_transition_regressions(
        X_nat_tr_base  = X_nat_tr_base,
        X_shf_tr_base  = X_shf_tr_base,
        X_nat_tr_haz   = X_nat_tr_haz,
        X_shf_tr_haz   = X_shf_tr_haz,
        Y_pseudo_ar    = Y_pseudo_ar,
        Y_init_ar      = Y_init[id_tr_ar],
        R_tr           = R_tr,
        D_tr           = D_tr,
        exit_idx       = exit_idx,
        rem_idx        = rem_idx,
        n_ar           = n_ar,
        cl_tr_ar       = cl_tr_ar,
        cl_exit_ar     = cl_exit_ar,
        cl_rem_ar      = cl_rem_ar,
        pool_g_death     = pool_g_death,
        fit_dex_pooled   = fit_dex_pooled,
        cn_pool          = cn_pool,
        cols_pool        = cols_pool,
        pool_q_exit      = pool_q_exit,
        fit_qexit_pooled = fit_qexit_pooled,
        cn_qexit_pool    = cn_qexit_pool,
        cols_pool_qexit  = cols_pool_qexit,
        pool_time_basis  = pool_time_basis,
        D                = D,
        rows_tr          = rows_tr,
        D_shifted_tt     = D_shifted[[tt]],
        a_names          = a_names,
        tt               = tt,
        k                = k,
        t_min            = t_min,
        pooled_pred_tr   = pooled_pred_cache$tr[[tt]],
        pooled_pred_vl   = pooled_pred_cache$vl[[tt]],
        rows_vl          = rows_vl,
        X_nat_vl_base    = X_nat_vl_base,
        X_shf_vl_base    = X_shf_vl_base,
        is_binom         = is_binom,
        tmax           = tmax,
        scale_info     = scale_info,
        sl_remain      = sl_remain,
        sl_death       = sl_death,
        sl_y           = sl_y,
        sl_recursive   = sl_recursive,
        sl_rec_early   = sl_rec_early,
        rec_transition = rec_transition,
        inner_v        = inner_v,
        seed           = seed,
        f_idx          = f_idx,
        parallel       = parallel,
        reg_workers    = reg_workers,
        sl_workers     = sl_workers
      )

      r_rem  <- reg_res$g_remain
      r_dex  <- reg_res$g_death_exit
      r_exit <- reg_res$Q_exit
      r_qrem <- reg_res$Q_rem
      r_vl   <- reg_res$valid

      sl_rem      <- r_rem$sl_rem;    fit_rem     <- r_rem$fit_rem
      p_rem_const <- r_rem$p_rem_const
      p_rem_nat   <- r_rem$p_rem_nat; p_rem_shf   <- r_rem$p_rem_shf

      sl_dex      <- r_dex$sl_dex;    fit_dex     <- r_dex$fit_dex
      p_dex_const <- r_dex$p_dex_const
      p_dex_nat   <- r_dex$p_dex_nat; p_dex_shf   <- r_dex$p_dex_shf

      sl_qexit    <- r_exit$sl_qexit; fit_qexit   <- r_exit$fit_qexit
      cn_qexit    <- r_exit$cn_qexit; muY         <- r_exit$muY
      q_death_nat <- r_exit$q_death_nat; q_dc_nat  <- r_exit$q_dc_nat
      q_death_shf <- r_exit$q_death_shf; q_dc_shf  <- r_exit$q_dc_shf

      sl_qrem     <- r_qrem$sl_qrem;  fit_qrem    <- r_qrem$fit_qrem
      cn_qrem     <- r_qrem$cn_qrem;  muP         <- r_qrem$muP
      q_rem_nat   <- r_qrem$q_rem_nat; q_rem_shf  <- r_qrem$q_rem_shf

      rm(reg_res, r_rem, r_dex, r_exit, r_qrem); gc()

      q_death_nat_raw <- q_death_nat; q_dc_nat_raw <- q_dc_nat
      q_death_shf_raw <- q_death_shf; q_dc_shf_raw <- q_dc_shf

      if (length(absorb)) {
        D_block_tr  <- D[rows_tr]
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
      } else if (pooled_dex_used && tt == tmax) {
        meta <- sl_meta(sl_dex_pooled, f_idx, 0L, "g_death_exit_pooled",
                        length(all_Y_pool), n_val_tt)
      } else {
        meta <- NULL
      }
      if (!is.null(meta)) sl_chunks_here[[length(sl_chunks_here) + 1L]] <- meta

      if (!is.null(sl_qexit)) {
        meta <- sl_meta(sl_qexit, f_idx, tt, "Q_exit", length(exit_idx), n_val_tt)
      } else if (pooled_qexit_used && tt == tmax) {
        meta <- sl_meta(sl_qexit_pooled, f_idx, 0L, "Q_exit_pooled",
                        length(all_Y_qexit_pool), n_val_tt)
      } else {
        meta <- NULL
      }
      if (!is.null(meta)) sl_chunks_here[[length(sl_chunks_here) + 1L]] <- meta
      meta <- sl_meta(sl_qrem,  f_idx, tt, "Q_rem",  length(rem_idx),  n_val_tt)
      if (!is.null(meta)) sl_chunks_here[[length(sl_chunks_here) + 1L]] <- meta

      rm(fit_rem, fit_dex, fit_qexit, fit_qrem); gc(FALSE)

      if (!is.null(r_vl)) {
        p_rem_nat_vl   <- r_vl$p_rem_nat;   p_rem_shf_vl   <- r_vl$p_rem_shf
        p_dex_nat_vl   <- r_vl$p_dex_nat;   p_dex_shf_vl   <- r_vl$p_dex_shf
        q_death_nat_vl <- r_vl$q_death_nat; q_dc_nat_vl    <- r_vl$q_dc_nat
        q_death_shf_vl <- r_vl$q_death_shf; q_dc_shf_vl    <- r_vl$q_dc_shf
        q_rem_nat_vl   <- r_vl$q_rem_nat;   q_rem_shf_vl   <- r_vl$q_rem_shf

        if (length(absorb)) {
          D_block_vl     <- D[rows_vl]
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

        alive_end_vl  <- as.integer(D[[alive]][rows_vl])
        icu_end_vl    <- as.integer(D[[in_state]][rows_vl])
        D_vl          <- as.integer(alive_end_vl == 0L)
        R_vl          <- as.integer(alive_end_vl == 1L & icu_end_vl == 1L)
        rem_idx_vl    <- which(R_vl == 1L)
        exit_idx_vl   <- which(R_vl == 0L)
        q_exit_obs_vl_raw <- ifelse(D_vl == 1L, q_death_nat_vl, q_dc_nat_vl)
        q_exit_obs_tr_raw <- ifelse(D_tr == 1L, q_death_nat_raw, q_dc_nat_raw)

        m_grem_tr  <- binom_metric(R_tr, p_rem_nat)
        m_grem_vl  <- binom_metric(R_vl, p_rem_nat_vl)
        m_gdeath_tr <- binom_metric(D_tr[exit_idx],    p_dex_nat[exit_idx])
        m_gdeath_vl <- binom_metric(D_vl[exit_idx_vl], p_dex_nat_vl[exit_idx_vl])
        m_qrem_tr  <- gauss_metric(Y_init[id_tr_ar[rem_idx]], q_rem_nat[rem_idx])
        m_qrem_vl  <- gauss_metric(Y_init[id_vl_ar[rem_idx_vl]], q_rem_nat_vl[rem_idx_vl])

        m_qexit_tr <- if (isTRUE(is_binom))
          binom_metric(Y_init[id_tr_ar[exit_idx]], q_exit_obs_tr_raw[exit_idx])
        else gauss_metric(Y_init[id_tr_ar[exit_idx]], q_exit_obs_tr_raw[exit_idx])

        m_qexit_vl <- if (isTRUE(is_binom))
          binom_metric(Y_init[id_vl_ar[exit_idx_vl]], q_exit_obs_vl_raw[exit_idx_vl])
        else gauss_metric(Y_init[id_vl_ar[exit_idx_vl]], q_exit_obs_vl_raw[exit_idx_vl])

        branch_diag_here[[length(branch_diag_here) + 1L]] <- data.table::data.table(
          fold = f_idx, t = tt,
          n_tr_ar = length(id_tr_ar), n_vl_ar = length(id_vl_ar),
          n_tr_rem = length(rem_idx),  n_vl_rem = length(rem_idx_vl),
          n_tr_exit = length(exit_idx), n_vl_exit = length(exit_idx_vl),

          target_grem_tr = m_grem_tr$target,  pred_grem_tr = m_grem_tr$pred,
          target_grem_vl = m_grem_vl$target,  pred_grem_vl = m_grem_vl$pred,

          target_gdeath_tr = m_gdeath_tr$target, pred_gdeath_tr = m_gdeath_tr$pred,
          target_gdeath_vl = m_gdeath_vl$target, pred_gdeath_vl = m_gdeath_vl$pred,

          target_qrem_tr = m_qrem_tr$target,  pred_qrem_tr = m_qrem_tr$pred,
          target_qrem_vl = m_qrem_vl$target,  pred_qrem_vl = m_qrem_vl$pred,

          target_qexit_tr = m_qexit_tr$target, pred_qexit_tr = m_qexit_tr$pred,
          target_qexit_vl = m_qexit_vl$target, pred_qexit_vl = m_qexit_vl$pred,

          qexit_type = if (isTRUE(is_binom)) "binomial" else "gaussian",

          grem_tr_brier = m_grem_tr$brier,     grem_vl_brier = m_grem_vl$brier,
          grem_tr_auc   = m_grem_tr$auc,       grem_vl_auc   = m_grem_vl$auc,
          grem_tr_cal_int   = m_grem_tr$cal_int,   grem_vl_cal_int   = m_grem_vl$cal_int,
          grem_tr_cal_slope = m_grem_tr$cal_slope, grem_vl_cal_slope = m_grem_vl$cal_slope,

          gdeath_tr_brier = m_gdeath_tr$brier,   gdeath_vl_brier = m_gdeath_vl$brier,
          gdeath_tr_auc   = m_gdeath_tr$auc,     gdeath_vl_auc   = m_gdeath_vl$auc,
          gdeath_tr_cal_int   = m_gdeath_tr$cal_int,   gdeath_vl_cal_int   = m_gdeath_vl$cal_int,
          gdeath_tr_cal_slope = m_gdeath_tr$cal_slope, gdeath_vl_cal_slope = m_gdeath_vl$cal_slope,

          qrem_tr_mse  = m_qrem_tr$mse,   qrem_vl_mse  = m_qrem_vl$mse,
          qrem_tr_rmse = m_qrem_tr$rmse,  qrem_vl_rmse = m_qrem_vl$rmse,
          qrem_tr_mae  = m_qrem_tr$mae,   qrem_vl_mae  = m_qrem_vl$mae,
          qrem_tr_cor  = m_qrem_tr$cor,   qrem_vl_cor  = m_qrem_vl$cor,
          qrem_tr_cal_int   = m_qrem_tr$cal_int,   qrem_vl_cal_int   = m_qrem_vl$cal_int,
          qrem_tr_cal_slope = m_qrem_tr$cal_slope, qrem_vl_cal_slope = m_qrem_vl$cal_slope,

          qexit_tr_brier = if (isTRUE(is_binom)) m_qexit_tr$brier else NA_real_,
          qexit_vl_brier = if (isTRUE(is_binom)) m_qexit_vl$brier else NA_real_,
          qexit_tr_auc   = if (isTRUE(is_binom)) m_qexit_tr$auc   else NA_real_,
          qexit_vl_auc   = if (isTRUE(is_binom)) m_qexit_vl$auc   else NA_real_,

          qexit_tr_mse  = if (!isTRUE(is_binom)) m_qexit_tr$mse  else NA_real_,
          qexit_vl_mse  = if (!isTRUE(is_binom)) m_qexit_vl$mse  else NA_real_,
          qexit_tr_rmse = if (!isTRUE(is_binom)) m_qexit_tr$rmse else NA_real_,
          qexit_vl_rmse = if (!isTRUE(is_binom)) m_qexit_vl$rmse else NA_real_,
          qexit_tr_mae  = if (!isTRUE(is_binom)) m_qexit_tr$mae  else NA_real_,
          qexit_vl_mae  = if (!isTRUE(is_binom)) m_qexit_vl$mae  else NA_real_,
          qexit_tr_cor  = if (!isTRUE(is_binom)) m_qexit_tr$cor  else NA_real_,
          qexit_vl_cor  = if (!isTRUE(is_binom)) m_qexit_vl$cor  else NA_real_,

          qexit_tr_cal_int   = m_qexit_tr$cal_int,   qexit_vl_cal_int   = m_qexit_vl$cal_int,
          qexit_tr_cal_slope = m_qexit_tr$cal_slope, qexit_vl_cal_slope = m_qexit_vl$cal_slope
        )
      }

      if (any(!at_risk_tr)) {
        nat_train[!at_risk_tr, tt] <- nat_train[!at_risk_tr, tt + 1L]
        shf_train[!at_risk_tr, tt] <- shf_train[!at_risk_tr, tt + 1L]
      }
      if (any(!at_risk_vl)) {
        nat_valid[!at_risk_vl, tt] <- nat_valid[!at_risk_vl, tt + 1L]
        shf_valid[!at_risk_vl, tt] <- shf_valid[!at_risk_vl, tt + 1L]
      }

      exit_branch_nat_tr <- p_dex_nat * q_death_nat + (1 - p_dex_nat) * q_dc_nat
      exit_branch_shf_tr <- p_dex_shf * q_death_shf + (1 - p_dex_shf) * q_dc_shf

      fold_diag_here[[length(fold_diag_here) + 1L]] <- data.table::data.table(
        fold = f_idx, t = tt,

        n_train_at_risk = sum(at_risk_tr),
        n_valid_at_risk = sum(at_risk_vl),
        n_train_death   = sum(D_tr == 1L),
        n_train_dc      = sum(C_tr == 1L),
        n_train_remain  = sum(R_tr == 1L),

        used_const_rem   = if (!is.null(sl_rem))   isTRUE(sl_rem$used_const)   else NA,
        used_const_dex   = if (!is.null(sl_dex)) {
          isTRUE(sl_dex$used_const)
        } else if (pooled_dex_used) FALSE else if (!is.na(p_dex_const)) TRUE else NA,
        used_const_qexit = if (!is.null(sl_qexit)) {
          isTRUE(sl_qexit$used_const)
        } else if (pooled_qexit_used) FALSE else NA,
        used_const_qrem  = if (!is.null(sl_qrem))  isTRUE(sl_qrem$used_const)  else NA,

        p_rem_nat_mean = mean_or_na(p_rem_nat), p_rem_nat_sd = sd_or_na(p_rem_nat),
        p_rem_shf_mean = mean_or_na(p_rem_shf), p_rem_shf_sd = sd_or_na(p_rem_shf),

        p_dex_nat_mean = mean_or_na(p_dex_nat), p_dex_nat_sd = sd_or_na(p_dex_nat),
        p_dex_shf_mean = mean_or_na(p_dex_shf), p_dex_shf_sd = sd_or_na(p_dex_shf),

        q_rem_nat_mean = mean_or_na(q_rem_nat), q_rem_nat_sd = sd_or_na(q_rem_nat),
        q_rem_shf_mean = mean_or_na(q_rem_shf), q_rem_shf_sd = sd_or_na(q_rem_shf),

        q_exit_nat_mean = mean_or_na(exit_branch_nat_tr), q_exit_nat_sd = sd_or_na(exit_branch_nat_tr),
        q_exit_shf_mean = mean_or_na(exit_branch_shf_tr), q_exit_shf_sd = sd_or_na(exit_branch_shf_tr),

        Q_nat_pre_mean = mean_or_na(Q_nat_ar), Q_nat_pre_sd = sd_or_na(Q_nat_ar),
        Q_nat_pre_min  = min_or_na(Q_nat_ar),  Q_nat_pre_max = max_or_na(Q_nat_ar),

        Q_shf_pre_mean = mean_or_na(Q_shf_ar), Q_shf_pre_sd = sd_or_na(Q_shf_ar),
        Q_shf_pre_min  = min_or_na(Q_shf_ar),  Q_shf_pre_max = max_or_na(Q_shf_ar),

        Y_target_mean = mean_or_na(Y_pseudo_ar), Y_target_sd = sd_or_na(Y_pseudo_ar),
        Y_target_min  = min_or_na(Y_pseudo_ar),  Y_target_max = max_or_na(Y_pseudo_ar),

        Q_nat_vl_mean = if (!is.null(r_vl)) mean_or_na(Q_nat_vl) else NA_real_,
        Q_shf_vl_mean = if (!is.null(r_vl)) mean_or_na(Q_shf_vl) else NA_real_
      )

      rm(sl_rem, sl_dex, sl_qexit, sl_qrem)
      gc()

    } # end tt loop

    list(
      fold        = f_idx,
      valid_ids   = vl_ids,
      nat_valid   = nat_valid,
      shf_valid   = shf_valid,
      sl_meta     = if (length(sl_chunks_here))
        data.table::rbindlist(sl_chunks_here, use.names = TRUE, fill = TRUE) else NULL,
      fold_diag   = if (length(fold_diag_here))
        data.table::rbindlist(fold_diag_here, use.names = TRUE, fill = TRUE) else NULL,
      branch_diag = if (length(branch_diag_here))
        data.table::rbindlist(branch_diag_here, use.names = TRUE, fill = TRUE) else NULL
    )
  } # end fold_worker

  res_by_fold <- par_lapply(
    X        = seq_along(folds),
    FUN      = function(f_idx) fold_worker(f_idx, folds, Y_init, row_index, cluster_by_id),
    workers  = fold_workers,
    parallel = parallel,
    seed     = TRUE
  )

  fold_diag_chunks   <- list()
  branch_diag_chunks <- list()

  failed <- sapply(res_by_fold, function(x) inherits(x, "error") || is.character(x))
  if (any(failed)) {
    for (i in which(failed)) message(sprintf("fold %d: %s", i, as.character(res_by_fold[[i]])))
    stop("qreg: fold worker(s) failed -- see messages above")
  }

  for (fr in res_by_fold) {
    vl_ids <- fr$valid_ids
    pred_nat_all[vl_ids, ] <- fr$nat_valid
    pred_shf_all[vl_ids, ] <- fr$shf_valid
    if (!is.null(fr$sl_meta))     sl_chunks[[length(sl_chunks) + 1L]]               <- fr$sl_meta
    if (!is.null(fr$fold_diag))   fold_diag_chunks[[length(fold_diag_chunks) + 1L]] <- fr$fold_diag
    if (!is.null(fr$branch_diag)) branch_diag_chunks[[length(branch_diag_chunks) + 1L]] <- fr$branch_diag
  }

  fold_diag   <- if (length(fold_diag_chunks))
    data.table::rbindlist(fold_diag_chunks, use.names = TRUE, fill = TRUE) else NULL
  branch_diag <- if (length(branch_diag_chunks))
    data.table::rbindlist(branch_diag_chunks, use.names = TRUE, fill = TRUE) else NULL

  psi_plugin_nat_scaled <- mean(pred_nat_all[, 1L], na.rm = TRUE)
  psi_plugin_shf_scaled <- mean(pred_shf_all[, 1L], na.rm = TRUE)
  psi_scaled            <- psi_plugin_shf_scaled

  q1_vec     <- pred_shf_all[, 1L]
  ok_q1      <- is.finite(q1_vec)
  n_eff      <- sum(ok_q1)
  if (n_eff < 2L) stop("qreg: fewer than 2 finite Q predictions at t=1.", call. = FALSE)
  se_naive_scaled <- stats::sd(q1_vec[ok_q1]) / sqrt(n_eff)

  se_scaled      <- se_naive_scaled
  se_boot_scaled <- NULL
  boot_psi_scaled <- NULL

  if (n_boot > 0L) {
    unique_cl        <- unique(cluster_by_id)
    n_cl             <- length(unique_cl)
    use_cluster_boot <- n_cl < N

    boot_psi_scaled <- numeric(n_boot)
    for (b in seq_len(n_boot)) {
      if (use_cluster_boot) {
        boot_cl <- sample(unique_cl, n_cl, replace = TRUE)
        bidx    <- unlist(lapply(boot_cl, function(cl) which(cluster_by_id == cl)),
                          use.names = FALSE)
      } else {
        bidx <- sample(N, N, replace = TRUE)
      }
      N_b             <- length(bidx)
      Y_init_b        <- Y_init[bidx]
      row_index_b     <- row_index[bidx, , drop = FALSE]
      cluster_by_id_b <- cluster_by_id[bidx]
      folds_b         <- make_subject_folds_qreg(seq_len(N_b), v = v,
                                                  cluster_by_id = cluster_by_id_b,
                                                  seed = seed + b)
      res_b <- par_lapply(
        X        = seq_along(folds_b),
        FUN      = function(f_idx) fold_worker(f_idx, folds_b, Y_init_b,
                                               row_index_b, cluster_by_id_b,
                                               verbose_fold = FALSE),
        workers  = fold_workers,
        parallel = parallel,
        seed     = TRUE
      )
      failed_b <- sapply(res_b, function(x) inherits(x, "error") || is.character(x))
      if (any(failed_b)) {
        warning(sprintf("qreg bootstrap rep %d: fold worker failed, skipping", b))
        boot_psi_scaled[b] <- NA_real_
        next
      }
      pred_shf_b <- matrix(NA_real_, nrow = N_b, ncol = tmax + 1L)
      pred_shf_b[, tmax + 1L] <- Y_init_b
      for (fr in res_b) pred_shf_b[fr$valid_ids, ] <- fr$shf_valid
      boot_psi_scaled[b] <- mean(pred_shf_b[, 1L], na.rm = TRUE)
    }
    ok_b           <- is.finite(boot_psi_scaled)
    se_boot_scaled <- if (sum(ok_b) >= 2L) stats::sd(boot_psi_scaled[ok_b]) else NA_real_
    se_scaled      <- se_boot_scaled
  }

  Y_obs_scaled <- mean(Y_init, na.rm = TRUE)

  t_vec <- seq_len(tmax)
  diag_table <- data.frame(
    t          = t_vec,
    n_at_risk  = vapply(t_vec, function(tt) sum(!is.na(row_index[, tt])), integer(1L)),
    mean_Q_nat = vapply(t_vec, function(tt) {
      ar <- !is.na(row_index[, tt]); mean(pred_nat_all[ar, tt], na.rm = TRUE)
    }, numeric(1L)),
    mean_Q_shf = vapply(t_vec, function(tt) {
      ar <- !is.na(row_index[, tt]); mean(pred_shf_all[ar, tt], na.rm = TRUE)
    }, numeric(1L))
  )
  diag_table$delta      <- diag_table$mean_Q_shf - diag_table$mean_Q_nat
  diag_table            <- diag_table[diag_table$n_at_risk > 0L, , drop = FALSE]

  if (scale_info$bounded) {
    psi         <- scale_info$from_unit(psi_scaled)
    psi_natural <- scale_info$from_unit(psi_plugin_nat_scaled)
    psi_shifted <- scale_info$from_unit(psi_plugin_shf_scaled)
    se_naive    <- scale_info$y_rng * se_naive_scaled
    se_boot     <- if (!is.null(se_boot_scaled)) scale_info$y_rng * se_boot_scaled else NULL
    se          <- scale_info$y_rng * se_scaled
    Y_obs       <- scale_info$from_unit(Y_obs_scaled)
    diag_table$mean_Q_nat <- scale_info$from_unit(diag_table$mean_Q_nat)
    diag_table$mean_Q_shf <- scale_info$from_unit(diag_table$mean_Q_shf)
    diag_table$delta      <- scale_info$y_rng * diag_table$delta
  } else {
    psi         <- psi_scaled
    psi_natural <- psi_plugin_nat_scaled
    psi_shifted <- psi_plugin_shf_scaled
    se_naive    <- se_naive_scaled
    se_boot     <- se_boot_scaled
    se          <- se_scaled
    Y_obs       <- Y_obs_scaled
  }

  out <- list(
    psi   = psi,
    se    = se,
    ci    = c(psi - 1.96 * se, psi + 1.96 * se),
    Y_obs = Y_obs,

    predictions = list(
      natural = pred_nat_all,
      shifted = pred_shf_all
    ),

    diagnostics = list(
      recursion_diag = fold_diag,
      branch_cal     = branch_diag,
      sl_summary     = if (length(sl_chunks))
        data.table::rbindlist(sl_chunks, use.names = TRUE, fill = TRUE) else NULL,
      diag_table     = diag_table
    ),

    decomposition = list(
      psi_plugin_nat  = psi_natural,
      psi_plugin_shf  = psi_shifted,
      se_naive        = se_naive,
      se_boot         = se_boot,
      boot_psi_dist   = if (!is.null(boot_psi_scaled)) {
        if (scale_info$bounded) scale_info$from_unit(boot_psi_scaled) else boot_psi_scaled
      } else NULL
    ),

    settings = list(
      n              = length(ids),
      outcome_family = outcome_family,
      pool_g_death   = pool_g_death,
      pool_q_exit    = pool_q_exit,
      pool_time      = pool_time,
      variable_info  = list(
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
    )
  )
  class(out) <- c("qreg_fit", "list")
  out
}
