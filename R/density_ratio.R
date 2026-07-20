
make_group_folds <- function(groups, v = 5L, seed = 1L) {
  g <- as.character(groups)
  n <- length(g)
  if (n == 0L) return(integer())

  g2 <- g
  g2[is.na(g2)] <- "__NA_GROUP__"

  set.seed(as.integer(seed))
  v_eff <- max(1L, min(as.integer(v), length(unique(g2))))

  folds <- origami::make_folds(cluster_ids = g2, V = v_eff)
  as.integer(origami::folds2foldvec(folds))
}


build_stack <- function(
    D, t, a_names, policy_spec_fun,
    base_ok, lag_names_by_depth, k,
    id, time,
    L_NOW = NULL,
    tol_equal = 0,
    t_min = 1L
) {
  DT <- data.table::as.data.table(D)
  data.table::setorderv(DT, c(id, time))

  rows_t <- which(DT[[time]] == t)
  if (!length(rows_t)) return(NULL)
  at_risk <- rows_t

  kk <- max(0L, min(as.integer(k), as.integer(t - t_min)))
  hist_cols <- if (kk > 0L) unlist(lag_names_by_depth[seq_len(kk)], use.names = FALSE) else character(0)
  cols_hist <- intersect(unique(c(base_ok, hist_cols)), names(DT))

  X_hist <- if (length(cols_hist)) {
    as.data.frame(DT[at_risk, ..cols_hist])
  } else {
    data.frame(`(Intercept)` = 1)
  }

  L_now_cols <- intersect(if (is.null(L_NOW)) character(0) else L_NOW, names(DT))
  if (length(L_now_cols)) {
    X_hist <- cbind(X_hist, as.data.frame(DT[at_risk, ..L_now_cols]))
  }

  A_obs <- as.data.frame(DT[at_risk, ..a_names, with = FALSE])
  A_obs[] <- lapply(A_obs, function(z) suppressWarnings(as.numeric(z)))

  block_t <- DT[rows_t]
  spec <- try(policy_spec_fun(block_t, t, a_names), silent = TRUE)

  A_shift <- A_obs
  if (inherits(spec, "data.frame") && nrow(spec) == nrow(block_t)) {
    for (a in intersect(names(spec), a_names)) {
      A_shift[[a]] <- suppressWarnings(as.numeric(spec[[a]]))
    }
  }
  for (a in names(A_shift)) {
    v  <- suppressWarnings(as.numeric(A_shift[[a]]))
    vo <- A_obs[[a]]
    bad <- !is.finite(v); if (any(bad)) v[bad] <- vo[bad]
    A_shift[[a]] <- v
  }

  X_obs   <- cbind(X_hist, A_obs)
  X_shift <- cbind(X_hist, A_shift)
  names(X_shift) <- names(X_obs)
  X <- rbind(X_obs, X_shift)

  S <- c(rep.int(0L, nrow(X_obs)), rep.int(1L, nrow(X_shift)))
  id_in <- as.character(DT[[id]][at_risk])

  list(
    X        = X,
    S        = as.numeric(S),
    X_obs    = X_obs,
    id_pair  = c(id_in, id_in),
    rows_idx = at_risk
  )
}



dr_from_prob <- function(p_mat, beta, bounds, denom_cap) {
  p_mat    <- as.matrix(p_mat)
  p_mat[]  <- pmin(1 - bounds, pmax(bounds, p_mat))
  psi      <- p_mat
  psi[]    <- p_mat / (1 - pmin(p_mat, denom_cap))
  r        <- as.vector(psi %*% beta)
  r[!is.finite(r)] <- 1
  r
}


#' Estimate density ratios for a modified treatment policy
#'
#' Fits per-time-point treatment models and computes the instantaneous density
#' ratio \eqn{r_t = d\tilde{P}(A_t | H_t) / dP(A_t | H_t)} comparing the
#' modified treatment policy (MTP) to the natural course. The output is passed
#' directly to [sdr()] or [itmle()] as the `weight_object` argument.
#'
#' Cross-fitting is used throughout: treatment models are trained on one fold
#' and density ratios predicted on the held-out fold, so the same fold
#' structure is inherited by the downstream estimator.
#'
#' @param df A `data.frame` in long format (one row per subject per time point).
#' @param a_names Character vector of treatment variable names. Can contain
#'   one or more variables for joint multi-treatment policies (e.g.
#'   `c("A1", "A2")`). The joint density ratio for the full treatment vector
#'   is estimated at each time point. Pass the same `a_names` to the
#'   downstream estimators (`sdr()`, `itmle()`, `qreg()`).
#' @param tmax Integer. Maximum follow-up time (number of time points).
#' @param baseline Character vector of baseline covariate names.
#'   Default `NULL` (no baseline covariates).
#' @param tv_names Character vector of time-varying covariate names.
#'   Default `character(0)`.
#' @param no_lag_vars Character vector of variables that should not be lagged.
#' @param policy_names Character vector of column names holding shifted
#'   treatment values under the MTP.
#' @param sl_g SuperLearner library for treatment models. Required.
#' @param k Integer. Number of lags to include. Default `2`.
#' @param inner_v Integer. Inner cross-validation folds for SuperLearner.
#'   Default `10L`.
#' @param cluster Character. Cluster variable name for clustered fold
#'   assignment. `NULL` uses independent subject folds.
#' @param v Integer. Number of outer cross-fitting folds. Default `5L`.
#' @param seed Integer random seed. Default `1L`.
#' @param policy_spec_fun A function `(D_block, t, a_names)` returning a
#'   `data.table` with shifted treatment values for time `t`. Used when the
#'   MTP cannot be pre-computed as static columns.
#' @param id Character. Subject identifier column name. Default `"id"`.
#' @param time Character. Time column name. Default `"trial_time"`.
#' @param bounds Numeric. Probability floor for clipping predicted probabilities
#'   before conversion to density ratios (standard pathway only). Default `1e-5`.
#' @param dr_sl Logical. Selects the metalearner used to combine base learners.
#'   `FALSE` (default): standard SuperLearner with NNLS metalearner -- each base
#'   learner outputs a propensity score (probability scale) and the metalearner
#'   combines them with non-negative least squares; the combined probability is
#'   then converted to a density ratio via `p / (1 - p)`.
#'   `TRUE`: the Wu-Benkeser metalearner ([method.WB_dr]) -- combines base
#'   learners directly in density-ratio space by minimising a log density-ratio
#'   loss rather than a squared-error loss on the probability scale.  This can
#'   give better-calibrated density ratios when the propensity is far from 0.5
#'   and extreme ratios are a concern.  When `dr_sl = TRUE` the SuperLearner
#'   object returns density ratios directly rather than probabilities.
#'   **Important:** the Wu-Benkeser pathway requires base learners whose
#'   `predict` method returns density ratios on the positive real line, not
#'   probabilities in `[0, 1]`.  Standard SuperLearner wrappers (e.g.
#'   `SL.glm`, `SL.xgboost`) are not compatible -- you must supply custom
#'   wrappers that internally fit a classifier and convert predictions to
#'   density ratios before returning them.  This package does not currently
#'   include such wrappers; built-in DR-returning wrappers compatible with
#'   the WB pathway are planned for version 1.0.
#' @param drop_small_cluster_splits Logical. Drop time points where a fold
#'   has too few clusters to fit a model. Default `TRUE`.
#' @param parallel_t Logical. Parallelise across time points. Default `FALSE`.
#' @param t_workers Integer. Number of workers for time-point parallelism.
#'   `NULL` auto-detects.
#' @param fold_workers Integer. Workers for parallelising across outer
#'   cross-fitting folds within each time point via [parallel::mclapply()]
#'   (fork-based; Linux/Mac only). `NULL` disables. When combined with
#'   `parallel_t`, the nested `mclapply` scheme limits use of multi-threaded
#'   or GPU-based learners. Default `NULL`.
#' @param sl_workers Integer. Workers for parallel learner evaluation within
#'   each SuperLearner call via [SuperLearner::mcSuperLearner()] (fork-based;
#'   Linux/Mac only). `NULL` disables. Default `NULL`. **Note:** incompatible
#'   with `dbarts`/`SL.dbarts` learners -- remove dbarts from `sl_g` when
#'   `sl_workers` is set.
#'
#' @section Parallelism:
#'   `t_workers`, `fold_workers`, and `sl_workers` can be used independently
#'   or together. Using a single level is robust with all learners. Combining
#'   two or more creates nested `mclapply` calls; in that case multi-threaded
#'   or GPU-based learners (e.g. xgboost with CUDA, OpenMP-based methods) may
#'   crash in the child processes and should be avoided or limited to one thread.
#'
#' @return A list with components:
#'   \describe{
#'     \item{`weights_dt`}{A `data.table` with one row per subject-time
#'       containing `Rt_t` (instantaneous density ratio) and `global_fold`.
#'       Pass this to [sdr()] or [itmle()] as `weight_object`.}
#'     \item{`sl_summary`}{`data.table` of SuperLearner learner weights.
#'       One row per learner per (fold, time-point). Useful for checking
#'       which treatment models dominate across time points.}
#'     \item{`fold_diag`}{`data.table` of per-fold, per-time-point
#'       diagnostics: mean predicted probabilities, effective sample sizes
#'       (ESS) for the density ratios, and calibration summaries for
#'       each treatment model. ESS collapse (ESS much smaller than n)
#'       indicates extreme density ratios and warrants caution.}
#'   }
#'
#' @references
#' Diaz I, Williams N, Hoffman KL, Schenck EJ (2021). Nonparametric Causal
#' Effects Based on Longitudinal Modified Treatment Policies. *JASA*
#' 118(542):846-857.
#'
#' Williams NT, Diaz I (2023). lmtp: An R package for estimating the causal
#' effects of modified treatment policies. *Observational Studies*.
#'
#' @seealso [sdr()], [itmle()]
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
#' # Policy: increase treatment probability by 0.3
#' policy_fn <- function(D_block, t, a_names) {
#'   out <- D_block[, ..a_names, drop = FALSE]
#'   out[[a_names[1]]] <- pmin(D_block[[a_names[1]]] + 0.3, 1)
#'   out
#' }
#'
#' wr <- density_ratio(
#'   df              = df,
#'   a_names         = "A",
#'   tmax            = 5L,
#'   baseline        = "age",
#'   tv_names        = "L1",
#'   sl_g            = c("SL.mean", "SL.glm"),
#'   k               = 1L,
#'   inner_v         = 3L,
#'   v               = 3L,
#'   seed            = 1L,
#'   id              = "id",
#'   time            = "time",
#'   policy_spec_fun = policy_fn
#' )
#' head(wr$weights_dt)
#' }
#'
#' # Multi-treatment: joint density ratio for two binary treatments (A1, A2)
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
#' # Joint policy: shift both treatments upward by 0.3
#' policy_fn2 <- function(D_block, t, a_names) {
#'   out <- D_block[, ..a_names, drop = FALSE]
#'   out[[a_names[1]]] <- pmin(D_block[[a_names[1]]] + 0.3, 1)
#'   out[[a_names[2]]] <- pmin(D_block[[a_names[2]]] + 0.3, 1)
#'   out
#' }
#'
#' wr2 <- density_ratio(
#'   df              = df2,
#'   a_names         = c("A1", "A2"),
#'   tmax            = 5L,
#'   baseline        = "age",
#'   tv_names        = "L1",
#'   sl_g            = c("SL.mean", "SL.glm"),
#'   k               = 1L,
#'   inner_v         = 3L,
#'   v               = 3L,
#'   seed            = 1L,
#'   id              = "id",
#'   time            = "time",
#'   policy_spec_fun = policy_fn2
#' )
#' head(wr2$weights_dt)
#' }
#'
#' @export
density_ratio <- function(
    df, a_names, tmax,
    baseline = NULL, tv_names = character(0),
    no_lag_vars = character(0),
    policy_names = character(0),
    sl_g, k = 2, inner_v = 10L,
    cluster = NULL, v = 5L,
    seed = 1L,
    policy_spec_fun,
    id = "id", time = "trial_time",
    bounds = 1e-5,
    dr_sl = FALSE,
    drop_small_cluster_splits = TRUE,
    parallel_t = FALSE,
    t_workers = NULL,
    fold_workers = NULL,
    sl_workers = NULL
) {

  `%||%` <- function(x, y) if (!is.null(x)) x else y

  keep_cols <- unique(c(
    id, time, a_names,
    baseline, tv_names,
    no_lag_vars, policy_names, cluster
  ))
  keep_cols <- intersect(keep_cols, names(df))
  df <- as.data.frame(df)[, keep_cols, drop = FALSE]

  data_check(data = df, id = id, time = time,
             a_names = a_names, tv_names = tv_names,
             bs = baseline, tmin = 1, tmax = tmax)

  DT <- data.table::as.data.table(df)
  data.table::setorderv(DT, c(id, time))
  rm(df); gc()

  t_min <- DT[, .SD[1L], by = id][[time]] |> max()

  hist_out <- prepare_history_lags(
    DT          = DT,
    a_names     = a_names,
    tv_names     = tv_names,
    k           = k,
    id      = id,
    no_lag_vars = no_lag_vars
  )
  DT                 <- hist_out$DT
  lag_names_by_depth <- hist_out$lag_names_by_depth

  CL <- if (!is.null(cluster) && cluster %in% names(DT)) cluster else id
  DT[["global_fold"]] <- make_group_folds(DT[[CL]], v = v, seed = seed)

  base_ok <- intersect(baseline %||% character(0), names(DT))

  Rt_t <- rep(NA_real_, nrow(DT))

  fit_one_t <- function(t) {

    sl_rows  <- vector("list", 0L)
    diag_rows <- vector("list", 0L)

    stk <- build_stack(
      DT, t, a_names, policy_spec_fun,
      base_ok, lag_names_by_depth, k,
      id = id, time = time,
      L_NOW  = tv_names,
      tol_equal = 0,
      t_min = t_min
    )
    if (is.null(stk)) return(NULL)

    n_nat    <- nrow(stk$X_obs)
    A_cols   <- intersect(a_names, names(stk$X_obs))
    A_shft_t <- stk$X[seq_len(n_nat) + n_nat, A_cols, drop = FALSE]
    no_shift <- all(stk$X_obs[, A_cols, drop = FALSE] == A_shft_t, na.rm = TRUE)

    changed <- rowSums(
      abs(as.matrix(stk$X_obs[, A_cols, drop = FALSE]) -
            as.matrix(A_shft_t[, A_cols, drop = FALSE])),
      na.rm = TRUE
    ) > 0

    A_nat_mat <- as.matrix(stk$X_obs[, A_cols, drop = FALSE])
    A_shf_mat <- as.matrix(A_shft_t[, A_cols, drop = FALSE])
    A_delta   <- A_shf_mat - A_nat_mat

    change_row <- data.table::data.table(
      t           = as.integer(t),
      n_nat       = n_nat,
      n_changed   = sum(changed, na.rm = TRUE),
      pct_changed = 100 * mean(changed, na.rm = TRUE)
    )

    if (length(A_cols) <= 3L) {
      for (ac in A_cols) {
        # pre/post means over whole population; delta quantiles among changed only
        # (unchanged have delta = 0 by construction and dilute the shift distribution)
        change_row[[paste0("mean_A_pre_",    ac)]] <- mean(A_nat_mat[, ac], na.rm = TRUE)
        change_row[[paste0("mean_A_post_",   ac)]] <- mean(A_shf_mat[, ac], na.rm = TRUE)
        if (any(changed)) {
          d <- A_delta[changed, ac]
          change_row[[paste0("mean_A_delta_",   ac)]] <- mean(d, na.rm = TRUE)
          change_row[[paste0("median_A_delta_", ac)]] <- median(d, na.rm = TRUE)
          change_row[[paste0("q05_A_delta_",    ac)]] <- as.numeric(stats::quantile(d, 0.05, na.rm = TRUE))
          change_row[[paste0("q95_A_delta_",    ac)]] <- as.numeric(stats::quantile(d, 0.95, na.rm = TRUE))
        }
      }
    } else {
      change_row[, `:=`(mean_A_pre = NA_real_, mean_A_post = NA_real_,
                        mean_A_delta = NA_real_, median_A_delta = NA_real_,
                        q05_A_delta  = NA_real_, q95_A_delta    = NA_real_)]
    }

    if (no_shift) {
      message(sprintf("t=%d: no shift -- ratios set to 1 (n=%d)", t, n_nat))
      return(list(
        idx    = stk$rows_idx,
        r_t    = rep(1, n_nat),
        sl     = NULL,
        diag   = NULL,
        change = change_row
      ))
    }

    fvec_t     <- DT$global_fold[stk$rows_idx]
    folds_here <- sort(unique(fvec_t))

    n_nat   <- nrow(stk$X_obs)
    r_hat   <- rep(NA_real_, n_nat)

    method_sl <- if (isTRUE(dr_sl)) {
      method.WB_dr()
    } else {
      SuperLearner::method.NNLS
    }

    cl_t <- if (!is.null(cluster) && cluster %in% names(DT)) {
      DT[[cluster]][stk$rows_idx]
    } else NULL

    fit_one_fold <- function(kk) {
      message(sprintf("[t=%d][fold %d][pid=%d] %s", t, kk, Sys.getpid(), format(Sys.time(), "%H:%M:%S")))
      idx_valid <- which(fvec_t == kk)
      idx_train <- which(fvec_t != kk)

      if (length(idx_valid) == 0L || length(idx_train) < 2L) {
        return(list(idx_valid = idx_valid, r = rep(1, length(idx_valid)),
                    sl_tab = NULL, diag_tab = NULL))
      }

      rows_train_nat   <- idx_train
      rows_train_shift <- idx_train + n_nat
      rows_train       <- c(rows_train_nat, rows_train_shift)

      id_nat <- DT[[id]][stk$rows_idx][rows_train_nat]

      Y_train <- stk$S[rows_train]
      X_train <- stk$X[rows_train, , drop = FALSE]

      cluster_inner <- if (!is.null(cl_t)) {
        c(cl_t[rows_train_nat], cl_t[rows_train_nat])
      } else {
        NULL
      }
      id_inner <- if (is.null(cluster_inner)) c(id_nat, id_nat) else NULL

      sl_res <- sl_block_fit(
        Y             = Y_train,
        X             = X_train,
        family        = stats::binomial(),
        sl_lib        = sl_g,
        v_inner       = inner_v,
        id_inner      = id_inner,
        cluster_inner = cluster_inner,
        seed_inner    = est_seed_for(seed, fold = kk, t = t, component = "density_ratio", phase = "fit"),
        seed_cv_rows  = est_seed_for(seed, fold = kk, t = t, component = "density_ratio", phase = "cv_rows"),
        seed_sl_fit   = est_seed_for(seed, fold = kk, t = t, component = "density_ratio", phase = "sl_fit"),
        method        = method_sl,
        sl_workers    = sl_workers
      )
      fit <- sl_res$fit

      sl_tab <- if (!is.null(sl_res$fit$coef)) {
        data.table::data.table(
          t          = as.integer(t),
          outer_fold = as.integer(kk),
          learner    = names(sl_res$fit$coef),
          weight     = as.numeric(sl_res$fit$coef),
          cvRisk     = as.numeric(sl_res$fit$cvRisk)
        )
      } else NULL

      newdata_obs <- as.data.frame(stk$X_obs[idx_valid, , drop = FALSE])

      rows_valid_nat   <- idx_valid
      rows_valid_shift <- idx_valid + n_nat
      rows_valid_stack <- c(rows_valid_nat, rows_valid_shift)

      pred_cal <- as.numeric(sl_predict(
        fit,
        newdata = as.data.frame(stk$X[rows_valid_stack, , drop = FALSE])
      ))

      diag_tab <- data.table::data.table(
        t          = as.integer(t),
        outer_fold = as.integer(kk),
        row_type   = rep(c("natural", "shifted"), each = length(idx_valid)),
        pair_id    = rep(idx_valid, times = 2),
        S          = as.integer(stk$S[rows_valid_stack]),
        p_raw      = if (!isTRUE(dr_sl)) as.numeric(pred_cal) else NA_real_,
        r_raw      = if (isTRUE(dr_sl))  as.numeric(pred_cal) else NA_real_,
        changed    = c(changed[idx_valid], changed[idx_valid])
      )

      pred_val <- as.numeric(sl_predict(fit, newdata = newdata_obs))
      r <- if (isTRUE(dr_sl)) pred_val else dr_from_prob(pred_val, 1, bounds, 0.999)

      list(idx_valid = idx_valid, r = r, sl_tab = sl_tab, diag_tab = diag_tab)
    }

    fold_res <- if (!is.null(fold_workers)) {
      parallel::mclapply(folds_here, fit_one_fold,
                         mc.cores = as.integer(fold_workers), mc.set.seed = TRUE)
    } else {
      lapply(folds_here, fit_one_fold)
    }

    for (fr in fold_res) {
      if (length(fr$idx_valid))  r_hat[fr$idx_valid]                       <- fr$r
      if (!is.null(fr$sl_tab))   sl_rows[[length(sl_rows) + 1L]]           <- fr$sl_tab
      if (!is.null(fr$diag_tab)) diag_rows[[length(diag_rows) + 1L]]       <- fr$diag_tab
    }

    r_hat[!is.finite(r_hat) | is.na(r_hat)] <- 1

    sl_tab <- if (length(sl_rows)) data.table::rbindlist(sl_rows, use.names = TRUE, fill = TRUE) else NULL

    diag_tab <- if (length(diag_rows)) {
      data.table::rbindlist(diag_rows, use.names = TRUE, fill = TRUE)
    } else {
      NULL
    }

    return(list(
      idx    = stk$rows_idx,
      r_t    = r_hat,
      sl     = sl_tab,
      diag   = diag_tab,
      change = change_row
    ))
  }

  if (isTRUE(parallel_t)) {
    if (is.null(t_workers)) {
      t_workers <- min(as.integer(getOption("mc.cores", 1L)), as.integer(tmax))
    }

    if (is.na(t_workers) || t_workers < 1L) {
      t_workers <- 1L
    }

    t_workers <- min(as.integer(t_workers), as.integer(tmax))

    res_by_t <- parallel::mclapply(
      X           = seq_len(tmax),
      FUN         = fit_one_t,
      mc.cores    = t_workers,
      mc.set.seed = TRUE
    )
  } else {
    res_by_t <- lapply(seq_len(tmax), fit_one_t)
  }

  sl_chunks     <- list()
  diag_chunks   <- list()
  change_chunks <- list()

  for (res in res_by_t) {
    if (is.null(res)) next
    if (inherits(res, "try-error")) {
      stop(sprintf("density_ratio: worker failed -- %s",
                   conditionMessage(attr(res, "condition"))),
           call. = FALSE)
    }
    Rt_t[res$idx] <- res$r_t
    if (!is.null(res$sl))     sl_chunks[[length(sl_chunks) + 1L]]         <- res$sl
    if (!is.null(res$diag))   diag_chunks[[length(diag_chunks) + 1L]]     <- res$diag
    if (!is.null(res$change)) change_chunks[[length(change_chunks) + 1L]] <- res$change
  }

  sl_summary <- if (length(sl_chunks)) data.table::rbindlist(sl_chunks, use.names = TRUE, fill = TRUE) else NULL

  policy_change_summary <- if (length(change_chunks)) {
    data.table::rbindlist(change_chunks, use.names = TRUE, fill = TRUE)[order(t)]
  } else {
    NULL
  }

  pred_dt <- if (length(diag_chunks)) {
    data.table::rbindlist(diag_chunks, use.names = TRUE, fill = TRUE)
  } else {
    NULL
  }

  is_wb_path <- isTRUE(dr_sl) && !is.null(pred_dt) &&
                "r_raw" %in% names(pred_dt) && !all(is.na(pred_dt[["r_raw"]]))

  density_calibration <- if (!is.null(pred_dt)) {
    if (is_wb_path) {
      pred_dt[
        row_type == "natural",
        .(
          n      = .N,
          mean_r = mean(r_raw, na.rm = TRUE),
          r_min  = min(r_raw, na.rm = TRUE),
          r_q01  = as.numeric(stats::quantile(r_raw, 0.01, na.rm = TRUE)),
          r_q05  = as.numeric(stats::quantile(r_raw, 0.05, na.rm = TRUE)),
          r_q50  = as.numeric(stats::quantile(r_raw, 0.50, na.rm = TRUE)),
          r_q95  = as.numeric(stats::quantile(r_raw, 0.95, na.rm = TRUE)),
          r_q99  = as.numeric(stats::quantile(r_raw, 0.99, na.rm = TRUE)),
          r_max  = max(r_raw, na.rm = TRUE)
        ),
        by = .(t)
      ][order(t)]
    } else {
      pred_dt[
        ,
        .(
          n             = .N,
          shifted_rate  = mean(S, na.rm = TRUE),
          mean_pred     = mean(p_raw, na.rm = TRUE),
          brier         = mean((S - p_raw)^2, na.rm = TRUE),
          auc           = auc(S, p_raw),
          cal_intercept = calibration_intercept(S, p_raw),
          cal_slope     = calibration_slope(S, p_raw),
          p_min = min(p_raw, na.rm = TRUE),
          p_q01 = as.numeric(stats::quantile(p_raw, 0.01, na.rm = TRUE)),
          p_q05 = as.numeric(stats::quantile(p_raw, 0.05, na.rm = TRUE)),
          p_q50 = as.numeric(stats::quantile(p_raw, 0.50, na.rm = TRUE)),
          p_q95 = as.numeric(stats::quantile(p_raw, 0.95, na.rm = TRUE)),
          p_q99 = as.numeric(stats::quantile(p_raw, 0.99, na.rm = TRUE)),
          p_max = max(p_raw, na.rm = TRUE)
        ),
        by = .(t)
      ][order(t)]
    }
  } else {
    NULL
  }

  density_calibration_by_change <- if (!is.null(pred_dt) && "changed" %in% names(pred_dt)) {
    if (is_wb_path) {
      pred_dt[
        row_type == "natural",
        .(
          n      = .N,
          mean_r = mean(r_raw, na.rm = TRUE),
          r_min  = min(r_raw, na.rm = TRUE),
          r_q01  = as.numeric(stats::quantile(r_raw, 0.01, na.rm = TRUE)),
          r_q05  = as.numeric(stats::quantile(r_raw, 0.05, na.rm = TRUE)),
          r_q50  = as.numeric(stats::quantile(r_raw, 0.50, na.rm = TRUE)),
          r_q95  = as.numeric(stats::quantile(r_raw, 0.95, na.rm = TRUE)),
          r_q99  = as.numeric(stats::quantile(r_raw, 0.99, na.rm = TRUE)),
          r_max  = max(r_raw, na.rm = TRUE)
        ),
        by = .(t, changed)
      ][order(t, changed)]
    } else {
      pred_dt[
        ,
        .(
          n             = .N,
          shifted_rate  = mean(S, na.rm = TRUE),
          mean_pred     = mean(p_raw, na.rm = TRUE),
          brier         = mean((S - p_raw)^2, na.rm = TRUE),
          auc           = auc(S, p_raw),
          cal_intercept = calibration_intercept(S, p_raw),
          cal_slope     = calibration_slope(S, p_raw),
          p_min = min(p_raw, na.rm = TRUE),
          p_q01 = as.numeric(stats::quantile(p_raw, 0.01, na.rm = TRUE)),
          p_q05 = as.numeric(stats::quantile(p_raw, 0.05, na.rm = TRUE)),
          p_q50 = as.numeric(stats::quantile(p_raw, 0.50, na.rm = TRUE)),
          p_q95 = as.numeric(stats::quantile(p_raw, 0.95, na.rm = TRUE)),
          p_q99 = as.numeric(stats::quantile(p_raw, 0.99, na.rm = TRUE)),
          p_max = max(p_raw, na.rm = TRUE)
        ),
        by = .(t, changed)
      ][order(t, changed)]
    }
  } else {
    NULL
  }

  density_pair_summary <- if (!is.null(pred_dt) &&
                              all(c("row_type", "pair_id", "changed") %in% names(pred_dt))) {
    val_var <- if (is_wb_path) "r_raw" else "p_raw"
    pred_pair <- data.table::dcast(
      pred_dt[changed == TRUE],
      t + outer_fold + pair_id ~ row_type,
      value.var = val_var
    )

    if (all(c("natural", "shifted") %in% names(pred_pair))) {
      pred_pair[, delta := shifted - natural]
      pred_pair[
        ,
        .(
          n                 = .N,
          mean_delta        = mean(delta, na.rm = TRUE),
          median_delta      = median(delta, na.rm = TRUE),
          q05_delta         = as.numeric(stats::quantile(delta, 0.05, na.rm = TRUE)),
          q95_delta         = as.numeric(stats::quantile(delta, 0.95, na.rm = TRUE)),
          prop_shift_gt_nat = mean(delta > 0, na.rm = TRUE)
        ),
        by = t
      ][order(t)]
    } else {
      NULL
    }
  } else {
    NULL
  }

  Rt_t[!is.finite(Rt_t)] <- 1

  idx_w <- DT[[time]] >= 1L & DT[[time]] <= tmax

  weights_dt <- data.table::data.table(
    ..id        = DT[[id]][idx_w],
    ..time      = DT[[time]][idx_w],
    Rt_t        = Rt_t[idx_w],
    global_fold = as.integer(DT$global_fold[idx_w])
  )

  data.table::setnames(
    weights_dt,
    c("..id", "..time"),
    c(id, time)
  )

  weights_dt[, global_fold := as.integer(global_fold)]

  out <- list(
    weights_dt                    = weights_dt,
    sl_summary                    = sl_summary,
    pred_dt                       = pred_dt,
    density_calibration           = density_calibration,
    density_calibration_by_change = density_calibration_by_change,
    density_pair_summary          = density_pair_summary,
    policy_change_summary         = policy_change_summary,
    settings = list(
      id         = id,
      time       = time,
      tmax       = tmax,
      k          = k,
      dr_sl      = dr_sl,
      bounds     = bounds,
      v          = v,
      inner_v    = inner_v,
      seed       = seed,
      a_names    = as.character(a_names),
      baseline   = as.character(baseline  %||% character(0)),
      tv_names   = as.character(tv_names  %||% character(0)),
      cluster    = cluster,
      created_at = Sys.time()
    )
  )
  class(out) <- c("density_ratio_fit", "list")
  out
}

