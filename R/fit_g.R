
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

  block_t <- DT[rows_t]
  spec    <- policy_spec_fun(block_t, t, a_names)

  A_shift <- g_shift_values(A_obs, spec, t)

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


#' Estimate treatment mechanism weights for a modified treatment policy
#'
#' Fits per-time-point treatment models and computes the instantaneous density
#' ratio \eqn{r_t = d\tilde{P}(A_t | H_t) / dP(A_t | H_t)} comparing the
#' modified treatment policy (MTP) to the natural course. The output is passed
#' directly to [sdr()] or [itmle()] as the `weight_object` argument.
#'
#' Two estimation methods are available via the `method` argument:
#'
#' \describe{
#'   \item{`"classification"` (default)}{Stacks observed and shifted treatment
#'     rows and fits a binary classifier to distinguish the two worlds. The
#'     density ratio is recovered as the classifier odds ratio. Works for
#'     binary, continuous, categorical, and multi-treatment interventions.}
#'   \item{`"propensity"`}{Directly estimates the treatment mechanism
#'     \eqn{g(A|H)} and computes the policy density ratio analytically.
#'     Supports one numeric binary treatment (0/1), or one ordinal treatment
#'     supplied as an ordered factor. Binary treatment uses one binomial
#'     SuperLearner; ordinal treatment uses a sequence of binomial
#'     SuperLearners to estimate continuation probabilities
#'     \eqn{P(A > a_j | A \ge a_j, H)}. The resulting dose probabilities
#'     are non-negative and sum to one without a proportional-odds assumption.
#'     Binary policies may return probabilities in \eqn{[0,1]}.
#'     Ordinal policies must return deterministic category labels or factors
#'     from the treatment's declared levels.}
#' }
#'
#' Cross-fitting is used throughout: treatment models are trained on one fold
#' and density ratios predicted on the held-out fold, so the same fold
#' structure is inherited by the downstream estimator.
#'
#' For ordinal treatment, the policy is evaluated at every possible source
#' dose, holding history fixed. Its induced probability is
#' \deqn{g^*(a|h) = \sum_b I\{d(b,h)=a\} g(b|h).}
#' Summing all preimages supports many-to-one policies such as dose caps and
#' minimum-dose rules. The policy must depend only on the current treatment,
#' time, and pre-treatment history; precomputed shifted-value columns cannot
#' describe how the policy acts at hypothetical doses. Zero numerator
#' probabilities produce exact zero ratios. `bounds` clips each continuation
#' probability and floors the observed-dose denominator.
#'
#' Every declared ordinal dose must have training observations in each outer
#' fold at each time where the policy changes treatment. Sparse upper-dose
#' risk sets must also meet the usual SuperLearner sample-size and event
#' requirements; fitting errors propagate rather than substituting unit ratios.
#'
#' Factor treatment columns are encoded internally as numeric reference
#' indicators, omitting the first dose, whenever treatment is a predictor.
#' The same fixed columns are used for current treatment, lagged treatment,
#' policy-shifted prediction, and downstream transition, outcome, and targeting
#' models. Levels absent from a fitting subset retain their indicator columns.
#' Policies receive the original factor columns, with their levels and ordered
#' status preserved after shifts. A numeric
#' dose term in a GLM would otherwise impose equal changes in the linear
#' predictor per dose unit; reference indicators allow separate dose effects.
#' This encoding does not impose a monotone dose response.
#'
#' Factor levels and ordered status are saved in the `g_fit` settings.
#' [sdr()] and [itmle()] check their input factor columns against this schema
#' when given the complete object. [qreg()] infers the same coding directly
#' from its input factor columns. There is no separate levels argument.
#'
#' @param df A `data.frame` in long format (one row per subject per time point).
#' @param a_names Character vector of treatment variable names. Categorical
#'   treatments must be factors with explicitly declared levels; ordinal
#'   propensity estimation requires an ordered factor. Level order determines
#'   the continuation sequence, and the first level is the predictor reference.
#'   Numeric columns retain numeric coding for binary or continuous treatment.
#' @param tmax Integer. Maximum follow-up time (number of time points).
#' @param baseline Character vector of baseline covariate names.
#' @param tv_names Character vector of time-varying covariate names.
#' @param no_lag_vars Character vector of variables that should not be lagged.
#' @param policy_names Character vector of column names holding shifted
#'   treatment values under the MTP.
#' @param sl_g SuperLearner library for treatment models. Required.
#' @param k Integer. Number of lags to include. Default `2`.
#' @param inner_v Integer. Inner cross-validation folds for SuperLearner.
#'   Default `10L`.
#' @param cluster Character. Cluster variable name for clustered fold
#'   assignment.
#' @param v Integer. Number of outer cross-fitting folds. Default `5L`.
#' @param seed Integer random seed. Default `1L`.
#' @param policy_spec_fun A function `(D_block, t, a_names)` returning a
#'   `data.table` with shifted treatment values for time `t`. For binary
#'   treatment under `method = "propensity"`, the shifted value is interpreted
#'   as the target treatment probability (values in \eqn{[0,1]}) or a
#'   deterministic assignment (0 or 1). For ordinal propensity estimation,
#'   receive an ordered factor and return deterministic category labels or a
#'   factor using its declared levels. Numeric outputs for a factor are matched
#'   to category labels, never interpreted as internal integer codes. The
#'   function is also called on hypothetical treatments at every declared level.
#' @param id Character. Subject identifier column name. Default `"id"`.
#' @param time Character. Time column name. Default `"trial_time"`.
#' @param bounds Numeric. Probability floor for clipping. Default `1e-5`.
#' @param method Character. Estimation method: `"classification"` (default) or
#'   `"propensity"`. Propensity supports one numeric binary 0/1 treatment
#'   or one ordered factor. Use classification for unordered factors,
#'   continuous treatment, or multiple treatments.
#' @param dr_sl Logical. When `method = "classification"`, use the
#'   Wu-Benkeser metalearner in density-ratio space. Default `FALSE`.
#' @param drop_small_cluster_splits Logical. Default `TRUE`.
#' @param parallel Logical. Parallelise across time points. Default `FALSE`.
#' @param t_workers Integer. Workers for time-point parallelism.
#' @param fold_workers Integer. Workers for fold parallelism.
#' @param sl_workers Integer. Workers for SuperLearner parallelism.
#' @param verbose Logical. Print per-time progress messages. Default `TRUE`.
#'
#' @return A list of class `"g_fit"` with components:
#'   \describe{
#'     \item{`weights_dt`}{A `data.table` with `Rt_t` (instantaneous density
#'       ratio) and `global_fold` per subject-time. Pass to [sdr()] or
#'       [itmle()] as `weight_object`.}
#'     \item{`sl_summary`}{SuperLearner learner weights per fold and time;
#'       ordinal models also include the continuation `threshold`.}
#'     \item{`density_calibration`}{Per-time calibration diagnostics (AUC,
#'       Brier, calibration slope). Only populated for
#'       `method = "classification"`; `NULL` for `method = "propensity"`.}
#'     \item{`settings`}{Named list of all call settings including `method`
#'       and `is_binary`, plus factor `treatment_levels`, `treatment_ordered`,
#'       `treatment_encoding`, and `propensity_model`.}
#'   }
#'
#' @seealso [sdr()], [itmle()], [balance_check()], [weight_diagnostics()]
#'
#' @references
#' Diaz I, Williams N, Hoffman KL, Schenck EJ (2021). Nonparametric Causal
#' Effects Based on Longitudinal Modified Treatment Policies. *JASA*
#' 118(542):846-857.
#'
#' @export
fit_g <- function(
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
    method = "classification",
    dr_sl  = FALSE,
    drop_small_cluster_splits = TRUE,
    parallel = FALSE,
    t_workers = NULL,
    fold_workers = NULL,
    sl_workers = NULL,
    verbose = TRUE
) {

  `%||%` <- function(x, y) if (!is.null(x)) x else y

  method <- match.arg(method, c("propensity", "classification"))

  is_propensity <- method == "propensity"

  if (!isTRUE(parallel)) {
    t_workers    <- NULL
    fold_workers <- NULL
    sl_workers   <- NULL
  }

  keep_cols <- unique(c(
    id, time, a_names,
    baseline, tv_names,
    no_lag_vars, policy_names, cluster
  ))
  keep_cols <- intersect(keep_cols, names(df))
  df <- as.data.frame(df)[, keep_cols, drop = FALSE]

  data_check(data = df, id = id, time = time,
             a_names = a_names, tv_names = tv_names,
             bs = baseline, tmin = 1, tmax = tmax,
             verbose = verbose)

  DT <- data.table::as.data.table(df)
  data.table::setorderv(DT, c(id, time))
  rm(df); gc()

  t_min <- DT[, .SD[1L], by = id][[time]] |> max()

  # Explicit support controls treatment coding without guessing from unique values.
  treatment_spec <- factor_treatment_levels(DT, a_names)
  ordinal_levels <- if (length(a_names) == 1L) treatment_spec[[a_names[1L]]] else NULL
  learner_design <- function(X) encode_treatment_design(X, treatment_spec)

  is_binary_a <- is_propensity && length(a_names) == 1L &&
    !is.factor(DT[[a_names[1L]]]) && all(na.omit(DT[[a_names[1L]]]) %in% c(0, 1))
  is_ordinal_a <- is_propensity && length(a_names) == 1L && is.ordered(DT[[a_names[1L]]])

  if (is_propensity && !is_binary_a && !is_ordinal_a) {
    stop(
      'method = "propensity" requires a single numeric binary treatment (0/1) ',
      'or an ordered factor. Declare ordinal treatment with factor(..., ordered = TRUE); ',
      'use method = "classification" for unordered factors or other treatments.',
      call. = FALSE
    )
  }
  if (is_propensity && (length(bounds) != 1L || !is.finite(bounds) || bounds <= 0 || bounds >= 0.5))
    stop("bounds must be a finite number strictly between 0 and 0.5.", call. = FALSE)

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

  method_sl <- if (isTRUE(dr_sl)) method.WB_dr() else SuperLearner::method.NNLS

  fit_one_t <- function(t) {

    sl_rows   <- vector("list", 0L)
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
    rows_idx <- stk$rows_idx
    A_cols   <- intersect(a_names, names(stk$X_obs))
    A_shft_t <- stk$X[seq_len(n_nat) + n_nat, A_cols, drop = FALSE]
    changed <- Reduce(`|`, lapply(A_cols, function(a) stk$X_obs[[a]] != A_shft_t[[a]]),
                      init = rep(FALSE, n_nat))
    no_shift <- !any(changed)

    change_row <- data.table::data.table(
      t           = as.integer(t),
      n_nat       = n_nat,
      n_changed   = sum(changed, na.rm = TRUE),
      pct_changed = 100 * mean(changed, na.rm = TRUE)
    )

    if (length(A_cols) <= 3L) {
      for (ac in A_cols) {
        # Means and numeric dose differences have no meaning for category labels.
        if (is.factor(stk$X_obs[[ac]])) next
        change_row[[paste0("mean_A_pre_",    ac)]] <- mean(stk$X_obs[[ac]], na.rm = TRUE)
        change_row[[paste0("mean_A_post_",   ac)]] <- mean(A_shft_t[[ac]], na.rm = TRUE)
        if (any(changed)) {
          d <- A_shft_t[[ac]][changed] - stk$X_obs[[ac]][changed]
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

    policy_map <- NULL
    if (is_ordinal_a) {
      policy_map <- ordinal_policy_map(DT[rows_idx], t, a_names[1L],
                                       ordinal_levels, policy_spec_fun)
      # Check identity on every possible dose, not only observed assignments.
      no_shift <- all(policy_map == rep(seq_along(ordinal_levels), each = n_nat))
    }

    if (no_shift) {
      message(sprintf("t=%d: no shift -- ratios set to 1 (n=%d)", t, n_nat))
      return(list(
        idx    = rows_idx,
        r_t    = rep(1, n_nat),
        sl     = NULL,
        diag   = NULL,
        change = change_row
      ))
    }

    # Evaluate policy at A=0 and A=1 to get marginal policy distribution
    # (binary propensity path only — deterministic or stochastic policies handled uniformly)
    d_at_1 <- NULL
    d_at_0 <- NULL
    if (is_propensity && is_binary_a) {
      a_col_b <- a_names[1L]
      bt  <- DT[rows_idx]
      ba1 <- data.table::copy(bt); ba1[[a_col_b]] <- 1L
      ba0 <- data.table::copy(bt); ba0[[a_col_b]] <- 0L
      sp1 <- policy_spec_fun(ba1, t, a_names)
      sp0 <- policy_spec_fun(ba0, t, a_names)
      d_at_1 <- pmin(1, pmax(0, as.numeric(sp1[[a_col_b]])))
      d_at_0 <- pmin(1, pmax(0, as.numeric(sp0[[a_col_b]])))
      rm(bt, ba1, ba0, sp1, sp0)
    }

    fvec_t     <- DT$global_fold[stk$rows_idx]
    folds_here <- sort(unique(fvec_t))

    r_hat <- rep(NA_real_, n_nat)

    cl_t <- if (!is.null(cluster) && cluster %in% names(DT)) {
      DT[[cluster]][stk$rows_idx]
    } else NULL

    fit_one_fold <- function(kk) {
      message(sprintf("[t=%d][fold %d][pid=%d] %s", t, kk, Sys.getpid(), format(Sys.time(), "%H:%M:%S")))
      idx_valid <- which(fvec_t == kk)
      idx_train <- which(fvec_t != kk)

      if (length(idx_valid) == 0L || length(idx_train) < 2L) {
        if (is_ordinal_a && length(idx_valid))
          stop(sprintf("Ordinal propensity at t=%d, fold=%d: fewer than two training observations.",
                       t, kk), call. = FALSE)
        return(list(idx_valid = idx_valid, r = rep(1, length(idx_valid)),
                    sl_tab = NULL, diag_tab = NULL))
      }

      rows_train_nat   <- idx_train
      rows_train_shift <- idx_train + n_nat
      rows_train       <- c(rows_train_nat, rows_train_shift)

      id_nat <- DT[[id]][stk$rows_idx][rows_train_nat]

      cl_fold <- if (!is.null(cl_t)) cl_t[rows_train_nat] else NULL

      # ---- propensity path ----
      if (is_propensity) {
        hist_cols <- setdiff(names(stk$X_obs), A_cols)
        X_hist    <- learner_design(stk$X_obs[, hist_cols, drop = FALSE])
        X_tr      <- as.data.frame(X_hist[rows_train_nat, , drop = FALSE])
        X_vl      <- as.data.frame(X_hist[idx_valid,      , drop = FALSE])

        A_obs_all <- if (is_ordinal_a) as.integer(stk$X_obs[[A_cols[1L]]]) else
          as.numeric(stk$X_obs[[A_cols[1L]]])

        if (is_ordinal_a) {
          ordinal_fit <- ordinal_propensity_fit(
            A = A_obs_all[rows_train_nat], X = X_tr, newX = X_vl,
            levels = ordinal_levels, sl_g = sl_g, inner_v = inner_v,
            ids = id_nat, cluster = cl_fold, seed = seed, fold = kk, t = t,
            bounds = bounds, sl_workers = sl_workers
          )
          r <- ordinal_policy_ratio(
            ordinal_fit$prob, policy_map[idx_valid, , drop = FALSE],
            A_obs_all[idx_valid], bounds
          )
          return(list(idx_valid = idx_valid, r = r,
                      sl_tab = ordinal_fit$sl_tab, diag_tab = NULL))
        }

        # Binary: one binomial SL on P(A=1|W)

        sl_res <- sl_block_fit(
          Y             = A_obs_all[rows_train_nat],
          X             = X_tr,
          family        = stats::binomial(),
          sl_lib        = sl_g,
          v_inner       = inner_v,
          id_inner      = if (is.null(cl_fold)) id_nat else NULL,
          cluster_inner = cl_fold,
          seed_inner    = est_seed_for(seed, fold = kk, t = t, component = "g_prop", phase = "fit"),
          seed_cv_rows  = est_seed_for(seed, fold = kk, t = t, component = "g_prop", phase = "cv_rows"),
          seed_sl_fit   = est_seed_for(seed, fold = kk, t = t, component = "g_prop", phase = "sl_fit"),
          method        = SuperLearner::method.NNLS,
          sl_workers    = sl_workers
        )

        p_hat <- pmin(1 - bounds, pmax(bounds,
          as.numeric(sl_predict(sl_res$fit, newdata = X_vl))
        ))

        a_obs_vl <- A_obs_all[idx_valid]

        # g*(A|W) via marginal policy distribution: P_d(A=1|W) = d1*p + d0*(1-p)
        d1       <- d_at_1[idx_valid]
        d0       <- d_at_0[idx_valid]
        g_star_1 <- d1 * p_hat + d0 * (1 - p_hat)
        g_star_0 <- 1 - g_star_1
        g_obs    <- a_obs_vl * p_hat + (1 - a_obs_vl) * (1 - p_hat)
        g_star   <- a_obs_vl * g_star_1 + (1 - a_obs_vl) * g_star_0
        r        <- g_star / pmax(g_obs, bounds)
        r[!is.finite(r)] <- 1

        sl_tab <- if (!is.null(sl_res$fit$coef)) data.table::data.table(
          t = as.integer(t), outer_fold = as.integer(kk),
          learner = names(sl_res$fit$coef),
          weight  = as.numeric(sl_res$fit$coef),
          cvRisk  = as.numeric(sl_res$fit$cvRisk)
        ) else NULL

        rm(sl_res); gc(FALSE)
        return(list(idx_valid = idx_valid, r = r, sl_tab = sl_tab, diag_tab = NULL))
      }

      # ---- classification path (original) ----
      Y_train <- stk$S[rows_train]
      X_train <- learner_design(stk$X[rows_train, , drop = FALSE])

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

      newdata_obs <- learner_design(stk$X_obs[idx_valid, , drop = FALSE])

      rows_valid_nat   <- idx_valid
      rows_valid_shift <- idx_valid + n_nat
      rows_valid_stack <- c(rows_valid_nat, rows_valid_shift)

      pred_cal <- as.numeric(sl_predict(
        fit,
        newdata = learner_design(stk$X[rows_valid_stack, , drop = FALSE])
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

      rm(sl_res, fit); gc(FALSE)
      list(idx_valid = idx_valid, r = r, sl_tab = sl_tab, diag_tab = diag_tab)
    }

    fold_res <- par_lapply(folds_here, fit_one_fold, workers = fold_workers, parallel = parallel)

    for (fr in fold_res) {
      if (inherits(fr, "try-error"))
        stop(sprintf("fit_g: fold worker crashed at t=%d -- %s",
                     t, conditionMessage(attr(fr, "condition"))),
             call. = FALSE)
      if (length(fr$idx_valid))  r_hat[fr$idx_valid]                       <- fr$r
      if (!is.null(fr$sl_tab))   sl_rows[[length(sl_rows) + 1L]]           <- fr$sl_tab
      if (!is.null(fr$diag_tab)) diag_rows[[length(diag_rows) + 1L]]       <- fr$diag_tab
    }
    rm(stk, fold_res); gc(FALSE)

    if (anyNA(r_hat))
      stop(sprintf(
        "fit_g: %d subject(s) at t=%d have no density ratio after all folds",
        sum(is.na(r_hat)), t), call. = FALSE)

    r_hat[!is.finite(r_hat)] <- 1

    sl_tab <- if (length(sl_rows)) data.table::rbindlist(sl_rows, use.names = TRUE, fill = TRUE) else NULL

    diag_tab <- if (length(diag_rows)) {
      data.table::rbindlist(diag_rows, use.names = TRUE, fill = TRUE)
    } else {
      NULL
    }

    return(list(
      idx    = rows_idx,
      r_t    = r_hat,
      sl     = sl_tab,
      diag   = diag_tab,
      change = change_row
    ))
  }

  res_by_t <- par_lapply(seq_len(tmax), fit_one_t, workers = t_workers, parallel = parallel)

  sl_chunks     <- list()
  diag_chunks   <- list()
  change_chunks <- list()

  for (res in res_by_t) {
    if (is.null(res)) next
    if (inherits(res, "try-error")) {
      stop(sprintf("fit_g: worker failed -- %s",
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

  # Calibration diagnostics — classification path only; propensity path has no pred_dt
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
      method     = method,
      is_binary  = is_binary_a,
      treatment_levels = treatment_spec,
      treatment_ordered = factor_treatment_ordered(DT, treatment_spec),
      treatment_encoding = if (length(treatment_spec)) "reference" else NULL,
      propensity_model = if (is_ordinal_a) "continuation_ratio" else if (is_propensity) "binary" else NULL,
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
  class(out) <- c("g_fit", "list")
  out
}
