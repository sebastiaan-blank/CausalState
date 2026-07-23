data_check <- function(data, id, time, a_names, tv_names,
                       bs = NULL, y_col = NULL, alive_col = NULL, in_state_col = NULL,
                       tmin = 1, tmax = NULL) {

  `%||%` <- function(x, y) if (is.null(x)) y else x


  bs <- bs %||% character(0)
  y_col <- y_col %||% character(0)
  alive_col <- alive_col %||% character(0)
  in_state_col <- in_state_col %||% character(0)

  role_list <- list(
    id             = id,
    time           = time,
    treatment      = a_names,
    `time-varying` = tv_names,
    baseline       = bs,
    outcome        = y_col,
    alive          = alive_col,
    in_state            = in_state_col
  )
  role_list <- role_list[lengths(role_list) > 0]

  rn <- names(role_list)
  overlaps <- character(0)
  for (i in seq_along(role_list)) {
    for (j in seq_len(i - 1L)) {
      common <- intersect(role_list[[i]], role_list[[j]])
      if (length(common))
        overlaps <- c(overlaps, sprintf("%s & %s: %s",
                                        rn[j], rn[i], paste(common, collapse = ", ")))
    }
  }
  if (length(overlaps)) {
    stop(sprintf("Columns assigned to more than one role:\n  %s",
                 paste(overlaps, collapse = "\n  ")))
  }


  required_cols <- unique(c(id, time, a_names, tv_names,
                            bs, y_col, alive_col, in_state_col))

  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Missing columns in dataset:\n  %s",
      paste(missing_cols, collapse = "\n  ")
    ))
  }

  trial_data <- data %>% dplyr::filter(.data[[time]] >= tmin)
  if (!is.null(tmax)) {
    trial_data <- trial_data %>% dplyr::filter(.data[[time]] <= tmax)
  }

  if (nrow(trial_data) == 0L) {
    stop("No rows in requested trial period.")
  }

  if (any(trial_data[[time]] != as.integer(trial_data[[time]]), na.rm = TRUE)) {
    stop("Time column must be integer-valued in the trial period.")
  }

  dup_n <- trial_data %>%
    dplyr::count(.data[[id]], .data[[time]]) %>%
    dplyr::filter(n > 1L) %>%
    nrow()

  if (dup_n > 0L) {
    stop(sprintf("Duplicated id-time rows in trial period: %d duplicated id-time cells", dup_n))
  }

  check_vars <- unique(c(a_names, tv_names, bs, y_col, alive_col, in_state_col))

  non_finite <- trial_data %>%
    dplyr::summarise(dplyr::across(dplyr::all_of(check_vars), ~ {
      if (is.numeric(.x) || is.integer(.x)) {
        sum(!is.finite(.x))
      } else {
        sum(is.na(.x))
      }
    })) %>%
    tidyr::pivot_longer(dplyr::everything()) %>%
    dplyr::filter(value > 0)

  if (nrow(non_finite) > 0) {
    stop(sprintf(
      "Non-finite or missing values in trial period [t >= %d%s]:\n  %s",
      tmin,
      if (!is.null(tmax)) sprintf(", t <= %d", tmax) else "",
      paste(sprintf("%-30s: %d", non_finite$name, non_finite$value),
            collapse = "\n  ")
    ))
  }

  cat("data_check passed\n")
  cat(sprintf("  Patients      : %d\n", dplyr::n_distinct(data[[id]])))
  cat(sprintf("  Time range    : %d to %d\n",
              min(data[[time]], na.rm = TRUE),
              max(data[[time]], na.rm = TRUE)))
  cat(sprintf("  Trial rows    : %d  (t >= %d%s)\n",
              nrow(trial_data), tmin,
              if (!is.null(tmax)) sprintf(", t <= %d", tmax) else ""))
  cat(sprintf("  Columns checked: %d  (%d baseline, %d time-varying, %d treatment)\n",
              length(check_vars), length(bs), length(tv_names), length(a_names)))
  cat("Resolved columns:\n")
  cat(sprintf("  id           : %s\n", id))
  cat(sprintf("  time         : %s\n", time))
  cat(sprintf("  treatment    : %s\n", paste(a_names, collapse = ", ")))
  cat(sprintf("  time-varying : %s\n", paste(tv_names, collapse = ", ")))
  cat(sprintf("  baseline     : %s\n", paste(bs, collapse = ", ")))
  cat(sprintf("  outcome      : %s\n", paste(y_col, collapse = ", ")))
  cat(sprintf("  alive        : %s\n", paste(alive_col, collapse = ", ")))
  cat(sprintf("  in_state          : %s\n", paste(in_state_col, collapse = ", ")))

  invisible(TRUE)
}

make_lags <- function(DT, vars, k, id, time = NULL, overwrite = FALSE) {
  if (k <= 0L || !length(vars)) return(DT)

  stopifnot(id %in% names(DT))
  if (!is.null(time)) {
    stopifnot(time %in% names(DT))
    data.table::setorderv(DT, c(id, time))
  }

  for (nm in intersect(vars, names(DT))) {
    for (jj in seq_len(k)) {
      lnm <- paste0(nm, "_lag", jj)
      if (overwrite || !lnm %in% names(DT)) {
        DT[, (lnm) := data.table::shift(get(nm), jj, type = "lag"),
           by = id]
      }
    }
  }
  DT
}


prepare_history_lags <- function(DT, a_names, tv_names, k, id,
                                  time = NULL, no_lag_vars = character(0)) {
  tv_hist <- intersect(c(a_names, tv_names), names(DT))

  tv_hist <- setdiff(tv_hist, no_lag_vars)

  if (k > 0L && length(tv_hist)) {
    DT <- make_lags(DT, tv_hist, k = k, id = id, time = time, overwrite = FALSE)
    lag_names_by_depth <- lapply(
      seq_len(k),
      function(d) intersect(paste0(tv_hist, "_lag", d), names(DT))
    )
  } else {
    lag_names_by_depth <- vector("list", length = max(0L, k))
  }

  list(DT = DT, lag_names_by_depth = lag_names_by_depth)
}


est_seed_for <- function(seed, fold = 0L, t = 0L, component, phase = "fit") {
  component_code <- switch(
    as.character(component),
    "density_ratio"       = 10L,
    "g_remain"            = 20L,
    "g_death_exit"        = 30L,
    "g_death_exit_pooled" = 31L,
    "Q_exit"              = 40L,
    "Q_exit_pooled"       = 41L,
    "Q_remain"            = 50L,
    "target"              = 60L,
    stop("Unknown seed component: ", component, call. = FALSE)
  )

  phase_code <- switch(
    as.character(phase),
    "fit"     = 0L,
    "cv_rows" = 1L,
    "sl_fit"  = 2L,
    stop("Unknown seed phase: ", phase, call. = FALSE)
  )

  out <- as.double(seed) +
    1000000 * as.integer(fold) +
    10000   * as.integer(t) +
    100     * component_code +
    phase_code

  as.integer(out %% .Machine$integer.max)
}

auc <- function(y, p) {
  ok <- is.finite(y) & is.finite(p) & !is.na(y) & !is.na(p)
  y <- y[ok]; p <- p[ok]
  if (length(unique(y)) < 2L) return(NA_real_)
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(p, ties.method = "average")
  as.numeric((sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0))
}


logit <- function(p, eps = 1e-6) {
  p <- pmin(1 - eps, pmax(eps, p))
  stats::qlogis(p)
}


calibration_intercept <- function(y, p) {
  ok <- is.finite(y) & is.finite(p) & !is.na(y) & !is.na(p)
  y <- y[ok]
  p <- p[ok]

  if (length(unique(y)) < 2L) return(NA_real_)

  lp <- logit(p)

  fit <- tryCatch(
    stats::glm(y ~ 1, family = stats::binomial(), offset = lp),
    error = function(e) NULL
  )

  if (is.null(fit)) return(NA_real_)

  as.numeric(stats::coef(fit)[1])
}


calibration_slope <- function(y, p) {
  ok <- is.finite(y) & is.finite(p) & !is.na(y) & !is.na(p)
  y <- y[ok]
  p <- p[ok]

  if (length(unique(y)) < 2L) return(NA_real_)

  lp <- logit(p)

  fit <- tryCatch(
    stats::glm(y ~ lp, family = stats::binomial()),
    error = function(e) NULL
  )

  if (is.null(fit)) return(NA_real_)

  as.numeric(stats::coef(fit)[2])
}

binom_metric <- function(y, p) {
  ok <- is.finite(y) & is.finite(p) & !is.na(y) & !is.na(p)
  y <- y[ok]
  p <- p[ok]
  
  data.table::data.table(
    n         = length(y),
    target    = if (length(y)) mean(y) else NA_real_,
    pred      = if (length(y)) mean(p) else NA_real_,
    error     = if (length(y)) mean(p) - mean(y) else NA_real_,
    brier     = if (length(y)) mean((y - p)^2) else NA_real_,
    auc       = auc(y, p),
    cal_int   = calibration_intercept(y, p),
    cal_slope = calibration_slope(y, p)
  )
}

gauss_metric <- function(y, p) {
  ok <- is.finite(y) & is.finite(p) & !is.na(y) & !is.na(p)
  y <- y[ok]
  p <- p[ok]
  
  if (!length(y)) {
    return(data.table::data.table(
      n = 0L, target = NA_real_, pred = NA_real_, error = NA_real_,
      rmse = NA_real_, mae = NA_real_, cor = NA_real_,
      cal_int = NA_real_, cal_slope = NA_real_
    ))
  }
  
  fit <- tryCatch(stats::lm(y ~ p), error = function(e) NULL)
  
  data.table::data.table(
    n      = length(y),
    target = mean(y),
    pred   = mean(p),
    error  = mean(p) - mean(y),
    mse  = mean((y - p)^2),
    rmse = sqrt(mean((y - p)^2)),
    mae  = mean(abs(y - p)),
    cor    = if (length(y) > 2L && stats::sd(y) > 0 && stats::sd(p) > 0) {
      stats::cor(y, p)
    } else NA_real_,
    cal_int = if (!is.null(fit)) unname(stats::coef(fit)[1]) else NA_real_,
    cal_slope = if (!is.null(fit) && length(stats::coef(fit)) >= 2L) {
      unname(stats::coef(fit)[2])
    } else NA_real_
  )
}

cluster_se <- function(ic_by_id, id, cluster_col) {
  stopifnot(requireNamespace("data.table", quietly = TRUE))
  DT <- data.table::as.data.table(ic_by_id)

  stopifnot(all(c(id, "ic", cluster_col) %in% names(DT)))

  DT <- DT[!is.na(get(cluster_col))]
  N  <- nrow(DT)
  if (N <= 1L) return(NA_real_)

  ic0 <- DT[["ic"]] - mean(DT[["ic"]], na.rm = TRUE)

  Sdt <- DT[, .(S = sum(ic0, na.rm = TRUE)), by = .(cl = get(cluster_col))]
  G   <- nrow(Sdt)

  if (G <= 1L) return(stats::sd(ic0, na.rm = TRUE) / sqrt(N))

  S <- Sdt$S
  sqrt((G / (G - 1)) * sum((S - mean(S))^2) / (N^2))
}

get_folds_from_weights <- function(weights_dt, id, ids,
                                       cluster_by_id = NULL, cluster_se_only = FALSE) {
  W <- data.table::as.data.table(weights_dt)
  if (!(id %in% names(W))) stop("weights_dt must contain id column.", call. = FALSE)
  if (!("global_fold" %in% names(W))) {
    stop("weights_dt must contain 'global_fold'. No fallback folding is implemented.", call. = FALSE)
  }

  data.table::set(W, j = ".__id", value = as.character(W[[id]]))

  ff <- W[, .(global_fold = data.table::first(na.omit(as.integer(global_fold)))), by = .__id]
  if (anyNA(ff$global_fold)) stop("Some ids have missing global_fold in weights_dt.", call. = FALSE)

  fold_vec <- ff$global_fold[match(ids, ff$.__id)]
  if (anyNA(fold_vec)) stop("Some ids in 'ids' not found in weights_dt (cannot map global_fold).", call. = FALSE)

  fold_ids <- sort(unique(fold_vec))
  folds <- lapply(fold_ids, function(k) {
    list(training_set = which(fold_vec != k),
         validation_set = which(fold_vec == k))
  })
  names(folds) <- paste0("fold_", fold_ids)

  list(fold_vec = fold_vec, folds = folds)
}


scale_y <- function(D, y, outcome_family, y_bounds = NULL, bounds = 1e-5) {
  fam_req  <- match.arg(outcome_family, c("binomial", "gaussian"))
  is_binom <- (fam_req == "binomial")

  bounded <- !is_binom
  if (bounded) {
    if (is.null(y_bounds)) {
      y_min <- suppressWarnings(min(D[[y]], na.rm = TRUE))
      y_max <- suppressWarnings(max(D[[y]], na.rm = TRUE))
    } else {
      stopifnot(
        is.numeric(y_bounds), length(y_bounds) == 2L,
        is.finite(y_bounds[1]), is.finite(y_bounds[2]),
        y_bounds[2] > y_bounds[1]
      )
      y_min <- as.numeric(y_bounds[1])
      y_max <- as.numeric(y_bounds[2])
    }
    y_rng <- max(1e-12, y_max - y_min)
  } else {
    y_min <- 0
    y_max <- 1
    y_rng <- 1
  }

  to_unit <- if (bounded) {
    function(x) (as.numeric(x) - y_min) / y_rng
  } else {
    as.numeric
  }

  from_unit <- if (bounded) {
    function(u) y_min + y_rng * as.numeric(u)
  } else {
    as.numeric
  }

  clip <- if (is.null(bounds) || bounds <= 0) {
    function(x) {
      x <- as.numeric(x)
      pmin(pmax(x, 0), 1)
    }
  } else {
    b <- as.numeric(bounds)
    stopifnot(is.finite(b), b >= 0, b <= 0.5)
    function(x) {
      x <- as.numeric(x)
      pmax(pmin(x, 1 - b), b)
    }
  }

  list(
    is_binom  = is_binom,
    bounded   = bounded,
    y_min     = y_min,
    y_max     = y_max,
    y_rng     = y_rng,
    to_unit   = to_unit,
    from_unit = from_unit,
    clip      = clip
  )
}

prepare_long <- function(DT, id, time, y, cluster = NULL) {
  stopifnot(requireNamespace("data.table", quietly = TRUE))
  D <- data.table::as.data.table(DT)

  id   <- as.character(id);   stopifnot(length(id) == 1L, nzchar(id))
  time <- as.character(time); stopifnot(length(time) == 1L, nzchar(time))
  y    <- as.character(y);    stopifnot(length(y) == 1L, nzchar(y))
  if (!is.null(cluster)) cluster <- as.character(cluster)

  miss <- setdiff(c(id, time, y), names(D))
  if (length(miss)) stop(sprintf("Missing column(s) in DT: %s", paste(miss, collapse = ", ")), call. = FALSE)

  data.table::setorderv(D, c(id, time))

  D[, `:=`(
    .__id   = .SD[[1L]],
    .__time = as.integer(.SD[[2L]]),
    .__y    = suppressWarnings(as.numeric(.SD[[3L]]))
  ), .SDcols = c(id, time, y)]

  if (!is.null(cluster)) D[, .__cl := get(cluster)]

  start_t <- suppressWarnings(min(D$.__time, na.rm = TRUE))
  if (!is.finite(start_t)) stop("No actionable time found (no rows).", call. = FALSE)

  list(D = D, start_t = start_t)
}


build_Y_by_id <- function(D, scale_info) {
  Y_by_id <- D[!is.na(.__y), .(Y = .__y[.N]), by = .__id]   # last non-NA
  ids  <- unique(D$.__id)
  Y_vec <- Y_by_id$Y[match(ids, Y_by_id$.__id)]

  if (any(!is.finite(Y_vec))) stop("Non-finite Y for some ids.", call. = FALSE)
  if (scale_info$bounded) Y_vec <- scale_info$to_unit(Y_vec)

  list(ids = ids, Y = as.numeric(Y_vec))
}


build_id_time_index <- function(D, ids, tmax) {
  N <- length(ids)
  row_index <- matrix(NA_integer_, nrow = N, ncol = tmax)

  m  <- match(D$.__id, ids)
  tt <- as.integer(D$.__time)
  keep <- !is.na(m) & tt >= 1L & tt <= tmax

  row_index[cbind(m[keep], tt[keep])] <- which(keep)
  row_index
}


build_density_ratios <- function(weights_dt, id, time, ids, tmax, row_index = NULL, trim = 0.99) {
  stopifnot(requireNamespace("data.table", quietly = TRUE))
  id   <- as.character(id)
  time <- as.character(time)

  WW <- data.table::copy(data.table::as.data.table(weights_dt))

  if (!"Rt_t" %in% names(WW)) {
    stop("weights_dt must contain column 'Rt_t' with per-time *instantaneous* density ratios.")
  }

  if (!is.null(trim) && trim < 1) {
    trim_val <- stats::quantile(WW$Rt_t[is.finite(WW$Rt_t)], trim, na.rm = TRUE)
    WW[, Rt_t := pmin(Rt_t, trim_val)]
  }

  data.table::set(WW, j = "id_std",   value = as.character(WW[[id]]))
  data.table::set(WW, j = "time_std", value = as.integer(WW[[time]]))

  WW <- WW[time_std >= 1L & time_std <= tmax]

  WW <- WW[, .(Rt_t = data.table::first(Rt_t)), by = .(id_std, time_std)]

  wide <- data.table::dcast(
    WW,
    formula   = id_std ~ time_std,
    value.var = "Rt_t",
    fill      = NA_real_
  )

  id_tab <- data.table::data.table(id_std = as.character(ids))
  wide   <- merge(id_tab, wide, by = "id_std", all.x = TRUE, sort = FALSE)

  dens <- matrix(NA_real_, nrow = nrow(wide), ncol = tmax)

  time_cols <- setdiff(names(wide), "id_std")
  for (tt_name in time_cols) {
    tt <- suppressWarnings(as.integer(tt_name))
    if (!is.na(tt) && tt >= 1L && tt <= tmax) {
      dens[, tt] <- as.numeric(wide[[tt_name]])
    }
  }

  bad <- !is.finite(dens) | dens <= 0
  dens[bad] <- NA_real_

  if (!is.null(row_index)) {
    stopifnot(nrow(row_index) == nrow(dens), ncol(row_index) >= tmax)
    observed_mask <- !is.na(row_index[, seq_len(tmax), drop = FALSE])

    miss_obs <- observed_mask & is.na(dens)
    if (any(miss_obs)) {
      warning(sprintf(
        "Missing/invalid Rt_t for %d observed idxtime cells; setting those to 1.",
        sum(miss_obs)
      ), call. = FALSE)
      dens[miss_obs] <- 1
    }
  }

  dens[is.na(dens)] <- 1

  rownames(dens) <- as.character(ids)
  colnames(dens) <- paste0("t", seq_len(tmax))

  dens
}

make_cols <- function(tt, baseline, tv_names, a_names, k, all_names, t_min = 1L) {
  cols <- intersect(baseline, all_names)

  now  <- intersect(unique(c(tv_names, a_names)), all_names)
  cols <- unique(c(cols, now))

  use <- unique(c(tv_names, a_names))
  kk  <- max(0L, min(as.integer(k), as.integer(tt - t_min)))
  if (kk > 0L && length(use)) {
    lag_cols <- unlist(lapply(seq_len(kk), function(j) paste0(use, "_lag", j)))
    cols <- c(cols, lag_cols)
  }

  unique(intersect(cols, all_names))
}

make_design <- function(D, rows, cols) {
  stopifnot(requireNamespace("data.table", quietly = TRUE))
  n <- length(rows)
  if (!n) return(data.frame())

  if (!length(cols)) {
    out <- data.table::data.table(`(Intercept)` = rep(1, n))
    data.table::setDF(out)
    return(out)
  }

  out <- D[rows, ..cols]
  data.table::setDF(out)  # by reference; avoids as.data.frame copy
  out
}


sl_meta <- function(sl_res, fold, t, component, n_train, n_val) {
  if (is.null(sl_res) || isTRUE(sl_res$used_const) ||
      inherits(sl_res$fit, "sdr_const_fit") || is.null(sl_res$fit$coef)) {
    Veff <- if (!is.null(sl_res) && !is.null(sl_res$v_eff)) sl_res$v_eff else NA_integer_
    return(data.table::data.table(
      fold = fold, t = t, component = component,
      learner = "const", weight = 1,
      cv_risk = NA_real_,
      v_eff = as.integer(Veff),
      n_train = as.integer(n_train),
      n_val   = as.integer(n_val),
      used_const = TRUE
    ))
  }

  coef <- sl_res$fit$coef
  learners <- names(coef)

  cv_r <- rep(NA_real_, length(learners))
  if (!is.null(sl_res$fit$cvRisk)) {
    cv_vec <- as.numeric(sl_res$fit$cvRisk)
    names(cv_vec) <- names(sl_res$fit$cvRisk)
    cv_r <- cv_vec[learners]
  }

  data.table::data.table(
    fold = fold, t = t, component = component,
    learner = learners,
    weight  = as.numeric(coef),
    cv_risk = as.numeric(cv_r),
    v_eff   = as.integer(sl_res$v_eff),
    n_train = as.integer(n_train),
    n_val   = as.integer(n_val),
    used_const = FALSE
  )
}



apply_absorb_branch <- function(q_vec, D_block, t, branch, scale_info, absorb_rules) {
  if (!length(absorb_rules)) return(as.numeric(q_vec))

  env <- list2env(as.list(D_block), parent = parent.frame())
  assign("t", t, envir = env)
  assign("branch", branch, envir = env)

  out <- as.numeric(q_vec)

  for (rule in absorb_rules) {
    br <- if (!is.null(rule$branch)) rule$branch else "any"
    if (!(br %in% c("any", branch))) next

    cond <- eval(rule$cond, envir = env)
    if (length(cond) == 1L) cond <- rep(isTRUE(cond), nrow(D_block))
    cond <- as.logical(cond)
    cond[is.na(cond)] <- FALSE
    if (!any(cond)) next

    val <- eval(rule$value, envir = env)
    if (length(val) == 1L) val <- rep(as.numeric(val), nrow(D_block))
    val <- as.numeric(val)

    if (isTRUE(scale_info$bounded)) val <- scale_info$to_unit(val)
    val <- scale_info$clip(val)

    out[cond] <- val[cond]
  }

  out
}

#' Define an absorbing-state override rule
#'
#' Constructs a rule that overrides Q-model predictions for subjects in a
#' specific absorbing state. Pass one or more rules to the `absorb` argument
#' of [sdr()] or [itmle()].
#'
#' @param cond An unquoted expression (evaluated row-wise in the long-format
#'   data) that identifies rows to which the rule applies. Variables from the
#'   dataset and the special symbols `t` (current time) and `branch`
#'   (`"death"` or `"dc"`) are in scope.
#' @param value An unquoted expression giving the outcome value to assign when
#'   `cond` is `TRUE`. May reference data variables or constants. Scaled and
#'   clipped automatically via `scale_info`.
#' @param branch Which exit branch the rule applies to: `"any"` (both death
#'   and discharge, the default), `"death"`, or `"dc"` (discharge).
#'
#' @return A named list with elements `branch`, `cond`, and `value` (the latter
#'   two as unevaluated expressions via [base::substitute()]).
#'
#' @examples
#' # Force outcome to 0 for patients who die at any time
#' absorb_rule(alive == 0, value = 0, branch = "death")
#'
#' # Force outcome to 1 for discharged patients when t >= 5
#' absorb_rule(in_state == 0 & t >= 5, value = 1, branch = "dc")
#'
#' @export
absorb_rule <- function(cond, value, branch = c("any", "death", "dc")) {
  list(
    branch = match.arg(branch),
    cond   = substitute(cond),
    value  = substitute(value)
  )
}

align_cols <- function(X, cn) {
  X <- as.data.frame(X)
  miss <- setdiff(cn, names(X))
  if (length(miss)) {
    for (m in miss) {
      fill <- if (m == "(Intercept)") 1 else 0
      X[[m]] <- rep_len(fill, nrow(X))
    }
  }
  X <- X[, cn, drop = FALSE]
  for (j in seq_along(X)) X[[j]] <- as.numeric(X[[j]])
  X
}

can_fit_bin <- function(y, min_n = 30L, min_events = 5L) {
  y <- y[is.finite(y)]
  n <- length(y)
  if (n == 0L) return(list(ok = FALSE, p = NA_real_))
  if (n == 1L) return(list(ok = FALSE, p = as.numeric(y[1L] == 1L)))
  p1 <- sum(y == 1L)
  p0 <- n - p1
  if (n < min_n) {
    if (p1 >= min_events && p0 >= min_events) return(list(ok = TRUE, p = NA_real_))
    if (p1 == 0L) return(list(ok = FALSE, p = 0))
    if (p0 == 0L) return(list(ok = FALSE, p = 1))
    return(list(ok = FALSE, p = p1 / n))
  }
  if (p1 < min_events) return(list(ok = FALSE, p = p1 / n))  # near-zero event rate
  if (p0 < min_events) return(list(ok = FALSE, p = p1 / n))  # near-one event rate
  list(ok = TRUE, p = NA_real_)
}

rep_const <- function(p, n) rep(as.numeric(p), n)

patch_shifted_design <- function(X_nat, D_shifted_tt, rows, a_names, tt, k, t_min) {
  X_shf <- X_nat
  
  affected <- intersect(a_names, names(X_nat))
  
  if (length(affected)) {
    for (nm in affected) {
      X_shf[[nm]] <- D_shifted_tt[[nm]][rows]
    }
  }
  
  X_shf
}

#' Cluster bootstrap standard error from influence curve values
#'
#' Estimates a standard error by resampling clusters (with replacement) from
#' the influence curve contributions. This is a post-estimation utility -- it
#' operates on the IC values already stored in the estimator output, not on
#' the full dataset.
#'
#' @param ic Numeric vector of per-subject influence curve values.
#' @param cl Vector (same length as `ic`) of cluster identifiers. Subjects
#'   with `NA` cluster are dropped.
#' @param B Integer number of bootstrap replications. Default `2000L`.
#'
#' @return A single numeric value: the bootstrap standard error of the
#'   cluster-mean influence curve.
#'
#' @noRd
cluster_boot_se <- function(ic, cl, B = 2000L) {
  ok <- is.finite(ic) & !is.na(cl)
  ic <- ic[ok]; cl <- cl[ok]

  dt <- data.table::data.table(cl = cl, ic = ic)
  S  <- dt[, .(S = sum(ic), n = .N), by = cl]
  G  <- nrow(S)
  if (G < 2L) stop("Need >=2 clusters")

  n_tot <- sum(S$n)

  boot <- replicate(B, {
    idx <- sample.int(G, size = G, replace = TRUE)
    sum(S$S[idx]) / n_tot
  })

  stats::sd(boot)
}

extract_weights_dt <- function(weight_object) {
  if (is.data.frame(weight_object)) {
    return(weight_object)
  }
  
  if (is.list(weight_object) && "weights_dt" %in% names(weight_object)) {
    return(weight_object$weights_dt)
  }
  
  stop(
    "`weight_object` must be either a weights data.frame/data.table ",
    "or a list/object containing `weights_dt`.",
    call. = FALSE
  )
}





#

overwrite_policy_history_for_Q <- function(
    D, t, a_names, policy_spec_fun,
    id, time,
    k = 2,
    t_min = 1L
){
  stopifnot(requireNamespace("data.table", quietly = TRUE))
  id   <- as.character(id); time <- as.character(time)
  
  DT <- data.table::as.data.table(data.table::copy(D))
  
  if (!".__id" %in% names(DT))   data.table::set(DT, j = ".__id",   value = as.character(DT[[id]]))
  if (!".__time" %in% names(DT)) data.table::set(DT, j = ".__time", value = as.integer(DT[[time]]))
  
  data.table::setorderv(DT, c(".__id", ".__time"))
  
  if (!length(a_names)) return(DT)
  a_names <- intersect(a_names, names(DT))
  if (!length(a_names)) return(DT)

  rows_tt <- DT[.__time == t, which = TRUE]
  if (!length(rows_tt)) return(DT)
  
  block <- DT[rows_tt]
  spec  <- policy_spec_fun(block, t, a_names)
  if (is.null(spec) || !length(spec)) return(DT)
  
  if (is.data.frame(spec)) {
    for (a in intersect(a_names, names(spec))) {
      vals <- suppressWarnings(as.numeric(spec[[a]]))
      DT[rows_tt, (a) := vals]
    }
  } else {
    for (a in intersect(names(spec), a_names)) {
      s <- spec[[a]]
      if (is.atomic(s) && is.numeric(s)) {
        vals <- if (length(s) == 1L) rep(as.numeric(s), nrow(block)) else as.numeric(s)
        DT[rows_tt, (a) := vals]
      } else if (is.list(s)) {
        m <- if (!is.null(s$idx) && length(s$idx) == nrow(block)) as.logical(s$idx) else rep(TRUE, nrow(block))
        m[is.na(m)] <- FALSE
        if (!any(m)) next
        val_any <- if (!is.null(s$value)) s$value else if (!is.null(s$mean)) s$mean else s$setto
        vals <- if (length(val_any) == 1L) rep(as.numeric(val_any), nrow(block)) else as.numeric(val_any)
        idx  <- rows_tt[which(m)]
        DT[idx, (a) := vals[m]]
      }
    }
  }

  DT
}


par_n_workers <- function(workers, n, parallel = FALSE) {
  if (!isTRUE(parallel)) return(NULL)
  if (is.null(workers) || is.na(workers) || workers <= 1L) return(NULL)
  min(as.integer(workers), as.integer(n))
}

par_lapply <- function(X, FUN, workers = NULL, parallel = FALSE, seed = TRUE) {
  ww <- par_n_workers(workers, length(X), parallel)
  if (is.null(ww)) return(lapply(X, FUN))
  parallel::mclapply(
    X, FUN,
    mc.cores          = ww,
    mc.set.seed       = seed,
    mc.preschedule    = FALSE,
    mc.allow.recursive = TRUE
  )
}


mean_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) mean(x) else NA_real_
}

sd_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) >= 2L) stats::sd(x) else NA_real_
}

q_or_na <- function(x, p) {
  x <- x[is.finite(x)]
  if (length(x)) stats::quantile(x, probs = p, na.rm = TRUE, names = FALSE) else NA_real_
}

min_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) min(x) else NA_real_
}

max_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) max(x) else NA_real_
}

rmse_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) sqrt(mean(x^2)) else NA_real_
}

mae_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) mean(abs(x)) else NA_real_
}


expand_to_horizon <- function(DT, id, time, alive, in_state, tmax) {
  stopifnot(requireNamespace("data.table", quietly = TRUE))
  DT <- data.table::as.data.table(data.table::copy(DT))
  data.table::setorderv(DT, c(id, time))

  time_nm     <- time
  alive_nm    <- alive
  in_state_nm <- in_state

  last_dt <- DT[, {
    lt <- max(get(time_nm), na.rm = TRUE)
    rr <- .SD[get(time_nm) == lt][.N]
    list(last_time    = lt,
         alive_last   = as.integer(rr[[alive_nm]]),
         in_state_last = as.integer(rr[[in_state_nm]]))
  }, by = c(id)]

  grid <- data.table::CJ(
    tmp_id = unique(DT[[id]]),
    tmp_t  = seq_len(tmax),
    unique = TRUE
  )
  data.table::setnames(grid, c("tmp_id", "tmp_t"), c(id, time))

  full <- merge(grid, DT,      by = c(id, time), all.x = TRUE, sort = FALSE)
  full <- merge(full, last_dt, by = id,           all.x = TRUE, sort = FALSE)

  t_vec      <- full[[time]]
  alive_vec  <- full[[alive]]
  in_state_vec <- full[[in_state]]
  lt_vec     <- full[["last_time"]]

  bad_gap <- full[(t_vec <= lt_vec) & (is.na(alive_vec) | is.na(in_state_vec))]
  if (nrow(bad_gap)) {
    stop("Found missing alive/in_state BEFORE last_time for some ids (gaps in long data). Fix upstream.",
         call. = FALSE)
  }

  after <- t_vec > lt_vec
  full[after, (alive)    := alive_last]
  full[after, (in_state) := in_state_last]
  full[, c("alive_last", "in_state_last") := NULL]

  data.table::setorderv(full, c(id, time))
  full
}


# ------------------------------------------------------------------
# contrast -- risk difference, risk ratio, and odds ratio from two fitted estimators
# ------------------------------------------------------------------

#' Compute contrasts between two fitted estimators
#'
#' Takes an intervention-arm fit and a reference (from [sdr()] or [itmle()], or the
#' crude observed mean) and returns the risk difference (RD), risk ratio (RR), and
#' odds ratio (OR) with standard errors derived from the per-subject influence curves.
#'
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
#'     \item{\code{RD}, \code{se_RD}, \code{ci_RD}}{Risk difference and 95\% Wald CI.}
#'     \item{\code{RR}, \code{se_log_RR}, \code{ci_RR}}{Risk ratio, SE on log scale,
#'       95\% CI (exponentiated from log scale). \code{NULL} for Gaussian outcomes.}
#'     \item{\code{OR}, \code{se_log_OR}, \code{ci_OR}}{Odds ratio, SE on log scale,
#'       95\% CI (exponentiated). \code{NULL} for Gaussian outcomes.}
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
#' sim_ex <- function(n = 800L, tmax = 5L) {
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
    psi0 <- mean(dt_y$y)
    ic0  <- dt_y[, .(id, ic0 = y - psi0)]
  } else {
    psi0 <- fit0$psi
    ic0  <- data.table::as.data.table(fit0$ic_df)[, .(id, ic0 = ic)]
  }

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
    cl_ok <- cl[ok]
    dt    <- data.table::data.table(cl = cl_ok, ic = ic_vec[ok])
    S     <- dt[, .(S = sum(ic)), by = cl]
    G     <- nrow(S)
    n_ok  <- sum(ok)
    if (G < 2L) return(stats::sd(ic_vec[ok]) / sqrt(n_ok))
    Sbar  <- mean(S$S)
    sqrt((G / (G - 1L)) * sum((S$S - Sbar)^2) / n_ok^2)
  }

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
  se0     <- if (obs_ref) .se(merged$ic0) else fit0$se

  if (gaussian) {
    tbl <- data.frame(
      Estimate   = c(psi1, psi0, rd),
      `CI Lower` = c(psi1 - 1.96 * fit1$se, psi0 - 1.96 * se0, rd - 1.96 * se_rd),
      `CI Upper` = c(psi1 + 1.96 * fit1$se, psi0 + 1.96 * se0, rd + 1.96 * se_rd),
      check.names = FALSE,
      row.names   = c("Intervention", ref_lbl, "Mean Difference")
    )
  } else {
    tbl <- data.frame(
      Estimate   = c(psi1, psi0, rd, rr, or),
      `CI Lower` = c(
        psi1 - 1.96 * fit1$se,
        psi0 - 1.96 * se0,
        rd   - 1.96 * se_rd,
        exp(log(rr) - 1.96 * se_log_rr),
        exp(log(or) - 1.96 * se_log_or)
      ),
      `CI Upper` = c(
        psi1 + 1.96 * fit1$se,
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
