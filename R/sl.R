
make_validRows_cluster_binom <- function(Y01, cluster, v) {
  g <- as.character(cluster)
  dt <- data.table::data.table(g = g, y = Y01)
  G  <- dt[, .(n = .N, y1 = sum(y == 1L)), by = g]
  data.table::setorder(G, -y1, -n, g)  # deterministic order

  fold_y1 <- integer(v)
  fold_n  <- integer(v)
  fold_groups <- vector("list", v)

  for (i in seq_len(nrow(G))) {
    fi <- order(fold_y1, fold_n)[1]  # deterministic tie-break
    fold_groups[[fi]] <- c(fold_groups[[fi]], G$g[i])
    fold_y1[fi] <- fold_y1[fi] + G$y1[i]
    fold_n[fi]  <- fold_n[fi]  + G$n[i]
  }

  lapply(seq_len(v), function(fi) which(g %in% fold_groups[[fi]]))
}

make_validRows_cluster_gauss <- function(cluster, v) {
  g <- as.character(cluster)
  G <- sort(table(g), decreasing = TRUE)   # deterministic
  ug <- names(G)

  fold_n <- integer(v)
  fold_groups <- vector("list", v)

  for (grp in ug) {
    fi <- which.min(fold_n)
    fold_groups[[fi]] <- c(fold_groups[[fi]], grp)
    fold_n[fi] <- fold_n[fi] + as.integer(G[[grp]])
  }

  lapply(seq_len(v), function(fi) which(g %in% fold_groups[[fi]]))
}

with_seed <- function(seed, expr) {
  old <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (is.null(old)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
    } else {
      assign(".Random.seed", old, envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(as.integer(seed))
  force(expr)
}

make_validRows_iid_strat_binom <- function(Y01, v, seed = 1L) {

  idx0 <- which(Y01 == 0L)
  idx1 <- which(Y01 == 1L)

  idx0 <- with_seed(seed + 1000L, sample(idx0))
  idx1 <- with_seed(seed + 2000L, sample(idx1))

  folds <- vector("list", v)

  if (length(idx1)) for (k in seq_along(idx1)) folds[[ (k - 1L) %% v + 1L ]] <- c(folds[[ (k - 1L) %% v + 1L ]], idx1[k])
  if (length(idx0)) for (k in seq_along(idx0)) folds[[ (k - 1L) %% v + 1L ]] <- c(folds[[ (k - 1L) %% v + 1L ]], idx0[k])

  lapply(folds, sort)
}

make_validRows_iid <- function(n, v, seed = 1L) {
  idx <- with_seed(seed + 3000L, sample.int(n))
  cuts <- cut(seq_len(n), breaks = v, labels = FALSE)
  split(idx, cuts)
}


cap_v_binom_cluster_events <- function(Y01, cluster, v) {
  g <- as.character(cluster)
  dt <- data.table::data.table(g = g, y = Y01)
  G  <- dt[, .(y1 = any(y == 1L), y0 = any(y == 0L)), by = g]
  v2 <- min(as.integer(v), nrow(G), sum(G$y1), sum(G$y0))
  as.integer(v2)
}

train_splits_ok_binom <- function(Y01, validRows, min_events) {
  n <- length(Y01)
  for (vr in validRows) {
    tr <- setdiff(seq_len(n), vr)
    tab <- table(Y01[tr])
    if (length(tab) < 2L) return(FALSE)
    if (any(tab < min_events)) return(FALSE)
  }
  TRUE
}

method.NNloglik_safe <- SuperLearner::method.NNloglik()

local({
  original_computeCoef <- method.NNloglik_safe$computeCoef
  method.NNloglik_safe$computeCoef <<- function(Z, Y, libraryNames, verbose, obsWeights, control,
                                                ...,
                                                .orig = original_computeCoef) {
    if (ncol(Z) == 1L) {
      return(list(coef = setNames(1, libraryNames), cvRisk = NA, optimizer = "single_learner"))
    }
    Z <- pmin(pmax(Z, 1e-4), 1 - 1e-4)
    col_vars <- apply(Z, 2, function(x) var(x, na.rm = TRUE))
    row_vars <- apply(Z, 1, function(x) var(x, na.rm = TRUE))
    col_vars[!is.finite(col_vars)] <- 0
    row_vars[!is.finite(row_vars)] <- 0
    if (all(col_vars < 1e-6) || all(row_vars < 1e-6)) {
      message(sprintf("[method.NNloglik_safe][pid=%d] degenerate Z -> NNLS fallback", Sys.getpid()))
      return(SuperLearner::method.NNLS()$computeCoef(Z, Y, libraryNames, verbose, obsWeights, control, ...))
    }
    tryCatch(
      .orig(Z, Y, libraryNames, verbose, obsWeights, control, ...),
      error = function(e) {
        msg <- sprintf("[method.NNloglik_safe][pid=%d] optim failed -> NNLS fallback: %s",
                       Sys.getpid(), conditionMessage(e))
        message(msg)
        try(cat(paste0(format(Sys.time(), "%F %T"), " ", msg, "\n"),
                file = .SL_WARN_FILE,
                append = TRUE), silent = TRUE)
        SuperLearner::method.NNLS()$computeCoef(Z, Y, libraryNames, verbose, obsWeights, control, ...)
      }
    )
  }
})

#
#
#
method.WB_dr <- function(dr_floor = 1e-10) {
  list(
    require = NULL,
    computeCoef = function(Z, Y, libraryNames, verbose, obsWeights, control, ...) {
      Z <- as.matrix(Z)
      K <- ncol(Z)
      if (K == 1L) {
        return(list(
          cvRisk = setNames(NA_real_, libraryNames),
          coef   = setNames(1, libraryNames)
        ))
      }
      Z[] <- pmax(Z, dr_floor)
      s2  <- 1 - 2 * as.numeric(Y)
      cvRisk <- vapply(seq_len(K), function(k) {
        mean(s2 * log(Z[, k]))
      }, numeric(1L))
      names(cvRisk) <- libraryNames
      loss_fn <- function(alpha) {
        b       <- exp(alpha - max(alpha))
        b       <- b / sum(b)
        psi_bar <- pmax(as.vector(Z %*% b), dr_floor)
        sum(s2 * log(psi_bar))
      }
      opt <- tryCatch(
        stats::optim(rep(0, K), loss_fn, method = "BFGS",
                     control = list(maxit = 500L, reltol = 1e-8)),
        error = function(e) {
          warning(sprintf("[method.WB_dr] BFGS error: %s -- using equal weights",
                          conditionMessage(e)))
          NULL
        }
      )
      if (is.null(opt) || opt$convergence > 1L) {
        if (!is.null(opt)) {
          warning(sprintf(
            "[method.WB_dr] BFGS did not converge (code %d) -- using equal weights",
            opt$convergence))
        }
        beta <- rep(1 / K, K)
      } else {
        a    <- opt$par
        beta <- exp(a - max(a))
        beta <- beta / sum(beta)
      }
      names(beta) <- libraryNames
      list(cvRisk = cvRisk, coef = beta)
    },
    computePred = function(predY, coef, control, ...) {
      predY    <- as.matrix(predY)
      predY[]  <- pmax(predY, dr_floor)
      as.vector(predY %*% coef)
    }
  )
}

sl_block_fit <- function(
    Y, X, family, sl_lib,
    v_inner = 5L,
    id_inner = NULL,
    cluster_inner = NULL,
    min_per_fold = 20L,
    min_events = 2L,
    seed_inner = 1L,
    seed_cv_rows = NULL,
    seed_sl_fit = NULL,
    method = SuperLearner::method.NNLS,
    sl_workers = NULL
) {
  stopifnot(requireNamespace("SuperLearner", quietly = TRUE))

  X <- as.data.frame(X)
  n <- length(Y)

  if (is.null(seed_cv_rows)) seed_cv_rows <- seed_inner
  if (is.null(seed_sl_fit))  seed_sl_fit  <- seed_inner + 4000L

  seed_inner   <- as.integer(seed_inner)
  seed_cv_rows <- as.integer(seed_cv_rows)
  seed_sl_fit  <- as.integer(seed_sl_fit)

  if (is.null(cluster_inner) && !is.null(id_inner)) {
    cluster_inner <- id_inner
  }

  is_binom <- identical(family$family, "binomial")
  if (is_binom) {
    Y01 <- as.integer(Y > 0)
  }

  if (!is.null(cluster_inner)) {
    G <- length(unique(cluster_inner))
    v_eff <- min(
      as.integer(v_inner),
      as.integer(G),
      as.integer(floor(n / min_per_fold))
    )
  } else {
    v_eff <- min(
      as.integer(v_inner),
      as.integer(n),
      as.integer(floor(n / min_per_fold))
    )
  }

  if (is.na(v_eff) || v_eff < 2L) {
    stop(sprintf("sl_block_fit: v_eff = %s < 2 (n=%d)", v_eff, n), call. = FALSE)
  }

  validRows <- NULL
  v_use     <- v_eff

  if (!is_binom) {
    if (!is.null(cluster_inner)) {
      validRows <- make_validRows_cluster_gauss(cluster_inner, v_use)
    } else {
      validRows <- make_validRows_iid(n, v_use, seed = seed_cv_rows)
    }

    Y_use <- Y

  } else {
    tab <- table(Y01)

    if (length(tab) < 2L || any(tab < min_events)) {
      stop(
        sprintf(
          "sl_block_fit: too few events (table: %s)",
          paste(names(tab), tab, sep = "=", collapse = ", ")
        ),
        call. = FALSE
      )
    }

    if (!is.null(cluster_inner)) {
      v_use <- cap_v_binom_cluster_events(Y01, cluster_inner, v_use)

      while (v_use >= 2L) {
        vr <- make_validRows_cluster_binom(Y01, cluster_inner, v_use)

        if (train_splits_ok_binom(Y01, vr, min_events)) {
          validRows <- vr
          break
        }

        v_use <- v_use - 1L
      }
    }

    if (is.null(validRows)) {
      v_use <- min(
        as.integer(v_inner),
        as.integer(n),
        as.integer(floor(n / min_per_fold))
      )

      if (is.na(v_use) || v_use < 2L) {
        stop(sprintf("sl_block_fit: v < 2 in iid fallback (n=%d)", n), call. = FALSE)
      }

      while (v_use >= 2L) {
        vr <- make_validRows_iid_strat_binom(
          Y01,
          v_use,
          seed = seed_cv_rows
        )

        if (train_splits_ok_binom(Y01, vr, min_events)) {
          validRows <- vr
          break
        }

        v_use <- v_use - 1L
      }

      if (is.null(validRows)) {
        stop(
          sprintf("sl_block_fit: cannot form valid binomial folds (n=%d, v=%d)", n, v_use),
          call. = FALSE
        )
      }
    }

    Y_use <- Y01
  }

  cv_control <- SuperLearner::SuperLearner.CV.control(
    V          = as.integer(v_use),
    shuffle    = FALSE,
    stratifyCV = FALSE,
    validRows  = validRows
  )

  use_mc <- !is.null(sl_workers)
  if (use_mc) {
    old_cores <- getOption("mc.cores")
    options(mc.cores = as.integer(sl_workers))
    on.exit(options(mc.cores = old_cores), add = TRUE)
  }
  sl_fn <- if (use_mc) SuperLearner::mcSuperLearner else SuperLearner::SuperLearner

  fit_sl <- with_seed(
    seed_sl_fit,
    sl_fn(
      Y          = Y_use,
      X          = X,
      family     = family,
      SL.library = sl_lib,
      method     = method,
      cvControl  = cv_control,
      control    = SuperLearner::SuperLearner.control(saveFitLibrary = TRUE),
      id         = if (!is.null(cluster_inner)) cluster_inner else seq_len(length(Y_use)),
      verbose    = FALSE
    )
  )

  sl_tab <- if (!is.null(fit_sl$coef)) {
    data.table::data.table(
      learner = names(fit_sl$coef),
      weight  = as.numeric(fit_sl$coef),
      v_eff   = as.integer(v_use)
    )
  } else {
    NULL
  }

  fit_sl$Z                <- NULL
  fit_sl$library.predict  <- NULL
  fit_sl$SL.predict       <- NULL
  fit_sl$Y                <- NULL
  fit_sl$X                <- NULL
  fit_sl$times            <- NULL

  list(
    fit        = fit_sl,
    sl_tab     = sl_tab,
    v_eff      = as.integer(v_use),
    used_const = FALSE
  )
}

sl_predict <- function(fit, newdata) {
  n <- NROW(newdata)

  if (is.list(fit) && !inherits(fit, "SuperLearner") && !is.null(fit$fit)) {
    fit <- fit$fit
  }

  if (inherits(fit, "sdr_const_fit") || (is.list(fit) && !is.null(fit$const))) {
    return(rep_len(as.numeric(fit$const), n))
  }

  if (!inherits(fit, "SuperLearner")) {
    stop("sl_predict: 'fit' is not a SuperLearner (class = ",
         paste(class(fit), collapse = ", "), ").")
  }

  #
  if (is.null(fit$control)) fit$control <- list()
  sfl <- fit$control$saveFitLibrary
  if (!is.logical(sfl) || length(sfl) != 1L) {
    fit$control$saveFitLibrary <- TRUE  # prevent `!` on non-logical
  }

  pr <- SuperLearner::predict.SuperLearner(
    fit,
    newdata = as.data.frame(newdata),
    onlySL  = TRUE
  )
  p <- pr$pred
  if (is.vector(p)) return(as.numeric(p))
  if (is.matrix(p) || is.data.frame(p)) return(as.numeric(p[, 1L]))
  stop("Unexpected structure of SL predictions.")
}
