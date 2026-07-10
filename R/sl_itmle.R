extract_offset <- function(X, newX) {
  off_col <- "._sl_offset"
  if (off_col %in% colnames(X)) {
    offset_tr <- X[[off_col]]
    X         <- X[, colnames(X) != off_col, drop = FALSE]
  } else {
    offset_tr <- rep(0, nrow(X))
  }
  if (off_col %in% colnames(newX)) {
    offset_pr <- newX[[off_col]]
    newX      <- newX[, colnames(newX) != off_col, drop = FALSE]
  } else {
    offset_pr <- rep(0, nrow(newX))
  }
  list(X = X, newX = newX, offset_tr = offset_tr, offset_pr = offset_pr)
}

tgt_clip <- function(p, bounds = 1e-5) {
  pmin(pmax(p, bounds), 1 - bounds)
}

#' Custom SuperLearner wrappers for the iTMLE targeting step
#'
#' @description
#' A family of SuperLearner-compatible learner functions designed for the
#' iTMLE targeting step. They differ from standard wrappers in one critical
#' way: the logit offset is passed as a **column** in `X` (named
#' `._sl_offset`) rather than via the `offset` argument. This is necessary
#' because SuperLearner's internal cross-validation subsetting drops the
#' `offset` vector, whereas a column in `X` is correctly subset.
#'
#' Each learner:
#' 1. Extracts `._sl_offset` from `X`/`newX` via `extract_offset()`.
#' 2. Fits a fluctuation model with that offset.
#' 3. Returns predictions on the probability scale, clipped to `[bounds, 1-bounds]`.
#'
#' Learners available:
#' \describe{
#'   \item{`SL.tgt.empty`}{No fluctuation: returns `expit(offset)` unchanged.
#'     Acts as the "no update" option.}
#'   \item{`SL.tgt.intercept`}{One-parameter intercept fluctuation (standard
#'     TMLE update).}
#'   \item{`SL.tgt.glm`}{Main-terms logistic GLM fluctuation with offset.}
#'   \item{`SL.tgt.glmnet`}{Penalised logistic regression (elastic net) with
#'     offset via `glmnet::cv.glmnet()`.}
#'   \item{`SL.tgt.xgboost`}{Gradient boosted trees with offset via
#'     `base_margin`. CUDA-capable.}
#' }
#'
#' Pre-configured variants (`SL.tmle_*`) expose specific hyperparameter
#' choices and are the recommended building blocks for `sl_tmle`:
#' `SL.tmle_empty`, `SL.tmle_intercept`, `SL.tmle_glm`,
#' `SL.tmle_glmnet_ridge`, `SL.tmle_glmnet_enet`, `SL.tmle_glmnet_lasso`,
#' `SL.tmle_xgb_d1`, `SL.tmle_xgb_d3`, `SL.tmle_xgb_d6`.
#'
#' The default `sl_tmle` character vector contains a recommended set of these.
#'
#' @param Y Numeric outcome vector (on `[0, 1]` after scaling).
#' @param X `data.frame` of covariates, **must** include a column named
#'   `._sl_offset` containing the logit of the current Q estimate.
#' @param newX `data.frame` for prediction, same structure as `X`.
#' @param family Passed by SuperLearner; should be `binomial()`.
#' @param obsWeights Numeric vector of observation weights (the iTMLE clever
#'   covariate / density ratio product).
#' @param id Subject identifiers (passed by SuperLearner, not used directly).
#' @param bounds Numeric. Clipping bound for predictions. Default `1e-5`.
#' @param ... Additional arguments (ignored).
#'
#' @return A list with elements `pred` (numeric predictions on the probability
#'   scale) and `fit` (a fitted object with a `predict` method).
#'
#' @name sl_itmle
#' @aliases SL.tgt.empty SL.tgt.intercept SL.tgt.glm SL.tgt.glmnet
#'   SL.tgt.xgboost SL.tmle_empty SL.tmle_intercept SL.tmle_glm
#'   SL.tmle_glmnet_ridge SL.tmle_glmnet_enet SL.tmle_glmnet_lasso
#'   SL.tmle_xgb_d1 SL.tmle_xgb_d3 SL.tmle_xgb_d6
#'
#' @seealso [itmle()]
#'
#' @export SL.tgt.empty
#' @export SL.tgt.intercept
#' @export SL.tgt.glm
#' @export SL.tgt.glmnet
#' @export SL.tgt.xgboost
#' @export SL.tmle_glm
#' @export SL.tmle_glmnet_ridge
#' @export SL.tmle_glmnet_enet
#' @export SL.tmle_glmnet_lasso
#' @export SL.tmle_xgb_d1
#' @export SL.tmle_xgb_d3
#' @export SL.tmle_xgb_d6
SL.tgt.empty <- function(Y, X, newX, family, obsWeights, id,
                         bounds = 1e-5, ...) {
  pcs <- extract_offset(X, newX)
  pred_Q <- tgt_clip(plogis(pcs$offset_pr), bounds)

  list(
    pred = pred_Q,
    fit  = structure(list(bounds = bounds), class = "SL.tgt.empty")
  )
}

predict.SL.tgt.empty <- function(object, newdata, ...) {
  pcs <- extract_offset(newdata, newdata)
  bounds <- if (!is.null(object$bounds)) object$bounds else 1e-5
  tgt_clip(plogis(pcs$offset_pr), bounds)
}

SL.tgt.intercept <- function(Y, X, newX, family, obsWeights, id,
                             bounds = 1e-5, ...) {
  pcs <- extract_offset(X, newX)
  Y_b <- tgt_clip(Y, bounds)

  fit <- tryCatch(
    suppressWarnings(stats::glm(
      Y_b ~ 1,
      family  = stats::quasibinomial(),
      offset  = pcs$offset_tr,
      weights = obsWeights
    )),
    error = function(e) NULL
  )

  eps <- if (!is.null(fit)) as.numeric(stats::coef(fit)[1L]) else 0
  if (!is.finite(eps)) eps <- 0

  pred_Q <- tgt_clip(plogis(pcs$offset_pr + eps), bounds)

  list(
    pred = pred_Q,
    fit  = structure(
      list(eps = eps, bounds = bounds),
      class = "SL.tgt.intercept"
    )
  )
}

predict.SL.tgt.intercept <- function(object, newdata, ...) {
  pcs <- extract_offset(newdata, newdata)
  bounds <- if (!is.null(object$bounds)) object$bounds else 1e-5
  tgt_clip(plogis(pcs$offset_pr + object$eps), bounds)
}

SL.tgt.glm <- function(Y, X, newX, family, obsWeights, id,
                       bounds = 1e-5, ...) {
  pcs <- extract_offset(X, newX)

  if (ncol(pcs$X) == 0L) {
    return(SL.tgt.intercept(
      Y, X, newX, family, obsWeights, id,
      bounds = bounds, ...
    ))
  }

  Y_b  <- tgt_clip(Y, bounds)
  df_t <- as.data.frame(pcs$X)
  df_t$.__y <- Y_b

  fmla <- stats::reformulate(colnames(pcs$X), response = ".__y")

  fit <- tryCatch(
    suppressWarnings(stats::glm(
      fmla,
      data    = df_t,
      family  = stats::quasibinomial(),
      offset  = pcs$offset_tr,
      weights = obsWeights
    )),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    pred_Q <- tgt_clip(plogis(pcs$offset_pr), bounds)
    return(list(
      pred = pred_Q,
      fit  = structure(
        list(bounds = bounds),
        class = "SL.tgt.glm.failed"
      )
    ))
  }

  X_pred_mat <- stats::model.matrix(
    fmla[-2L],
    data = as.data.frame(pcs$newX)
  )

  eps    <- as.numeric(X_pred_mat %*% stats::coef(fit))
  pred_Q <- tgt_clip(plogis(pcs$offset_pr + eps), bounds)

  list(
    pred = pred_Q,
    fit  = structure(
      list(
        coef    = stats::coef(fit),
        formula = fmla,
        bounds  = bounds
      ),
      class = "SL.tgt.glm"
    )
  )
}

predict.SL.tgt.glm <- function(object, newdata, ...) {
  pcs <- extract_offset(newdata, newdata)
  bounds <- if (!is.null(object$bounds)) object$bounds else 1e-5

  if (is.null(object$coef)) {
    return(tgt_clip(plogis(pcs$offset_pr), bounds))
  }

  X_pred_mat <- stats::model.matrix(
    object$formula[-2L],
    data = as.data.frame(pcs$newX)
  )

  eps <- as.numeric(X_pred_mat %*% object$coef)
  tgt_clip(plogis(pcs$offset_pr + eps), bounds)
}

SL.tgt.glmnet <- function(Y, X, newX, family, obsWeights, id,
                          alpha = 0.5, nfolds = 3L,
                          s_select = "lambda.min",
                          bounds = 1e-5, ...) {
  pcs <- extract_offset(X, newX)

  if (ncol(pcs$X) == 0L) {
    return(SL.tgt.intercept(
      Y, X, newX, family, obsWeights, id,
      bounds = bounds, ...
    ))
  }

  X_mat    <- as.matrix(pcs$X)
  newX_mat <- as.matrix(pcs$newX)
  Y_b      <- tgt_clip(Y, bounds)
  Y_mat    <- cbind(1 - Y_b, Y_b)

  fit <- tryCatch(
    suppressWarnings(glmnet::cv.glmnet(
      x       = X_mat,
      y       = Y_mat,
      family  = "binomial",
      offset  = pcs$offset_tr,
      weights = obsWeights,
      alpha   = alpha,
      nfolds  = nfolds
    )),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    pred_Q <- tgt_clip(plogis(pcs$offset_pr), bounds)
    return(list(
      pred = pred_Q,
      fit  = structure(
        list(bounds = bounds),
        class = "SL.tgt.glmnet.failed"
      )
    ))
  }

  pred_Q <- as.numeric(stats::predict(
    fit,
    newx      = newX_mat,
    newoffset = pcs$offset_pr,
    s         = s_select,
    type      = "response"
  ))

  pred_Q <- tgt_clip(pred_Q, bounds)

  list(
    pred = pred_Q,
    fit  = structure(
      list(
        object   = fit,
        s_select = s_select,
        alpha    = alpha,
        bounds   = bounds
      ),
      class = "SL.tgt.glmnet"
    )
  )
}

predict.SL.tgt.glmnet <- function(object, newdata, ...) {
  pcs <- extract_offset(newdata, newdata)
  bounds <- if (!is.null(object$bounds)) object$bounds else 1e-5

  if (is.null(object$object)) {
    return(tgt_clip(plogis(pcs$offset_pr), bounds))
  }

  newX_mat <- as.matrix(pcs$newX)

  pred_Q <- as.numeric(stats::predict(
    object$object,
    newx      = newX_mat,
    newoffset = pcs$offset_pr,
    s         = object$s_select,
    type      = "response"
  ))

  tgt_clip(pred_Q, bounds)
}

SL.tgt.xgboost <- function(Y, X, newX, family, obsWeights, id,
                           nrounds = 100L, max_depth = 2L, eta = 0.1,
                           subsample = 0.8, colsample_bytree = 0.8,
                           max_delta_step = 0,
                           nthread = 1L, use_cuda = TRUE,
                           bounds = 1e-5, ...) {
  pcs <- extract_offset(X, newX)

  if (ncol(pcs$X) == 0L) {
    return(SL.tgt.intercept(
      Y, X, newX, family, obsWeights, id,
      bounds = bounds, ...
    ))
  }

  X_mat    <- as.matrix(pcs$X)
  newX_mat <- as.matrix(pcs$newX)
  Y_b      <- tgt_clip(Y, bounds)

  dtrain <- xgboost::xgb.DMatrix(
    data        = X_mat,
    label       = Y_b,
    weight      = obsWeights,
    base_margin = pcs$offset_tr
  )

  params <- list(
    objective        = "binary:logistic",
    eval_metric      = "logloss",
    eta              = eta,
    max_depth        = max_depth,
    max_delta_step   = max_delta_step,
    subsample        = subsample,
    colsample_bytree = colsample_bytree,
    nthread          = nthread,
    tree_method      = "hist",
    device           = if (isTRUE(use_cuda)) "cuda" else "cpu"
  )

  fit <- tryCatch(
    xgboost::xgb.train(
      params  = params,
      data    = dtrain,
      nrounds = nrounds,
      verbose = 0
    ),
    error = function(e) {
      message(sprintf("[SL.tgt.xgboost] xgb error: %s", conditionMessage(e)))
      NULL
    }
  )

  if (is.null(fit)) {
    pred_Q <- tgt_clip(plogis(pcs$offset_pr), bounds)
    return(list(
      pred = pred_Q,
      fit  = structure(
        list(bounds = bounds),
        class = "SL.tgt.xgboost.failed"
      )
    ))
  }

  dtest <- xgboost::xgb.DMatrix(
    data        = newX_mat,
    base_margin = pcs$offset_pr
  )

  pred_Q <- predict(fit, dtest)
  pred_Q <- tgt_clip(pred_Q, bounds)

  list(
    pred = pred_Q,
    fit  = structure(
      list(
        object         = fit,
        use_cuda       = use_cuda,
        max_delta_step = max_delta_step,
        bounds         = bounds
      ),
      class = "SL.tgt.xgboost"
    )
  )
}

predict.SL.tgt.xgboost <- function(object, newdata, ...) {
  pcs <- extract_offset(newdata, newdata)
  bounds <- if (!is.null(object$bounds)) object$bounds else 1e-5

  if (is.null(object$object)) {
    return(tgt_clip(plogis(pcs$offset_pr), bounds))
  }

  newX_mat <- as.matrix(pcs$newX)

  dtest <- xgboost::xgb.DMatrix(
    data        = newX_mat,
    base_margin = pcs$offset_pr
  )

  pred_Q <- predict(object$object, dtest)
  tgt_clip(pred_Q, bounds)
}

predict.SL.tgt.glm.failed <- function(object, newdata, ...) {
  pcs <- extract_offset(newdata, newdata)
  bounds <- if (!is.null(object$bounds)) object$bounds else 1e-5
  tgt_clip(plogis(pcs$offset_pr), bounds)
}

predict.SL.tgt.glmnet.failed <- function(object, newdata, ...) {
  pcs <- extract_offset(newdata, newdata)
  bounds <- if (!is.null(object$bounds)) object$bounds else 1e-5
  tgt_clip(plogis(pcs$offset_pr), bounds)
}

predict.SL.tgt.xgboost.failed <- function(object, newdata, ...) {
  pcs <- extract_offset(newdata, newdata)
  bounds <- if (!is.null(object$bounds)) object$bounds else 1e-5
  tgt_clip(plogis(pcs$offset_pr), bounds)
}


SL.tmle_empty     <- SL.tgt.empty
SL.tmle_intercept <- SL.tgt.intercept

predict.SL.tmle_empty     <- predict.SL.tgt.empty
predict.SL.tmle_intercept <- predict.SL.tgt.intercept

SL.tmle_glm <- function(Y, X, newX, family, obsWeights, id, ...) {
  SL.tgt.glm(Y, X, newX, family, obsWeights, id, ...)
}
predict.SL.tmle_glm <- predict.SL.tgt.glm

SL.tmle_glmnet_ridge <- function(Y, X, newX, family, obsWeights, id, ...) {
  SL.tgt.glmnet(
    Y, X, newX, family, obsWeights, id,
    alpha = 0,
    nfolds = 3L,
    s_select = "lambda.1se",
    ...
  )
}
predict.SL.tmle_glmnet_ridge <- predict.SL.tgt.glmnet

SL.tmle_glmnet_enet <- function(Y, X, newX, family, obsWeights, id, ...) {
  SL.tgt.glmnet(Y, X, newX, family, obsWeights, id,
                alpha = 0.5, nfolds = 3L, s_select = "lambda.min", ...)
}
predict.SL.tmle_glmnet_enet <- predict.SL.tgt.glmnet

SL.tmle_glmnet_lasso <- function(Y, X, newX, family, obsWeights, id, ...) {
  SL.tgt.glmnet(Y, X, newX, family, obsWeights, id,
                alpha = 1, nfolds = 3L, s_select = "lambda.min", ...)
}
predict.SL.tmle_glmnet_lasso <- predict.SL.tgt.glmnet

SL.tmle_xgb_d1 <- function(Y, X, newX, family, obsWeights, id, ...) {
  SL.tgt.xgboost(
    Y, X, newX, family, obsWeights, id,
    nrounds = 50L,
    eta = 0.1,
    max_delta_step = 1,
    nthread = 1L,
    use_cuda = F,
    ...
  )
}
predict.SL.tmle_xgb_d1 <- predict.SL.tgt.xgboost

SL.tmle_xgb_d3 <- function(Y, X, newX, family, obsWeights, id, ...) {
  SL.tgt.xgboost(
    Y, X, newX, family, obsWeights, id,
    nrounds = 50L,
    eta = 0.1,
    max_delta_step = 3,
    nthread = 1L,
    use_cuda = F,
    ...
  )
}
predict.SL.tmle_xgb_d3 <- predict.SL.tgt.xgboost

SL.tmle_xgb_d6 <- function(Y, X, newX, family, obsWeights, id, ...) {
  SL.tgt.xgboost(
    Y, X, newX, family, obsWeights, id,
    nrounds = 50L,
    eta = 0.1,
    max_delta_step = 6,
    nthread = 1L,
    use_cuda = F,
    ...
  )
}
predict.SL.tmle_xgb_d6 <- predict.SL.tgt.xgboost

sl_tmle <- c(
  "SL.tmle_empty",
  "SL.tmle_intercept",
  "SL.tmle_glm",
  "SL.tmle_glmnet_ridge",
  "SL.tmle_xgb_d1",
  "SL.tmle_xgb_d3",
  "SL.tmle_xgb_d6"
)
