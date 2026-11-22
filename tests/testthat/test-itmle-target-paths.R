library(SuperLearner)

# -------------------------------------------------------------------
# Helpers to build small deterministic fixtures for itmle_target_fit()
# -------------------------------------------------------------------

# Build a fixture suitable for direct calls to itmle_target_fit.
# All rows are "at risk" at this inner step.  A subset is given weight = 0
# so we can test the pred_ok_tr branch (zero-weight training rows).
make_target_fit_fixture <- function(
    n_tr = 60L,
    n_vl = 30L,
    n_zero_weight_tr = 10L,
    seed = 42L
) {
  set.seed(seed)

  # Training design: two covariates + offset column added by target_fit.
  X_tr <- data.frame(
    x1 = stats::rnorm(n_tr),
    x2 = stats::rnorm(n_tr)
  )
  X_vl <- data.frame(
    x1 = stats::rnorm(n_vl),
    x2 = stats::rnorm(n_vl)
  )

  # Offsets within clip range so the identity fluctuation is non-degenerate.
  offset_tr <- stats::qlogis(stats::runif(n_tr, 0.2, 0.8))
  offset_vl <- stats::qlogis(stats::runif(n_vl, 0.2, 0.8))

  # Outcomes depend on x1 (so a non-trivial learner has something to fit).
  Y_tr <- stats::plogis(0.5 + 0.8 * X_tr$x1 + offset_tr)
  Y_tr <- pmin(pmax(Y_tr + stats::rnorm(n_tr, sd = 0.02), 1e-3), 1 - 1e-3)

  # Weights: positive for the first (n_tr - n_zero_weight_tr); zero for the
  # last block.  This simulates at-risk training patients with cum density
  # ratio == 0 at the current (t, s).
  weights <- rep(1.0, n_tr)
  if (n_zero_weight_tr > 0L) {
    weights[(n_tr - n_zero_weight_tr + 1L):n_tr] <- 0.0
  }

  list(
    Y_train        = Y_tr,
    X_train        = X_tr,
    X_pred_tr      = X_tr,
    X_pred_vl      = X_vl,
    offset_train   = offset_tr,
    offset_pred_tr = offset_tr,
    offset_pred_vl = offset_vl,
    weights        = weights,
    zero_tr_idx    = if (n_zero_weight_tr > 0L)
      (n_tr - n_zero_weight_tr + 1L):n_tr else integer(0)
  )
}


# -------------------------------------------------------------------
# T1 - Null-learner-only sanity: SL.tmle_empty returns plogis(offset),
# so eps should be numerically zero everywhere.
# -------------------------------------------------------------------
test_that("T1: null learner library gives eps ~= 0 everywhere", {
  fx <- make_target_fit_fixture(n_tr = 60L, n_vl = 30L, n_zero_weight_tr = 10L)

  res <- suppressMessages(suppressWarnings(
    CausalState:::itmle_target_fit(
      Y_train        = fx$Y_train,
      X_train        = fx$X_train,
      X_pred_tr      = fx$X_pred_tr,
      X_pred_vl      = fx$X_pred_vl,
      offset_train   = fx$offset_train,
      offset_pred_tr = fx$offset_pred_tr,
      offset_pred_vl = fx$offset_pred_vl,
      weights        = fx$weights,
      sl_lib         = "SL.tmle_empty",
      v_sl_inner     = 3L,
      parallel       = FALSE,
      seed           = 1L
    )
  ))

  # Tolerance 1e-10, not exact: the offset passes through
  # qlogis(clip(plogis(offset))) which can carry a float residual.
  expect_lt(max(abs(res$eps_tr)), 1e-10)
  expect_lt(max(abs(res$eps_vl)), 1e-10)

  # sl_meta has one row per candidate; empty is the only one so weight = 1.
  expect_equal(nrow(res$sl_meta), 1L)
  expect_equal(res$sl_meta$weight[1L], 1, tolerance = 1e-10)
})


# -------------------------------------------------------------------
# T2 - Fit-set row uses Z (not predict(fit)):
# For a fit-set row r (indexing fit$Z), the realised p_tr at the
# corresponding training row i = fit_ok_idx[r] should equal
# method$computePred(fit$Z[r, ], fit$coef) and differ from
# predict(fit, X[i, ])$pred.
#
# Because itmle_target_fit does not expose fit directly, we exercise the
# code by calling itmle_target_fit twice: once with the real SL fit, then
# we independently refit SL using the same seed/data and compare.
# -------------------------------------------------------------------
test_that("T2: fit-set rows use Z%*%coef, not the full-refit predict", {
  skip_if_not_installed("glmnet")

  fx <- make_target_fit_fixture(n_tr = 80L, n_vl = 20L,
                                n_zero_weight_tr = 0L, seed = 12L)

  # A library that actually fits so Z and library.predict differ.
  lib <- c("SL.tmle_glm", "SL.tmle_empty")

  res <- suppressMessages(suppressWarnings(
    CausalState:::itmle_target_fit(
      Y_train        = fx$Y_train,
      X_train        = fx$X_train,
      X_pred_tr      = fx$X_pred_tr,
      X_pred_vl      = fx$X_pred_vl,
      offset_train   = fx$offset_train,
      offset_pred_tr = fx$offset_pred_tr,
      offset_pred_vl = fx$offset_pred_vl,
      weights        = fx$weights,
      sl_lib         = lib,
      v_sl_inner     = 3L,
      parallel       = FALSE,
      seed           = 7L
    )
  ))

  # Independently refit SuperLearner with the same recipe as
  # itmle_target_fit, so we can read off fit$Z and compare.
  fit_ok <- is.finite(fx$weights) & fx$weights > 0 &
    is.finite(fx$Y_train) & is.finite(fx$offset_train) &
    stats::complete.cases(fx$X_pred_tr)
  fit_ok_idx <- which(fit_ok)
  n_fit <- length(fit_ok_idx)

  X_train_aug <- cbind(as.data.frame(fx$X_train),
                      "._sl_offset" = fx$offset_train)
  X_pred_tr_aug <- cbind(as.data.frame(fx$X_pred_tr),
                        "._sl_offset" = fx$offset_pred_tr)

  V_eff <- min(3L, n_fit)
  inner_folds <- CausalState:::with_seed(7L + 1000L, {
    origami::make_folds(n = n_fit, V = V_eff,
                       fold_fun = origami::folds_vfold)
  })
  sl_validRows <- lapply(inner_folds, function(f) f$validation_set)

  clip_fun <- function(x) pmin(pmax(x, 1e-5), 1 - 1e-5)

  fit <- CausalState:::with_seed(7L + 2000L, suppressWarnings(
    SuperLearner::SuperLearner(
      Y          = clip_fun(fx$Y_train[fit_ok_idx]),
      X          = X_train_aug[fit_ok_idx, , drop = FALSE],
      SL.library = lib,
      family     = stats::binomial(),
      obsWeights = fx$weights[fit_ok_idx],
      id         = seq_along(fit_ok_idx),
      cvControl  = list(V = V_eff, validRows = sl_validRows),
      env        = asNamespace("CausalState")
    )
  ))

  # Skip if the ensemble collapsed to the null learner (Z == plogis(offset)).
  skip_if(all(fit$coef[fit$libraryNames != "SL.tmle_empty_All"] == 0),
          "Ensemble collapsed to the null learner; Z == refit on identity.")

  r <- 1L
  i <- fit_ok_idx[r]

  manual_z <- as.numeric(
    fit$method$computePred(
      predY   = fit$Z[r, , drop = FALSE],
      coef    = fit$coef,
      control = fit$control
    )
  )
  manual_pr <- as.numeric(
    stats::predict(fit, X_pred_tr_aug[i, , drop = FALSE])$pred
  )

  # realised p_tr at the original row position
  p_tr_i <- clip_fun(stats::plogis(res$eps_tr[i] + fx$offset_pred_tr[i]))

  expect_lt(abs(p_tr_i - clip_fun(manual_z)), 1e-10)
  expect_gt(abs(p_tr_i - clip_fun(manual_pr)), 1e-6)
})


# -------------------------------------------------------------------
# T3 - Symmetric correction across tr (zero weight) and vl:
# Construct a training row and a validation row with the same
# predictors and offset.  The training row has weight 0 so it is NOT
# in the fit.  Expect the same eps for both.
# -------------------------------------------------------------------
test_that("T3: zero-weight tr and vl rows with identical H receive identical eps", {
  fx <- make_target_fit_fixture(n_tr = 60L, n_vl = 29L,
                                n_zero_weight_tr = 10L, seed = 3L)

  # Append one extra validation row that matches a zero-weight training row.
  tr_match <- fx$zero_tr_idx[1L]
  X_pred_vl_extra <- rbind(
    fx$X_pred_vl,
    fx$X_pred_tr[tr_match, , drop = FALSE]
  )
  offset_vl_extra <- c(fx$offset_pred_vl, fx$offset_pred_tr[tr_match])

  res <- suppressMessages(suppressWarnings(
    CausalState:::itmle_target_fit(
      Y_train        = fx$Y_train,
      X_train        = fx$X_train,
      X_pred_tr      = fx$X_pred_tr,
      X_pred_vl      = X_pred_vl_extra,
      offset_train   = fx$offset_train,
      offset_pred_tr = fx$offset_pred_tr,
      offset_pred_vl = offset_vl_extra,
      weights        = fx$weights,
      sl_lib         = c("SL.tmle_glm", "SL.tmle_empty"),
      v_sl_inner     = 3L,
      parallel       = FALSE,
      seed           = 5L
    )
  ))

  # eps on the matched zero-weight training row vs the matched validation row.
  eps_tr_i <- res$eps_tr[tr_match]
  eps_vl_j <- res$eps_vl[nrow(X_pred_vl_extra)]

  expect_lt(abs(eps_tr_i - eps_vl_j), 1e-10)
})


# -------------------------------------------------------------------
# T4 - Zero-weight patient outcome invariance:
# Perturbing a zero-weight patient's Y must not change fit$coef,
# fit$Z structure, or that patient's predicted correction.
# -------------------------------------------------------------------
test_that("T4: perturbing a zero-weight patient's Y leaves fit and prediction unchanged", {
  fx_a <- make_target_fit_fixture(n_tr = 60L, n_vl = 20L,
                                  n_zero_weight_tr = 10L, seed = 11L)

  # Perturb only the zero-weight patient's Y value in a second fixture.
  fx_b <- fx_a
  tr_match <- fx_a$zero_tr_idx[1L]
  fx_b$Y_train[tr_match] <- 1 - fx_a$Y_train[tr_match]

  res_a <- suppressMessages(suppressWarnings(
    CausalState:::itmle_target_fit(
      Y_train        = fx_a$Y_train,
      X_train        = fx_a$X_train,
      X_pred_tr      = fx_a$X_pred_tr,
      X_pred_vl      = fx_a$X_pred_vl,
      offset_train   = fx_a$offset_train,
      offset_pred_tr = fx_a$offset_pred_tr,
      offset_pred_vl = fx_a$offset_pred_vl,
      weights        = fx_a$weights,
      sl_lib         = c("SL.tmle_glm", "SL.tmle_empty"),
      v_sl_inner     = 3L,
      parallel       = FALSE,
      seed           = 101L
    )
  ))
  res_b <- suppressMessages(suppressWarnings(
    CausalState:::itmle_target_fit(
      Y_train        = fx_b$Y_train,
      X_train        = fx_b$X_train,
      X_pred_tr      = fx_b$X_pred_tr,
      X_pred_vl      = fx_b$X_pred_vl,
      offset_train   = fx_b$offset_train,
      offset_pred_tr = fx_b$offset_pred_tr,
      offset_pred_vl = fx_b$offset_pred_vl,
      weights        = fx_b$weights,
      sl_lib         = c("SL.tmle_glm", "SL.tmle_empty"),
      v_sl_inner     = 3L,
      parallel       = FALSE,
      seed           = 101L
    )
  ))

  # Fitted ensemble weights (via sl_meta) must be identical.
  expect_equal(res_a$sl_meta$weight, res_b$sl_meta$weight, tolerance = 1e-10)

  # The zero-weight patient's predicted correction must be identical.
  expect_equal(res_a$eps_tr[tr_match], res_b$eps_tr[tr_match],
               tolerance = 1e-10)
  # And every outer-validation row is unchanged too.
  expect_equal(res_a$eps_vl, res_b$eps_vl, tolerance = 1e-10)
})


# -------------------------------------------------------------------
# T5 - Zero-weight tr gets a non-zero eps when the ensemble is
# non-trivial (i.e., learner coef is not concentrated on empty).
# -------------------------------------------------------------------
test_that("T5: zero-weight tr rows get non-zero eps under a non-trivial ensemble", {
  fx <- make_target_fit_fixture(n_tr = 80L, n_vl = 20L,
                                n_zero_weight_tr = 10L, seed = 23L)

  res <- suppressMessages(suppressWarnings(
    CausalState:::itmle_target_fit(
      Y_train        = fx$Y_train,
      X_train        = fx$X_train,
      X_pred_tr      = fx$X_pred_tr,
      X_pred_vl      = fx$X_pred_vl,
      offset_train   = fx$offset_train,
      offset_pred_tr = fx$offset_pred_tr,
      offset_pred_vl = fx$offset_pred_vl,
      weights        = fx$weights,
      sl_lib         = c("SL.tmle_glm", "SL.tmle_empty"),
      v_sl_inner     = 3L,
      parallel       = FALSE,
      seed           = 31L
    )
  ))

  nonempty_rows <- res$sl_meta$learner != "SL.tmle_empty_All"
  skip_if(!any(nonempty_rows) || all(res$sl_meta$weight[nonempty_rows] == 0),
          "Ensemble collapsed to the null learner.")

  eps_tr_zw <- res$eps_tr[fx$zero_tr_idx]
  expect_true(max(abs(eps_tr_zw)) > 1e-8,
              info = "Zero-weight tr patients should carry a non-zero correction when the ensemble is non-trivial.")
})


# -------------------------------------------------------------------
# T6 - All-zero-coef guard.  There is no public API for injecting an
# all-zero coef vector mid-function, so we exercise the EQUIVALENT
# branch via a learner that always produces pred == offset (SL.tmle_empty).
# In that case, even though coef may be 1 on the null learner, the eps
# must still be numerically zero.
#
# For a direct all-zero coef test we would need to monkey-patch
# SuperLearner internals, which is brittle.  The null-learner case
# tests the output invariant we care about.
# -------------------------------------------------------------------
test_that("T6: null-learner outcome matches the all-zero-coef output invariant", {
  fx <- make_target_fit_fixture(n_tr = 60L, n_vl = 30L,
                                n_zero_weight_tr = 10L, seed = 77L)

  res <- suppressMessages(suppressWarnings(
    CausalState:::itmle_target_fit(
      Y_train        = fx$Y_train,
      X_train        = fx$X_train,
      X_pred_tr      = fx$X_pred_tr,
      X_pred_vl      = fx$X_pred_vl,
      offset_train   = fx$offset_train,
      offset_pred_tr = fx$offset_pred_tr,
      offset_pred_vl = fx$offset_pred_vl,
      weights        = fx$weights,
      sl_lib         = "SL.tmle_empty",
      v_sl_inner     = 3L,
      parallel       = FALSE,
      seed           = 203L
    )
  ))

  expect_true(all(abs(res$eps_tr) < 1e-10))
  expect_true(all(abs(res$eps_vl) < 1e-10))
})
