library(SuperLearner)

test_that("factor treatment coding permits unrestricted dose effects", {
  dose <- ordered(rep(0:3, each = 100L), levels = 0:3)
  probabilities <- c(0.1, 0.9, 0.2, 0.8)
  outcome <- unlist(lapply(probabilities, function(p)
    c(rep(1, round(100 * p)), rep(0, 100 - round(100 * p)))))
  D <- data.frame(A = dose)
  spec <- CausalState:::factor_treatment_levels(D, "A")
  X <- CausalState:::encode_treatment_design(D, spec)
  unordered <- factor(as.character(dose), levels = levels(dose), ordered = FALSE)
  reference <- stats::model.matrix(~ unordered, contrasts.arg = list(unordered = "contr.treatment"))[, -1L]
  expect_equal(unname(as.matrix(X)), unname(reference))
  fit <- stats::glm(outcome ~ ., data = X, family = stats::binomial())
  expect_equal(as.numeric(stats::predict(fit, type = "response")),
               probabilities[as.integer(dose)], tolerance = 1e-8)
  relabelled <- ordered(c("0", "2", "5", "9")[as.integer(dose)], levels = c("0", "2", "5", "9"))
  other <- CausalState:::encode_treatment_design(data.frame(A = relabelled),
                                                list(A = levels(relabelled)))
  expect_identical(X, other)
  previous <- options(contrasts = c("contr.sum", "contr.poly"))
  modified <- tryCatch(CausalState:::encode_treatment_design(D, spec),
                       finally = options(previous))
  expect_identical(modified, X)
})

test_that("fixed factor levels cover absent categories and every treatment history form", {
  D <- data.table::data.table(
    A = ordered(c(0, 3), levels = 0:3), A_lag1 = ordered(c(3, 0), levels = 0:3),
    A_t2 = ordered(c(1, 2), levels = 0:3), continuous = c(0.2, 0.8),
    B = factor(c("control", "treated"), levels = c("control", "treated")))
  spec <- CausalState:::factor_treatment_levels(D, c("A", "B", "continuous"))
  data.table::setattr(D, "cs_treatment_levels", spec)
  X <- CausalState:::make_design(D, 1:2, names(D))
  expect_setequal(names(X), c(paste0("A__level_", 2:4),
                             paste0("A_lag1__level_", 2:4),
                             paste0("A_t2__level_", 2:4), "continuous", "B__level_2"))
  expect_equal(X$A__level_2, c(0, 0))
  expect_equal(X$A__level_4, c(0, 1))
  expect_identical(X$continuous, D$continuous)
  expect_equal(X$B__level_2, c(0, 1))
  subset <- CausalState:::make_design(D, 1L, names(D))
  expect_identical(names(subset), names(X))
  expect_equal(subset, X[1L, , drop = FALSE])
})

test_that("policy evaluation preserves factor semantics and updates the correct indicators", {
  D <- data.table::data.table(id = 1:4, time = 2L,
                              A = ordered(0:3, levels = 0:3),
                              A_lag1 = ordered(3:0, levels = 0:3))
  spec <- CausalState:::factor_treatment_levels(D, "A")
  data.table::setattr(D, "cs_treatment_levels", spec)
  original <- CausalState:::make_design(D, 1:4, c("A", "A_lag1"))
  policy <- function(D, t, a_names) {
    stopifnot(is.ordered(D$A), identical(levels(D$A), as.character(0:3)))
    # The returned factor deliberately has only the observed destination levels.
    data.frame(A = factor(levels(D$A)[pmax(as.integer(D$A), 3L)]))
  }
  shifted <- CausalState:::overwrite_policy_history_for_Q(
    D, 2L, "A", policy, "id", "time")
  expect_true(is.ordered(shifted$A))
  expect_identical(levels(shifted$A), as.character(0:3))
  expect_identical(as.character(shifted$A), c("2", "2", "2", "3"))
  patched <- CausalState:::patch_shifted_design(original, shifted, 1:4, "A", 2L, 1L, 1L)
  expect_equal(patched$A__level_2, c(0, 0, 0, 0))
  expect_equal(patched$A__level_3, c(1, 1, 1, 0))
  expect_equal(patched$A__level_4, c(0, 0, 0, 1))
  history <- paste0("A_lag1__level_", 2:4)
  expect_identical(patched[history], original[history])
  expect_identical(as.character(D$A), as.character(0:3))
  expect_error(CausalState:::overwrite_policy_history_for_Q(
    D, 2L, "A", function(D,t,a_names) list(A = "invalid"), "id", "time"),
    "outside factor levels")
  masked <- CausalState:::g_shift_values(data.frame(A = D$A),
    list(A = list(idx = c(TRUE,FALSE,TRUE,FALSE), value = "2")), 2L)
  expect_identical(as.character(masked$A), c("2","1","2","3"))
  expect_true(is.ordered(masked$A))
})

test_that("factor schemas are inferred and checked against fit_g metadata", {
  D <- data.frame(A = ordered(0:3, levels = 0:3))
  expected <- list(A = levels(D$A))
  weights <- list(settings = list(treatment_levels = expected, treatment_ordered = c(A = TRUE)))
  expect_identical(CausalState:::resolve_treatment_levels(D, "A", weights), expected)
  expect_identical(CausalState:::resolve_treatment_levels(D, "A"), expected)
  reordered <- data.frame(A = ordered(D$A, levels = rev(levels(D$A))))
  expect_error(CausalState:::resolve_treatment_levels(reordered, "A", weights), "factor levels differ")
  nominal <- data.frame(A = factor(D$A, ordered = FALSE))
  expect_error(CausalState:::resolve_treatment_levels(nominal, "A", weights), "ordering differs")
  expect_error(CausalState:::resolve_treatment_levels(data.frame(A = 0:3), "A", weights),
               "factor levels differ")
  expect_error(CausalState:::encode_treatment_design(data.frame(A = D$A, A__level_2 = 1),
                                                     expected), "conflict")
  binary <- data.frame(A = c(0, 1), A_lag1 = c(1, 0))
  expect_identical(CausalState:::factor_treatment_levels(binary, "A"), list())
  expect_identical(CausalState:::encode_treatment_design(binary, list()), binary)
})

test_that("public estimators use factor declarations without a treatment_levels argument", {
  for (fun in list(fit_g, sdr, itmle, qreg))
    expect_false("treatment_levels" %in% names(formals(fun)))
})

test_that("all estimator paths and targeting receive numeric treatment indicators", {
  seen <- new.env(parent = emptyenv())
  seen$models <- list()
  seen$target <- list()
  block_fit <- CausalState:::sl_block_fit
  target_fit <- CausalState:::itmle_target_fit
  testthat::local_mocked_bindings(
    sl_block_fit = function(Y, X, ...) {
      seen$models[[length(seen$models) + 1L]] <- names(X)
      block_fit(Y = Y, X = X, ...)
    },
    itmle_target_fit = function(...) {
      args <- list(...)
      seen$target[[length(seen$target) + 1L]] <- names(args$X_train)
      do.call(target_fit, args)
    },
    .package = "CausalState"
  )
  D <- sim_cat(n = 1000L, tmax = 3L, seed = 42L)
  policy <- function(D, t, a_names) {
    stopifnot(is.ordered(D$A))
    rank <- as.integer(D$A)
    rank[D$L1 > 0.5] <- pmax(rank[D$L1 > 0.5], 3L)
    list(A = levels(D$A)[rank])
  }
  g_args <- list(df = D, a_names = "A", tmax = 3L, baseline = "age", tv_names = "L1",
                 sl_g = "SL.mean", k = 1L, inner_v = 2L, v = 2L, time = "time",
                 policy_spec_fun = policy, verbose = FALSE)
  w <- suppressMessages(do.call(fit_g, c(g_args, list(method = "propensity"))))
  expect_true(any(vapply(seen$models, function(nm)
    all(paste0("A_lag1__level_", 2:4) %in% nm), logical(1L))))
  classified <- suppressMessages(do.call(fit_g, c(g_args, list(method = "classification"))))
  expect_identical(classified$settings$treatment_encoding, "reference")
  expect_true(any(vapply(seen$models, function(nm)
    all(paste0("A__level_", 2:4) %in% nm), logical(1L))))

  args <- list(df = D, weight_object = w, tmax = 3L, id = "id", time = "time",
               alive = "alive", in_state = "in_state", y = "Y", baseline = "age",
               tv_names = "L1", a_names = "A", sl_remain = "SL.mean", sl_death = "SL.mean",
               sl_recursive = "SL.mean", sl_y = "SL.mean", k = 1L, inner_v = 2L,
               policy_spec_fun = policy, verbose = FALSE)
  for (pooled in c(FALSE, TRUE)) {
    for (estimator in list(sdr, qreg, itmle)) {
      est_args <- c(args, list(pool_g_death = pooled, pool_q_exit = pooled))
      if (identical(estimator, qreg)) {
        est_args$weight_object <- NULL
        est_args$v <- 2L
      }
      if (identical(estimator, itmle)) {
        est_args$sl_tmle <- "SL.tgt.intercept"
        est_args$v_target_itmle <- 2L
        est_args$v_sl_inner_itmle <- 2L
      }
      result <- suppressMessages(do.call(estimator, est_args))
      expect_true(is.finite(result$psi))
      expect_true(is.finite(result$se))
      expect_identical(result$settings$treatment_levels, list(A = levels(D$A)))
      expect_identical(result$settings$treatment_encoding, "reference")
    }
  }
  raw <- c("A", "A_lag1", "A_t1")
  expect_false(any(unlist(seen$models) %in% raw))
  expect_false(any(unlist(seen$target) %in% raw))
  expect_true(any(vapply(seen$target, function(nm)
    all(paste0("A_t1__level_", 2:4) %in% nm), logical(1L))))
}) 

test_that("joint factor and continuous classification keeps both treatment types", {
  D <- data.frame(id = 1:400, time = 1L,
    A = factor(rep(c("none","low","medium","high"),100),
               levels = c("none","low","medium","high")),
    B = seq(0,1,length.out=400))
  original <- D
  policy <- function(D,t,a_names) {
    stopifnot(is.factor(D$A), is.numeric(D$B))
    list(A = "medium", B = D$B + 0.25)
  }
  w <- suppressMessages(fit_g(
    D, c("A","B"), 1L, sl_g = "SL.mean", k = 0L, inner_v = 2L, v = 2L,
    time = "time", policy_spec_fun = policy, verbose = FALSE
  ))
  expect_identical(w$settings$treatment_levels, list(A = levels(D$A)))
  expect_identical(w$settings$treatment_ordered, c(A = FALSE))
  expect_true(all(is.finite(w$weights_dt$Rt_t)))
  expect_equal(w$policy_change_summary$mean_A_delta_B, 0.25)
  expect_false("mean_A_delta_A" %in% names(w$policy_change_summary))
  expect_identical(D, original)
})
