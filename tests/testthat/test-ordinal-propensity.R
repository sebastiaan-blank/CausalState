library(SuperLearner)

ordinal_data <- function(n = 800L, labels = c("none", "low", "medium", "high")) {
  data.frame(id = seq_len(n), time = 1L,
             A = ordered(rep(labels, length.out = n), levels = labels),
             W = rep(c(-1, 1), each = length(labels), length.out = n))
}

fit_ordinal_test <- function(d, policy, ...) {
  suppressMessages(fit_g(
    df = d, a_names = "A", tmax = max(d$time), baseline = "W",
    sl_g = "SL.mean", k = 0L, inner_v = 2L, v = 2L, seed = 19L,
    id = "id", time = "time", policy_spec_fun = policy,
    method = "propensity", verbose = FALSE, ...
  ))
}

test_that("ordinal factor ratios agree with empirical training probabilities", {
  d <- ordinal_data()
  policy <- function(D, t, a_names) {
    value <- as.character(D$A)
    value[D$W > 0 & as.integer(D$A) < 3L] <- "medium"
    data.frame(A = factor(value, levels = levels(D$A), ordered = TRUE))
  }
  w <- fit_ordinal_test(d, policy)
  expect_s3_class(w, "g_fit")
  expect_identical(w$settings$propensity_model, "continuation_ratio")
  expect_identical(w$settings$treatment_levels, list(A = levels(d$A)))
  expect_identical(w$settings$treatment_ordered, c(A = TRUE))
  expect_false(w$settings$is_binary)
  expect_setequal(w$sl_summary$threshold, levels(d$A)[1:3])
  expect_false(any(grepl("^mean_A_", names(w$policy_change_summary))))
  a <- as.integer(d$A)
  for (fold in unique(w$weights_dt$global_fold)) {
    validation <- w$weights_dt$global_fold == fold
    p <- tabulate(a[!validation], nbins = 4L) / sum(!validation)
    expected <- ifelse(d$W <= 0, 1,
                       ifelse(a < 3, 0, ifelse(a == 3, sum(p[1:3]) / p[3], 1)))
    expect_equal(w$weights_dt$Rt_t[validation], expected[validation], tolerance = 1e-12)
  }
  expect_true(is.ordered(d$A))
})

test_that("saturated ordinal shifts sum both upper-boundary preimages", {
  d <- ordinal_data()
  w <- fit_ordinal_test(d, function(D, t, a_names)
    list(A = levels(D$A)[pmin(as.integer(D$A) + 1L, 4L)]))
  a <- as.integer(d$A)
  for (fold in unique(w$weights_dt$global_fold)) {
    validation <- w$weights_dt$global_fold == fold
    p <- tabulate(a[!validation], nbins = 4L) / sum(!validation)
    expected <- (c(0, p[1], p[2], sum(p[3:4])) / p)[a[validation]]
    expect_equal(w$weights_dt$Rt_t[validation], expected, tolerance = 1e-12)
  }
})

test_that("numeric factor labels are never confused with internal factor codes", {
  d <- ordinal_data(labels = c("0", "2", "5", "9"))
  w <- fit_ordinal_test(d, function(D, t, a_names)
    list(A = levels(D$A)[pmax(as.integer(D$A), 3L)]))
  a <- as.integer(d$A)
  for (fold in unique(w$weights_dt$global_fold)) {
    validation <- w$weights_dt$global_fold == fold
    p <- tabulate(a[!validation], nbins = 4L) / sum(!validation)
    expected <- c(0, 0, sum(p[1:3]) / p[3], 1)[a[validation]]
    expect_equal(w$weights_dt$Rt_t[validation], expected, tolerance = 1e-12)
  }
  zero <- fit_ordinal_test(d, function(D, t, a_names) list(A = 0))
  expect_true(all(zero$weights_dt$Rt_t[as.character(d$A) != "0"] == 0))
  expect_true(all(zero$weights_dt$Rt_t[as.character(d$A) == "0"] > 0))
})

test_that("identity ordinal policies produce exactly one without fitting", {
  d <- ordinal_data(n = 8L)
  for (policy in list(function(D, t, a_names) NULL,
                      function(D, t, a_names) data.frame(A = D$A))) {
    w <- fit_ordinal_test(d, policy)
    expect_identical(w$weights_dt$Rt_t, rep(1, nrow(d)))
    expect_null(w$sl_summary)
  }
})

test_that("categorical declarations and policy labels are validated", {
  d <- ordinal_data()
  policy <- function(D, t, a_names) list(A = "medium")
  expect_error(fit_ordinal_test(transform(d, A = as.integer(A)), policy), "ordered factor")
  expect_error(fit_ordinal_test(transform(d, A = factor(A, ordered = FALSE)), policy),
               "ordered factor")
  expect_error(fit_ordinal_test(transform(d, A = as.character(A)), policy), "numeric or a factor")
  one <- d
  one$A <- ordered(rep("none", nrow(d)), levels = "none")
  expect_error(fit_ordinal_test(one, policy), "at least two")
  expect_error(fit_ordinal_test(d, function(D, t, a_names) list(A = "invalid")),
               "outside factor levels")
  expect_error(fit_ordinal_test(d, function(D, t, a_names) list(A = 1)),
               "outside factor levels")
  expect_error(fit_ordinal_test(d, function(D, t, a_names) list(A = NA_character_)),
               "missing")
  expect_error(fit_ordinal_test(d, function(D, t, a_names) list(A = c("none", "low"))),
               "one value or 800 values")
  expect_error(fit_g(transform(d, B = A), c("A", "B"), 1L,
                    sl_g = "SL.mean", time = "time", policy_spec_fun = policy,
                    method = "propensity", verbose = FALSE), "single")
  nominal <- d
  nominal$A <- factor(nominal$A, levels = levels(d$A), ordered = FALSE)
  classified <- suppressMessages(fit_g(
    nominal, "A", 1L, sl_g = "SL.mean", time = "time", k = 0L, inner_v = 2L,
    policy_spec_fun = policy, verbose = FALSE))
  expect_identical(classified$settings$treatment_encoding, "reference")
  expect_identical(classified$settings$treatment_ordered, c(A = FALSE))
})

test_that("hypothetical doses are checked even if observed doses do not change", {
  d <- ordinal_data()
  lev <- levels(d$A)
  d$A <- ordered(rep(lev[1:3], length.out = nrow(d)), levels = lev)
  expect_error(fit_ordinal_test(d, function(D, t, a_names)
    list(A = ifelse(D$A == "high", "invalid", as.character(D$A)))), "outside factor levels")
  expect_error(fit_ordinal_test(d, function(D, t, a_names)
    list(A = ifelse(D$A == "high", "medium", as.character(D$A)))),
    "no training observations.*high")
})

test_that("ordinal factor fitting preserves cluster folds over time", {
  d <- ordinal_data()
  d$site <- rep(seq_len(80L), each = 10L)
  d <- rbind(d, transform(d, time = 2L))
  d <- d[order(d$id, d$time), ]
  w <- fit_ordinal_test(d, function(D, t, a_names)
    list(A = levels(D$A)[pmax(as.integer(D$A), 3L)]), cluster = "site")
  joined <- merge(d, as.data.frame(w$weights_dt), by = c("id", "time"))
  expect_true(all(tapply(joined$global_fold, joined$site,
                         function(x) length(unique(x))) == 1L))
  expect_true(all(is.finite(joined$Rt_t)))
  expect_true(all(joined$Rt_t[as.integer(joined$A) < 3L] == 0))
})

test_that("numeric binary stochastic propensity retains its existing semantics", {
  d <- ordinal_data()
  d$A <- rep(0:1, length.out = nrow(d))
  policy <- function(D, t, a_names) data.frame(A = 0.2 + 0.6 * D$A)
  w <- fit_ordinal_test(d, policy)
  expect_identical(w$settings$propensity_model, "binary")
  expect_identical(w$settings$treatment_levels, list())
  baseline <- transform(d, time = 0L, A = NA_real_)
  with_baseline <- fit_ordinal_test(rbind(baseline, d), policy)
  expect_identical(with_baseline$settings$propensity_model, "binary")
  expect_true(all(is.finite(with_baseline$weights_dt$Rt_t)))
  for (fold in unique(w$weights_dt$global_fold)) {
    validation <- w$weights_dt$global_fold == fold
    p <- mean(d$A[!validation])
    q <- 0.8 * p + 0.2 * (1 - p)
    a <- d$A[validation]
    expected <- ifelse(a == 1, q / p, (1 - q) / (1 - p))
    expect_equal(w$weights_dt$Rt_t[validation], expected, tolerance = 1e-12)
  }
})

test_that("ordinal factor weights work with every downstream estimator", {
  d <- sim_cat(n = 1000L, tmax = 2L, seed = 42L)
  expect_true(is.ordered(d$A))
  policy <- function(D, t, a_names) {
    rank <- as.integer(D$A)
    rank[D$L1 > 0.5] <- pmax(rank[D$L1 > 0.5], 3L)
    list(A = levels(D$A)[rank])
  }
  w <- suppressMessages(fit_g(
    d, "A", 2L, baseline = "age", tv_names = "L1",
    sl_g = c("SL.mean", "SL.glm"), k = 1L, inner_v = 2L, v = 2L,
    time = "time", policy_spec_fun = policy, method = "propensity", verbose = FALSE
  ))
  args <- list(
    df = d, weight_object = w, tmax = 2L, id = "id", time = "time",
    alive = "alive", in_state = "in_state", y = "Y", baseline = "age",
    tv_names = "L1", a_names = "A", sl_remain = "SL.mean", sl_death = "SL.mean",
    sl_recursive = "SL.mean", sl_y = "SL.mean", k = 1L, inner_v = 2L,
    policy_spec_fun = policy
  )
  for (estimator in list(sdr, qreg)) {
    estimator_args <- args
    if (identical(estimator, qreg)) estimator_args$weight_object <- NULL
    result <- suppressMessages(do.call(estimator, estimator_args))
    expect_true(is.finite(result$psi))
    expect_true(is.finite(result$se))
    expect_identical(result$settings$treatment_levels, list(A = levels(d$A)))
  }
  result <- suppressMessages(do.call(itmle, c(args, list(
    sl_tmle = "SL.tgt.glm", v_target_itmle = 2L, v_sl_inner_itmle = 2L
  ))))
  expect_true(is.finite(result$psi))
  expect_true(is.finite(result$se))
  expect_identical(result$settings$treatment_levels, list(A = levels(d$A)))
})

test_that("zero policy ratios survive weight assembly and cumulative multiplication", {
  weights <- data.frame(id = c(1, 1, 2, 2), time = c(1, 2, 1, 2),
                        Rt_t = c(0, 2, 1, 0))
  expect_no_warning(dens <- CausalState:::build_density_ratios(
    weights, "id", "time", ids = c("2", "1"), tmax = 2L,
    row_index = matrix(1L, 2L, 2L), trim = 1
  ))
  expect_equal(unname(dens), matrix(c(1, 0, 0, 2), nrow = 2L))
  expect_equal(apply(dens, 1L, prod), c("2" = 0, "1" = 0))
})
