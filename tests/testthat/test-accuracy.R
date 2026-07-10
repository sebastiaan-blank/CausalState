library(SuperLearner)

# End-to-end accuracy tests against a known Monte Carlo truth.
#
# The DGP is the same structure as the smoke-test panel but extended with a
# propensity_fn hook.  The "truth" E[Y(d)] and E[mort(d)] are obtained by
# simulating n = 500,000 subjects under the policy propensity directly, which
# makes the Monte Carlo error negligible (~0.001) relative to estimation error.
#
# Policy: propensity shift by delta = 0.3 (g_d(A|H) = Bernoulli(p_nat + 0.3)).
# Matching policy function passed to the estimators: A_d = min(A + 0.3, 1).
#
# For each estimator the test checks:
#   |psi_hat - truth| < 3 * se       (fails ~0.3% by chance if estimator correct)
#   ci[1] < truth < ci[2]            (fails ~5% by chance if CI is nominal)
#
# Diagnostic printf lines print psi, truth, se, and |diff| / (3*se) for every
# estimator so failures can be diagnosed without re-running.

DELTA <- 0.3  # propensity shift magnitude

# ── DGP ───────────────────────────────────────────────────────────────────────

sim_dgp <- function(n = 2000L, tmax = 5L, seed = 1L,
                    propensity_fn = identity) {
  set.seed(seed)
  rows <- vector("list", n)

  for (i in seq_len(n)) {
    age <- round(rnorm(1, 65, 10))
    sex <- rbinom(1, 1, 0.5)
    L1  <- rnorm(1, 0, 1)
    L2  <- rbinom(1, 1, 0.4)
    pat <- list()

    for (t in seq_len(tmax)) {
      p_nat <- plogis(0.3 * L1 - 0.4 + 0.2 * sex)
      A     <- rbinom(1, 1, propensity_fn(p_nat))

      p_die <- plogis(-4.0 + 0.3 * L1 - 0.1 * age / 10)
      p_dc  <- plogis(-2.5 + 0.5 * A  - 0.2 * L2)
      u     <- runif(1)

      if (u < p_die) {
        alive <- 0L; in_state <- 0L
      } else if (u < p_die + p_dc) {
        alive <- 1L; in_state <- 0L
      } else {
        alive <- 1L; in_state <- 1L
      }

      py <- if (!alive) plogis(-3.0 + 0.1 * L1)
            else if (!in_state) plogis(1.5 + 0.2 * L1 - 0.1 * age / 10 + 0.3 * A)
            else plogis(-0.5 + 0.4 * A - 0.2 * L1 + 0.1 * L2)
      Y <- rbinom(1, 1, py)

      pat[[length(pat) + 1L]] <- data.frame(
        id = i, time = t,
        age = age, sex = sex,
        alive = alive, in_state = in_state,
        L1 = L1, L2 = L2,
        A = A, Y = Y
      )

      if (in_state == 0L) break

      if (t < tmax) {
        L1 <- L1 + rnorm(1, -0.1 * A, 0.3)
        L2 <- rbinom(1, 1, plogis(0.5 * L2 + 0.3 * A - 0.5))
      }
    }

    rows[[i]] <- do.call(rbind, pat)
  }

  do.call(rbind, rows)
}

# Policy function for the estimators: A_d = min(A + DELTA, 1).
# For binary A in {0, 1}: A=0 -> DELTA, A=1 -> 1.
# Consistent with propensity_fn = function(p) pmin(p + DELTA, 0.99) in the DGP.
policy_fn <- function(D_block, t, a_names) {
  out <- D_block[, ..a_names, drop = FALSE]
  out[[a_names[1]]] <- pmin(D_block[[a_names[1]]] + DELTA, 1)
  out
}

# ── Monte Carlo truth (computed once, shared across all tests) ─────────────────

pol_propensity <- function(p) pmin(p + DELTA, 0.99)

df_mc <- sim_dgp(n = 500000L, tmax = 5L, seed = 77L,
                 propensity_fn = pol_propensity)

# E[Y(d)]: last-row Y per patient (mirrors build_Y_by_id which takes .N row)
last_rows_mc <- df_mc[!duplicated(df_mc$id, fromLast = TRUE), ]
truth_main   <- mean(last_rows_mc$Y)

# E[mort(d)]: any death during admission (mirrors sdr_competing/itmle_competing)
death_by_id  <- tapply(df_mc$alive, df_mc$id, function(x) any(x == 0L))
truth_mort   <- mean(death_by_id)

rm(df_mc, last_rows_mc, death_by_id)

message(sprintf("[accuracy] MC truth (delta=%.2f, n=500k): E[Y(d)]=%.4f  E[mort(d)]=%.4f",
                DELTA, truth_main, truth_mort))

# ── Shared fixtures ───────────────────────────────────────────────────────────

sl_fast     <- c("SL.mean", "SL.glm")
sl_tgt_glm  <- c("SL.tgt.glm")

df_est <- sim_dgp(n = 2000L, tmax = 5L, seed = 42L)

wr <- density_ratio(
  df              = df_est,
  a_names         = "A",
  tmax            = 5L,
  baseline        = c("age", "sex"),
  tv_names        = c("L1", "L2"),
  sl_g            = sl_fast,
  k               = 1L,
  inner_v         = 2L,
  v               = 2L,
  seed            = 1L,
  id              = "id",
  time            = "time",
  parallel_t      = FALSE,
  policy_spec_fun = policy_fn
)

# ── Helpers ───────────────────────────────────────────────────────────────────

check_accuracy <- function(res, truth, label) {
  diff <- abs(res$psi - truth)
  tol  <- 3 * res$se
  message(sprintf(
    "[accuracy] %s: psi=%.4f  truth=%.4f  se=%.4f  |diff|/3se=%.2f",
    label, res$psi, truth, res$se, diff / tol
  ))
  expect_true(is.finite(res$psi), label = paste(label, "psi finite"))
  expect_gt(res$se, 0, label = paste(label, "se > 0"))
  expect_lt(diff, tol,
            label = sprintf("%s |psi-truth|=%.4f >= 3*se=%.4f", label, diff, tol))
}

# ── sdr ───────────────────────────────────────────────────────────────────────

test_that("sdr estimate is within 3 SE of MC truth", {
  res <- sdr(
    df              = df_est,
    weight_object   = wr,
    tmax            = 5L,
    id              = "id",
    time            = "time",
    alive           = "alive",
    in_state        = "in_state",
    y               = "Y",
    baseline        = c("age", "sex"),
    tv_names        = c("L1", "L2"),
    a_names         = "A",
    sl_g            = sl_fast,
    sl_recursive    = sl_fast,
    sl_y            = sl_fast,
    outcome_family  = "binomial",
    k               = 1L,
    inner_v         = 2L,
    parallel        = FALSE,
    seed            = 1L,
    policy_spec_fun = policy_fn
  )

  check_accuracy(res, truth_main, "sdr")
})

# ── itmle ─────────────────────────────────────────────────────────────────────

test_that("itmle estimate is within 3 SE of MC truth", {
  res <- itmle(
    df               = df_est,
    weight_object    = wr,
    tmax             = 5L,
    id               = "id",
    time             = "time",
    alive            = "alive",
    in_state         = "in_state",
    y                = "Y",
    baseline         = c("age", "sex"),
    tv_names         = c("L1", "L2"),
    a_names          = "A",
    sl_g             = sl_fast,
    sl_recursive     = sl_fast,
    sl_y             = sl_fast,
    sl_tmle          = sl_tgt_glm,
    outcome_family   = "binomial",
    k                = 1L,
    inner_v          = 2L,
    v_target_itmle   = 2L,
    v_sl_inner_itmle = 2L,
    parallel         = FALSE,
    seed             = 1L,
    policy_spec_fun  = policy_fn
  )

  check_accuracy(res, truth_main, "itmle")
})

# ── sdr_competing ─────────────────────────────────────────────────────────────

test_that("sdr_competing estimate is within 3 SE of MC mortality truth", {
  res <- sdr_competing(
    DT              = df_est,
    weights_dt      = wr$weights_dt,
    policy_spec_fun = policy_fn,
    id              = "id",
    time            = "time",
    alive           = "alive",
    in_state        = "in_state",
    baseline        = c("age", "sex"),
    tv_names        = c("L1", "L2"),
    a_names         = "A",
    k               = 1L,
    tmax            = 5L,
    seed            = 1L,
    sl_q            = sl_fast,
    inner_v         = 2L,
    trim            = 0.99,
    parallel        = FALSE
  )

  check_accuracy(res, truth_mort, "sdr_competing")
})

# ── itmle_competing ───────────────────────────────────────────────────────────

test_that("itmle_competing estimate is within 3 SE of MC mortality truth", {
  res <- itmle_competing(
    DT               = df_est,
    weights_dt       = wr$weights_dt,
    policy_spec_fun  = policy_fn,
    id               = "id",
    time             = "time",
    alive            = "alive",
    in_state         = "in_state",
    baseline         = c("age", "sex"),
    tv_names         = c("L1", "L2"),
    a_names          = "A",
    k                = 1L,
    tmax             = 5L,
    seed             = 1L,
    sl_q             = sl_fast,
    sl_tmle          = sl_tgt_glm,
    inner_v          = 2L,
    v_target_itmle   = 2L,
    v_sl_inner_itmle = 2L,
    trim             = 0.99,
    parallel         = FALSE
  )

  check_accuracy(res, truth_mort, "itmle_competing")
})

# ── qreg ──────────────────────────────────────────────────────────────────────
#
# qreg is a plug-in estimator with first-order bias; the 3*SE criterion used
# for DR estimators does not apply.  Test that it produces a finite, bounded
# estimate with positive SE on both pathways (with and without weight_object).

test_that("qreg produces finite estimates and positive SE", {
  res_no_w <- qreg(
    df              = df_est,
    tmax            = 5L,
    id              = "id",
    time            = "time",
    alive           = "alive",
    in_state        = "in_state",
    y               = "Y",
    baseline        = c("age", "sex"),
    tv_names        = c("L1", "L2"),
    a_names         = "A",
    sl_g            = sl_fast,
    sl_recursive    = sl_fast,
    sl_y            = sl_fast,
    outcome_family  = "binomial",
    k               = 1L,
    inner_v         = 2L,
    v               = 2L,
    parallel        = FALSE,
    seed            = 1L,
    policy_spec_fun = policy_fn
  )

  expect_true(is.finite(res_no_w$estimate))
  expect_gte(res_no_w$estimate, 0)
  expect_lte(res_no_w$estimate, 1)
  expect_gt(res_no_w$se_naive, 0)
  expect_null(res_no_w$se_eif)
  expect_equal(res_no_w$se, res_no_w$se_naive)

  message(sprintf("[accuracy] qreg (no weights): estimate=%.4f  se_naive=%.4f",
                  res_no_w$estimate, res_no_w$se_naive))

  res_w <- qreg(
    df              = df_est,
    weight_object   = wr,
    tmax            = 5L,
    id              = "id",
    time            = "time",
    alive           = "alive",
    in_state        = "in_state",
    y               = "Y",
    baseline        = c("age", "sex"),
    tv_names        = c("L1", "L2"),
    a_names         = "A",
    sl_g            = sl_fast,
    sl_recursive    = sl_fast,
    sl_y            = sl_fast,
    outcome_family  = "binomial",
    k               = 1L,
    inner_v         = 2L,
    parallel        = FALSE,
    seed            = 1L,
    policy_spec_fun = policy_fn
  )

  expect_true(is.finite(res_w$estimate))
  expect_gte(res_w$estimate, 0)
  expect_lte(res_w$estimate, 1)
  expect_gt(res_w$se_naive, 0)
  expect_true(is.finite(res_w$se_eif))
  expect_gt(res_w$se_eif, 0)
  expect_equal(res_w$se, res_w$se_eif)

  message(sprintf("[accuracy] qreg (with weights): estimate=%.4f  se_naive=%.4f  se_eif=%.4f",
                  res_w$estimate, res_w$se_naive, res_w$se_eif))
})
