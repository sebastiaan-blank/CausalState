library(SuperLearner)

# End-to-end smoke test using a small simulated ICU panel.
#
# Data structure:
#   - One row per patient per time point WHILE IN ICU.
#   - Trajectory ends when the patient exits: the last row carries the exit
#     state (alive=1/in_state=0 for discharge, alive=0/in_state=0 for death).
#   - Patients still in ICU at tmax have alive=1, in_state=1 on their last row.
#   - A is treatment given in that ICU period (0 only makes clinical sense
#     when in_state=1; the last row for discharged/dead patients still records A
#     because treatment was given before the exit transition was observed).
#   - Y is the outcome of interest (e.g. 28-day survival); it is defined on
#     every row and is modelled by Q_exit at the time of ICU departure and by
#     Q_rem for patients still present at tmax.

sim_panel <- function(n = 150L, tmax = 3L, seed = 42L) {
  set.seed(seed)
  rows <- vector("list", n)

  for (i in seq_len(n)) {
    age <- round(rnorm(1, 65, 10))
    sex <- rbinom(1, 1, 0.5)
    L1  <- rnorm(1, 0, 1)
    L2  <- rbinom(1, 1, 0.4)
    pat <- list()

    for (t in seq_len(tmax)) {
      # Treatment given at start of this ICU period
      A <- rbinom(1, 1, plogis(0.3 * L1 - 0.4 + 0.2 * sex))

      # Exit determination: die / discharge / remain
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

      # Y: outcome determined at exit (Q_exit) or at tmax (Q_rem).
      # Dead → poor outcome; discharged alive → generally good; still in ICU
      # at tmax → intermediate.
      py <- if (!alive)  plogis(-3.0 + 0.1 * L1)
            else if (!in_state) plogis(1.5  + 0.2 * L1 - 0.1 * age / 10 + 0.3 * A)
            else           plogis(-0.5 + 0.4 * A   - 0.2 * L1 + 0.1 * L2)
      Y <- rbinom(1, 1, py)

      pat[[length(pat) + 1L]] <- data.frame(
        id = i, time = t,
        age = age, sex = sex,
        alive = alive, in_state = in_state,
        L1 = L1, L2 = L2,
        A = A, Y = Y
      )

      # Stop adding rows once patient has exited ICU
      if (in_state == 0L) break

      # Update time-varying covariates for next period (only while in ICU)
      if (t < tmax) {
        L1 <- L1 + rnorm(1, -0.1 * A, 0.3)
        L2 <- rbinom(1, 1, plogis(0.5 * L2 + 0.3 * A - 0.5))
      }
    }

    rows[[i]] <- do.call(rbind, pat)
  }

  do.call(rbind, rows)
}

# ── Shared fixtures ───────────────────────────────────────────────────────

df_smoke <- sim_panel(n = 1000L, tmax = 3L, seed = 42L)

# Verify structure before running estimators
stopifnot(all(c("id","time","alive","in_state","age","sex","L1","L2","A","Y") %in% names(df_smoke)))
last_rows <- df_smoke[!duplicated(df_smoke$id, fromLast = TRUE), ]
stopifnot(all(last_rows$in_state == 0L | last_rows$time == 3L))

# Soft upward shift: A → min(A + 0.5, 1)
policy_up <- function(D_block, t, a_names) {
  out <- D_block[, ..a_names, drop = FALSE]
  out[[a_names[1]]] <- pmin(D_block[[a_names[1]]] + 0.5, 1)
  out
}

sl_fast        <- c("SL.mean", "SL.glm")
sl_tmle_simple <- c("SL.tgt.empty", "SL.tgt.intercept")

wr_smoke <- density_ratio(
  df              = df_smoke,
  a_names         = "A",
  tmax            = 3L,
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
  policy_spec_fun = policy_up
)

# ── density_ratio ─────────────────────────────────────────────────────────

test_that("density_ratio returns valid weights", {
  wdt <- wr_smoke$weights_dt
  expect_s3_class(wdt, "data.frame")
  expect_gt(nrow(wdt), 0L)
  expect_true("Rt_t" %in% names(wdt))
  expect_true(all(is.finite(wdt$Rt_t)))
  expect_true(all(wdt$Rt_t >= 0))
})

# ── sdr ───────────────────────────────────────────────────────────────────

test_that("sdr returns a finite psi with positive se", {
  res <- sdr(
    df              = df_smoke,
    weight_object   = wr_smoke,
    tmax            = 3L,
    id              = "id",
    time            = "time",
    alive           = "alive",
    in_state        = "in_state",
    y               = "Y",
    baseline        = c("age", "sex"),
    tv_names        = c("L1", "L2"),
    a_names         = "A",
    sl_remain       = sl_fast,
    sl_death        = sl_fast,
    sl_recursive    = sl_fast,
    sl_y            = sl_fast,
    outcome_family  = "binomial",
    k               = 1L,
    inner_v         = 2L,
    parallel        = FALSE,
    seed            = 1L,
    policy_spec_fun = policy_up
  )

  expect_true(is.list(res))
  expect_true(is.finite(res$psi))
  expect_gte(res$psi, 0)
  expect_lte(res$psi, 1)
  expect_true(is.finite(res$se))
  expect_gt(res$se, 0)
})

test_that("sdr with pooled g_death + Q_exit runs under each pool_time mode", {
  for (mode in c("spline", "linear", "factor")) {
    res <- sdr(
      df              = df_smoke,
      weight_object   = wr_smoke,
      tmax            = 3L,
      id              = "id", time = "time",
      alive = "alive", in_state = "in_state", y = "Y",
      baseline = c("age", "sex"),
      tv_names = c("L1", "L2"),
      a_names  = "A",
      sl_remain = sl_fast, sl_death = sl_fast,
      sl_recursive = sl_fast, sl_y = sl_fast,
      outcome_family = "binomial",
      k = 1L, inner_v = 2L, parallel = FALSE, seed = 1L,
      policy_spec_fun = policy_up,
      pool_g_death = TRUE, pool_q_exit = TRUE, pool_time = mode
    )
    expect_true(is.finite(res$psi), info = paste("mode =", mode))
    expect_gt(res$se, 0, label = paste("se, mode =", mode))
    expect_identical(res$settings$pool_time, mode)
  }
})

# ── itmle ─────────────────────────────────────────────────────────────────

test_that("itmle returns finite psi, positive se, and ordered ci", {
  res <- itmle(
    df               = df_smoke,
    weight_object    = wr_smoke,
    tmax             = 3L,
    id               = "id",
    time             = "time",
    alive            = "alive",
    in_state         = "in_state",
    y                = "Y",
    baseline         = c("age", "sex"),
    tv_names         = c("L1", "L2"),
    a_names          = "A",
    sl_remain        = sl_fast,
    sl_death         = sl_fast,
    sl_recursive     = sl_fast,
    sl_y             = sl_fast,
    sl_tmle          = sl_tmle_simple,
    outcome_family   = "binomial",
    k                = 1L,
    inner_v          = 2L,
    v_target_itmle   = 2L,
    v_sl_inner_itmle = 2L,
    parallel         = FALSE,
    seed             = 1L,
    policy_spec_fun  = policy_up
  )

  expect_true(is.list(res))
  expect_true(is.finite(res$psi))
  expect_gte(res$psi, 0)
  expect_lte(res$psi, 1)
  expect_true(is.finite(res$se))
  expect_gt(res$se, 0)
  expect_length(res$ci, 2)
  expect_true(all(is.finite(res$ci)))
  expect_lt(res$ci[1], res$psi)
  expect_gt(res$ci[2], res$psi)
})

# ── qreg ──────────────────────────────────────────────────────────────────────

test_that("qreg without weight_object returns finite estimate and naive SE only", {
  res <- qreg(
    df              = df_smoke,
    tmax            = 3L,
    id              = "id",
    time            = "time",
    alive           = "alive",
    in_state        = "in_state",
    y               = "Y",
    baseline        = c("age", "sex"),
    tv_names        = c("L1", "L2"),
    a_names         = "A",
    sl_remain       = sl_fast,
    sl_death        = sl_fast,
    sl_recursive    = sl_fast,
    sl_y            = sl_fast,
    outcome_family  = "binomial",
    k               = 1L,
    inner_v         = 2L,
    v               = 2L,
    parallel        = FALSE,
    seed            = 1L,
    policy_spec_fun = policy_up
  )

  expect_true(is.list(res))
  expect_true(is.finite(res$psi))
  expect_gte(res$psi, 0)
  expect_lte(res$psi, 1)
  expect_true(is.finite(res$decomposition$se_naive))
  expect_gt(res$decomposition$se_naive, 0)
  expect_null(res$decomposition$se_eif)
  expect_equal(res$se, res$decomposition$se_naive)
  expect_length(res$ci, 2)
  expect_true(all(is.finite(res$ci)))
  expect_lt(res$ci[1], res$psi)
  expect_gt(res$ci[2], res$psi)
})

test_that("qreg with weight_object computes EIF-based SE and se == se_eif", {
  res <- qreg(
    df              = df_smoke,
    weight_object   = wr_smoke,
    tmax            = 3L,
    id              = "id",
    time            = "time",
    alive           = "alive",
    in_state        = "in_state",
    y               = "Y",
    baseline        = c("age", "sex"),
    tv_names        = c("L1", "L2"),
    a_names         = "A",
    sl_remain       = sl_fast,
    sl_death        = sl_fast,
    sl_recursive    = sl_fast,
    sl_y            = sl_fast,
    outcome_family  = "binomial",
    k               = 1L,
    inner_v         = 2L,
    parallel        = FALSE,
    seed            = 1L,
    policy_spec_fun = policy_up
  )

  expect_true(is.list(res))
  expect_true(is.finite(res$psi))
  expect_gte(res$psi, 0)
  expect_lte(res$psi, 1)
  expect_true(is.finite(res$decomposition$se_naive))
  expect_gt(res$decomposition$se_naive, 0)
  expect_true(is.finite(res$decomposition$se_eif))
  expect_gt(res$decomposition$se_eif, 0)
  expect_equal(res$se, res$decomposition$se_eif)
  expect_length(res$ci, 2)
  expect_true(all(is.finite(res$ci)))
  expect_lt(res$ci[1], res$psi)
  expect_gt(res$ci[2], res$psi)
})
