library(SuperLearner)

# End-to-end accuracy tests against known Monte Carlo truths.
#
# Three clinically motivated DGPs:
#
#   1. Binary treatment: ICU care bundle (A ∈ {0,1}).
#      Policy: mandate bundle for patients with elevated biomarker (L1 > 1).
#      Returns 0L / 1L — no fractional treatment values.
#
#   2. Continuous treatment: vasopressor dose (A ≥ 0, mcg/kg/min).
#      Policy: increase dose by 0.3 units, capped at 2.0.
#
#   3. Two treatments: antibiotic decision (A1 ∈ {0,1}) + steroid dose (A2 ≥ 0).
#      Policy: mandate antibiotics for L1 > 1; increase steroid dose by 0.2.
#
# For each DGP, E[Y(d)] is obtained by simulating n = 500,000 subjects under
# the policy (Monte Carlo error < 0.002 for all scenarios).
#
# SDR and iTMLE: test |psi_hat - truth| < 3 * se.
# qreg: calibration check only — plug-in estimate ≈ observed mean under the
#   natural course.  Accuracy against truth is not expected and not tested.

sl_fast    <- c("SL.mean", "SL.glm")
sl_tgt_glm <- c("SL.tgt.glm")

# ── Shared helper ─────────────────────────────────────────────────────────────

check_accuracy <- function(res, truth, label) {
  diff <- abs(res$psi - truth)
  tol  <- 3 * res$se
  message(sprintf(
    "[accuracy] %s: psi=%.4f  truth=%.4f  se=%.4f  |diff|/3se=%.2f",
    label, res$psi, truth, res$se, diff / tol
  ))
  expect_true(is.finite(res$psi), label = paste(label, "psi finite"))
  expect_gt(res$se, 0,            label = paste(label, "se > 0"))
  expect_lt(diff, tol,
            label = sprintf("%s |psi-truth|=%.4f >= 3*se=%.4f", label, diff, tol))
}

# ══════════════════════════════════════════════════════════════════════════════
# DGP 1 — Binary treatment: ICU care bundle
# ══════════════════════════════════════════════════════════════════════════════
#
# Patients are ICU admissions.  A ∈ {0,1}: 1 = evidence-based care bundle
# implemented.  L1 is a continuous biomarker (higher = more severe).  L2 is a
# binary comorbidity flag.
#
# Natural propensity: bundle more likely for higher L1 and male sex.
# Policy: mandate bundle for any patient with L1 > 1.0 (moderate–severe
#   elevation).  Patients already receiving the bundle are unaffected.
#   d(A, H) = max(A, I(L1 > 1))  →  always 0L or 1L.

sim_bin <- function(n = 2000L, tmax = 5L, seed = 1L,
                    apply_policy = FALSE) {
  set.seed(seed)
  rows <- vector("list", n)

  for (i in seq_len(n)) {
    age <- round(rnorm(1, 65, 10))
    sex <- rbinom(1, 1, 0.5)
    L1  <- rnorm(1, 0, 1)
    L2  <- rbinom(1, 1, 0.4)
    pat <- list()

    for (t in seq_len(tmax)) {
      p_A   <- plogis(0.5 * L1 - 0.3 + 0.2 * sex)
      A_nat <- rbinom(1, 1, p_A)
      A     <- if (apply_policy && L1 > 1.0) 1L else A_nat

      p_die <- plogis(-4.0 + 0.3 * L1 - 0.1 * age / 10)
      p_dc  <- plogis(-2.5 + 0.8 * A  - 0.2 * L2)
      u     <- runif(1)

      if (u < p_die) {
        alive <- 0L; in_state <- 0L
      } else if (u < p_die + p_dc) {
        alive <- 1L; in_state <- 0L
      } else {
        alive <- 1L; in_state <- 1L
      }

      py <- if (!alive)     plogis(-3.0 + 0.1 * L1)
            else if (!in_state) plogis(1.5 + 0.4 * A - 0.1 * age / 10 + 0.2 * L1)
            else            plogis(-0.5 + 0.5 * A - 0.3 * L1 + 0.1 * L2)
      Y <- rbinom(1, 1, py)

      pat[[length(pat) + 1L]] <- data.frame(
        id = i, time = t,
        age = age, sex = sex,
        L1 = L1, L2 = L2,
        A = A, alive = alive, in_state = in_state, Y = Y
      )

      if (in_state == 0L) break
      if (t < tmax) {
        L1 <- L1 + rnorm(1, -0.2 * A, 0.3)
        L2 <- rbinom(1, 1, plogis(0.5 * L2 + 0.3 * A - 0.5))
      }
    }

    rows[[i]] <- do.call(rbind, pat)
  }

  do.call(rbind, rows)
}

# Policy: mandate bundle for all patients with L1 > 1.0.
# Returns 0L or 1L — binary, no fractional treatment values.
policy_bin <- function(D_block, t, a_names) {
  out <- D_block[, ..a_names, drop = FALSE]
  out[[a_names[1]]] <- pmax(D_block[[a_names[1]]],
                             as.integer(D_block[["L1"]] > 1.0))
  out
}

df_mc_bin   <- sim_bin(n = 500000L, tmax = 5L, seed = 77L, apply_policy = TRUE)
last_bin    <- df_mc_bin[!duplicated(df_mc_bin$id, fromLast = TRUE), ]
truth_bin   <- mean(last_bin$Y)
rm(df_mc_bin, last_bin)
message(sprintf("[accuracy] truth_bin (n=500k): E[Y(d)] = %.4f", truth_bin))

df_est_bin <- sim_bin(n = 2000L, tmax = 5L, seed = 42L)
obs_mean_bin <- mean(df_est_bin[!duplicated(df_est_bin$id, fromLast = TRUE), "Y"])

wr_bin <- density_ratio(
  df              = df_est_bin,
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
  policy_spec_fun = policy_bin
)

# ── sdr [binary A] ────────────────────────────────────────────────────────────

test_that("sdr [binary A]: estimate within 3 SE of MC truth", {
  res <- sdr(
    df              = df_est_bin,
    weight_object   = wr_bin,
    tmax            = 5L,
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
    policy_spec_fun = policy_bin
  )
  check_accuracy(res, truth_bin, "sdr [binary A]")
})

# ── itmle [binary A] ──────────────────────────────────────────────────────────

test_that("itmle [binary A]: estimate within 3 SE of MC truth", {
  res <- itmle(
    df               = df_est_bin,
    weight_object    = wr_bin,
    tmax             = 5L,
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
    sl_tmle          = sl_tgt_glm,
    outcome_family   = "binomial",
    k                = 1L,
    inner_v          = 2L,
    v_target_itmle   = 2L,
    v_sl_inner_itmle = 2L,
    parallel         = FALSE,
    seed             = 1L,
    policy_spec_fun  = policy_bin
  )
  check_accuracy(res, truth_bin, "itmle [binary A]")
})

# ── qreg [binary A]: calibration check ───────────────────────────────────────
#
# qreg is a plug-in with first-order bias and should not be compared to the
# causal truth.  The correct diagnostic is: natural-course plug-in ≈ observed
# mean of Y.  A poor fit here indicates Q-model misspecification.

test_that("qreg [binary A]: natural-course plug-in close to observed mean", {
  policy_nat <- function(D_block, t, a_names) D_block[, ..a_names, drop = FALSE]

  res <- qreg(
    df              = df_est_bin,
    tmax            = 5L,
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
    policy_spec_fun = policy_nat
  )

  diff_cal <- abs(res$estimate - obs_mean_bin)
  message(sprintf(
    "[accuracy] qreg [binary A] nat: plug-in=%.4f  obs_mean=%.4f  diff=%.4f",
    res$estimate, obs_mean_bin, diff_cal
  ))
  expect_true(is.finite(res$estimate))
  expect_gt(res$se_naive, 0)
  expect_lt(diff_cal, 3 * res$se_naive,
            label = "qreg nat-course plug-in outside 3*se_naive of observed mean")
})


# ══════════════════════════════════════════════════════════════════════════════
# DGP 2 — Continuous treatment: vasopressor dose
# ══════════════════════════════════════════════════════════════════════════════
#
# Patients in septic shock.  A ≥ 0: norepinephrine dose (mcg/kg/min).
# Natural dose is severity-driven (higher L1 → higher dose) with noise.
# Policy: increase dose by 0.3 mcg/kg/min, capped at 2.0.
# Additive shift of a continuous dose is the canonical LMTP use case.

DOSE_SHIFT <- 0.3

sim_cont <- function(n = 2000L, tmax = 5L, seed = 1L, dose_shift = 0) {
  set.seed(seed)
  rows <- vector("list", n)

  for (i in seq_len(n)) {
    age <- round(rnorm(1, 65, 10))
    L1  <- rnorm(1, 0, 1)    # hemodynamic instability (higher = worse)
    pat <- list()

    for (t in seq_len(tmax)) {
      A_nat <- pmax(0, rnorm(1, 0.6 + 0.5 * L1, 0.3))
      A     <- pmin(A_nat + dose_shift, 2.0)

      p_die <- plogis(-4.5 + 0.5 * L1 - 0.08 * age / 10 - 0.5 * A)
      p_dc  <- plogis(-2.0 + 0.5 * A  - 0.3 * L1)
      u     <- runif(1)

      if (u < p_die) {
        alive <- 0L; in_state <- 0L
      } else if (u < p_die + p_dc) {
        alive <- 1L; in_state <- 0L
      } else {
        alive <- 1L; in_state <- 1L
      }

      py <- if (!alive)     plogis(-3.0 + 0.1 * L1)
            else            plogis(0.5  + 0.9 * A  - 0.4 * L1 - 0.05 * age / 10)
      Y <- rbinom(1, 1, py)

      pat[[length(pat) + 1L]] <- data.frame(
        id = i, time = t,
        age = age, L1 = L1,
        A = A, alive = alive, in_state = in_state, Y = Y
      )

      if (in_state == 0L) break
      if (t < tmax)
        L1 <- L1 + rnorm(1, -0.4 * A, 0.4)
    }

    rows[[i]] <- do.call(rbind, pat)
  }

  do.call(rbind, rows)
}

policy_cont <- function(D_block, t, a_names) {
  out <- D_block[, ..a_names, drop = FALSE]
  out[[a_names[1]]] <- pmin(D_block[[a_names[1]]] + DOSE_SHIFT, 2.0)
  out
}

df_mc_cont  <- sim_cont(n = 500000L, tmax = 5L, seed = 77L, dose_shift = DOSE_SHIFT)
last_cont   <- df_mc_cont[!duplicated(df_mc_cont$id, fromLast = TRUE), ]
truth_cont  <- mean(last_cont$Y)
rm(df_mc_cont, last_cont)
message(sprintf("[accuracy] truth_cont (n=500k): E[Y(d)] = %.4f", truth_cont))

df_est_cont <- sim_cont(n = 2000L, tmax = 5L, seed = 42L)

wr_cont <- density_ratio(
  df              = df_est_cont,
  a_names         = "A",
  tmax            = 5L,
  baseline        = "age",
  tv_names        = "L1",
  sl_g            = sl_fast,
  k               = 1L,
  inner_v         = 2L,
  v               = 2L,
  seed            = 1L,
  id              = "id",
  time            = "time",
  parallel_t      = FALSE,
  policy_spec_fun = policy_cont
)

# ── sdr [continuous A] ────────────────────────────────────────────────────────

test_that("sdr [continuous A]: estimate within 3 SE of MC truth", {
  res <- sdr(
    df              = df_est_cont,
    weight_object   = wr_cont,
    tmax            = 5L,
    id              = "id",
    time            = "time",
    alive           = "alive",
    in_state        = "in_state",
    y               = "Y",
    baseline        = "age",
    tv_names        = "L1",
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
    policy_spec_fun = policy_cont
  )
  check_accuracy(res, truth_cont, "sdr [continuous A]")
})

# ── itmle [continuous A] ──────────────────────────────────────────────────────

test_that("itmle [continuous A]: estimate within 3 SE of MC truth", {
  res <- itmle(
    df               = df_est_cont,
    weight_object    = wr_cont,
    tmax             = 5L,
    id               = "id",
    time             = "time",
    alive            = "alive",
    in_state         = "in_state",
    y                = "Y",
    baseline         = "age",
    tv_names         = "L1",
    a_names          = "A",
    sl_remain        = sl_fast,
    sl_death         = sl_fast,
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
    policy_spec_fun  = policy_cont
  )
  check_accuracy(res, truth_cont, "itmle [continuous A]")
})

# ── qreg [continuous A]: calibration + finite shifted estimate ────────────────

test_that("qreg [continuous A]: natural-course calibration and finite shifted estimate", {
  policy_nat <- function(D_block, t, a_names) D_block[, ..a_names, drop = FALSE]
  obs_mean_cont <- mean(df_est_cont[!duplicated(df_est_cont$id, fromLast = TRUE), "Y"])

  res_nat <- qreg(
    df              = df_est_cont,
    tmax            = 5L,
    id              = "id",
    time            = "time",
    alive           = "alive",
    in_state        = "in_state",
    y               = "Y",
    baseline        = "age",
    tv_names        = "L1",
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
    policy_spec_fun = policy_nat
  )

  diff_cal <- abs(res_nat$estimate - obs_mean_cont)
  message(sprintf(
    "[accuracy] qreg [cont A] nat: plug-in=%.4f  obs_mean=%.4f  diff=%.4f",
    res_nat$estimate, obs_mean_cont, diff_cal
  ))
  expect_true(is.finite(res_nat$estimate))
  expect_gt(res_nat$se_naive, 0)
  expect_lt(diff_cal, 3 * res_nat$se_naive,
            label = "qreg cont nat-course plug-in outside 3*se_naive of observed mean")

  res_shf <- qreg(
    df              = df_est_cont,
    tmax            = 5L,
    id              = "id",
    time            = "time",
    alive           = "alive",
    in_state        = "in_state",
    y               = "Y",
    baseline        = "age",
    tv_names        = "L1",
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
    policy_spec_fun = policy_cont
  )

  message(sprintf(
    "[accuracy] qreg [cont A] shifted: estimate=%.4f  se_naive=%.4f",
    res_shf$estimate, res_shf$se_naive
  ))
  expect_true(is.finite(res_shf$estimate))
  expect_gt(res_shf$se_naive, 0)
  expect_true(res_shf$psi_shifted != res_shf$psi_natural,
              label = "qreg cont shifted estimate equals natural (policy had no effect)")
})


# ══════════════════════════════════════════════════════════════════════════════
# DGP 3 — Two treatments: antibiotic (binary) + steroid dose (continuous)
# ══════════════════════════════════════════════════════════════════════════════
#
# ICU sepsis patients with two independent treatment decisions:
#   A1 ∈ {0,1}: antibiotic prescription (binary).
#   A2 ≥ 0:     corticosteroid dose (mcg/kg/h, continuous).
#
# Policy:
#   A1: mandate antibiotics for patients with L1 > 1.0 (severe sepsis).
#   A2: increase steroid dose by 0.2 units, capped at 2.0.

sim_two <- function(n = 2000L, tmax = 5L, seed = 1L, apply_policy = FALSE) {
  set.seed(seed)
  rows <- vector("list", n)

  for (i in seq_len(n)) {
    age <- round(rnorm(1, 65, 10))
    sex <- rbinom(1, 1, 0.5)
    L1  <- rnorm(1, 0, 1)     # severity (infection + hemodynamics)
    pat <- list()

    for (t in seq_len(tmax)) {
      # A1: antibiotic
      p_A1   <- plogis(0.5 * L1 + 0.2 * sex - 0.3)
      A1_nat <- rbinom(1, 1, p_A1)
      A1     <- if (apply_policy && L1 > 1.0) 1L else A1_nat

      # A2: steroid dose (independent of antibiotic decision)
      A2_nat <- pmax(0, rnorm(1, 0.4 + 0.3 * L1, 0.3))
      A2     <- if (apply_policy) pmin(A2_nat + 0.2, 2.0) else A2_nat

      p_die <- plogis(-4.0 + 0.4 * L1 - 0.1 * age / 10 - 0.6 * A1 - 0.3 * A2)
      p_dc  <- plogis(-2.5 + 0.7 * A1 + 0.3 * A2 - 0.2 * L1)
      u     <- runif(1)

      if (u < p_die) {
        alive <- 0L; in_state <- 0L
      } else if (u < p_die + p_dc) {
        alive <- 1L; in_state <- 0L
      } else {
        alive <- 1L; in_state <- 1L
      }

      py <- if (!alive) plogis(-3.0 + 0.1 * L1)
            else        plogis(0.5 + 0.6 * A1 + 0.5 * A2 - 0.3 * L1 - 0.05 * age / 10)
      Y <- rbinom(1, 1, py)

      pat[[length(pat) + 1L]] <- data.frame(
        id = i, time = t,
        age = age, sex = sex, L1 = L1,
        A1 = A1, A2 = A2,
        alive = alive, in_state = in_state, Y = Y
      )

      if (in_state == 0L) break
      if (t < tmax)
        L1 <- L1 + rnorm(1, -0.2 * A1 - 0.1 * A2, 0.3)
    }

    rows[[i]] <- do.call(rbind, pat)
  }

  do.call(rbind, rows)
}

# Policy: A1 → mandate for L1 > 1 (binary output); A2 → +0.2, capped at 2.0
policy_two <- function(D_block, t, a_names) {
  out <- D_block[, ..a_names, drop = FALSE]
  out[[a_names[1]]] <- pmax(D_block[[a_names[1]]],
                             as.integer(D_block[["L1"]] > 1.0))
  out[[a_names[2]]] <- pmin(D_block[[a_names[2]]] + 0.2, 2.0)
  out
}

df_mc_two  <- sim_two(n = 500000L, tmax = 5L, seed = 77L, apply_policy = TRUE)
last_two   <- df_mc_two[!duplicated(df_mc_two$id, fromLast = TRUE), ]
truth_two  <- mean(last_two$Y)
rm(df_mc_two, last_two)
message(sprintf("[accuracy] truth_two (n=500k): E[Y(d)] = %.4f", truth_two))

df_est_two <- sim_two(n = 2000L, tmax = 5L, seed = 42L)

wr_two <- density_ratio(
  df              = df_est_two,
  a_names         = c("A1", "A2"),
  tmax            = 5L,
  baseline        = c("age", "sex"),
  tv_names        = "L1",
  sl_g            = sl_fast,
  k               = 1L,
  inner_v         = 2L,
  v               = 2L,
  seed            = 1L,
  id              = "id",
  time            = "time",
  parallel_t      = FALSE,
  policy_spec_fun = policy_two
)

# ── sdr [binary + continuous A] ───────────────────────────────────────────────

test_that("sdr [binary + continuous A]: estimate within 3 SE of MC truth", {
  res <- sdr(
    df              = df_est_two,
    weight_object   = wr_two,
    tmax            = 5L,
    id              = "id",
    time            = "time",
    alive           = "alive",
    in_state        = "in_state",
    y               = "Y",
    baseline        = c("age", "sex"),
    tv_names        = "L1",
    a_names         = c("A1", "A2"),
    sl_remain       = sl_fast,
    sl_death        = sl_fast,
    sl_recursive    = sl_fast,
    sl_y            = sl_fast,
    outcome_family  = "binomial",
    k               = 1L,
    inner_v         = 2L,
    parallel        = FALSE,
    seed            = 1L,
    policy_spec_fun = policy_two
  )
  check_accuracy(res, truth_two, "sdr [binary + continuous A]")
})

# ── itmle [binary + continuous A] ─────────────────────────────────────────────

test_that("itmle [binary + continuous A]: estimate within 3 SE of MC truth", {
  res <- itmle(
    df               = df_est_two,
    weight_object    = wr_two,
    tmax             = 5L,
    id               = "id",
    time             = "time",
    alive            = "alive",
    in_state         = "in_state",
    y                = "Y",
    baseline         = c("age", "sex"),
    tv_names         = "L1",
    a_names          = c("A1", "A2"),
    sl_remain        = sl_fast,
    sl_death         = sl_fast,
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
    policy_spec_fun  = policy_two
  )
  check_accuracy(res, truth_two, "itmle [binary + continuous A]")
})
