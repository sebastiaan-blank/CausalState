# sim_data.R
#
# Toy ICU-style panel generators used in vignettes, examples, and internal
# accuracy tests. Each one produces a long-format data.frame with one row per
# patient per time-step while in the active state; the trajectory ends when
# the patient transitions to an absorbing state (death or discharge).
#
# Column contract (shared across all three):
#   id, time, alive (0/1), in_state (0/1), Y (binary outcome)
# Covariates and treatments vary per DGP -- see individual @return sections.

#' Simulated ICU panel with a single binary treatment
#'
#' Generates a longitudinal panel of ICU-style trajectories with one binary
#' treatment (`A`), competing exit events (death vs discharge), and a binary
#' outcome (`Y`). Intended for vignettes, examples, and quick smoke tests.
#'
#' @param n Number of subjects.
#' @param tmax Maximum follow-up (time-steps while in-state).
#' @param seed Integer seed.
#' @param apply_policy Logical. When `TRUE`, forces `A = 1` for subjects with
#'   `L1 > 1.0` at each time-step -- used internally to generate Monte-Carlo
#'   truth for accuracy tests. Leave at `FALSE` for regular observed data.
#'
#' @return A `data.frame` with one row per subject-time. Columns:
#'   `id`, `time`, `age`, `sex`, `L1`, `L2`, `A`, `alive`, `in_state`, `Y`.
#'
#' @seealso [sim_cont()], [sim_multi()], [sdr()], [itmle()]
#'
#' @export
sim_bin <- function(n = 2000L, tmax = 5L, seed = 1L,
                    apply_policy = FALSE) {
  set.seed(seed)
  rows <- vector("list", n)

  for (i in seq_len(n)) {
    age <- round(stats::rnorm(1, 65, 10))
    sex <- stats::rbinom(1, 1, 0.5)
    L1  <- stats::rnorm(1, 0, 1)
    L2  <- stats::rbinom(1, 1, 0.4)
    pat <- list()

    for (t in seq_len(tmax)) {
      p_A   <- stats::plogis(0.5 * L1 - 0.3 + 0.2 * sex)
      A_nat <- stats::rbinom(1, 1, p_A)
      A     <- if (apply_policy && L1 > 1.0) 1L else A_nat

      p_die <- stats::plogis(-4.0 + 0.3 * L1 - 0.1 * age / 10)
      p_dc  <- stats::plogis(-2.5 + 0.8 * A  - 0.2 * L2)
      u     <- stats::runif(1)

      if (u < p_die) {
        alive <- 0L; in_state <- 0L
      } else if (u < p_die + p_dc) {
        alive <- 1L; in_state <- 0L
      } else {
        alive <- 1L; in_state <- 1L
      }

      py <- if (!alive)     stats::plogis(-3.0 + 0.1 * L1)
            else if (!in_state) stats::plogis(1.5 + 0.4 * A - 0.1 * age / 10 + 0.2 * L1)
            else            stats::plogis(-0.5 + 0.5 * A - 0.3 * L1 + 0.1 * L2)
      Y <- stats::rbinom(1, 1, py)

      pat[[length(pat) + 1L]] <- data.frame(
        id = i, time = t,
        age = age, sex = sex,
        L1 = L1, L2 = L2,
        A = A, alive = alive, in_state = in_state, Y = Y
      )

      if (in_state == 0L) break
      if (t < tmax) {
        L1 <- L1 + stats::rnorm(1, -0.2 * A, 0.3)
        L2 <- stats::rbinom(1, 1, stats::plogis(0.5 * L2 + 0.3 * A - 0.5))
      }
    }

    rows[[i]] <- do.call(rbind, pat)
  }

  do.call(rbind, rows)
}

#' Simulated ICU panel with a single continuous treatment
#'
#' Longitudinal ICU-style trajectories with one continuous treatment (`A`,
#' bounded roughly in `[0, 2]`), competing exit events, and a binary outcome.
#' Continuous-treatment analogue of [sim_bin()].
#'
#' @param n Number of subjects.
#' @param tmax Maximum follow-up.
#' @param seed Integer seed.
#' @param dose_shift Numeric dose increment applied to every time-step's
#'   observed `A`, capped at `2.0`. Used internally to generate Monte-Carlo
#'   truth. Leave at `0` for regular observed data.
#'
#' @return A `data.frame` with columns
#'   `id`, `time`, `age`, `L1`, `A`, `alive`, `in_state`, `Y`.
#'
#' @seealso [sim_bin()], [sim_multi()], [sdr()], [itmle()]
#'
#' @export
sim_cont <- function(n = 2000L, tmax = 5L, seed = 1L, dose_shift = 0) {
  set.seed(seed)
  rows <- vector("list", n)

  for (i in seq_len(n)) {
    age <- round(stats::rnorm(1, 65, 10))
    L1  <- stats::rnorm(1, 0, 1)
    pat <- list()

    for (t in seq_len(tmax)) {
      A_nat <- pmax(0, stats::rnorm(1, 0.6 + 0.5 * L1, 0.3))
      A     <- pmin(A_nat + dose_shift, 2.0)

      p_die <- stats::plogis(-4.5 + 0.5 * L1 - 0.08 * age / 10 - 0.5 * A)
      p_dc  <- stats::plogis(-2.0 + 0.5 * A  - 0.3 * L1)
      u     <- stats::runif(1)

      if (u < p_die) {
        alive <- 0L; in_state <- 0L
      } else if (u < p_die + p_dc) {
        alive <- 1L; in_state <- 0L
      } else {
        alive <- 1L; in_state <- 1L
      }

      py <- if (!alive)     stats::plogis(-3.0 + 0.1 * L1)
            else            stats::plogis(0.5  + 0.9 * A  - 0.4 * L1 - 0.05 * age / 10)
      Y <- stats::rbinom(1, 1, py)

      pat[[length(pat) + 1L]] <- data.frame(
        id = i, time = t,
        age = age, L1 = L1,
        A = A, alive = alive, in_state = in_state, Y = Y
      )

      if (in_state == 0L) break
      if (t < tmax)
        L1 <- L1 + stats::rnorm(1, -0.4 * A, 0.4)
    }

    rows[[i]] <- do.call(rbind, pat)
  }

  do.call(rbind, rows)
}

#' Simulated ICU panel with multiple treatments (any mix of binary + continuous)
#'
#' Longitudinal ICU-style trajectories with an arbitrary number of binary
#' and/or continuous treatments, competing exit events, and a binary outcome.
#' Generalises [sim_bin()] / [sim_cont()] to multi-treatment MTPs.
#'
#' Binary treatments are named `A_b1, A_b2, ...` (first `n_binary` columns);
#' continuous treatments are named `A_c1, A_c2, ...`. Each treatment enters
#' the exit-event and outcome models with a coefficient scaled by
#' `1 / sqrt(k)` where `k` is the treatment index within its type, so adding
#' more treatments does not blow up effect sizes.
#'
#' @param n Number of subjects.
#' @param tmax Maximum follow-up.
#' @param seed Integer seed.
#' @param n_binary Number of binary treatments (>= 0).
#' @param n_continuous Number of continuous treatments (>= 0).
#' @param apply_policy Logical. When `TRUE`, forces every binary treatment to
#'   `1` for subjects with `L1 > 1.0` and adds `+0.2` (capped at `2.0`) to
#'   every continuous treatment at each time-step -- used internally to
#'   generate Monte-Carlo truth for accuracy tests. Leave at `FALSE` for
#'   regular observed data.
#'
#' @return A `data.frame` with columns
#'   `id`, `time`, `age`, `sex`, `L1`, `A_b1..A_b{n_binary}`,
#'   `A_c1..A_c{n_continuous}`, `alive`, `in_state`, `Y`.
#'
#' @seealso [sim_bin()], [sim_cont()], [sdr()], [itmle()]
#'
#' @export
sim_multi <- function(n = 2000L, tmax = 5L, seed = 1L,
                      n_binary = 1L, n_continuous = 1L,
                      apply_policy = FALSE) {
  n_binary     <- as.integer(n_binary)
  n_continuous <- as.integer(n_continuous)
  if (n_binary < 0L || n_continuous < 0L)
    stop("sim_multi: n_binary and n_continuous must be >= 0", call. = FALSE)
  if (n_binary + n_continuous == 0L)
    stop("sim_multi: at least one treatment required", call. = FALSE)

  a_bin_names  <- if (n_binary     > 0L) paste0("A_b", seq_len(n_binary))     else character(0)
  a_cont_names <- if (n_continuous > 0L) paste0("A_c", seq_len(n_continuous)) else character(0)
  # per-treatment coefficient scaling: 1st gets full weight, later ones down-weighted
  w_bin  <- if (n_binary     > 0L) 1 / sqrt(seq_len(n_binary))     else numeric(0)
  w_cont <- if (n_continuous > 0L) 1 / sqrt(seq_len(n_continuous)) else numeric(0)

  set.seed(seed)
  rows <- vector("list", n)

  for (i in seq_len(n)) {
    age <- round(stats::rnorm(1, 65, 10))
    sex <- stats::rbinom(1, 1, 0.5)
    L1  <- stats::rnorm(1, 0, 1)
    pat <- list()

    for (t in seq_len(tmax)) {
      A_bin  <- numeric(n_binary)
      A_cont <- numeric(n_continuous)

      for (j in seq_len(n_binary)) {
        p_A <- stats::plogis(0.5 * L1 + 0.2 * sex - 0.3 - 0.1 * (j - 1))
        A_nat <- stats::rbinom(1, 1, p_A)
        A_bin[j] <- if (apply_policy && L1 > 1.0) 1L else A_nat
      }
      for (j in seq_len(n_continuous)) {
        A_nat <- pmax(0, stats::rnorm(1, 0.4 + 0.3 * L1 - 0.1 * (j - 1), 0.3))
        A_cont[j] <- if (apply_policy) pmin(A_nat + 0.2, 2.0) else A_nat
      }

      eff_bin  <- if (n_binary     > 0L) sum(w_bin  * A_bin)  else 0
      eff_cont <- if (n_continuous > 0L) sum(w_cont * A_cont) else 0

      p_die <- stats::plogis(-4.0 + 0.4 * L1 - 0.1 * age / 10 - 0.6 * eff_bin - 0.3 * eff_cont)
      p_dc  <- stats::plogis(-2.5 + 0.7 * eff_bin + 0.3 * eff_cont - 0.2 * L1)
      u     <- stats::runif(1)

      if (u < p_die) {
        alive <- 0L; in_state <- 0L
      } else if (u < p_die + p_dc) {
        alive <- 1L; in_state <- 0L
      } else {
        alive <- 1L; in_state <- 1L
      }

      py <- if (!alive) stats::plogis(-3.0 + 0.1 * L1)
            else        stats::plogis(0.5 + 0.6 * eff_bin + 0.5 * eff_cont - 0.3 * L1 - 0.05 * age / 10)
      Y <- stats::rbinom(1, 1, py)

      row <- data.frame(
        id = i, time = t,
        age = age, sex = sex, L1 = L1,
        alive = alive, in_state = in_state, Y = Y
      )
      for (j in seq_len(n_binary))     row[[a_bin_names[j]]]  <- A_bin[j]
      for (j in seq_len(n_continuous)) row[[a_cont_names[j]]] <- A_cont[j]
      # Reorder columns so treatments sit between L1 and alive
      lead  <- c("id", "time", "age", "sex", "L1")
      trail <- c("alive", "in_state", "Y")
      row <- row[, c(lead, a_bin_names, a_cont_names, trail)]

      pat[[length(pat) + 1L]] <- row

      if (in_state == 0L) break
      if (t < tmax)
        L1 <- L1 + stats::rnorm(1, -0.2 * eff_bin - 0.1 * eff_cont, 0.3)
    }

    rows[[i]] <- do.call(rbind, pat)
  }

  do.call(rbind, rows)
}
