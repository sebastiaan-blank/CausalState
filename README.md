
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CausalState

<!-- badges: start -->

[![License: AGPL
v3](https://img.shields.io/badge/License-AGPL_v3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0.html)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

**CausalState** provides Sequential Doubly Robust (SDR) and
infinite-dimensional Targeted Maximum Likelihood (iTMLE) estimators for
longitudinal modified treatment policies (MTPs) in care-episode settings
where patients can transition irreversibly out of an active state – for
example, ICU discharge or in-hospital death.

The distinctive feature of the package is that the MTP can **shift the
transition dynamics themselves**, not only the terminal outcome. A
policy that shortens ICU stay changes both when patients leave and what
outcomes they experience after leaving; CausalState handles both
simultaneously.

## Problem setting

The estimators apply when:

- Data are in **long format**, one row per patient per time point, with
  rows present only while the patient is in the active state.
- At each time point a patient can remain, discharge (alive exit), or
  die (dead exit); no return to the active state is possible.
- The treatment policy of interest is a **modified treatment policy**
  (MTP): a data-adaptive shift of the observed treatment rather than a
  static intervention.
- Interest is in the **mean counterfactual outcome** under the MTP,
  estimated from observational or trial data.

The package is designed for **medium-horizon** episodes of 5–10 time
points. Beyond roughly ten points, cumulative density-ratio products
tend to collapse positivity regardless of estimator choice.

## Estimators

Both estimators are sequentially doubly robust (2^K-robust, Luedtke et
al. 2017): consistent whenever, at each time point, either the treatment
model or the outcome model is correctly specified.

|                    | SDR (`sdr()`)                                             | iTMLE (`itmle()`)                                                                      |
|--------------------|-----------------------------------------------------------|----------------------------------------------------------------------------------------|
| Update step        | EIF pseudo-outcome (Diaz et al. 2021)                     | Infinite-dimensional TMLE fluctuation (Luedtke et al. 2017)                            |
| SE formula         | Centered: `sd(IC) / sqrt(n)`, E\[IC\] = 0 by construction | Second-moment: `sqrt(mean(IC^2) / n)`, conservative when targeting is near-convergence |
| Natural-course run | Collapses to `mean(Y)` – not a model check                | Does not collapse – genuine Q-model calibration check                                  |
| Extra inputs       | None beyond Q/g libraries                                 | Targeting SL library (`sl_tmle`)                                                       |

`qreg()` is a pure Q-recursion plug-in (no update, no DR guarantees)
included as a weight-independent sensitivity check.

## Workflow

    fit_g()  ->  sdr() / itmle() / qreg()

`fit_g()` runs first for SDR and iTMLE. It fits treatment models at each
time point (classification by default), computes the instantaneous
density ratios r_t = dP~(A_t\|H_t) / dP(A_t\|H_t), and packages them
with fold assignments that the downstream estimators inherit. Running it
once and passing the result to multiple estimators is the intended
pattern.

Categorical treatment is declared in the data as a factor. For ordinal
treatment, set the dose order explicitly:

``` r
patient_data$A <- factor(patient_data$A, levels = 0:3, ordered = TRUE)
```

`fit_g(method = "propensity", ...)` then estimates the ordinal dose
probabilities using sequential binomial SuperLearners from `sl_g`.
Classification remains the default and also supports ordinary factors.

CausalState reads the factor levels automatically and gives every
learner the same numeric reference indicators for treatment and its
history, including after policy shifts. GLMs can estimate separate
category effects. Policies receive the original factor column, with
levels and ordering preserved. The levels and ordered status are saved
in the `fit_g()` result and checked by `sdr()` and `itmle()`. `qreg()`
reads its input factors directly. There is no separate levels argument;
numeric binary and continuous treatments retain numeric coding.

## Quick example

``` r
library(CausalState)
library(SuperLearner)

# --- 1. Define the policy -------------------------------------------------
# Threshold policy: treat everyone with L1 > 1.0
policy_bin <- function(D_block, t, a_names) {
  out <- D_block[, ..a_names, drop = FALSE]
  out[[a_names[1]]] <- pmax(D_block[[a_names[1]]],
                             as.integer(D_block[["L1"]] > 1.0))
  out
}

# --- 2. Density ratios ----------------------------------------------------
sl_lib <- c("SL.mean", "SL.glm")   # replace with richer library in practice

wr <- fit_g(
  df              = patient_data,   # long-format data frame
  a_names         = "A",
  tmax            = 7L,
  baseline        = c("age", "sex"),
  tv_names        = c("L1", "L2"),
  sl_g            = sl_lib,
  k               = 5L,
  inner_v         = 5L,
  v               = 5L,
  seed            = 1L,
  id              = "id",
  time            = "time",
  policy_spec_fun = policy_bin
)

# --- 3a. SDR estimate -----------------------------------------------------
res_sdr <- sdr(
  df              = patient_data,
  weight_object   = wr,
  tmax            = 7L,
  id              = "id", time = "time",
  alive           = "alive", in_state = "in_state",
  y               = "Y",
  baseline        = c("age", "sex"),
  tv_names        = c("L1", "L2"),
  a_names         = "A",
  sl_remain       = sl_lib, sl_death = sl_lib,
  sl_recursive    = sl_lib, sl_y = sl_lib,
  outcome_family  = "binomial",
  k               = 5L, inner_v = 5L,
  seed            = 1L,
  policy_spec_fun = policy_bin
)

cat(sprintf("SDR:   psi = %.3f  SE = %.3f  95%% CI [%.3f, %.3f]\n",
            res_sdr$psi, res_sdr$se, res_sdr$ci[1], res_sdr$ci[2]))

# --- 3b. iTMLE estimate ---------------------------------------------------
res_itmle <- itmle(
  df              = patient_data,
  weight_object   = wr,
  tmax            = 7L,
  id              = "id", time = "time",
  alive           = "alive", in_state = "in_state",
  y               = "Y",
  baseline        = c("age", "sex"),
  tv_names        = c("L1", "L2"),
  a_names         = "A",
  sl_remain       = sl_lib, sl_death = sl_lib,
  sl_recursive    = sl_lib, sl_y = sl_lib,
  sl_target       = sl_tmle,
  outcome_family  = "binomial",
  k               = 5L, inner_v = 5L,
  seed            = 1L,
  policy_spec_fun = policy_bin
)

cat(sprintf("iTMLE: psi = %.3f  SE = %.3f  95%% CI [%.3f, %.3f]\n",
            res_itmle$psi, res_itmle$se, res_itmle$ci[1], res_itmle$ci[2]))

# --- 4. Risk difference ---------------------------------------------------
ctr <- contrast(res_sdr, res_sdr_nat)   # intervention vs natural course
ctr$RD; ctr$ci_RD
```

## Key design choices

**Asymmetric g/Q regularisation.** Treatment models (g) are tuned
sharper for accurate density-ratio discrimination; outcome models (Q)
are more heavily regularised to suppress recursion noise. This asymmetry
is the primary tuning lever.

**Uniform clipping.** Every SL prediction across the entire pipeline (g
and Q models at every time point) is clipped to \[bounds, 1-bounds\]
using the same `bounds` parameter (default 1e-5). The only exception is
the Wu-Benkeser direct density-ratio metalearner, which clips in
density-ratio space via `dr_floor` – see `vignette("wb-metalearner")`.

**Weight reuse.** `fit_g()` is designed to be run once and shared across
`sdr()`, `itmle()`, and competing-event variants. Trimming is applied
globally at consumption time, so all estimators that share a weight
object operate on identically trimmed weights and produce directly
comparable estimates.

**Parallelism via `mclapply`.** Process-level parallelism is available
at the fold level (`parallel = TRUE`) and within-fold regression level
(`reg_workers`). Does not work on Windows. Set BLAS and learner thread
counts to 1 when enabling process-level parallelism to avoid
oversubscription.

## Installation

``` r
# install.packages("remotes")
remotes::install_github("sebastiaan-blank/CausalState")
```

## Getting started

See `vignette("getting-started")` for a complete worked example with a
simulated ICU dataset, including data structure, policy definition,
SuperLearner library choices, diagnostics, and a competing-event
analysis.

## Citation

    Blank S (2026). CausalState: SDR and iTMLE for State-Aware Longitudinal
    Modified Treatment Policies. R package version 0.10.2.
    https://github.com/sebastiaan-blank/CausalState

## References

Diaz I, Williams N, Hoffman KL, Schenck EJ (2023). Nonparametric Causal
Effects Based on Longitudinal Modified Treatment Policies. *JASA*
118(542):846-857.

Luedtke AR, Sofrygin O, van der Laan MJ, Carone M (2017). Sequential
Double Robustness in Right-Censored Longitudinal Models.
arXiv:1705.02459.

Williams NT, Diaz I (2023). lmtp: An R package for estimating the causal
effects of modified treatment policies. *Observational Studies*.

## License

AGPL-3. See `LICENSE.md` for details.
