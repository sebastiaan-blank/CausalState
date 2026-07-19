# CLAUDE.md — CausalState

R package for longitudinal causal inference under state transitions.

---

## Working agreement (read first, applies always)

- **Statistical correctness is the top priority** — above style, speed, and
  tidiness. A change that risks altering numerical or statistical behaviour is
  not an acceptable change, however clean it looks.
- **Discuss before changing.** Propose edits and get explicit sign-off before
  touching anything that affects results: estimator logic, the EIF / targeting
  steps, density-ratio handling, SuperLearner configuration, or family/link
  choices. No silent "improvements" or refactors-in-passing.
- **Assume the source has changed since the last session.** On every startup,
  re-read the relevant files before reasoning about or editing them. Do not
  trust earlier conversation context about what the code currently contains.
- **The current work is reorganisation, and it must be behaviour-preserving.**
  Splitting files and renaming must not change outputs. Where feasible, confirm
  results are unchanged before and after.

---

## What CausalState is

A general-purpose package for **longitudinal causal analysis under state
transitions** — settings such as ICU → discharge or ward → discharge, where the
discharge state and discharge time are themselves flexibly modelled so that the
treatment policy can shift them, not only the terminal outcome Y. "The
intervention can alter the transition dynamics, not just Y" is the distinctive
idea the package is built around.

Two main estimators, both implementations of **published** methods (not novel
methodology — see References):

- **SDR** (sequentially doubly robust) — the LMTP-SDR estimator of Díaz et al.,
  which generalises the SDR construction of Luedtke et al. (2017) to the
  modified-treatment-policy / density-ratio setting. 2τ-robust.
- **iTMLE** (infinite-dimensional TMLE) — from the same Luedtke et al. (2017)
  construction, Section 5 / Algorithm 4. Also sequentially doubly robust
  (2^τ-robust): achieves SDR through an infinite-dimensional targeting step
  rather than the backward regression approach used by SDR.

Both are heavily based on the **LMTP** framework and the `lmtp` package.

The package also includes a **simpler competing-event variant** of both
estimators for the strictly-survival case (ICU / in-hospital mortality). It
lives in the same package — rather than standalone — specifically to **reuse the
density-ratio weights** computed for the main estimators.

> NOTE: applications (e.g. any MIMIC-IV diuretics analysis) are *not* part of
> this package; they are downstream users of it.

---

## Current status & priorities

- Main estimator functions: **done**; core estimator logic is unchanged
  recently — the main recent additions are an expanded set of diagnostics.
- The **survival / competing-event** estimator: **needs to be rewritten**, but
  is **not** a current priority.
- **Active task: reorganisation** of the codebase (see below).

---

## Active task — reorganisation

Goal is structure, with **zero behavioural change**. Split the two monolithic
scripts into files organised by concern. Note: R packages share a single flat
namespace, so the `R/` file layout is for human navigation, not enforced
scoping — the discipline is that a shared helper *physically lives in a shared
file*, never in an estimator file.

Proposed `R/` layout (layered: generic → fit engine → estimator-shared →
estimator-specific → diagnostics):

```
R/
  general_helpers.R  # generic utilities (validation, small math/data helpers); no estimator logic
  sl.R               # SuperLearner FITTING machinery: stratification, fold orchestration, fit driver.
                     #   Standard learner wrappers are external and NOT included here.
                     #   Used by density_ratio, sdr, and itmle.
  density_ratio.R    # g-side density-ratio estimation; consumes sl.R
  shared_helpers.R   # helpers shared by sdr & itmle (e.g. EIF assembly, pseudo-outcome scaffolding)
  sdr.R              # SDR estimator + SDR-specific helpers + SDR competing-event / survival variant
  itmle.R            # iTMLE estimator + iTMLE-specific helpers + iTMLE survival variant
  sl_itmle.R         # custom iTMLE SL wrappers — custom because the targeting offset
                     #   (offset-as-column hack) isn't handled by standard wrappers. Own file
                     #   so it also serves as an easy, tunable worked example users can adapt.
  diagnostics.R      # diagnostics: fold-level / ESS / EIF
```

Two helper tiers, deliberately distinct: `general_helpers.R` is estimator-
agnostic plumbing; `shared_helpers.R` is estimation logic common to both SDR and
iTMLE. Estimator-specific helpers stay inline in `sdr.R` / `itmle.R`.

On SL setup: the package does not hand-hold general SuperLearner configuration —
choosing learners for the g-side and SDR Q-models is the user's job (if you
can't set up an SL, you're not ready for this package yet). The iTMLE wrappers in
`sl_itmle.R` are the exception: they must be custom because of the targeting
offset, and they live in their own file both for that reason and so they double
as a clear, tunable example users can adapt (swap learners, change xgboost
settings, etc.).

Naming, now that the file carries the scope:

- Estimator-only functions keep `sdr_` / `itmle_` and live in their own file.
- Helpers used by **both** estimators move into `shared_helpers.R` and **shed
  the `sdr_` prefix** — the current `sdr_` on shared helpers is a misnomer and is
  the main thing to fix in the rename pass.
- Keep names short; a tiny module tag on shared helpers is optional since the
  file already signals scope.
- Reconcile divergent variable names for the same quantity across functions.

---

## Architecture & design decisions (settled — do not relitigate)

- **Flat functional design.** Plain functions; no S4 / R6 / OOP layers. Keep the
  call graph shallow and explicit.
- **Parallelism via `mclapply`, not `future`.** Do not introduce a
  `future` / `furrr` backend.
- **Separate g and Q modules.** Treatment/censoring (g) and outcome (Q)
  machinery stay in distinct modules. The density-ratio / estimator components
  are separated both to enable weight reuse across estimators and to improve
  speed.
- **Hand-maintained NAMESPACE.** Do not regenerate it from roxygen; edit by hand
  and keep exports deliberate.
- **License: AGPL-3.**

---

## SuperLearner conventions

The user supplies their own SuperLearner stack — these are **design principles
from development**, not configs the package imposes.

- **Asymmetric g/Q regularisation:** g-models tuned sharper (better
  density-ratio discrimination); Q-models more regularised (stability). This
  asymmetry is the guiding principle across all learner choices.
- **Q must cover the intervention trajectory, not just the natural arm.** A
  model that fits E[Y|H, A_obs] well may extrapolate poorly under the shift.
  The Q SL library needs enough flexibility to capture the counterfactual arm
  without chasing recursion noise — this is the primary design tension for the
  Q stack.
- **Recursive Q-regression is the hardest fitting problem** in these estimators.
  Noise accumulates backward through time, so the SL library for intermediate
  time steps should be more parsimonious than the terminal-Y library: favour
  shallow trees and strong regularisation over complex learners.
- **Density-ratio metalearner:** `density_ratio()` accepts `method = :wb_dr`
  to use the Wu-Benkeser DR metalearner (log-DR loss, analytical gradient) in
  place of the default NNLS. The WB metalearner combines base learners in DR
  space; `sl_predict` on a `:wb_dr` fit returns density ratios directly (not
  probabilities). Pair with `ExpTiltLearner` (parametric KLIEP) for a
  complementary direct-DR component alongside tree-based classifiers.
- **All learners are CPU-only when SL folds run in parallel** — concurrent GPU
  fits contend for VRAM. Set `nthread = 1` on tree learners and rely on
  Julia/R's process-level parallelism for the outer folds, not within-learner
  threading.
- The iTMLE targeting step uses the **offset-as-column hack** to survive
  SuperLearner's CV fold subsetting (offset carried as a data column, not via
  the `offset` argument, so it gets subset correctly inside folds). This is
  exactly why the iTMLE wrappers are custom and package-owned (`sl_itmle.R`). Do
  not "clean this up" into a normal offset — it breaks under cross-validation.

---

## Environment & build

- **OS:** Ubuntu 24.04 · **GPU:** NVIDIA RTX · **RAM:** 196 GB DDR5 (frequently heavily used by concurrent R analysis).
- **Threads:** `mclapply` controls process-level parallelism; set BLAS and
  learner thread counts **explicitly** to avoid oversubscription.
- Dev loop (adjust to the actual project setup):
  - `devtools::load_all()` — load during development
  - `devtools::test()` — run tests
  - `devtools::document()` — man pages only; **must not overwrite the
    hand-maintained NAMESPACE**
  - `R CMD check` / `devtools::check()` — full check before pushing

---

## Gotchas

- **xgboost + CUDA:** version 3.2.1.1 breaks GPU support — stay on **3.1.1.1**.
  If GPU training suddenly fails after an xgboost update, check this first.
- **Density-ratio ESS collapse:** cumulative density ratios compound across time
  and can collapse effective sample size. Watch ESS diagnostics rather than
  assuming the weights are well-behaved.

---

## Statistical invariants (confirmed current — must be preserved)

These are correct as of the current code and must **not** be altered without
discussion. If the code ever appears to violate one of these, treat it as a bug
to flag, not a style choice to "tidy." (Per the working agreement, still re-read
the source at startup rather than assuming these are where the logic lives.)

- `Q_rem` at the terminal block (`tt == TMAX`) uses **`sl_y`** (the outcome library) and the outcome-appropriate family (binomial for binary Y, gaussian for continuous Y). At earlier time steps it uses `sl_recursive` with gaussian family.
- In iTMLE, pseudo-outcome propagation uses `shf_train`, not the `sdr_eif`
  output.
- `Q_exit` and `Q_rem` predictions are kept **separate until after the EIF
  update** (mixing them earlier was a structural bug).

---

## References (core set for this file; fuller list belongs in the manual)

Methods implemented here, with attribution:

- **Díaz I, Williams N, Hoffman KL, Schenck EJ.** Nonparametric Causal Effects
  Based on Longitudinal Modified Treatment Policies. *JASA* 118(542):846–857
  (online 2021 / print 2023). doi:10.1080/01621459.2021.1955691.
  — the LMTP estimators implemented here (SDR + TMLE).
- **Luedtke AR, Sofrygin O, van der Laan MJ, Carone M.** Sequential Double
  Robustness in Right-Censored Longitudinal Models. arXiv:1705.02459 (2017).
  — the SDR / iTMLE construction the above generalises.
- **Rotnitzky A, Robins J, Babino L.** On the multiply robust estimation of the
  mean of the g-functional. arXiv:1705.08582 (2017).
  — multiple-robustness companion to SDR.
- **Williams NT, Díaz I.** lmtp: An R package for estimating the causal effects
  of modified treatment policies. *Observational Studies* (2023).
  — the package this work is based on.
- **Díaz I, Hoffman KL, Hejazi NS.** Causal survival analysis under competing
  risks using longitudinal modified treatment policies. *Lifetime Data Analysis*
  (2023). doi:10.1007/s10985-023-09606-7 (arXiv:2202.03513).
  — basis for the competing-event / survival variant.

Foundational MTP lineage (for the package manual / vignettes, not needed here):
Haneuse & Rotnitzky (2013, *Stat Med* 32(30):5260–5277); Díaz Muñoz & van der
Laan (2012, *Biometrics* 68(2):541–549); Díaz & van der Laan (2018, "Stochastic
Treatment Regimes").
