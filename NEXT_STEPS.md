# Next Steps

This file is for you, not the package — it's not shipped (excluded by
`.Rbuildignore`). Delete or rename when no longer needed.

## 1. First things to edit

Before pushing to GitHub, do a find-and-replace across the skeleton:

- `YOUR-GITHUB-HANDLE` → your actual GitHub handle
  (in `DESCRIPTION`, `README.Rmd`, `README.md`)
- `REPLACE_WITH_YOUR_EMAIL@example.com` → your real email
  (in `DESCRIPTION` only)

Then update `LICENSE` if 2026 isn't right or if the copyright holder name
needs adjustment.

## 2. Push to GitHub

```bash
cd /path/to/CausalState
git init
git add -A
git commit -m "Initial commit: skeleton + imported research pipeline"
gh repo create CausalState --public --source=. --remote=origin --push
```

If `gh` is not installed: create the empty repo on github.com first, then:

```bash
git remote add origin git@github.com:YOUR-HANDLE/CausalState.git
git branch -M main
git push -u origin main
```

## 3. Verify the package builds

In an R session at the package root:

```r
# install once
install.packages(c("devtools", "usethis", "available", "pkgdown"))

# verify parse + load
devtools::load_all()       # should source R/functions.R cleanly

# generate man pages (will be empty until you add roxygen comments)
devtools::document()

# full check — expect WARNINGs about missing docs, no ERRORs
devtools::check()
```

Expected initial check output:
- **ERRORS**: 0 (if anything errors, send me the message and I'll fix)
- **WARNINGS**: ~1 ("no visible global function definition" for a few base
  R functions; fixed by adding `importFrom(stats, ...)` for each one as
  check flags them)
- **NOTES**: several ("undocumented code objects"); resolved as you add
  roxygen documentation per function

## 4. Refactoring roadmap

Order I'd suggest, smallest commit per step:

### Phase 1: make it builds-cleanly
1. Run `devtools::check()`. For each "no visible global function" warning,
   add the appropriate `importFrom()` line to `NAMESPACE`.
2. Decide whether `data_check` stays in dplyr/tidyr or gets rewritten in
   base / data.table for one fewer dependency.

### Phase 2: split functions.R into thematic files
Suggested split, based on the function clusters in your script:

```
R/
├── data-utils.R          # data_check, sdr_prepare_long, sdr_scale_y,
│                         # sdr_build_*, sdr_get_folds_from_weights,
│                         # expand_to_horizon
├── superlearner-utils.R  # .mc_cores, sl_block_fit, sl_predict
├── density-ratio.R       # density_ratio + sdr_build_density_ratios,
│                         # sdr_compute_weights
├── design.R              # sdr_mk_cols, sdr_make_design, sdr_sl_meta,
│                         # sdr_patch_shifted_design, itmle_mk_cols,
│                         # overwrite_policy_history_for_Q
├── absorb.R              # absorb_rule, sdr_apply_absorb_branch
├── eif.R                 # sdr_eif, itmle_target_diag_one,
│                         # itmle_weight, itmle_update
├── inference.R           # cluster_boot_se
├── learners-tgt.R        # SL.tgt.* + predict.SL.tgt.*
├── learners-tmle.R       # SL.tmle_*
├── sdr.R                 # sdr() + sdr_icu_mortality_competing()
└── itmle.R               # itmle() + itmle_target_fit, itmle_inner_target,
                          # target_lapply
```

Move one file at a time; run `devtools::test()` between moves; commit
each successful move separately.

### Phase 3: add roxygen documentation
For each exported function, add a roxygen header:

```r
#' One-line title (sentence case, no period)
#'
#' Optional paragraph describing what it does.
#'
#' @param data A long-format data.table with one row per id × t.
#' @param id Character. Name of the id column.
#' ...
#' @return Description of the return value.
#' @export
#' @examples
#' \dontrun{
#'   # example call
#' }
sdr <- function(data, id, time, ...) {
  ...
}
```

For internal helpers, use `@keywords internal` and `@noRd` (no man page
generated).

### Phase 4: enable roxygen-driven NAMESPACE
Once every exported function has `#' @export`, edit DESCRIPTION:

```
Roxygen: list(markdown = TRUE, roclets = c("collate", "namespace", "rd"))
```

Run `devtools::document()` — NAMESPACE is now auto-generated. Delete the
hand-written exports section but keep any `import(data.table)` / S3method
lines that don't come from roxygen, or move them to a dedicated file
(e.g. an `@importFrom data.table :=` tag somewhere).

### Phase 5: write the real getting-started vignette
The simulated dataset only needs to be small (~100 patients, ~7 days,
4-state). Once the vignette runs end-to-end on `devtools::build_vignettes()`,
you have a real package.

### Phase 6: tag v0.1.0
When the vignette runs, all tests pass on R CMD check with 0 errors / 0
warnings, and you have at least one published preprint / paper using it:

```r
usethis::use_version("minor")  # bumps to 0.1.0
```

Then on GitHub: Releases → Draft new release → tag v0.1.0. If you connect
the repo to Zenodo first, each release auto-generates a citable DOI.

## 5. Separate applied analysis repos

Keep downstream analyses in separate repos that install the package via
`remotes::install_github()`. This separation forces clean API decisions
because you become a user of your own package, not just its author.

## 6. Long-term: Julia port

### Why Julia

The R implementation has reached the practical ceiling of what fork-based
parallelism can safely provide. The nested `mclapply` scheme (fold × reg ×
sl workers) works well with one level at a time, but combining two or more
restricts learner choice — multi-threaded and GPU-based learners are unsafe in
grandchild processes after a fork, and R offers no clean alternative compatible
with the existing structure.

Julia removes these constraints entirely: true shared-memory `Threads.@spawn`
parallelism (no fork), CUDA works freely at any nesting depth, JIT-compiled to
native code, syntax close enough to R that the estimator logic ports directly.
Python is less suitable — the GIL limits true CPU-level threading.

### Ecosystem mapping

| R | Julia |
|---|-------|
| SuperLearner | MLJ.jl (`Stack`) |
| data.table | DataFrames.jl |
| glm / glmnet | GLM.jl / Lasso.jl |
| xgboost (GPU) | EvoTrees.jl (pure Julia, native CUDA) |
| earth / MARS | Approximated by EvoTrees stumps + Lasso |
| dbarts / BART | Dropped — limited value in practice; overfits or gets suppressed |
| mclapply (nested) | Not needed — Threads.@spawn shares memory natively |
| Neural nets | Flux.jl (native CUDA, fits in MLJ as a learner) |

BART is not a significant loss. In practice it either overfits the noisy
high-dimensional covariate space and gets suppressed by the metalearner, or
produces extreme out-of-fold predictions that contaminate the stack.

### Automated tuning — the key upgrade

In R, "tuning" means hand-crafting multiple wrappers with fixed hyperparameters
and letting the metalearner pick. In Julia, MLJ's `TunedModel` wraps any
learner and does CV-based hyperparameter search as a first-class model object
that sits directly in the stack:

```julia
xgb_tuned = TunedModel(
    model     = EvoTreeClassifier(),
    resampling = CV(nfolds = 5),
    tuning    = Grid(resolution = 8),
    range     = [range(model, :max_depth,  values = [3, 4, 5]),
                 range(model, :lambda,      values = [0.1, 0.5, 1.0, 5.0]),
                 range(model, :eta,         values = [0.05, 0.1, 0.3])]
)
```

Each stack component self-tunes before the metalearner combines them. The
tuning uses the same CV folds as the stacking — no data leakage.

### Learner strategy: asymmetric g / Q regularisation

**g-side (density ratio)**: moderate regularisation is fine. The WB DR
metalearner's log loss naturally suppresses extreme OOF predictions, and the
doubly-robust structure absorbs some g-misspecification. xgboost dominates in
practice; HAL/earth contribute little at the g-side.

**Q-side (SDR/iTMLE)**: regularisation-first. The backward pseudo-outcome
recursion is the sensitive path. Aggressive learners (deep trees, BART)
systematically over-predict at each recursive step; the bounding of
pseudo-outcomes provides partial self-correction but does not eliminate the
bias. The Q-side tuning grid should be constrained to regularised regimes —
shallow trees, high L2 penalty, conservative lambda — and should not include
unconstrained aggressive configurations.

The Julia `TunedModel` makes this structural: define separate tuning ranges for
g-side and Q-side stacks rather than hand-crafting conservative wrappers.

### Neural nets for Q-models

Untested in the SDR/iTMLE setting but theoretically motivated. Trees discretise
the covariate space; sharp leaf boundaries produce discontinuous OOF predictions
that overfit aggressively when pseudo-outcomes are noisy. Neural nets with smooth
activations (ReLU, GELU) and dropout produce continuous, smoother predictions.
Dropout is implicit model averaging over subnetworks — a natural regulariser
whose effect is gradual rather than discrete. Early stopping directly targets
held-out loss. Whether smoother OOF behaviour translates to better pseudo-outcome
recursion in practice is an open empirical question worth testing in Julia, where
Flux.jl + CUDA makes it practical in a way R does not.

### Implementation phases

1. **Core data utilities**: lag construction, fold helpers, data_check — ~1:1 port from sl.R / shared_helpers.R
2. **MLJ stacking machinery**: `sl_block_fit` equivalent with custom `validRows`; Wu-Benkeser DR metalearner as MLJ combiner; verify custom fold injection API
3. **Density ratio**: `build_stack`, cross-fitting loop with `Threads.@spawn` over time × fold freely
4. **SDR estimator**: transition regressions, EIF assembly, backward recursion
5. **iTMLE estimator**: targeting step — verify whether MLJ subsetting handles offset-as-column or whether the hack is still needed
6. **Competing events variant**: lowest priority, after main estimators validated
7. **Validation**: numerical agreement with R on simulated data via `RCall.jl`

### Key open question

MLJ's `Stack` uses its own CV strategy. Confirm early whether a fully custom
fold vector can be injected (equivalent to SuperLearner's `validRows`) or
whether a custom `ResamplingStrategy` must be implemented. This determines the
complexity of phase 2.
