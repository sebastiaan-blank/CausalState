<!-- README.md is generated from README.Rmd. Please edit that file -->

# CausalState

<!-- badges: start -->
[![R-CMD-check](https://github.com/YOUR-GITHUB-HANDLE/CausalState/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/YOUR-GITHUB-HANDLE/CausalState/actions/workflows/R-CMD-check.yaml)
[![License: AGPL v3](https://img.shields.io/badge/License-AGPL_v3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

**CausalState** provides Sequential Doubly Robust (SDR) and
infinite-dimensional Targeted Maximum Likelihood (iTMLE) estimators for
longitudinal modified treatment policies in settings with absorbing state
transitions — for example, ICU, ward, or emergency department care
episodes where treatment is applicable while in the active state and
becomes structurally inapplicable after discharge or death.

The package complements [`lmtp`](https://github.com/nt-williams/lmtp) by
handling settings where the timing of state transitions is itself
endogenous to the treatment policy.

> **Note:** This README is a placeholder. Edit `README.Rmd` and run
> `devtools::build_readme()` to regenerate.
