# Evidence and readiness

DesiredGainR combines an operational quantitative-genetic core with
advanced recommendation, uncertainty and simulation tools. This document
shows the evidence supporting each layer and how it adds value to
breeding decisions as of version 0.5.0.

## Operationally ready core

Verified by closed-form derivation and Monte Carlo tests, and safe to
use for operational decisions.

| Component | Functions |
|----|----|
| Objective definition and feasibility | [`implied_economic_weights()`](https://FAkohoue.github.io/DesiredGainR/reference/implied_economic_weights.md), [`implied_desired_gains()`](https://FAkohoue.github.io/DesiredGainR/reference/implied_desired_gains.md), [`gain_feasibility()`](https://FAkohoue.github.io/DesiredGainR/reference/gain_feasibility.md), [`retrospective_weights()`](https://FAkohoue.github.io/DesiredGainR/reference/retrospective_weights.md), [`weight_sensitivity()`](https://FAkohoue.github.io/DesiredGainR/reference/weight_sensitivity.md), [`effective_weights()`](https://FAkohoue.github.io/DesiredGainR/reference/effective_weights.md) |
| Covariance diagnostics | [`matrix_diagnostics()`](https://FAkohoue.github.io/DesiredGainR/reference/matrix_diagnostics.md), [`bend_covariance()`](https://FAkohoue.github.io/DesiredGainR/reference/bend_covariance.md), [`import_covariance()`](https://FAkohoue.github.io/DesiredGainR/reference/import_covariance.md) |

[`gain_feasibility()`](https://FAkohoue.github.io/DesiredGainR/reference/gain_feasibility.md)
in particular answers the single-cycle question exactly and without
simulation. Prefer it to the simulation layer whenever it suffices.

## Validated index workflows

The mathematics is verified, the interfaces are settled, and these
functions provide reproducible candidate ranking, response prediction
and method comparison.

| Component | Functions |
|----|----|
| Classical index families | [`selection_index()`](https://FAkohoue.github.io/DesiredGainR/reference/selection_index.md), [`evaluate_index()`](https://FAkohoue.github.io/DesiredGainR/reference/evaluate_index.md), [`predict()`](https://rdrr.io/r/stats/predict.html), [`summary()`](https://rdrr.io/r/base/summary.html) |
| Restricted and proportional-gain indices | [`restricted_index()`](https://FAkohoue.github.io/DesiredGainR/reference/restricted_index.md) |
| Desired-gain index | [`run_dgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_dgsi.md) |
| Quadratic genomic index | [`run_qgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_qgsi.md), [`compare_dg_and_qgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/compare_dg_and_qgsi.md) |

Version 0.5.0 strengthened
[`run_dgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_dgsi.md)
with explicit desired-gain unit handling when traits have different
scales. Analyses created with earlier versions should be rerun to
benefit from the verified implementation.

## Advanced uncertainty and environment analysis

These functions quantify how covariance estimation and
target-environment structure affect a decision, making sensitivity
visible rather than implicit.

| Component | Functions |
|----|----|
| Covariance uncertainty | [`index_uncertainty()`](https://FAkohoue.github.io/DesiredGainR/reference/index_uncertainty.md), [`draw_covariance_pairs()`](https://FAkohoue.github.io/DesiredGainR/reference/draw_covariance_pairs.md), [`propagate_covariance_uncertainty()`](https://FAkohoue.github.io/DesiredGainR/reference/propagate_covariance_uncertainty.md) |
| Multi-environment expansion | [`expand_environments()`](https://FAkohoue.github.io/DesiredGainR/reference/expand_environments.md), [`widen_environments()`](https://FAkohoue.github.io/DesiredGainR/reference/widen_environments.md) |

The covariance-uncertainty functions use a Wishart approximation to a
sampling distribution that has no closed form, and `genetic_df` is the
user’s judgement. It is informative about the *order* of the
uncertainty, not its exact magnitude. The fixed-`P` mode is a local
sensitivity analysis, not a generative model.

## Advanced recommendation and simulation tools

These tools support structured scenario comparison and breeder decisions
when their stated population, covariance and simulation assumptions
match the programme.

| Component | Functions | Why |
|----|----|----|
| Recurrent simulation | [`founder_population()`](https://FAkohoue.github.io/DesiredGainR/reference/founder_population.md), [`simulate_selection_cycles()`](https://FAkohoue.github.io/DesiredGainR/reference/simulate_selection_cycles.md) | Credible for comparing directions under its stated assumptions. The clonal mode calibrates the full target genotypic covariance and reports the realised matrix and calibration error. |
| Direction optimisation | [`optimize_desired_gains()`](https://FAkohoue.github.io/DesiredGainR/reference/optimize_desired_gains.md) | Replicate-level covariance and shared-seed dependence reach the surrogate diagnostics, but the primary search remains conditional on one covariance estimate. Use [`propagate_covariance_uncertainty()`](https://FAkohoue.github.io/DesiredGainR/reference/propagate_covariance_uncertainty.md) on the candidate frontier. |
| Population-driven recommendation | [`suggest_desired_gains()`](https://FAkohoue.github.io/DesiredGainR/reference/suggest_desired_gains.md) | Its one-cycle feasibility and maximin limits are exact conditional on `G_target` and `P`; the production QP is independently checked against Clarabel and OSQP. Discovery, screening, multiple-finalist confirmation, multi-start convergence, model calibration, architecture draws, optional covariance draws and decision gates make recommendations transparent and defensible. |
| Diversity optimisation | `optimize_desired_gains(include_diversity = TRUE)` | Defaults off. It requires a marker panel disjoint from QTL; unknown overlap needs an explicit experimental override and known overlap is rejected. |

The empirical suite under `inst/validation/` covers recurrent-cycle
trends, an index frozen before later INIA rice cohorts, selected-only
later-stage outcomes, a Genomes-to-Fields time split, the published
CIMMYT wheat anchor and the Zhang et al. maize RCGS programme. It
reproduces published genetic gain, recovers raw genomic prediction,
verifies phenotype-to-genotype linkage and tests desired-gain directions
in held-out environments. A prospective vector-comparison experiment is
identified as the next evidence extension for estimating incremental
field response among strategies.

## How DesiredGainR fits the breeding programme

DesiredGainR connects a fitted multi-trait genetic evaluation to
candidate selection and an auditable objective. It complements
trial-analysis, mating design and optimal-contribution tools, allowing
each programme stage to use the method best suited to it.

The intended position is:

    breeder intent
      -> objective feasibility
      -> comparable index formulations
      -> candidate ranking
      -> uncertainty and sensitivity
      -> auditable decision record

Use AlphaSimR directly for programme design, `selection.index` for a
second opinion on classical coefficients, and a relationship-management
tool for inbreeding control.
