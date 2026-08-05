# Package index

## Package overview

Main documentation entry point for DesiredGainR, and the standalone
Breeder’s Guide covering the objective-setting decision and the choice
among index families in narrative form, for readers who approve a
selection decision without running R themselves.

- [`DesiredGainR`](https://FAkohoue.github.io/DesiredGainR/reference/DesiredGainR-package.md)
  [`DesiredGainR-package`](https://FAkohoue.github.io/DesiredGainR/reference/DesiredGainR-package.md)
  : DesiredGainR: Auditable Multi-trait Selection Indices
- [`open_desiredgain_guide()`](https://FAkohoue.github.io/DesiredGainR/reference/open_desiredgain_guide.md)
  : Locate or Open the DesiredGainR Breeder's Guide

## Defining the breeding objective

The distinctive layer of the package. Economic weights define aggregate
merit. Desired gains define a response direction. An exact algebraic map
connects these statements under the stated covariance model, but it does
not give them the same biological meaning. This section covers that map,
a feasibility test derived from the achievable-response ellipsoid,
recovery of the objective implicit in past selection decisions, and a
sensitivity analysis reporting whether the decision depends on the
stated weights at all.

- [`implied_economic_weights()`](https://FAkohoue.github.io/DesiredGainR/reference/implied_economic_weights.md)
  : Map desired gains to the aggregate weights they imply
- [`implied_desired_gains()`](https://FAkohoue.github.io/DesiredGainR/reference/implied_desired_gains.md)
  : Map economic weights to the desired gains they imply
- [`define_desired_gain_intervals()`](https://FAkohoue.github.io/DesiredGainR/reference/define_desired_gain_intervals.md)
  : Define acceptable intervals for desired genetic gains
- [`suggest_desired_gains()`](https://FAkohoue.github.io/DesiredGainR/reference/suggest_desired_gains.md)
  : Suggest desired-gain directions from the current population
- [`gain_feasibility()`](https://FAkohoue.github.io/DesiredGainR/reference/gain_feasibility.md)
  : Test whether a desired-gain vector is attainable
- [`retrospective_weights()`](https://FAkohoue.github.io/DesiredGainR/reference/retrospective_weights.md)
  : Recover the weights implied by past selection decisions
- [`weight_sensitivity()`](https://FAkohoue.github.io/DesiredGainR/reference/weight_sensitivity.md)
  : Assess how sensitive a selection decision is to the stated objective
- [`effective_weights()`](https://FAkohoue.github.io/DesiredGainR/reference/effective_weights.md)
  : Report the effective contribution of each trait to an index

## Obtaining covariance matrices

Labelled approximations to the genetic covariance matrix for programmes
without a fitted multi-trait model, including the covariance of
across-environment adjusted means recommended for operational use by the
CGIAR implementation guideline. Each method records its estimand,
provenance and assumptions rather than relabelling a covariance of
predictions as a genetic covariance.

- [`estimate_genetic_covariance()`](https://FAkohoue.github.io/DesiredGainR/reference/estimate_genetic_covariance.md)
  : Estimate a working genetic covariance matrix
- [`matrix_diagnostics()`](https://FAkohoue.github.io/DesiredGainR/reference/matrix_diagnostics.md)
  : Report the numerical conditioning of a covariance matrix
- [`bend_covariance()`](https://FAkohoue.github.io/DesiredGainR/reference/bend_covariance.md)
  : Bend a covariance matrix to positive definiteness

## Building a selection index

The classical and restricted index families, the iterative search
applied to the desired-gain index, and the quadratic genomic selection
index. A common comparison objective and internal result adapters allow
these families to be evaluated on the same candidates, trait units,
target, utility and validation values.

- [`selection_index()`](https://FAkohoue.github.io/DesiredGainR/reference/selection_index.md)
  : Apply a classical multi-trait selection method
- [`selection_information()`](https://FAkohoue.github.io/DesiredGainR/reference/selection_information.md)
  : Define the information used by a selection index
- [`generalized_index()`](https://FAkohoue.github.io/DesiredGainR/reference/generalized_index.md)
  : Fit an index from a general information model
- [`restricted_index()`](https://FAkohoue.github.io/DesiredGainR/reference/restricted_index.md)
  : Restricted and proportional-gain selection indices
- [`restricted_breeding_values()`](https://FAkohoue.github.io/DesiredGainR/reference/restricted_breeding_values.md)
  : Project breeding values into a restricted response space
- [`run_dgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_dgsi.md)
  : Fit an optimised desired-gain selection index
- [`run_qgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_qgsi.md)
  [`run_qgsi_desired_gain()`](https://FAkohoue.github.io/DesiredGainR/reference/run_qgsi.md)
  : Fit and evaluate a quadratic genomic selection index
- [`run_dgsi_qgsi_pipeline()`](https://FAkohoue.github.io/DesiredGainR/reference/run_dgsi_qgsi_pipeline.md)
  : Run DGSI, QGSI, or both workflows

## Evaluating and comparing indices

The four criteria used by Rahimi and Debnath (2023) and by the
established selection-index software: the correlation between the index
and net merit, the expected gain in aggregate merit, efficiency relative
to direct selection on a main trait, and the coefficient of variation of
the index, alongside the expected response for every trait and the
heritability of the index itself. Also the tools for reusing a fitted
index on the next cycle and for asking how much of the answer is
determined by the data rather than by sampling error in the covariance
matrices.

- [`evaluate_index()`](https://FAkohoue.github.io/DesiredGainR/reference/evaluate_index.md)
  : Evaluate a selection index against the standard criteria
- [`compare_selection_methods()`](https://FAkohoue.github.io/DesiredGainR/reference/compare_selection_methods.md)
  : Compare selection methods on one common decision problem
- [`comparison_objective()`](https://FAkohoue.github.io/DesiredGainR/reference/comparison_objective.md)
  : Define one objective for a cross-family comparison
- [`summary(`*`<desiredgainr_index>`*`)`](https://FAkohoue.github.io/DesiredGainR/reference/summary.desiredgainr_index.md)
  : Summarise a fitted selection index
- [`predict(`*`<desiredgainr_index>`*`)`](https://FAkohoue.github.io/DesiredGainR/reference/predict.desiredgainr_index.md)
  : Score a new candidate set with a fitted index
- [`index_uncertainty()`](https://FAkohoue.github.io/DesiredGainR/reference/index_uncertainty.md)
  : Propagate sampling error in the covariance matrices into an index
- [`candidate_score_uncertainty()`](https://FAkohoue.github.io/DesiredGainR/reference/candidate_score_uncertainty.md)
  : Candidate-specific uncertainty in index scores and ranks
- [`index_information_efficiency()`](https://FAkohoue.github.io/DesiredGainR/reference/index_information_efficiency.md)
  : Cunningham information-deletion efficiency
- [`evaluate_restricted_response()`](https://FAkohoue.github.io/DesiredGainR/reference/evaluate_restricted_response.md)
  : Evaluate an achieved response against a desired-gain direction
- [`compare_dg_and_qgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/compare_dg_and_qgsi.md)
  : Compare an iteratively optimised DGSI with a QGSI

## Multi-environment structure

A multi-environment index is a single-environment index on an expanded
trait set, so the whole index layer applies once the separable
covariance has been built. The interaction share reported by the
expansion says whether modelling the environments separately was worth
the complexity.

- [`expand_environments()`](https://FAkohoue.github.io/DesiredGainR/reference/expand_environments.md)
  : Expand a trait-by-environment structure into an index problem
- [`widen_environments()`](https://FAkohoue.github.io/DesiredGainR/reference/widen_environments.md)
  : Reshape long-format multi-environment trial data for indexing

## Uncertainty in the recommendation

Every reported gain conditions on covariance matrices that are
themselves estimates. These functions report how far the answer moves
when that is taken into account, and separate covariance-estimation
uncertainty from Monte Carlo simulation error, because only the second
falls with more computing.

- [`index_uncertainty()`](https://FAkohoue.github.io/DesiredGainR/reference/index_uncertainty.md)
  : Propagate sampling error in the covariance matrices into an index
- [`draw_covariance_pairs()`](https://FAkohoue.github.io/DesiredGainR/reference/draw_covariance_pairs.md)
  : Draw jointly admissible genetic and residual covariance matrices
- [`propagate_covariance_uncertainty()`](https://FAkohoue.github.io/DesiredGainR/reference/propagate_covariance_uncertainty.md)
  : Propagate covariance uncertainty into a desired-gain recommendation

## Published reference values

- [`rahimi_debnath_2023()`](https://FAkohoue.github.io/DesiredGainR/reference/rahimi_debnath_2023.md)
  : Published selection-index results from Rahimi and Debnath (2023)

## Importing from upstream models

Adapters for the mixed-model packages that produce the covariance
matrices this package consumes. They reorder traits to the analysis
order rather than assuming it, distinguish a correlation matrix from a
covariance, and check the genetic and phenotypic matrices against each
other, because a silently transposed trait order produces an index that
is wrong in a way nothing downstream can detect.

- [`import_covariance()`](https://FAkohoue.github.io/DesiredGainR/reference/import_covariance.md)
  : Validate and record an extracted genetic covariance matrix

## Multi-cycle simulation

Comparing desired-gain directions over several breeding cycles, because
a direction that maximises response in the first cycle can exhaust the
genetic variation that response depends on within a few more. Founders
are built from the breeder’s own phased marker data, so the linkage
disequilibrium and allele-frequency structure of the simulation are
those of the target programme. The search treats the simulation as the
objective function and never substitutes a cheaper approximation for it.

- [`breeding_scenario()`](https://FAkohoue.github.io/DesiredGainR/reference/breeding_scenario.md)
  : Describe a breeding simulation scenario
- [`simulation_calibration()`](https://FAkohoue.github.io/DesiredGainR/reference/simulation_calibration.md)
  : Audit the calibration of an AlphaSimR setup
- [`stress_test_desired_gains()`](https://FAkohoue.github.io/DesiredGainR/reference/stress_test_desired_gains.md)
  : Stress-test a desired-gain vector across simulation scenarios
- [`founder_haplotypes()`](https://FAkohoue.github.io/DesiredGainR/reference/founder_haplotypes.md)
  : Validate and package phased founder haplotypes
- [`haplotypes_from_inbred_dosage()`](https://FAkohoue.github.io/DesiredGainR/reference/haplotypes_from_inbred_dosage.md)
  : Derive phased haplotypes from inbred-line dosage
- [`dosage_diagnostics()`](https://FAkohoue.github.io/DesiredGainR/reference/dosage_diagnostics.md)
  : Summarise heterozygosity and missingness in a dosage matrix
- [`founder_population()`](https://FAkohoue.github.io/DesiredGainR/reference/founder_population.md)
  : Build an AlphaSimR founder population with a target genetic
  covariance
- [`simulate_selection_cycles()`](https://FAkohoue.github.io/DesiredGainR/reference/simulate_selection_cycles.md)
  : Simulate recurrent selection under a fixed desired-gain direction
- [`optimize_desired_gains()`](https://FAkohoue.github.io/DesiredGainR/reference/optimize_desired_gains.md)
  : Search for the desired-gain direction giving the best multi-cycle
  outcome

## Example datasets

One simulated tropical maize programme with six traits, carrying the
antagonistic genetic correlations, the two-order-of-magnitude spread of
trait scales, and the internal consistency between markers and trait
values that a multi-trait selection index needs in order to demonstrate
anything.

- [`dgr_traits`](https://FAkohoue.github.io/DesiredGainR/reference/DesiredGainR-data.md)
  [`dgr_G`](https://FAkohoue.github.io/DesiredGainR/reference/DesiredGainR-data.md)
  [`dgr_P`](https://FAkohoue.github.io/DesiredGainR/reference/DesiredGainR-data.md)
  [`dgr_candidates`](https://FAkohoue.github.io/DesiredGainR/reference/DesiredGainR-data.md)
  [`dgr_gebv`](https://FAkohoue.github.io/DesiredGainR/reference/DesiredGainR-data.md)
  [`dgr_history`](https://FAkohoue.github.io/DesiredGainR/reference/DesiredGainR-data.md)
  [`dgr_hap1`](https://FAkohoue.github.io/DesiredGainR/reference/DesiredGainR-data.md)
  [`dgr_hap2`](https://FAkohoue.github.io/DesiredGainR/reference/DesiredGainR-data.md)
  [`dgr_map`](https://FAkohoue.github.io/DesiredGainR/reference/DesiredGainR-data.md)
  : Example breeding programme shipped with DesiredGainR
