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

The distinctive layer of the package. Economic weights and desired gains
are two parameterisations of one linear index, so each implies the other
exactly; those translations, a feasibility test derived from the
achievable-response ellipsoid, recovery of the objective implicit in
past selection decisions, and a sensitivity analysis reporting whether
the decision depends on the stated weights at all.

- [`implied_economic_weights()`](https://FAkohoue.github.io/DesiredGainR/reference/implied_economic_weights.md)
  : Translate desired gains into the economic weights they imply
- [`implied_desired_gains()`](https://FAkohoue.github.io/DesiredGainR/reference/implied_desired_gains.md)
  : Translate economic weights into the desired gains they imply
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

## Building a selection index

The classical index families through one interface, so that alternative
objectives and alternative families can be compared on identical data,
together with the iterative desired-gain index and the quadratic genomic
selection index.

- [`selection_index()`](https://FAkohoue.github.io/DesiredGainR/reference/selection_index.md)
  : Construct a classical multi-trait selection index
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
the index, alongside the expected response for every trait.

- [`evaluate_index()`](https://FAkohoue.github.io/DesiredGainR/reference/evaluate_index.md)
  : Evaluate a selection index against the standard criteria
- [`compare_dg_and_qgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/compare_dg_and_qgsi.md)
  : Compare DGSI and QGSI candidate rankings

## Multi-cycle simulation

Comparing desired-gain directions over several breeding cycles, because
a direction that maximises response in the first cycle can exhaust the
genetic variation that response depends on within a few more. Founders
are built from the breeder’s own phased marker data, so the linkage
disequilibrium and allele-frequency structure of the simulation are
those of the target programme. The search treats the simulation as the
objective function and never substitutes a cheaper approximation for it.

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
