# Compare selection methods on one common decision problem

This function accepts results from
[`selection_index()`](https://FAkohoue.github.io/DesiredGainR/reference/selection_index.md),
[`restricted_index()`](https://FAkohoue.github.io/DesiredGainR/reference/restricted_index.md),
[`generalized_index()`](https://FAkohoue.github.io/DesiredGainR/reference/generalized_index.md),
[`run_dgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_dgsi.md),
and
[`run_qgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_qgsi.md).
Internal adapters put response vectors in favourable original trait
units. A
[`comparison_objective()`](https://FAkohoue.github.io/DesiredGainR/reference/comparison_objective.md)
can then place them in trait or genetic-standard- deviation units and
evaluate every method against one fixed target and utility.

## Usage

``` r
compare_selection_methods(
  models,
  target_gains = NULL,
  objective = NULL,
  validation_data = NULL
)
```

## Arguments

- models:

  Named list containing at least two supported fitted objects.

- target_gains:

  Optional legacy target on the common fitted analysis scale.

- objective:

  Optional object returned by
  [`comparison_objective()`](https://FAkohoue.github.io/DesiredGainR/reference/comparison_objective.md).

- validation_data:

  Optional data frame containing one common set of objective-trait
  values for every candidate. One non-trait identifier column, or the
  row names, must match the fitted candidate identifiers.

## Value

An object of class `desiredgainr_method_comparison`. It contains method
summaries, model-based responses, optional validation responses,
fairness checks, rank correlations, and selected-set agreement.

## Details

Coefficients are never compared across families. Candidate ranks,
selected sets, per-trait responses, target alignment, and common utility
are compared. QGSI expected gains retain their documented
linear-regression approximation when `W` contains curvature. The `W = 0`
special case has the exact linear response. A quadratic common utility
is evaluated exactly only from `validation_data`. A singular common `G`
still permits response and decision comparisons. Mahalanobis alignment
and Satoh projection remain missing because those criteria require an
inverse.

`target_gains` is retained for backward compatibility. It uses the
analysis scale of the first fitted `desiredgainr_index` and requires
every model to have the same transformation. New cross-family analyses
should use `objective`.

## See also

[`comparison_objective()`](https://FAkohoue.github.io/DesiredGainR/reference/comparison_objective.md),
[`compare_dg_and_qgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/compare_dg_and_qgsi.md),
[`evaluate_index()`](https://FAkohoue.github.io/DesiredGainR/reference/evaluate_index.md),
[`candidate_score_uncertainty()`](https://FAkohoue.github.io/DesiredGainR/reference/candidate_score_uncertainty.md)
