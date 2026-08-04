# Compare multi-trait selection methods on one decision problem

A fair comparison uses the same candidates, trait orientation,
measurement scale, selection intensity, and definition of aggregate
merit. This function checks those conditions. It then places biological
response, objective attainment, ranking agreement, and selected-set
agreement in one result.

## Usage

``` r
compare_selection_methods(models, target_gains = NULL)
```

## Arguments

- models:

  Named list containing at least two fitted `desiredgainr_index`
  objects.

- target_gains:

  Optional named vector of favourable desired gains on the common fitted
  analysis scale.

## Value

An object of class `desiredgainr_method_comparison`. The `summary`
component contains one row per method. The `responses` component
contains one row per method and trait. The three matrices are
`rank_correlation`, `selected_jaccard`, and `selected_overlap`. The
`fairness` component records the conditions required for interpretation.

## Details

The function compares fitted objects returned by
[`selection_index()`](https://FAkohoue.github.io/DesiredGainR/reference/selection_index.md).
It keeps expected response separate from the observed differential among
the supplied candidates. Expected response follows selection-index
theory. Observed differential describes the present candidate set.

When `target_gains` is supplied, it must use the analysis scale stored
in the fitted objects. For example, a model fitted to
genetic-standard-deviation values requires targets in genetic standard
deviations. Positive values mean improvement because trait direction has
already been applied by
[`selection_index()`](https://FAkohoue.github.io/DesiredGainR/reference/selection_index.md).

Methods based on ranks or thresholds lack a closed-form expected
response in this implementation. Their expected-response fields remain
missing. Their observed differentials, ranks, and selected sets remain
comparable.

## References

Beavis W, Lamkey K, Mahama AA, Suza W (2023). Multiple Trait Selection.
In Suza WP and Lamkey KR, editors, *Quantitative Genetics for Plant
Breeding*. Iowa State University Digital Press.

Cunningham EP (1969). The relative efficiencies of selection indexes.
*Acta Agriculturae Scandinavica* 19:45-48.

## See also

[`selection_index()`](https://FAkohoue.github.io/DesiredGainR/reference/selection_index.md),
[`evaluate_index()`](https://FAkohoue.github.io/DesiredGainR/reference/evaluate_index.md),
[`candidate_score_uncertainty()`](https://FAkohoue.github.io/DesiredGainR/reference/candidate_score_uncertainty.md),
[`index_information_efficiency()`](https://FAkohoue.github.io/DesiredGainR/reference/index_information_efficiency.md)

## Examples

``` r
traits <- c("yield", "disease")
values <- data.frame(
  id = paste0("L", 1:30),
  yield = seq(-1.5, 1.5, length.out = 30),
  disease = seq(1.5, -1.5, length.out = 30)
)
G <- matrix(c(1, -0.2, -0.2, 0.8), 2,
  dimnames = list(traits, traits)
)
P <- matrix(c(1.8, -0.2, -0.2, 1.4), 2,
  dimnames = list(traits, traits)
)
objective <- c(yield = 1, disease = 0.5)
fits <- list(
  Smith_Hazel = selection_index(
    values, traits, id_col = "id", method = "smith_hazel",
    G = G, P = P, economic_weights = objective,
    lower_is_better = "disease", scale_traits = FALSE, n_select = 6
  ),
  Base = selection_index(
    values, traits, id_col = "id", method = "base",
    G = G, P = P, economic_weights = objective,
    lower_is_better = "disease", scale_traits = FALSE, n_select = 6
  )
)
compare_selection_methods(fits)
#> <desiredgainr_method_comparison>
#>   Methods: Smith_Hazel, Base 
#>   All recorded comparison conditions are satisfied.
#>       Method N_selected      R_HI  Delta_H        RE Worst_expected_attainment
#>  Smith_Hazel          6 0.7736293 1.281344 0.9408723                        NA
#>         Base          6 0.7718450 1.278389 0.9627089                        NA
#>  Mahalanobis_alignment
#>                     NA
#>                     NA
```
