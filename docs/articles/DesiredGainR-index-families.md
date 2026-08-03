# Choosing an index

## 1. Begin with the objective, not the method

An index family is appropriate only relative to a stated breeding
objective. Economic value, a desired response, hard eligibility limits
and non-linear value judgements are different objectives. They should
not be made to look equivalent by applying the same function name to all
of them.

| Breeder’s statement | Appropriate starting point |
|----|----|
| “I can defend the relative economic value of each trait.” | Smith–Hazel |
| “The measurements themselves already form the score.” | Base index |
| “I can defend a desired response direction.” | Pesek–Baker or Yamada |
| “I know acceptable gain intervals, but not one exact vector.” | Population-driven desired-gain suggestions |
| “Every selected candidate must pass these limits.” | Independent culling or a restricted index |
| “I need a simple, weight-light comparator.” | Mulamba–Mock rank sum |
| “Trait value is genuinely curved or interactive.” | Quadratic genomic selection index |

Tandem selection is included as a comparator, but it postpones the other
traits rather than balancing them in the same decision.

------------------------------------------------------------------------

## 2. Fit comparable linear families

All traits are oriented so that larger values are favourable. The
objective below is expressed in genetic standard-deviation units.

``` r
desired_gains <- c(
  GY = 1.0, PHT = 0.4, AD = 0.6,
  ASI = 0.5, EPP = 0.4, GLS = 0.6
)
economic_weights <- c(
  GY = 1.0, PHT = 0.2, AD = 0.5,
  ASI = 0.4, EPP = 0.3, GLS = 0.5
)
```

The two vectors are separate statements. `desired_gains` defines a
response direction, whereas `economic_weights` defines one common
aggregate merit for the merit-based comparison below.

``` r
smith_hazel <- selection_index(
  dgr_candidates, traits,
  method = "smith_hazel", G = dgr_G, P = dgr_P,
  economic_weights = economic_weights,
  lower_is_better = lower_is_better,
  scale_traits = TRUE, scale_by = "phenotypic",
  n_select = 20L, main_trait = "GY"
)

pesek_baker <- selection_index(
  dgr_candidates, traits,
  method = "pesek_baker", G = dgr_G, P = dgr_P,
  desired_gains = desired_gains,
  aggregate_weights = economic_weights,
  lower_is_better = lower_is_better,
  scale_traits = TRUE, scale_by = "phenotypic",
  n_select = 20L, main_trait = "GY"
)

base_index <- selection_index(
  dgr_candidates, traits,
  method = "base", G = dgr_G, P = dgr_P,
  economic_weights = economic_weights,
  lower_is_better = lower_is_better,
  scale_traits = TRUE, scale_by = "phenotypic",
  n_select = 20L, main_trait = "GY"
)

rank_sum <- selection_index(
  dgr_candidates, traits,
  method = "mulamba_mock",
  G = dgr_G, P = dgr_P,
  lower_is_better = lower_is_better,
  scale_traits = TRUE, scale_by = "phenotypic",
  n_select = 20L, main_trait = "GY"
)
```

Smith–Hazel and Pesek–Baker are compared here under one common
definition of aggregate merit. Without that common definition, their
index–merit correlation and aggregate-merit response would not answer
the same question.

------------------------------------------------------------------------

## 3. Compare biological responses before summary criteria

``` r
response_table <- data.frame(
  Trait = traits,
  Smith_Hazel = unname(smith_hazel$evaluation$expected_response[traits]),
  Pesek_Baker = unname(pesek_baker$evaluation$expected_response[traits]),
  Base = unname(base_index$evaluation$expected_response[traits])
)
response_table[-1] <- round(response_table[-1], 3)
response_table
#>   Trait Smith_Hazel Pesek_Baker  Base
#> 1    GY       0.489       0.372 0.532
#> 2   PHT       0.009       0.149 0.033
#> 3    AD       0.345       0.223 0.280
#> 4   ASI       0.523       0.186 0.495
#> 5   EPP       0.445       0.149 0.413
#> 6   GLS       0.485       0.223 0.435
```

This table is the primary comparison. It shows what is expected to
change in every trait. A high aggregate statistic cannot rescue an index
that moves a critical trait in the wrong direction or misses its minimum
requirement.

Then inspect the summary criteria.

``` r
criteria <- function(x) {
  with(x$evaluation, c(
    index_merit_correlation = R_HI,
    aggregate_merit_response = delta_H,
    relative_efficiency_GY = RE,
    index_heritability = h2_index
  ))
}
round(rbind(
  Smith_Hazel = criteria(smith_hazel),
  Pesek_Baker = criteria(pesek_baker),
  Base = criteria(base_index)
), 3)
#>             index_merit_correlation aggregate_merit_response
#> Smith_Hazel                   0.706                    1.249
#> Pesek_Baker                   0.420                    0.743
#> Base                          0.689                    1.218
#>             relative_efficiency_GY index_heritability
#> Smith_Hazel                  0.796              0.510
#> Pesek_Baker                  0.605              0.251
#> Base                         0.867              0.474
```

Interpret these criteria with their limits: the index–merit correlation
and aggregate-merit response require one common merit definition;
relative efficiency describes only the declared main trait; index
heritability describes repeatability of the score, not attainment of the
complete breeding objective.

------------------------------------------------------------------------

## 4. Compare the actual decision

Two indices can have similar expected responses while selecting
different parents. Compare ranks and selected-set overlap, then repeat
under covariance or weight uncertainty.

``` r
selected_ids <- function(x) {
  x$ranking$id[x$ranking$selected]
}

sets <- list(
  Smith_Hazel = selected_ids(smith_hazel),
  Pesek_Baker = selected_ids(pesek_baker),
  Base = selected_ids(base_index),
  Rank_sum = selected_ids(rank_sum)
)

length(Reduce(intersect, sets))
#> [1] 9
```

Use
[`index_uncertainty()`](https://FAkohoue.github.io/DesiredGainR/reference/index_uncertainty.md)
for covariance resampling and
[`weight_sensitivity()`](https://FAkohoue.github.io/DesiredGainR/reference/weight_sensitivity.md)
when economic weights are uncertain. For minimum gain floors, inspect
the worst-trait and joint-attainment probabilities supplied by the
population-driven tools.

------------------------------------------------------------------------

## 5. Where DGSI and QGSI belong

The iterative desired-gain selection index (DGSI) searches for a
selected-group differential close to a desired vector. Read [Optimising
desired
gains](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-dgsi.md).

The quadratic genomic selection index (QGSI) ranks GEBVs under explicit
linear, squared and cross-product economic values. It does not accept
desired gains as quadratic weights. Read [Quadratic genomic
selection](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-qgsi.md).

Choose the simplest family that represents the programme’s real
objective and changes the decision in a stable, biologically defensible
way.
