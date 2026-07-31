# Evaluate a selection index against the standard criteria

A selection index cannot be judged from its coefficients. Rahimi and
Debnath (2023) evaluate every index on four criteria, and this function
reproduces them so that alternative index families and alternative
objectives can be compared on identical data.

## Usage

``` r
evaluate_index(
  coefficients,
  G,
  P,
  aggregate_weights = NULL,
  selection_intensity = NA_real_,
  main_trait = names(coefficients)[1L],
  scores = NULL
)
```

## Arguments

- coefficients:

  Named numeric vector of index coefficients, in the same trait space as
  `G` and `P`.

- G:

  Genetic variance-covariance matrix, named by trait.

- P:

  Phenotypic variance-covariance matrix, named by trait.

- aggregate_weights:

  Named weights defining net merit. When `NULL`, the criteria that
  depend on net merit are returned as `NA`.

- selection_intensity:

  Standardised selection intensity.

- main_trait:

  Trait against which relative efficiency is computed.

- scores:

  Optional vector of realised index scores, used for \\CV_I\\.

## Value

An object of class `desiredgainr_evaluation`.

## Details

Let \\b\\ be the index coefficients, \\a\\ the aggregate weights
defining net merit \\H = a^\mathsf{T}g\\, and \\k\\ the standardised
selection intensity. The criteria are:

- \\R\_{HI}\\:

  The correlation between the index and net merit, \\b^\mathsf{T}Ga /
  \sqrt{b^\mathsf{T}Pb \cdot a^\mathsf{T}Ga}\\. Maximising it maximises
  the response in aggregate merit.

- \\\Delta H\\:

  The expected gain in aggregate merit, \\k R\_{HI}
  \sqrt{a^\mathsf{T}Ga}\\.

- \\\Delta_j\\:

  The expected response for each trait, \\k (Gb)\_j /
  \sqrt{b^\mathsf{T}Pb}\\.

- RE:

  Efficiency relative to direct selection on a single main trait.

- \\CV_I\\:

  The coefficient of variation of the index scores.

## Interpreting relative efficiency

Relative efficiency below 1 is not a defect. Every value reported by
Rahimi and Debnath was below 1, because a multi-trait index deliberately
trades response in the main trait for response in the remainder of the
objective. Hence RE should be read alongside the per-trait responses
that direct selection would sacrifice, and not as a pass-or-fail test.

## References

Rahimi M, Debnath S (2023). *Scientific Reports* 13:18977.
[doi:10.1038/s41598-023-46368-6](https://doi.org/10.1038/s41598-023-46368-6)

## See also

[`selection_index()`](https://FAkohoue.github.io/DesiredGainR/reference/selection_index.md)
