# Evaluate a selection index against the standard criteria

A selection index cannot be judged from its coefficients. This function
reports classical summaries and the expected response of every trait.
However, the summaries answer different questions and are comparable
only under their stated conditions.

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

  The correlation between the selection index \\I = b^\mathsf{T}x\\ and
  aggregate genetic merit \\H = a^\mathsf{T}g\\, \\b^\mathsf{T}Ga /
  \sqrt{b^\mathsf{T}Pb \cdot a^\mathsf{T}Ga}\\. It is bounded from -1
  to 1. A negative value selects against the stated merit, zero
  indicates no linear association, and a value nearer one indicates
  closer agreement. Compare indices only when the same \\a\\ defines
  merit.

- \\\Delta H\\:

  The expected gain in aggregate merit, \\k R\_{HI}
  \sqrt{a^\mathsf{T}Ga}\\. The result is named `delta_H` and may appear
  as `dH` in compact print output. Positive values improve the stated
  merit. Compare magnitudes only under the same merit definition,
  population and selection intensity.

- \\\Delta_j\\:

  The expected response for each trait, \\k (Gb)\_j /
  \sqrt{b^\mathsf{T}Pb}\\. Its sign, magnitude and units provide the
  primary biological interpretation of an index.

- RE:

  Relative efficiency: response in the declared main trait divided by
  its response under direct selection at the same intensity. RE equal to
  one matches direct response, 0.8 retains 80 percent, a value above one
  exceeds direct response through correlated information, and a negative
  value moves the main trait unfavourably. It says nothing about
  response in the other traits.

- \\CV_I\\:

  Coefficient of variation of the index scores, \\100s_I/\|\bar I\|\\.
  It depends on the arbitrary zero of the score and is undefined for a
  centred index. It is a descriptive legacy statistic, not a robust
  criterion for ranking indices.

- \\h^2_I\\:

  The heritability of the index treated as a trait, \\b^\mathsf{T}Gb /
  b^\mathsf{T}Pb\\, whose square root is the correlation between the
  score and its own additive genetic component. It is bounded from zero
  to one for compatible covariance matrices and does not measure
  agreement with a desired-gain objective.

## Criteria required for a robust desired-gain recommendation

No scalar above is sufficient on its own. A defensible comparison also
requires: (i) every per-trait response, (ii) exact-ray feasibility or
the worst-trait and joint-attainment probabilities for minimum floors,
(iii) response-direction alignment under covariance uncertainty from
[`index_uncertainty()`](https://FAkohoue.github.io/DesiredGainR/reference/index_uncertainty.md),
(iv) rank and selected-set stability under perturbation, and (v)
diversity or coancestry diagnostics for repeated-cycle use. These
criteria separate biological attainment from statistical stability.

## Interpreting relative efficiency

Relative efficiency below 1 is not a defect. Every value reported by
Rahimi and Debnath was below 1, because a multi-trait index deliberately
trades response in the main trait for response in the remainder of the
objective. Hence RE should be read alongside the per-trait responses
that direct selection would sacrifice, and not as a pass-or-fail test.
Conversely, a high RE does not prove that the multi-trait objective was
attained.

## References

Rahimi M, Debnath S (2023). *Scientific Reports* 13:18977.
[doi:10.1038/s41598-023-46368-6](https://doi.org/10.1038/s41598-023-46368-6)

## See also

[`selection_index()`](https://FAkohoue.github.io/DesiredGainR/reference/selection_index.md)
