# Project breeding values into a restricted response space

Satoh (2024) showed that restricted breeding values are linear
projections of ordinary breeding values. For constraints \\C g=0\\, the
projection is

## Usage

``` r
restricted_breeding_values(
  breeding_values,
  G,
  direction = NULL,
  constraint_matrix = NULL,
  id_col = NULL,
  lower_is_better = NULL
)
```

## Arguments

- breeding_values:

  Candidate-by-trait matrix of breeding values or estimated breeding
  values.

- G:

  Genetic covariance matrix for the same traits.

- direction:

  Optional named proportional desired-gain direction.

- constraint_matrix:

  Optional restriction matrix with traits in columns.

- id_col:

  Optional identifier column in breeding_values.

- lower_is_better:

  Traits for which smaller original values are favourable. When
  supplied, direction entries are favourable magnitudes.

## Value

An object of class desiredgainr_restricted_bv.

## Details

\$\$g_R=\[I-GC^\mathsf{T}(CGC^\mathsf{T})^{-1}C\]g.\$\$

A proportional desired-gain direction \\d\\ defines a one-dimensional
space. Its projection is

\$\$g_R=d(d^\mathsf{T}G^{-1}d)^{-1}d^\mathsf{T}G^{-1}g =\beta d.\$\$

The scalar \\\beta\\ measures progress along the requested direction. It
combines response magnitude and proportional agreement in one quantity.

## References

Satoh M (2024). Characteristics of restricted selection indices and
geometrical interpretation of restricted breeding values. *Journal of
Animal Breeding and Genetics* 141:353-363.
[doi:10.1111/jbg.12845](https://doi.org/10.1111/jbg.12845)

## See also

[`restricted_index()`](https://FAkohoue.github.io/DesiredGainR/reference/restricted_index.md),
[`evaluate_restricted_response()`](https://FAkohoue.github.io/DesiredGainR/reference/evaluate_restricted_response.md)
