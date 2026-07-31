# Report the effective contribution of each trait to an index

Index coefficients are not comparable across traits measured on
different scales. Crosbie et al. (1980) observed that assigning equal
weights to unstandardised traits places most of the selection pressure
on whichever trait carries the largest genetic variance. Hence the
interpretable quantity is not \\b_j\\ but \\b_j\sigma\_{g_j}\\, the
contribution of trait \\j\\ per unit of the genetic variation actually
available.

## Usage

``` r
effective_weights(
  coefficients,
  G,
  P = NULL,
  warn = FALSE,
  dominance_threshold = NULL
)
```

## Arguments

- coefficients:

  Named numeric vector of index coefficients.

- G:

  Genetic variance-covariance matrix, named by trait.

- P:

  Optional phenotypic variance-covariance matrix. When supplied, the
  phenotypic effective weights are reported alongside the genetic ones.

- warn:

  Whether to warn when one trait exceeds `dominance_threshold`. The
  default is `FALSE`, because a concentrated effective weight can arise
  from a deliberately asymmetric objective as readily as from mismatched
  trait scales, and warning on the former would be noise.

- dominance_threshold:

  Share of the total absolute effective weight above which a warning is
  issued for a single trait. The default adapts to the number of traits,
  since with two traits a share of 0.75 is merely a three-to-one
  emphasis whereas with eight traits it is extreme.

## Value

A `data.table` with one row per trait giving the coefficient, the
genetic and phenotypic effective weights, and their shares. Shares are
`NA` when the effective weights sum to zero.

## Details

Use this diagnostic whenever traits are analysed on their original
scales, and treat a single trait dominating the effective weights as a
signal to standardise before proceeding.

## References

Crosbie TM, Mock JJ, Smith OS (1980), as discussed in Guimaraes PHR et
al. (2021) *Euphytica* 217:95.
