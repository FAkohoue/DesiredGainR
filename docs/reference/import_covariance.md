# Validate and record an extracted genetic covariance matrix

Mixed-model packages report multi-trait variance components in different
shapes. Except for simple named ASReml component vectors, this function
does not extract a fitted object: the caller supplies the intended
covariance term and this function validates it. The errors it prevents
are not subtle: a covariance in the wrong trait order, a correlation
matrix mistaken for a covariance, or a genetic-by-environment block read
as a genetic one. This function converts the common formats and
validates what it produces.

## Usage

``` r
import_covariance(
  x,
  trait_cols,
  source = c("matrix", "asreml", "sommer", "breedR", "bglr", "stagewise"),
  estimand = "additive genetic covariance",
  P = NULL,
  is_correlation = FALSE,
  variances = NULL
)
```

## Arguments

- x:

  A covariance matrix, or a named vector of variance components.

- trait_cols:

  Trait names, in the order the rest of the analysis uses.

- source:

  Which upstream package produced `x`; see Details.

- estimand:

  What `x` represents, recorded in the provenance.

- P:

  Optional phenotypic covariance, checked for admissibility with `x`.

- is_correlation:

  Whether `x` is a correlation matrix that must be rescaled by
  `variances`.

- variances:

  Named trait variances, required when `is_correlation` is `TRUE`.

## Value

A list of class `desiredgainr_imported_covariance` holding the ordered
matrix, its diagnostics and its provenance.

## What is checked

The trait names are matched to `trait_cols` and reordered, not assumed
to already be in that order. Symmetry, finiteness and positive
semidefiniteness are verified. Where `P` is supplied as well, the pair
is checked for admissibility, because a genetic covariance exceeding the
phenotypic one in any direction is what produces heritabilities above
one downstream.

## Source labels

- `"asreml"`:

  A named vector or matrix of variance components from
  `summary(fit)$varcomp`, or an explicit covariance matrix.

- `"sommer"`:

  A matrix such as `fit$sigma[[term]]`, extracted by the caller.

- `"breedR"`, `"bglr"`, `"stagewise"`:

  A covariance matrix extracted by the caller; this function validates
  and orders it.

- `"matrix"`:

  Any matrix, for the case where the extraction has already been done.

Because the upstream APIs change between versions, this function accepts
the extracted matrix rather than the fitted object wherever possible.
That keeps the contract stable: the caller is responsible for pulling
the right term out of their model, and this function is responsible for
checking it is usable.

## See also

[`estimate_genetic_covariance()`](https://FAkohoue.github.io/DesiredGainR/reference/estimate_genetic_covariance.md),
[`matrix_diagnostics()`](https://FAkohoue.github.io/DesiredGainR/reference/matrix_diagnostics.md),
[`bend_covariance()`](https://FAkohoue.github.io/DesiredGainR/reference/bend_covariance.md)

## Examples

``` r
traits <- c("yield", "protein")
# A correlation matrix and separate variances, which is how several
# packages report multi-trait results.
correlation <- matrix(c(1, 0.3, 0.3, 1), 2,
  dimnames = list(traits, traits)
)
import_covariance(
  correlation, traits,
  source = "matrix",
  is_correlation = TRUE, variances = c(yield = 1.4, protein = 0.6)
)
#> <desiredgainr_imported_covariance>
#>   Source: matrix 
#>   Estimand: additive genetic covariance 
#>   Traits: yield, protein 
#>   Condition number: 2.886   Positive definite: yes
#>   Rescaled from a correlation matrix using the supplied variances.
```
