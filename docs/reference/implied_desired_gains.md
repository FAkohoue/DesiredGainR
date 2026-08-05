# Map economic weights to the desired gains they imply

Inverse of
[`implied_economic_weights()`](https://FAkohoue.github.io/DesiredGainR/reference/implied_economic_weights.md).
For economic weights \\w\\, the Smith-Hazel index requests the response
\$\$d = GP^{-1}Gw,\$\$ up to the scalar set by selection intensity.
Reading these implied gains is the quickest way for a breeder to
discover that a proposed weight vector asks for something other than
what was intended.

## Usage

``` r
implied_desired_gains(
  economic_weights,
  G,
  P,
  lower_is_better = NULL,
  gain_units = c("trait", "genetic_sd", "phenotypic_sd")
)
```

## Arguments

- economic_weights:

  Named numeric vector of economic weights, in the same space as
  `lower_is_better` implies.

- G:

  Genetic variance-covariance matrix, named by trait, in original trait
  units.

- P:

  Phenotypic variance-covariance matrix, named by trait, in original
  trait units.

- lower_is_better:

  Traits for which smaller original values are favourable. See
  [`implied_economic_weights()`](https://FAkohoue.github.io/DesiredGainR/reference/implied_economic_weights.md).

- gain_units:

  Units in which the implied gains are returned. Set this to match the
  units in which the weights were derived, so that a round trip through
  [`implied_economic_weights()`](https://FAkohoue.github.io/DesiredGainR/reference/implied_economic_weights.md)
  returns the vector it started from.

## Value

A named numeric vector of implied desired gains, carrying a `provenance`
attribute recording the units.

## See also

[`implied_economic_weights()`](https://FAkohoue.github.io/DesiredGainR/reference/implied_economic_weights.md)

## Examples

``` r
traits <- c("yield", "disease")
G <- matrix(c(0.60, -0.15, -0.15, 0.40), 2, dimnames = list(traits, traits))
P <- matrix(c(1.10, -0.20, -0.20, 0.90), 2, dimnames = list(traits, traits))

# The algebraic mapping is invertible when the units and model match.
d <- c(yield = 0.5, disease = 0.3)
w <- implied_economic_weights(
  d, G, P,
  lower_is_better = "disease", gain_units = "genetic_sd"
)
implied_desired_gains(
  w, G, P,
  lower_is_better = "disease", gain_units = "genetic_sd"
)
#>   yield disease 
#>     0.5     0.3 
#> attr(,"provenance")
#> [1] "Implied by the supplied economic weights through d = G P^-1 G w, expressed in genetic_sd units and up to the scalar set by selection intensity."
```
