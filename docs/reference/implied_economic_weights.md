# Map desired gains to the aggregate weights they imply

The Smith-Hazel economic index and the Yamada desired-gain index both
use linear candidate scores. They express different breeding objectives.
Setting their coefficient vectors equal gives an algebraic mapping.
Thus, \\b = P^{-1}Gw\\ equal to \\b = P^{-1}G(GP^{-1}G)^{-1}d\\ gives
\$\$w = G^{-1}PG^{-1}d.\$\$ Hence every desired-gain vector corresponds
to one implied aggregate-weight vector under the supplied covariance
model, provided the required inverses exist.
[`implied_desired_gains()`](https://FAkohoue.github.io/DesiredGainR/reference/implied_desired_gains.md)
performs the reverse mapping. The implied vector reproduces the same
index direction. It is not an independently estimated price, profit or
biological value.

## Usage

``` r
implied_economic_weights(
  desired_gains,
  G,
  P,
  lower_is_better = NULL,
  gain_units = c("trait", "genetic_sd", "phenotypic_sd")
)
```

## Arguments

- desired_gains:

  Named numeric vector of desired gains. When `lower_is_better` is
  supplied these are magnitudes in the favourable-direction space, so a
  positive value always means improvement.

- G:

  Genetic variance-covariance matrix, named by trait, in original trait
  units.

- P:

  Phenotypic variance-covariance matrix, named by trait, in original
  trait units.

- lower_is_better:

  Traits for which smaller original values are favourable. Supplying
  this orients `G` and `P` internally, so that `desired_gains` and the
  returned weights are both expressed as improvements. Omitting it means
  every quantity is interpreted in the raw trait direction, where
  improving a trait such as disease severity requires a negative desired
  gain.

- gain_units:

  Units of `desired_gains`: `"trait"` for original trait units,
  `"genetic_sd"` for genetic standard deviations, or `"phenotypic_sd"`
  for phenotypic standard deviations.

## Value

A named numeric vector of implied economic weights, carrying a
`provenance` attribute.

## Details

This is the most direct answer to a breeder who cannot state economic
weights but can state relative desired gains, and, more usefully, to a
breeder who proposes weights and needs to see what response those
weights actually request.

## Negative implied weights are expected, not an error

An implied weight can be negative for a trait the breeder wants to
improve. This is correct and interpretable. Where a trait receives a
large favourable correlated response from the rest of the objective,
achieving only the modest gain requested for it requires the index to
hold it back, so its implied weight turns negative.

When `lower_is_better` is supplied the weights are returned in the
favourable-direction space, in which larger is better for every trait. A
negative implied weight can hold back favourable change that
correlations would otherwise make larger than requested. The
desired-gain vector remains the breeder's statement of direction.

## The desired gains may be standardised while the weights are not

`gain_units = "genetic_sd"` means that `desired_gains` is already
expressed in genetic standard deviations. However, when `G` and `P` are
supplied in original trait units, this function converts the target to
those units and returns a weight per original unit. A coefficient per
tonne cannot be compared numerically with a coefficient per day or
centimetre.

Multiply each returned weight by that trait's genetic standard deviation
to express the effect associated with one genetic standard deviation.
This places the weights on a common scale; it does not standardise the
desired gains a second time. If `G` and `P` have already been
transformed to genetic standard-deviation units, the returned weights
are already per genetic standard deviation and must not be multiplied by
the original standard deviations.
[`effective_weights()`](https://FAkohoue.github.io/DesiredGainR/reference/effective_weights.md)
performs the corresponding calculation for index coefficients.

Consequently these coefficients should not be read as statements of
biological or economic importance. Covarrubias-Pazaran (2021) puts the
point plainly: the weights lack meaning, especially under strong genetic
correlations, and the desired response is the only decision of interest.
Use the implied weights to drive an index or to compare objectives, and
use
[`effective_weights()`](https://FAkohoue.github.io/DesiredGainR/reference/effective_weights.md)
when the question is which traits the index is really acting on.

## References

Covarrubias-Pazaran G (2021). *Practical implementation of selection
indices.* CGIAR Excellence in Breeding.

## See also

[`implied_desired_gains()`](https://FAkohoue.github.io/DesiredGainR/reference/implied_desired_gains.md),
[`gain_feasibility()`](https://FAkohoue.github.io/DesiredGainR/reference/gain_feasibility.md)

## Examples

``` r
traits <- c("yield", "disease")
G <- matrix(c(0.60, -0.15, -0.15, 0.40), 2, dimnames = list(traits, traits))
P <- matrix(c(1.10, -0.20, -0.20, 0.90), 2, dimnames = list(traits, traits))

# Disease severity should fall, so it is declared rather than signed by hand.
implied_economic_weights(
  c(yield = 0.5, disease = 0.3), G, P,
  lower_is_better = "disease"
)
#>     yield   disease 
#> 1.2212974 0.9845422 
#> attr(,"provenance")
#> [1] "Implied by the supplied desired gains through w = G^-1 P G^-1 d; not an independently estimated economic value. Expressed in the favourable-direction space, so a positive weight favours movement in the breeder-defined direction."
```
