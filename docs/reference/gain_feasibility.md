# Test whether a desired-gain vector is attainable

For any linear index the response vector is \\R =
iGb/\sqrt{b^\mathsf{T}Pb}\\, so the attainable responses lie on the
ellipsoid \\R^\mathsf{T}G^{-1}PG^{-1}R = i^{2}\\. Therefore a target
response \\d\\ requires exactly \$\$i\_{\mathrm{required}} =
\sqrt{d^\mathsf{T}G^{-1}PG^{-1}d},\$\$ and the selected proportion
delivering that intensity follows from the normal-truncation
relationship.

## Usage

``` r
gain_feasibility(
  desired_gains,
  G,
  P,
  n_candidates,
  n_select = NULL,
  selection_proportion = NULL,
  lower_is_better = NULL,
  gain_units = c("trait", "genetic_sd", "phenotypic_sd")
)
```

## Arguments

- desired_gains:

  Named numeric vector of target responses. When `lower_is_better` is
  supplied these are magnitudes in the favourable-direction space, so a
  positive value always means improvement.

- G:

  Genetic variance-covariance matrix, named by trait, in original trait
  units.

- P:

  Phenotypic variance-covariance matrix, named by trait, in original
  trait units.

- n_candidates:

  Number of selection candidates available.

- n_select:

  Number of candidates to be selected. Supply either this or
  `selection_proportion`.

- selection_proportion:

  Proportion to be selected, in `(0, 1]`.

- lower_is_better:

  Traits for which smaller original values are favourable. See
  [`implied_economic_weights()`](https://FAkohoue.github.io/DesiredGainR/reference/implied_economic_weights.md).

- gain_units:

  Units of `desired_gains`, as in
  [`implied_economic_weights()`](https://FAkohoue.github.io/DesiredGainR/reference/implied_economic_weights.md).

## Value

An object of class `desiredgainr_feasibility` reporting: (i) the
required selection intensity and the proportion that delivers it, (ii)
whether the target is attainable in a population of `n_candidates`,
(iii) the attainable response in the requested direction at the planned
intensity, and (iv) the per-trait shortfall.

The attainable response and the shortfall are returned twice:
`attainable_response` and `shortfall` are in original trait units,
whereas `attainable_response_input_units` and `shortfall_input_units`
are in whatever `gain_units` the target was stated in. Comparing a
request made in standard deviations against an answer given in trait
units is an easy way to misread the result.

## Details

Two consequences deserve emphasis. First, the classical Pesek-Baker
index honours only the *direction* of the desired-gain vector; scaling
every element by a constant leaves the index unchanged, because the
attainable magnitude is fixed by selection intensity. Second, an
antagonistic correlation structure can make a modest-looking target
unreachable at any practical intensity.

This function replaces the external feasibility check that
Covarrubias-Pazaran (2021) recommends performing in separate software,
and should be run before an optimisation is attempted rather than after
it fails.

## References

Covarrubias-Pazaran G (2021). *Practical implementation of selection
indices.* CGIAR Excellence in Breeding.

## Examples

``` r
traits <- c("yield", "disease")
G <- matrix(c(0.60, -0.15, -0.15, 0.40), 2, dimnames = list(traits, traits))
P <- matrix(c(1.10, -0.20, -0.20, 0.90), 2, dimnames = list(traits, traits))
gain_feasibility(
  desired_gains = c(yield = 1.0, disease = -0.8),
  G = G, P = P, n_candidates = 500, n_select = 50,
  gain_units = "genetic_sd"
)
#> <desiredgainr_feasibility>
#>   Required selection intensity: 1.5203
#>   Requires selecting the top 16.014% (81 of 500 candidates)
#>   Planned intensity: 1.7550 (top 10.0%)
#>   Feasible at planned intensity: yes
#>   Feasible anywhere in this population: yes
#>   Attainable fraction of the requested gain: 115.4%
```
