# Define one objective for a cross-family comparison

A comparison objective fixes the yardstick before fitted methods are
inspected. Desired gains define a response target. `aggregate_weights`
and `W` define one common utility
\$\$U(g)=a^\mathsf{T}g+g^\mathsf{T}Wg.\$\$ Every method is evaluated
against this same objective.

## Usage

``` r
comparison_objective(
  desired_gains = NULL,
  aggregate_weights = NULL,
  W = NULL,
  G = NULL,
  gain_units = c("trait", "genetic_sd")
)
```

## Arguments

- desired_gains:

  Optional named favourable desired-gain vector.

- aggregate_weights:

  Optional named linear aggregate-weight vector.

- W:

  Optional symmetric matrix of squared and cross-product utility
  weights. Its names must match the objective traits.

- G:

  Optional genetic covariance matrix in original trait coordinates. It
  is required when `gain_units = "genetic_sd"`. It also supplies one
  common geometry for target alignment and common-merit evaluation.

- gain_units:

  Coordinate system for the objective. Use either original trait units
  or genetic standard deviations. `desired_gains` are expressed in these
  coordinates. `aggregate_weights` and `W` must act on gains in the same
  coordinates.

## Value

An object of class `desiredgainr_comparison_objective`.

## See also

[`compare_selection_methods()`](https://FAkohoue.github.io/DesiredGainR/reference/compare_selection_methods.md)

## Examples

``` r
traits <- c("yield", "disease")
G <- matrix(c(1, -0.2, -0.2, 0.8), 2,
  dimnames = list(traits, traits)
)
comparison_objective(
  desired_gains = c(yield = 1, disease = 0.5),
  G = G,
  gain_units = "genetic_sd"
)
#> $traits
#> [1] "yield"   "disease"
#> 
#> $desired_gains
#>   yield disease 
#>     1.0     0.5 
#> 
#> $aggregate_weights
#> NULL
#> 
#> $W
#> NULL
#> 
#> $G
#>         yield disease
#> yield     1.0    -0.2
#> disease  -0.2     0.8
#> 
#> $gain_units
#> [1] "genetic_sd"
#> 
#> attr(,"class")
#> [1] "desiredgainr_comparison_objective" "list"                             
```
