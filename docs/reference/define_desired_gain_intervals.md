# Define acceptable intervals for desired genetic gains

Breeders are often more confident about an acceptable range than about
one exact desired-gain vector. This function records those ranges in
familiar raw trait directions and converts them to the
favourable-direction convention used by DesiredGainR. Thus a requested
reduction of 0.25 to 1.0 genetic standard deviations in disease severity
can be entered as `lower = -1` and `upper = -0.25` while declaring the
disease trait in `lower_is_better`.

## Usage

``` r
define_desired_gain_intervals(
  lower,
  upper,
  lower_is_better = NULL,
  gain_units = c("genetic_sd", "phenotypic_sd", "trait"),
  horizon_cycles = NULL
)
```

## Arguments

- lower, upper:

  Named numeric vectors giving the lower and upper bounds of acceptable
  genetic change in the raw trait direction. Both vectors must contain
  the same traits. For a trait in `lower_is_better`, reductions are
  normally negative; for other traits, improvements are normally
  positive.

- lower_is_better:

  Traits for which a reduction is favourable.

- gain_units:

  Units of the bounds: genetic standard deviations, phenotypic standard
  deviations, or original trait units. Genetic standard deviations are
  recommended for eliciting intervals across traits.

- horizon_cycles:

  Optional number of selection cycles over which the bounds are intended
  to be achieved. Recording the horizon prevents a one-cycle target from
  being silently evaluated as a five-cycle target.

## Value

A `desiredgainr_gain_intervals` data frame containing favourable-
direction lower and upper bounds, with the raw bounds and units retained
as attributes.

## Details

The resulting intervals constrain the *relative desired-gain directions*
searched by
[`optimize_desired_gains()`](https://FAkohoue.github.io/DesiredGainR/reference/optimize_desired_gains.md).
They are not promises that every value in the interval is attainable.
Use
[`gain_feasibility()`](https://FAkohoue.github.io/DesiredGainR/reference/gain_feasibility.md)
to test an absolute point and multi-cycle simulation to compare
admissible directions.

## Examples

``` r
intervals <- define_desired_gain_intervals(
  lower = c(Yield = 0.5, Disease = -1.0, Quality = 0.1),
  upper = c(Yield = 1.5, Disease = -0.25, Quality = 0.8),
  lower_is_better = "Disease",
  gain_units = "genetic_sd"
)
intervals
#> <desiredgainr_gain_intervals>
#>   Units: genetic_sd 
#>   Bounds shown as favourable genetic change (larger is better):
#> [1] trait lower upper
#> <0 rows> (or 0-length row.names)
#>   Lower-is-better traits were entered as raw reductions: Disease 
```
