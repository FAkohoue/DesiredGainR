# Audit the calibration of an AlphaSimR setup

The report compares the target and realised covariance. It also checks
heritability, marker coverage, and marker overlap with quantitative
trait loci. Each check receives a clear status.

## Usage

``` r
simulation_calibration(
  setup,
  tolerances = list(variance = 0.05, correlation = 0.05, heritability = 0.05)
)
```

## Arguments

- setup:

  Object returned by
  [`founder_population()`](https://FAkohoue.github.io/DesiredGainR/reference/founder_population.md).

- tolerances:

  Named list with variance, correlation, and heritability limits.

## Value

An object of class desiredgainr_calibration.
