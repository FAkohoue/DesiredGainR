# Stress-test a desired-gain vector across simulation scenarios

The function applies matched replicate seeds across calibrated setups.
It reports mean gain, Monte Carlo standard error, tail risk, target
attainment, and regret. Replication can continue until the utility
standard error reaches a stated precision.

## Usage

``` r
stress_test_desired_gains(setups, desired_gains, options = list())
```

## Arguments

- setups:

  Named list of objects returned by
  [`founder_population()`](https://FAkohoue.github.io/DesiredGainR/reference/founder_population.md).

- desired_gains:

  Named desired-gain direction.

- options:

  Named list. Fields include simulation, minimum_gains, min_replicates,
  max_replicates, batch_size, utility_mcse, and seed.

## Value

An object of class desiredgainr_stress_test.
