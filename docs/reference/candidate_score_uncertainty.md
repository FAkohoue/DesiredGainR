# Candidate-specific uncertainty in index scores and ranks

For candidate \\i\\, the score standard error is
\\\sqrt{b^\mathsf{T}PEC_i b}\\. Monte Carlo draws then estimate
selection probabilities. Within-candidate correlations among variables
are retained. Prediction errors are treated as independent among
candidates.

## Usage

``` r
candidate_score_uncertainty(
  index,
  prediction_error_covariance,
  n_draws = 2000L,
  level = 0.95,
  seed = NULL
)
```

## Arguments

- index:

  A fitted DesiredGainR coefficient-based index.

- prediction_error_covariance:

  One shared covariance matrix, one matrix per candidate in a list, or a
  three-dimensional array.

- n_draws:

  Number of Monte Carlo draws.

- level:

  Coverage of the score interval.

- seed:

  Optional random seed.

## Value

An object of class desiredgainr_candidate_uncertainty.
