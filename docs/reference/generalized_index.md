# Fit an index from a general information model

The economic solution is \\b=P^{-1}Ca\\. The desired-gain solution is
\\b=P^{-1}C(C^\mathsf{T}P^{-1}C)^{-1}d\\. The expected response is
\\\Delta g=iC^\mathsf{T}b/\sqrt{b^\mathsf{T}Pb}\\. Hence the records and
objective traits may differ in number and meaning.

## Usage

``` r
generalized_index(
  model,
  objective,
  method = c("economic", "desired_gain"),
  n_select = NULL,
  selection_intensity = NULL,
  lower_is_better = NULL
)
```

## Arguments

- model:

  Object returned by
  [`selection_information()`](https://FAkohoue.github.io/DesiredGainR/reference/selection_information.md).

- objective:

  Named economic weights or desired gains. Desired gains must use the
  trait units of `model$G`. Convert genetic-standard-deviation gains by
  multiplying them by `sqrt(diag(model$G))`.

- method:

  Either economic or desired_gain.

- n_select:

  Optional number of candidates selected.

- selection_intensity:

  Optional standardised selection intensity.

- lower_is_better:

  Objective traits for which smaller original values are favourable.
  When supplied, objective entries are favourable magnitudes. The
  function gives these traits a negative sign in the original trait
  coordinates.

## Value

An object of class desiredgainr_generalized_index.

## See also

[`selection_information()`](https://FAkohoue.github.io/DesiredGainR/reference/selection_information.md),
[`selection_index()`](https://FAkohoue.github.io/DesiredGainR/reference/selection_index.md)
