# Evaluate an achieved response against a desired-gain direction

This function applies Satoh's one-dimensional restricted breeding-value
criterion to an achieved response. It reports progress along the
direction, proportional departure, and the residual Mahalanobis
distance.

## Usage

``` r
evaluate_restricted_response(response, direction, G, lower_is_better = NULL)
```

## Arguments

- response:

  Named achieved response vector.

- direction:

  Named proportional desired-gain direction.

- G:

  Genetic covariance matrix.

- lower_is_better:

  Traits for which smaller original values are favourable. When
  supplied, direction entries are favourable magnitudes.

## Value

An object of class desiredgainr_restricted_response.
