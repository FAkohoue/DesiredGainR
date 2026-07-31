# Recover the weights implied by past selection decisions

Where a programme has been selecting for years without a formal index,
the selection differentials it has achieved already encode its
objective. Covarrubias-Pazaran recovers the corresponding linear rule as
\$\$b = P^{-1}s,\$\$ where \\s\\ contains the differentials between the
selected group and the full population. This is how a selection index
was introduced at both the International Maize and Wheat Improvement
Center (CIMMYT) maize programmes and the International Rice Research
Institute (IRRI) programmes.

## Usage

``` r
retrospective_weights(selected_values, population_values, trait_cols, P = NULL)
```

## Arguments

- selected_values:

  Numeric matrix or data frame of trait values for the candidates that
  were historically selected, with one column per trait.

- population_values:

  Numeric matrix or data frame of trait values for the full population
  from which they were selected.

- trait_cols:

  Character vector naming and ordering the trait columns.

- P:

  Optional phenotypic variance-covariance matrix. When `NULL`, it is
  estimated from `population_values` and reported as such.

## Value

An object of class `desiredgainr_retrospective` containing the recovered
coefficients, the achieved selection differentials in both original and
standard-deviation units, and the provenance of `P`.

## Details

The recovered weights answer the question "what linear rule best
approximates the decisions already made?". They do not answer "what
should this programme select for next?". Therefore treat the result as a
starting point, inspect it with
[`weight_sensitivity()`](https://FAkohoue.github.io/DesiredGainR/reference/weight_sensitivity.md),
and adjust the emphasis deliberately before adopting it.

## References

Covarrubias-Pazaran G. *Bringing a selection index into the CIMMYT-Maize
programs* and *Bringing a selection index into the IRRI programs.* CGIAR
Excellence in Breeding.
