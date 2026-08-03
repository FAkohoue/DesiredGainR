# Draw jointly admissible genetic and residual covariance matrices

A pair \\(\mathbf{G}^{\*}, \mathbf{P}^{\*})\\ is only admissible when
\\\mathbf{P}^{\*} - \mathbf{G}^{\*}\\ is itself a covariance matrix.
Resampling the two independently produces pairs that violate this, and
an index built on such a pair reports a heritability above one. Drawing
the genetic and residual matrices separately and reassembling
\\\mathbf{P}^{\*} = \mathbf{G}^{\*} + \mathbf{E}^{\*}\\ makes every draw
admissible by construction.

## Usage

``` r
draw_covariance_pairs(
  G,
  P = NULL,
  genetic_df,
  residual_df = NULL,
  n_draws = 20L
)
```

## Arguments

- G:

  Genetic covariance matrix.

- P:

  Phenotypic covariance matrix. When `NULL`, only `G` is drawn and the
  result carries no residual.

- genetic_df, residual_df:

  Degrees of freedom. See
  [`index_uncertainty()`](https://FAkohoue.github.io/DesiredGainR/reference/index_uncertainty.md)
  for what these count; they are the number of independent genetic
  units, not the number of plots.

- n_draws:

  Number of draws.

## Value

A list of draws, each with `G`, `E` and `P`.

## See also

[`index_uncertainty()`](https://FAkohoue.github.io/DesiredGainR/reference/index_uncertainty.md),
[`optimize_desired_gains()`](https://FAkohoue.github.io/DesiredGainR/reference/optimize_desired_gains.md)

## Examples

``` r
traits <- c("yield", "protein")
G <- matrix(c(1.0, 0.2, 0.2, 0.5), 2, dimnames = list(traits, traits))
P <- matrix(c(2.5, 0.4, 0.4, 1.2), 2, dimnames = list(traits, traits))
draws <- draw_covariance_pairs(G, P,
  genetic_df = 40, residual_df = 200,
  n_draws = 5
)
length(draws)
#> [1] 5
# Every draw is admissible by construction.
all(vapply(draws, function(d) {
  min(eigen(d$P - d$G, symmetric = TRUE, only.values = TRUE)$values) > 0
}, logical(1)))
#> [1] TRUE
```
