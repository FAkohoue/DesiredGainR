# Score a new candidate set with a fitted index

An index fitted on one cycle is normally applied to the next. Doing that
by hand means reproducing the direction, centring and scaling
transformations exactly, and a transformation applied inconsistently
between fitting and scoring silently reorders the candidates.

## Usage

``` r
# S3 method for class 'desiredgainr_index'
predict(object, newdata, id_col = NULL, n_select = object$n_select, ...)
```

## Arguments

- object:

  A fitted `desiredgainr_index` using a coefficient-based method.

- newdata:

  Data frame or matrix carrying every trait in `object$trait_cols`, in
  the original trait units.

- id_col:

  Optional column of `newdata` holding candidate identifiers.

- n_select:

  Optional number of candidates to select. Defaults to the number used
  when fitting.

- ...:

  Unused.

## Value

A `data.table` with one row per candidate, giving the identifier, the
index score, the rank, and whether it was selected.

## Details

The centring constants and scaling factors are those computed when the
index was fitted, not recomputed from `newdata`. That is deliberate.
Recomputing them would score each cycle on its own mean, so a candidate
set that had advanced genetically would appear identical to one that had
not, and the scores would not be comparable across cycles.

The consequence is that `newdata` must be measured on the same scale and
in the same units as the data the index was fitted on. Where a trial has
been re-standardised or a trait recorded differently, refit rather than
predict.

## See also

[`selection_index()`](https://FAkohoue.github.io/DesiredGainR/reference/selection_index.md)

## Examples

``` r
set.seed(1)
traits <- c("yield", "protein")
G <- matrix(c(1.0, 0.2, 0.2, 0.5), 2, dimnames = list(traits, traits))
P <- matrix(c(2.5, 0.4, 0.4, 1.2), 2, dimnames = list(traits, traits))
cycle1 <- as.data.frame(matrix(
  stats::rnorm(80),
  ncol = 2, dimnames = list(paste0("A", 1:40), traits)
))
cycle2 <- as.data.frame(matrix(
  stats::rnorm(60),
  ncol = 2, dimnames = list(paste0("B", 1:30), traits)
))
fit <- selection_index(
  cycle1, traits,
  method = "smith_hazel", G = G, P = P,
  economic_weights = c(yield = 2, protein = 1), n_select = 8
)
head(predict(fit, cycle2, n_select = 5))
#>        id     score  rank selected
#>    <char>     <num> <int>   <lgcl>
#> 1:     B3 1.6706452     1     TRUE
#> 2:    B12 1.6507905     2     TRUE
#> 3:    B30 1.3458000     3     TRUE
#> 4:    B15 1.2368344     4     TRUE
#> 5:    B13 0.7910313     5     TRUE
#> 6:    B16 0.7326888     6    FALSE
```
