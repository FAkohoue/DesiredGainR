# Reshape long-format multi-environment trial data for indexing

Converts genotype-by-environment records in long format into the wide
trait-by-environment matrix that
[`expand_environments()`](https://FAkohoue.github.io/DesiredGainR/reference/expand_environments.md)
describes.

## Usage

``` r
widen_environments(
  data,
  id_col,
  environment_col,
  trait_cols,
  environments = NULL,
  missing_policy = c("error", "drop", "mean_impute")
)
```

## Arguments

- data:

  A data frame with one row per genotype, environment and trait
  combination.

- id_col, environment_col:

  Column names identifying the genotype and the environment.

- trait_cols:

  Trait column names.

- environments:

  Optional subset and ordering of environments.

- missing_policy:

  What to do with genotypes not present in every environment: `"error"`,
  `"drop"` to keep only complete genotypes, or `"mean_impute"`.

## Value

A data frame with one row per genotype and one column per
trait-environment combination, named to match
[`expand_environments()`](https://FAkohoue.github.io/DesiredGainR/reference/expand_environments.md).

## See also

[`expand_environments()`](https://FAkohoue.github.io/DesiredGainR/reference/expand_environments.md)

## Examples

``` r
long <- expand.grid(
  geno = paste0("g", 1:5),
  env = c("irrigated", "rainfed"),
  stringsAsFactors = FALSE
)
set.seed(1)
long$yield <- stats::rnorm(nrow(long))
long$protein <- stats::rnorm(nrow(long))
widen_environments(long, "geno", "env", c("yield", "protein"))
#>   geno yield_irrigated protein_irrigated yield_rainfed protein_rainfed
#> 1   g1      -0.6264538         1.5117812    -0.8204684     -0.04493361
#> 2   g2       0.1836433         0.3898432     0.4874291     -0.01619026
#> 3   g3      -0.8356286        -0.6212406     0.7383247      0.94383621
#> 4   g4       1.5952808        -2.2146999     0.5757814      0.82122120
#> 5   g5       0.3295078         1.1249309    -0.3053884      0.59390132
```
