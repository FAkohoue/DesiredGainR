# Expand a trait-by-environment structure into an index problem

Builds the covariance matrices and economic weights for an index over
trait-environment combinations, so that a multi-environment problem can
be handed to
[`selection_index()`](https://FAkohoue.github.io/DesiredGainR/reference/selection_index.md),
[`restricted_index()`](https://FAkohoue.github.io/DesiredGainR/reference/restricted_index.md)
or
[`evaluate_index()`](https://FAkohoue.github.io/DesiredGainR/reference/evaluate_index.md)
without further preparation.

## Usage

``` r
expand_environments(
  trait_cols,
  environments,
  G,
  P,
  environment_correlation,
  residual_correlation = 0,
  economic_weights = NULL,
  environment_weights = NULL,
  include_stability = FALSE
)
```

## Arguments

- trait_cols:

  Trait names.

- environments:

  Environment names.

- G:

  Genetic covariance between traits, named by trait.

- P:

  Phenotypic covariance between traits, named by trait.

- environment_correlation:

  Between-environment genetic correlation matrix, named by environment.
  A single number applies that correlation to every pair. Values near 1
  mean little interaction.

- residual_correlation:

  Between-environment residual correlation. Defaults to 0, which is
  right when environments are separate trials.

- economic_weights:

  Named economic weights per trait.

- environment_weights:

  Named relative value of each environment, normalised to sum to one.

- include_stability:

  Whether to return per-trait stability response contrasts and a
  `constraint_matrix` consumable by
  [`restricted_index()`](https://FAkohoue.github.io/DesiredGainR/reference/restricted_index.md).

## Value

A list of class `desiredgainr_multi_environment` holding the expanded
`G`, `P`, economic weights and labels.

## The model

For trait \\j\\ in environment \\e\\, the genetic value is

\$\$g\_{je} = a_j + (ge)\_{je}\$\$

with \\a_j\\ the environment-invariant main effect and \\(ge)\_{je}\\
the interaction. The covariance between the same trait in two
environments is then \\\mathbf{G}\_{jj}\\, and between different traits
and environments \\\mathbf{G}\_{jk}\\ scaled by the environmental
correlation. This gives the separable structure

\$\$\mathrm{Cov}(g\_{je}, g\_{kf}) = \mathbf{G}\_{jk} \\
\mathbf{C}\_{ef}\$\$

where \\\mathbf{C}\\ is the between-environment genetic correlation
matrix. A \\\mathbf{C}\\ with off-diagonals near one means the trait
ranks genotypes almost identically everywhere and a single-environment
index would have sufficed; values near zero mean the environments are
effectively different traits.

Separability is an assumption, and it is the one to question first if
the result looks wrong. Where a fitted unstructured multi-environment
model is available, pass its covariance directly through
[`import_covariance()`](https://FAkohoue.github.io/DesiredGainR/reference/import_covariance.md)
instead of building one here.

## Environment-specific economic weights

A trait need not be worth the same everywhere: yield in a marginal
environment serving a small area contributes less to net merit than
yield in the main target environment. `environment_weights` scales each
environment's contribution, and the economic weight of trait \\j\\ in
environment \\e\\ becomes \\w_j \times v_e\\.

## Stability as an objective

Setting `include_stability = TRUE` returns, for each trait, response
contrasts equal to the difference between the trait in each environment
and its mean across environments. Pass `stability$constraint_matrix` to
[`restricted_index()`](https://FAkohoue.github.io/DesiredGainR/reference/restricted_index.md)
to impose zero response in those contrasts. They are not appended as
duplicate index variables, which would make the covariance matrix
singular because every contrast is a linear combination of the original
trait-environment variables.

Stability is not free. Constraining it costs response in the mean, and
whether the trade is worth making is exactly the kind of question
[`gain_feasibility()`](https://FAkohoue.github.io/DesiredGainR/reference/gain_feasibility.md)
and
[`restricted_index()`](https://FAkohoue.github.io/DesiredGainR/reference/restricted_index.md)
exist to answer.

## See also

[`selection_index()`](https://FAkohoue.github.io/DesiredGainR/reference/selection_index.md),
[`import_covariance()`](https://FAkohoue.github.io/DesiredGainR/reference/import_covariance.md),
[`restricted_index()`](https://FAkohoue.github.io/DesiredGainR/reference/restricted_index.md)

## Examples

``` r
traits <- c("yield", "protein")
environments <- c("irrigated", "rainfed")
G <- matrix(c(1.0, 0.2, 0.2, 0.5), 2, dimnames = list(traits, traits))
P <- matrix(c(2.5, 0.4, 0.4, 1.2), 2, dimnames = list(traits, traits))

expanded <- expand_environments(
  trait_cols = traits, environments = environments,
  G = G, P = P,
  environment_correlation = 0.6,
  economic_weights = c(yield = 2, protein = 1),
  environment_weights = c(irrigated = 0.7, rainfed = 0.3)
)
dim(expanded$G)
#> [1] 4 4
expanded$economic_weights
#>   yield_irrigated protein_irrigated     yield_rainfed   protein_rainfed 
#>               1.4               0.7               0.6               0.3 
```
