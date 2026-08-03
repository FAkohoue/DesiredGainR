# Restricted and proportional-gain selection indices

Classical constrained index families. Each imposes a linear constraint
on the expected genetic response and maximises merit subject to it, in
closed form.

## Usage

``` r
restricted_index(
  values,
  trait_cols,
  method = c("kempthorne_nordskog", "tallis", "mallard", "harville",
    "restricted_smith_hazel"),
  G,
  P,
  economic_weights = NULL,
  restricted_traits = NULL,
  target_gains = NULL,
  constraint_matrix = NULL,
  penalty = Inf,
  lower_is_better = NULL,
  id_col = NULL,
  center_traits = TRUE,
  scale_traits = TRUE,
  scale_by = c("sample", "phenotypic"),
  n_select = NULL,
  selection_intensity = NULL,
  main_trait = trait_cols[1L]
)
```

## Arguments

- values:

  Data frame or matrix of candidate trait values.

- trait_cols:

  Trait column names.

- method:

  Index family; see Details.

- G, P:

  Genetic and phenotypic covariance matrices, named by trait.

- economic_weights:

  Named economic weights defining net merit. Required by
  zero-restriction families and used to choose among any free response
  directions left by a proportionality constraint.

- restricted_traits:

  Traits held to zero expected response, for `"kempthorne_nordskog"` and
  `"restricted_smith_hazel"`.

- target_gains:

  Named responses for the constrained traits. Proportions for
  `"tallis"`, `"mallard"` and `"harville"`.

- constraint_matrix:

  Optional numeric contrast matrix with one column per trait. For a
  zero-restriction method it enforces
  `constraint_matrix %*% expected_response = 0`; this supports, for
  example, stability contrasts returned by
  [`expand_environments()`](https://FAkohoue.github.io/DesiredGainR/reference/expand_environments.md).

- penalty:

  Deprecated. Finite soft penalties were never Harville's published
  method and are rejected; use an exact published constraint.

- lower_is_better:

  Traits for which smaller original values are favourable.

- id_col:

  Optional column of `values` holding candidate identifiers.

- center_traits, scale_traits, scale_by:

  Passed through to the same transformation
  [`selection_index()`](https://FAkohoue.github.io/DesiredGainR/reference/selection_index.md)
  applies.

- n_select:

  Number of candidates selected.

- selection_intensity:

  Optional standardised selection intensity.

- main_trait:

  Trait against which relative efficiency is computed.

## Value

An object of class `desiredgainr_index`, as returned by
[`selection_index()`](https://FAkohoue.github.io/DesiredGainR/reference/selection_index.md),
with an added `constraint` element recording the restriction and how
nearly it was met.

## Which family answers which question

- `"kempthorne_nordskog"`:

  *Hold these traits still.* Maximises the correlation with net merit
  subject to zero expected response in the restricted traits. Use when a
  trait must not move: a quality standard, a maturity window, a plant
  height ceiling.

- `"tallis"`, `"mallard"`, `"harville"`:

  *Move these traits in this ratio.* These published
  predetermined-proportional-gain formulations are algebraically
  equivalent up to an arbitrary scale of the coefficient vector. They
  constrain response proportions, not unattainable absolute response
  magnitudes. The separate method names retain citation and reporting
  provenance.

- `"restricted_smith_hazel"`:

  Smith-Hazel with the restriction applied by projection rather than by
  Lagrange multipliers. Algebraically equivalent to
  `"kempthorne_nordskog"` and retained because the projection form is
  what most textbooks present.

## The algebra

Let \\\mathbf{C}\\ be the \\k \times p\\ matrix picking out the
constrained traits. Expected response is proportional to
\\\mathbf{G}\mathbf{b}\\, so a zero-response restriction is
\\\mathbf{C}\mathbf{G}\mathbf{b} = \mathbf{0}\\. Maximising
\\\mathbf{b}^\mathsf{T}\mathbf{G}\mathbf{a}\\ subject to that and to a
scale convention gives

\$\$\mathbf{b} = \left\[\mathbf{I} - \mathbf{P}^{-1}\mathbf{G}
\mathbf{C}^\mathsf{T}(\mathbf{C}\mathbf{G}\mathbf{P}^{-1}\mathbf{G}
\mathbf{C}^\mathsf{T})^{-1}\mathbf{C}\mathbf{G}\right\]
\mathbf{P}^{-1}\mathbf{G}\mathbf{a}\$\$

which is Kempthorne and Nordskog (1959). The proportional-gain methods
replace trait selectors with contrasts that are zero precisely when the
responses have the requested proportions.

## Why these are worth having beside [`run_dgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_dgsi.md)

[`run_dgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_dgsi.md)
reaches a constrained solution by searching, which introduces a
dependence on the number of replicates and a selection optimism that has
to be corrected for. Where the constraint is linear, these families give
the same answer exactly. Comparing the two is a useful check: a large
divergence means the search has not converged, not that the constraint
is unusual.

## References

Kempthorne, O. and Nordskog, A.W. (1959) Restricted selection indices.
*Biometrics* 15, 10-19.

Tallis, G.M. (1962) A selection index for optimum genotype. *Biometrics*
18, 120-122.

Mallard, J. (1972) La theorie et le calcul des index de selection avec
restrictions: synthese critique. *Biometrics* 28, 713-735.

Kemp, C.D. and Harville, D.A. (1975) Index selection with
proportionality constraints. *Biometrics* 31, 223-225.
[doi:10.2307/2529722](https://doi.org/10.2307/2529722)

## See also

[`selection_index()`](https://FAkohoue.github.io/DesiredGainR/reference/selection_index.md),
[`run_dgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_dgsi.md),
[`gain_feasibility()`](https://FAkohoue.github.io/DesiredGainR/reference/gain_feasibility.md)

## Examples

``` r
set.seed(1)
traits <- c("yield", "height", "protein")
G <- matrix(
  c(
    1.0, 0.4, 0.2,
    0.4, 0.8, 0.1,
    0.2, 0.1, 0.5
  ), 3,
  dimnames = list(traits, traits)
)
P <- G + diag(c(1.2, 1.0, 0.9))
dimnames(P) <- list(traits, traits)
values <- as.data.frame(matrix(
  stats::rnorm(90),
  ncol = 3, dimnames = list(paste0("g", 1:30), traits)
))

# Improve yield and protein while holding height still.
fit <- restricted_index(
  values, traits,
  method = "kempthorne_nordskog",
  G = G, P = P,
  economic_weights = c(yield = 2, height = 0, protein = 1),
  restricted_traits = "height", n_select = 6
)
fit$constraint$achieved_response
#>       height 
#> -5.64808e-17 
```
