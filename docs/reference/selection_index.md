# Construct a classical multi-trait selection index

`selection_index()` builds any of the classical linear selection indices
from a common interface, so that alternative objectives and alternative
index families can be compared on identical data. All traits are first
oriented so that larger values are favourable.

## Usage

``` r
selection_index(
  values,
  trait_cols,
  id_col = NULL,
  method = c("smith_hazel", "base", "pesek_baker", "yamada", "mulamba_mock", "elston",
    "independent_culling", "tandem"),
  G = NULL,
  P = NULL,
  economic_weights = NULL,
  desired_gains = NULL,
  aggregate_weights = NULL,
  lower_is_better = NULL,
  center_traits = TRUE,
  scale_traits = TRUE,
  scale_by = c("sample", "phenotypic"),
  culling_thresholds = NULL,
  tandem_order = NULL,
  n_select = NULL,
  selection_intensity = NULL,
  main_trait = trait_cols[1L]
)
```

## Arguments

- values:

  Numeric matrix or data frame of candidate trait values, with one
  column per trait. Candidate identifiers are taken from `id_col` when
  supplied, and from the row names otherwise.

- trait_cols:

  Character vector naming and ordering the trait columns.

- id_col:

  Optional name of a column in `values` holding the candidate
  identifiers, matching the argument of the same name in
  [`run_dgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_dgsi.md)
  and
  [`run_qgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_qgsi.md).
  When `NULL`, row names are used. A data frame carrying identifiers in
  a column rather than in its row names will otherwise be labelled by
  row number, so a warning is issued when that appears to have happened.

- method:

  Index family. See Details.

- G:

  Genetic variance-covariance matrix, named by trait. Required by every
  method except `"base"`, `"mulamba_mock"`, `"elston"`,
  `"independent_culling"`, and `"tandem"`.

- P:

  Phenotypic variance-covariance matrix, named by trait. Required by
  `"smith_hazel"` and `"yamada"`.

- economic_weights:

  Named economic weights in the favourable-direction trait space.
  Required by `"smith_hazel"` and `"base"`, and optionally used to
  weight ranks in `"mulamba_mock"`. Negative values are valid for the
  economic indices. They can arise when correlated response would
  otherwise move a trait beyond its economically preferred level. Rank
  weights remain non-negative because they state relative emphasis after
  trait direction has been declared.

- desired_gains:

  Named non-negative desired gains in the favourable-direction trait
  space. Required by `"pesek_baker"` and `"yamada"`.

- aggregate_weights:

  Optional named weights defining one common net merit for evaluation.
  This is separate from `desired_gains`: when it is absent,
  merit-dependent criteria (`R_HI` and `Delta_H`) are `NA` for the
  desired-gain families. Implied weights can be explored explicitly with
  [`implied_economic_weights()`](https://FAkohoue.github.io/DesiredGainR/reference/implied_economic_weights.md),
  but are not silently substituted because a different desired-gain
  direction would then change the definition of merit being compared.

- lower_is_better:

  Traits for which smaller original values are favourable.

- center_traits:

  Whether to subtract the trait means before indexing. Centring does not
  change the ranking, because it shifts every score by the same
  constant, but it does place the index mean at zero and therefore makes
  the index coefficient of variation undefined. Set this to `FALSE` when
  reproducing published results that report `CV_I` on an index built
  from raw trait values.

- scale_traits:

  Whether to divide traits by their standard deviations before indexing.
  Covarrubias-Pazaran (2021) recommends this, and a desired gain of 1
  then means one standard deviation of progress. See `scale_by` for
  which standard deviation is used.

- scale_by:

  Source of the scaling factors when `scale_traits = TRUE`.

  `"sample"`, the default, divides by the standard deviations of the
  supplied candidates. This equals the population standard deviation
  only when the candidates are an unselected random sample of the
  population that `G` and `P` describe. Candidates at a late trial stage
  have already been selected, so their spread is narrower than the
  population's, and mixing that sample scale with population covariance
  matrices inflates the apparent heritability of every trait that
  selection has already narrowed.

  `"phenotypic"` divides by \\\sqrt{\operatorname{diag}(\mathbf{P})}\\,
  which comes from the same upstream model as `G` and `P` themselves.
  The scaled `P` then has a unit diagonal exactly, and the scaled `G`
  has the narrow-sense heritabilities on its diagonal. Prefer it
  whenever `P` is a genuine population estimate rather than a covariance
  computed from the candidates at hand.

- culling_thresholds:

  Named acceptance limits in the original trait units, required by
  `"elston"` and `"independent_culling"`. Supply a minimum for traits
  where larger values are favourable. Supply a maximum for traits named
  in `lower_is_better`. The function applies direction, centring, and
  scaling to these limits internally.

- tandem_order:

  Character vector giving the order in which traits are selected,
  required by `"tandem"`.

- n_select:

  Number of candidates selected.

- selection_intensity:

  Optional standardised selection intensity used for the
  expected-response calculations. When `NULL`, it is derived from
  `n_select` under normal truncation.

- main_trait:

  Trait against which relative efficiency is computed. Defaults to the
  first trait.

## Value

An object of class `desiredgainr_index` containing the coefficients, the
candidate scores and ranking, the evaluation criteria from
[`evaluate_index()`](https://FAkohoue.github.io/DesiredGainR/reference/evaluate_index.md),
the effective weights, and the provenance of every input.

## Available methods

- `"smith_hazel"`:

  The optimum economic index, \\b = P^{-1}Ga\\, maximising the
  correlation between the index and the aggregate genotype \\H =
  a^\mathsf{T}g\\. Requires `economic_weights`.

- `"base"`:

  The base index of Brim, Cockerham and Clark, \\b = a\\. It requires
  neither `P` nor `G` and is included because Rahimi and Debnath (2023)
  found it to match the optimum index almost exactly in their maize
  data. Where the two agree, the effort spent estimating covariance
  matrices has bought nothing, which is itself worth reporting.

- `"pesek_baker"`:

  The original desired-gain index of Pesek and Baker (1969), \\b =
  G^{-1}d\\, applied to genotypic values. Requires `desired_gains`.

- `"yamada"`:

  The desired-gain index of Yamada, Yokouchi and Nishida (1975) for
  phenotypic selection criteria, \\b = P^{-1}G(GP^{-1}G)^{-1}d\\. This
  is the formulation used by Joukhadar et al. (2024) and by
  [`run_dgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_dgsi.md).

- `"mulamba_mock"`:

  The rank-sum index of Mulamba and Mock (1978). Candidates are ranked
  for each trait in the favourable direction and the ranks are summed.
  It requires neither economic weights nor covariance matrices, and
  Guimaraes et al. (2021) found it to give the most balanced multi-trait
  response of the methods they compared.

- `"elston"`:

  The Elston multiplicative index. Candidates first meet every stated
  floor. Eligible candidates are ranked by the product of their margins
  above those floors. The logarithm of the product provides numerical
  stability.

- `"independent_culling"`:

  Not an index. Candidates must exceed a threshold for every trait.
  Included as a comparator because it is what most programmes actually
  do.

- `"tandem"`:

  Not an index. Candidates are selected sequentially, one trait at a
  time. Included as a comparator.

## Two formulations share the name Pesek-Baker, and they coincide

The literature applies the name to two expressions. Pesek and Baker
(1969) give \\b = G^{-1}d\\, which is `method = "pesek_baker"` here and
what Rahimi and Debnath (2023) implement. Yamada et al. (1975) give \\b
= P^{-1}G(GP^{-1}G)^{-1}d\\, which is `method = "yamada"` and what
Joukhadar et al. (2024) and
[`run_dgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_dgsi.md)
use.

For a square invertible \\G\\ these are the same index, because
\$\$P^{-1}G(GP^{-1}G)^{-1}d = P^{-1}G(G^{-1}PG^{-1})d = G^{-1}d.\$\$
Both are therefore retained as documented routes to one estimand rather
than as competing methods.

Two situations make the distinction real. First, when \\G\\ is singular
or rank deficient, the direct inverse does not exist and only the Yamada
route is available. Second, the two routes are not numerically
equivalent: one inverts \\G\\, the other inverts \\P\\ and then
\\GP^{-1}G\\, so an ill-conditioned matrix can send them apart. That is
the practical reason to state which route produced a published result,
and the reason
[`matrix_diagnostics()`](https://FAkohoue.github.io/DesiredGainR/reference/matrix_diagnostics.md)
reports conditioning.

## References

Brim CA, Cockerham HW, Clark C (1959). *Agronomy Journal* 51:42-46.

Guimaraes PHR, Melo PGS, Cordeiro ACC, Torga PP, Rangel PHN, de Castro
AP (2021). *Euphytica* 217:95.
[doi:10.1007/s10681-021-02819-7](https://doi.org/10.1007/s10681-021-02819-7)

Joukhadar R, Li Y, Thistlethwaite R, Forrest KL, Tibbits JF, Trethowan
R, Hayden MJ (2024). *Frontiers in Plant Science* 15:1337388.
[doi:10.3389/fpls.2024.1337388](https://doi.org/10.3389/fpls.2024.1337388)

Mulamba NN, Mock JJ (1978). *Egyptian Journal of Genetics and Cytology*
7:40-51.

Elston RC (1963). A weight-free index for ranking or selection with
respect to several traits at a time. *Biometrics* 19:85-97.

Pesek J, Baker RJ (1969). *Canadian Journal of Plant Science*
49:803-804. [doi:10.4141/cjps69-137](https://doi.org/10.4141/cjps69-137)

Rahimi M, Debnath S (2023). *Scientific Reports* 13:18977.
[doi:10.1038/s41598-023-46368-6](https://doi.org/10.1038/s41598-023-46368-6)

Smith HF (1936). *Annals of Eugenics* 7:240-250.

Yamada Y, Yokouchi K, Nishida A (1975). *Japanese Journal of Genetics*
50:33-41. [doi:10.1266/jjg.50.33](https://doi.org/10.1266/jjg.50.33)

## See also

[`evaluate_index()`](https://FAkohoue.github.io/DesiredGainR/reference/evaluate_index.md),
[`gain_feasibility()`](https://FAkohoue.github.io/DesiredGainR/reference/gain_feasibility.md),
[`run_dgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_dgsi.md)

## Examples

``` r
set.seed(1)
traits <- c("yield", "disease")
values <- matrix(
  c(stats::rnorm(40), stats::rnorm(40)),
  ncol = 2,
  dimnames = list(paste0("G", 1:40), traits)
)
G <- matrix(c(0.60, -0.15, -0.15, 0.40), 2, dimnames = list(traits, traits))
P <- matrix(c(1.10, -0.20, -0.20, 0.90), 2, dimnames = list(traits, traits))

fit <- selection_index(
  values, traits,
  method = "smith_hazel",
  G = G, P = P,
  economic_weights = c(yield = 1, disease = 0.5),
  lower_is_better = "disease", n_select = 4
)
fit
#> <desiredgainr_index>
#>   Method: smith_hazel 
#>   Candidates: 40  Traits: 2 
#>   Traits standardised: yes 
#>   Selected: 4 (10.0%), intensity 1.755
#>   Coefficients:
#>   yield disease 
#>  0.5646  0.2653 
#>   R_HI 0.7469  dH 1.3512  RE 0.9649
#>   CV_I undefined for a centred index
```
