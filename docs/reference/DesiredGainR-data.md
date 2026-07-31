# Example breeding programme shipped with DesiredGainR

A simulated tropical maize breeding programme used throughout the
examples, tests and vignettes. The material is patterned on a real
programme structure but represents no real germplasm and no real trial.

## Usage

``` r
dgr_traits

dgr_G

dgr_P

dgr_candidates

dgr_gebv

dgr_history

dgr_hap1

dgr_hap2

dgr_map
```

## Format

`dgr_traits` is a data frame with six rows and seven columns:

- trait:

  Trait abbreviation.

- description:

  Trait name in full.

- unit:

  Measurement unit.

- direction:

  Either `"increase"` or `"decrease"`, stating which direction
  constitutes improvement. Traits marked `"decrease"` are the ones to
  pass to `lower_is_better`.

- heritability:

  Narrow-sense heritability, reproduced exactly by
  `diag(dgr_G) / diag(dgr_P)`.

- mean:

  Population mean on the original scale.

- genetic_sd:

  Additive genetic standard deviation. These span a factor of roughly
  one hundred and twenty, from ears per plant to plant height, which is
  what makes standardisation and matrix conditioning matter.

`dgr_G` is a six by six additive genetic variance-covariance matrix in
original trait units, named by trait.

`dgr_P` is a six by six phenotypic variance-covariance matrix in
original trait units, named by trait. It equals `dgr_G` plus a residual
covariance whose correlations deliberately differ from the genetic ones,
so that the two matrices are genuinely different objects rather than one
rescaled into the other.

`dgr_candidates` is a data frame with two hundred rows and eight
columns: `GenoID`, a `Family` label identifying twenty full-sib families
of ten, and one adjusted mean per trait on the original scale. The
realised genetic covariance of the underlying genetic values equals
`dgr_G` to numerical precision.

`dgr_gebv` is a data frame with two hundred rows and seven columns:
`GenoID` and one genomic estimated breeding value per trait. These
follow the prediction identity, in which the predictor and the
prediction error are orthogonal, so their variance is the genetic
variance scaled by the trait reliability rather than the full genetic
variance.

`dgr_history` is a data frame with two hundred rows and two columns,
`GenoID` and a logical `selected`, recording which forty candidates a
previous cycle retained.

The weight vector that produced the decision is attached as the
`generating_weights` attribute, in the standardised favourable-direction
trait space described by `generating_scale`. It is supplied so that the
recovery achieved by
[`retrospective_weights()`](https://FAkohoue.github.io/DesiredGainR/reference/retrospective_weights.md)
can be checked rather than merely asserted; an example whose answer
cannot be verified demonstrates nothing. The recovery is close but not
exact, because `b = P^-1 s` returns the linear rule that best reproduces
the observed differentials rather than the weights that generated them.

`dgr_hap1` and `dgr_hap2` are marker-by-candidate integer matrices coded
0 or 1, giving the allele carried on the first and second homologue.
Markers are grouped into blocks of fifteen that share a small pool of
founder haplotypes, so markers within a block travel together while
markers in different blocks segregate independently. This is the shape
[`founder_haplotypes()`](https://FAkohoue.github.io/DesiredGainR/reference/founder_haplotypes.md)
requires.

See `dgr_hap1`.

`dgr_map` is a data frame with one row per marker and three columns:
`variant_id`, `chromosome`, and `position_bp`. The column names match
the defaults of
[`founder_haplotypes()`](https://FAkohoue.github.io/DesiredGainR/reference/founder_haplotypes.md),
so it can be passed without further argument.

## Details

Earlier releases demonstrated the package on uncorrelated standard
normal variates, which made the genetic covariance matrix effectively an
identity matrix. That cannot illustrate a multi-trait selection index,
because the purpose of an index is to resolve the tension between
correlated traits measured on different scales. These data therefore
carry three properties the earlier examples lacked: (i) genuinely
antagonistic genetic correlations, (ii) trait scales spanning two orders
of magnitude, and (iii) one internally consistent population, in which
the markers actually explain the trait values.

The awkward pair is grain yield and anthesis date. Their genetic
correlation is positive, yet the objective requires yield to rise while
the cycle shortens, so the two cannot both be pushed freely. Grain yield
is also strongly and negatively correlated with the anthesis-silking
interval, which is the conventional drought-stress indicator.
Consequently
[`gain_feasibility()`](https://FAkohoue.github.io/DesiredGainR/reference/gain_feasibility.md)
returns informative answers on these data rather than declaring every
target attainable.

The six traits are grain yield (GY), plant height (PHT), anthesis date
(AD), the anthesis-silking interval (ASI), ears per plant (EPP), and
grey leaf spot severity (GLS).
