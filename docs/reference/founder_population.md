# Build an AlphaSimR founder population with a target genetic covariance

The founder genomes come from the breeder's phased marker data through
[`founder_haplotypes()`](https://FAkohoue.github.io/DesiredGainR/reference/founder_haplotypes.md),
whereas the trait architecture is calibrated to the genetic covariance
matrix the breeder has already estimated. Therefore the simulation
reproduces both the observed germplasm structure and the observed trait
variances and correlations, rather than an assumed demography.

## Usage

``` r
founder_population(
  founders,
  G,
  h2,
  residual_covariance = NULL,
  n_qtl_per_chromosome = 100L,
  n_markers_per_chromosome = NULL,
  dominance_degree = NULL,
  dominance_variance = NULL,
  heritability = c("narrow", "broad"),
  seed = 42L
)
```

## Arguments

- founders:

  An object from
  [`founder_haplotypes()`](https://FAkohoue.github.io/DesiredGainR/reference/founder_haplotypes.md).

- G:

  Genetic variance-covariance matrix, named by trait. Without dominance
  this is the additive covariance. With dominance it is the total
  **genotypic** covariance, because that is what a clonal programme
  selects on; see Details.

- h2:

  Named heritabilities, or a single value applied to every trait.
  Narrow- or broad-sense according to `heritability`.

- residual_covariance:

  Optional named residual covariance matrix. When supplied it is passed
  directly to AlphaSimR and takes precedence over the marginal `h2`
  values. This is the route used when propagating a sampled residual
  covariance; `h2` remains recorded as the nominal setup input.

- n_qtl_per_chromosome:

  Number of quantitative trait loci simulated per chromosome.

- n_markers_per_chromosome:

  Optional number of neutral markers per chromosome, used for
  coancestry. When supplied, the markers are held disjoint from the
  quantitative trait loci where the installed AlphaSimR supports it.
  Strongly recommended: without a panel, diversity is measured on all
  segregating sites.

- dominance_degree:

  Optional named mean degree of dominance per trait. Supplying it adds
  dominance to the simulated traits, which the clonal mating system
  requires.

- dominance_variance:

  Optional named variance of the **degree of dominance** across loci.
  This is AlphaSimR's `varDD` and is not the dominance genetic variance;
  the two are different quantities.

- heritability:

  Whether `h2` is narrow-sense (the default, a statement about breeding
  values) or broad-sense (a statement about genotypic values). They
  differ exactly when dominance is simulated, which is when a clonal
  programme is being represented, so it must be stated rather than
  assumed.

- seed:

  Random seed. The caller's random number generator state is restored on
  exit.

## Value

A list of class `desiredgainr_sim_setup` holding the AlphaSimR
simulation parameters, the founder population, and the calibration
targets.

## Details

Trait variances are taken from the diagonal of `G` and trait
correlations from the corresponding correlation matrix.

## Additive and clonal programmes calibrate differently

For a self-pollinated or cross-pollinated programme, `G` is the additive
genetic covariance and is passed straight to AlphaSimR.

For a clonal programme, supply `dominance_degree`. `G` is then the total
**genotypic** covariance, because the unit of selection is the clone and
dominance is inherited intact. This needs a calibration step that
earlier versions omitted. AlphaSimR's `var` argument sets the *additive*
variance, so with dominance present the dominance variance is added on
top and the realised genotypic variance exceeds the supplied target.
Since scaling every quantitative trait locus effect by \\c\\ scales both
the additive and the dominance variance by \\c^2\\, the correction is
exact in one pass: the trait is built, the realised genotypic variance
measured, and the additive target rescaled by the ratio.

The correlations are calibrated too, by the same fixed point. Rescaling
the variances alone would leave the realised genotypic correlations
emergent, so the setup would accept a full covariance matrix as its
target and reproduce only its diagonal. Dominance perturbs the
correlation structure as well, and that perturbation is a smooth,
near-identity function of the additive correlations supplied, so
correcting both and projecting the correction back onto the set of valid
correlation matrices recovers the whole of `G`.

`setup$G_realised` records what was achieved and
`setup$calibration_error` the largest remaining deviation in each of the
variances and correlations. Convergence is not guaranteed for every
target: a strongly antagonistic correlation structure at a high
dominance degree may not be attainable by any additive-plus-dominance
architecture. Failure warns and is recorded in
`setup$calibration_converged` rather than passing silently.

## Which heritability

`heritability = "narrow"` sets the error variance so that \\V_A / V_P\\
equals `h2`; `"broad"` targets \\V_G / V_P\\. They coincide without
dominance and diverge with it, so a clonal programme must say which is
meant. Selection in a clonal programme acts on genotypic value, making
broad sense usually the relevant one.

## Markers for coancestry

`n_markers_per_chromosome` creates a neutral marker panel, held disjoint
from the quantitative trait loci where the installed AlphaSimR supports
it. Diversity should not be measured on the loci under selection:
changing their frequencies is what genetic gain is, so a relationship
matrix built from them rises whenever selection succeeds, and a
desired-gain direction would be penalised for working.

## See also

[`founder_haplotypes()`](https://FAkohoue.github.io/DesiredGainR/reference/founder_haplotypes.md),
[`simulate_selection_cycles()`](https://FAkohoue.github.io/DesiredGainR/reference/simulate_selection_cycles.md)
