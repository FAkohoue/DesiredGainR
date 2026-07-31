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
  n_qtl_per_chromosome = 100L,
  dominance_degree = NULL,
  dominance_variance = NULL,
  seed = 42L
)
```

## Arguments

- founders:

  An object from
  [`founder_haplotypes()`](https://FAkohoue.github.io/DesiredGainR/reference/founder_haplotypes.md).

- G:

  Genetic variance-covariance matrix, named by trait.

- h2:

  Named narrow-sense heritabilities, or a single value applied to every
  trait.

- n_qtl_per_chromosome:

  Number of quantitative trait loci simulated per chromosome.

- dominance_degree:

  Optional named mean degree of dominance per trait. Supplying it adds
  dominance to the simulated traits, which the clonal mating system
  requires.

- dominance_variance:

  Optional named variance of the degree of dominance.

- seed:

  Random seed. The caller's random number generator state is restored on
  exit.

## Value

A list of class `desiredgainr_sim_setup` holding the AlphaSimR
simulation parameters, the founder population, and the calibration
targets.

## Details

Trait variances are taken from the diagonal of `G` and trait
correlations from the corresponding correlation matrix. For a clonal
programme, supply `dominance_degree`, and `G` should then be the
**genotypic** covariance, because the unit of selection is the clone and
dominance is therefore inherited intact. For a self-pollinated or
cross-pollinated programme, `G` should be the additive genetic
covariance.

## See also

[`founder_haplotypes()`](https://FAkohoue.github.io/DesiredGainR/reference/founder_haplotypes.md),
[`simulate_selection_cycles()`](https://FAkohoue.github.io/DesiredGainR/reference/simulate_selection_cycles.md)
