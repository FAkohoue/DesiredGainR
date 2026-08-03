# Propagate covariance uncertainty into a desired-gain recommendation

[`optimize_desired_gains()`](https://FAkohoue.github.io/DesiredGainR/reference/optimize_desired_gains.md)
conditions on one genetic covariance estimate. Its frontier therefore
carries Monte Carlo error but treats the covariance as known, so a
direction that is only best because of how `G` happened to be estimated
is indistinguishable from one that is robustly best.

## Usage

``` r
propagate_covariance_uncertainty(
  setup,
  directions,
  genetic_df,
  P = NULL,
  residual_df = NULL,
  n_covariance_draws = 20L,
  n_replicates = 3L,
  n_cycles = 5L,
  include_diversity = FALSE,
  allow_unverified_diversity = FALSE,
  rank_weights = NULL,
  verbose = FALSE,
  seed = 42L,
  ...
)
```

## Arguments

- setup:

  A setup from
  [`founder_population()`](https://FAkohoue.github.io/DesiredGainR/reference/founder_population.md),
  used for the founder genomes, heritabilities and marker panel. Its
  `G_target` is the point estimate around which the draws are taken.

- directions:

  An object returned by
  [`optimize_desired_gains()`](https://FAkohoue.github.io/DesiredGainR/reference/optimize_desired_gains.md),
  or a matrix of desired-gain directions in original trait units, one
  per row, with columns named by trait. Passing the optimisation object
  is safer when its search used genetic- or phenotypic-SD units because
  the stored trait-unit simulation directions are used automatically.

- genetic_df:

  Degrees of freedom for the genetic covariance. This governs the width
  of everything reported and counts independent genetic units, not
  plots.

- P:

  Optional phenotypic covariance. When supplied, the residual is
  resampled too and each draw is admissible by construction.

- residual_df:

  Degrees of freedom for the residual covariance.

- n_covariance_draws:

  Number of covariance draws.

- n_replicates:

  Simulation replicates within each draw.

- n_cycles:

  Selection cycles per evaluation.

- include_diversity:

  Whether to append diversity as an objective.

- allow_unverified_diversity:

  Experimental override for a marker panel whose QTL disjointness could
  not be verified. Known overlap is rejected.

- rank_weights:

  Optional named scalarisation weights for every objective. Required to
  compute rank churn; when absent, rank churn is reported as not
  requested rather than averaging objectives in their raw units.

- verbose:

  Whether to report progress.

- seed:

  Random seed. The caller's stream is restored on exit.

- ...:

  Further arguments passed to
  [`simulate_selection_cycles()`](https://FAkohoue.github.io/DesiredGainR/reference/simulate_selection_cycles.md).

## Value

An object of class `desiredgainr_covariance_uncertainty`.

## Details

This function repeats the evaluation of a fixed set of directions across
draws from the sampling distribution of the covariance, rebuilding the
founder population for each draw, and reports how far the recommendation
moves.

## Why the founders are rebuilt

The genetic covariance is not a matrix the simulation multiplies by. It
is the architecture the trait is built from, so a different `G` means
different quantitative trait locus effects and a different founder
population. Each draw therefore calls
[`founder_population()`](https://FAkohoue.github.io/DesiredGainR/reference/founder_population.md)
again. That is what makes this expensive, and why the directions are
supplied rather than searched: running a full Bayesian optimisation
inside every covariance draw would multiply an already slow search by
`n_covariance_draws`.

The intended workflow is to run
[`optimize_desired_gains()`](https://FAkohoue.github.io/DesiredGainR/reference/optimize_desired_gains.md)
once and pass the complete optimisation object here. The trait-unit
directions actually used by the simulator are then selected
automatically; this avoids mistaking genetic-SD search coordinates for
original trait units.

## The two sources of uncertainty are separated

Within one covariance draw, replicate simulations differ by Monte Carlo
error alone. Across draws, they differ by both. The variance
decomposition reported as `variance_components` splits them, because
they have different remedies: Monte Carlo error falls with more
replicates, covariance uncertainty does not and can only be reduced by
estimating `G` better.

If the covariance component dominates, running more simulation
replicates is wasted effort.

## What to read

When `rank_weights` is supplied, `rank_churn` is the Spearman
correlation between the declared scalar ranking under each draw and
under the point estimate. Raw objectives are never averaged across
incompatible units. `frontier_membership` is the proportion of draws in
which each direction was non-dominated. A direction that is on the
frontier under the point estimate but in only a third of the draws is
not a robust recommendation.

## See also

[`optimize_desired_gains()`](https://FAkohoue.github.io/DesiredGainR/reference/optimize_desired_gains.md),
[`index_uncertainty()`](https://FAkohoue.github.io/DesiredGainR/reference/index_uncertainty.md),
[`draw_covariance_pairs()`](https://FAkohoue.github.io/DesiredGainR/reference/draw_covariance_pairs.md)

## Examples

``` r
# \donttest{
# See the "Uncertainty" section of the simulation vignette; this example
# needs AlphaSimR and takes several minutes.
# }
```
