# Simulate recurrent selection under a fixed desired-gain direction

Runs a recurrent-selection programme forward for a stated number of
cycles, selecting each cycle on a desired-gain index built from the
supplied direction, and records the cumulative genetic gain and the loss
of genetic variation. The purpose is to compare desired-gain directions,
because a direction that maximises response in the first cycle may
exhaust the variance that response depends on within a few more.

## Usage

``` r
simulate_selection_cycles(
  setup,
  desired_gains,
  n_cycles = 5L,
  mating_system = c("self", "outcross", "clonal"),
  n_parents = 20L,
  n_crosses = 50L,
  n_progeny_per_cross = 10L,
  n_selfing_generations = 3L,
  use_doubled_haploids = FALSE,
  lower_is_better = NULL,
  reestimate_index = TRUE,
  n_clonal_replicates = 1L,
  n_threads = 1L,
  seed = 42L,
  prediction = list(method = "phenotype")
)
```

## Arguments

- setup:

  An object from
  [`founder_population()`](https://FAkohoue.github.io/DesiredGainR/reference/founder_population.md).

- desired_gains:

  Named desired-gain direction. Only the direction matters, because the
  attainable magnitude is fixed by selection intensity.

- n_cycles:

  Number of selection cycles.

- mating_system:

  One of `"self"`, `"outcross"`, or `"clonal"`.

- n_parents:

  Number of parents recycled each cycle. Fewer parents increase
  short-term gain and reduce long-term potential.

- n_crosses:

  Number of crosses made each cycle.

- n_progeny_per_cross:

  Progeny evaluated per cross.

- n_selfing_generations:

  Selfing generations before evaluation, used by `"self"` when
  `use_doubled_haploids` is `FALSE`.

- use_doubled_haploids:

  Whether `"self"` produces doubled haploids instead of advancing by
  selfing.

- lower_is_better:

  Traits for which smaller values are favourable.

- reestimate_index:

  Whether to rebuild the index from each cycle's own simulated data. See
  Details.

- n_clonal_replicates:

  Ramets evaluated per genotype in a clonal programme. A clonal trial
  phenotypes several copies of each genotype and averages them, so its
  selection is more accurate than a single plot; leaving this at 1
  understates the response a clonal programme achieves. Only meaningful
  with `mating_system = "clonal"`.

- n_threads:

  Threads AlphaSimR may use. The default of 1 is deliberate: above one,
  the order in which random numbers are consumed is not guaranteed, so a
  run cannot be reproduced exactly from its seed. Values above 1 warn.

- seed:

  Random seed. The caller's random number generator state is restored on
  exit, and `setup$SP` is deep-cloned so that the caller's `SimParam` is
  not advanced by the simulation.

- prediction:

  A named list describing the selection criterion. The default
  `list(method = "phenotype")` preserves phenotypic selection.
  `method = "rrblup"` uses leakage-free cross-fitted RR-BLUP predictions
  and requires a marker panel; optional fields are `folds` (5),
  `max_iter` (10000), `update_training` (`TRUE`) and `max_training`
  (`NULL`).

## Value

An object of class `desiredgainr_simulation` containing a per-cycle
table of genetic means, genetic variances, mean relationship, effective
population size, prediction accuracy and prediction-calibration slope,
together with the cumulative gain per trait.

## What the simulation adds beyond the single-cycle formula

A single-cycle response prediction cannot distinguish desired-gain
directions beyond what the achievable-response ellipsoid already states,
and
[`gain_feasibility()`](https://FAkohoue.github.io/DesiredGainR/reference/gain_feasibility.md)
gives that answer exactly and without simulation. Simulation is
informative only because it represents the processes that make cycle
five differ from cycle one: (i) the reduction and reshaping of genetic
variance caused by truncation selection, (ii) drift and the accumulation
of relatedness in a finite population, (iii) the point at which an
antagonistic correlation drives a secondary trait past an unacceptable
level, and (iv) the compounding of estimation error when the index is
rebuilt each cycle.

## What a cycle means

This is the definition to check before interpreting any number below.

**Cycle 0** is the founder population, before any selection. It exists
so that the cycle 1 row has something to be a response *to*.

**Cycle \\t\\** records the population produced by selecting parents
from the cycle \\t-1\\ candidates and crossing them. Every row from
cycle 1 onward is therefore a transmitted selection response, and
`n_cycles = 1` gives exactly one such response.

Earlier versions measured the candidates rather than their selected
descendants. Because random mating does not shift a population mean,
that made the cycle 1 gain zero in expectation *for every desired-gain
direction*, so a one-cycle run could not distinguish directions at all.
Results from before this change should be re-run.

Genetic mean, variance, relatedness and inbreeding are all measured on
the same population, the cycle's response. `parent_inbreeding`
additionally reports the selected parents, which is a different and
smaller group.

## Mating systems

- `"self"`:

  Self-pollinated line development. Selected parents are intercrossed,
  and the resulting families are advanced by selfing or doubled haploidy
  before evaluation. Recombination therefore releases variation slowly.

- `"outcross"`:

  Recurrent selection in a random-mating population. Selected parents
  are intercrossed each cycle, so half the disequilibrium generated by
  selection decays per generation.

- `"clonal"`:

  Clonally propagated crops, with the sexual and clonal phases kept
  distinct. Recombination happens once per cycle, when selected parents
  are crossed to raise a seedling generation; every evaluation stage
  after that is a copy of a seedling, so no further meiosis occurs and
  dominance is transmitted intact. Selection therefore acts on total
  genotypic value rather than breeding value, and response is a
  broad-sense quantity. Requires a setup built with `dominance_degree`,
  and `G` should be the genotypic covariance. See
  [`founder_population()`](https://FAkohoue.github.io/DesiredGainR/reference/founder_population.md)
  for how the additive and dominance components are calibrated to that
  target.

## Index re-estimation

With `reestimate_index = TRUE`, the covariance of the observed selection
criterion is re-estimated each cycle. Phenotypic selection retains the
breeder-supplied `G_target`; it never estimates genetic covariance from
hidden simulated breeding values. RR-BLUP uses the calibrated-prediction
identity `Cov(A, Ahat) = Var(Ahat)` and reports held-out prediction
accuracy as a diagnostic. With `FALSE`, cycle-one coefficients are
reused. Thus no selection decision has access to simulation truth.

## See also

[`founder_population()`](https://FAkohoue.github.io/DesiredGainR/reference/founder_population.md),
[`gain_feasibility()`](https://FAkohoue.github.io/DesiredGainR/reference/gain_feasibility.md)
