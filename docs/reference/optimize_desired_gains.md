# Search for the desired-gain direction giving the best multi-cycle outcome

A single-cycle response calculation cannot distinguish desired-gain
directions beyond what the achievable-response ellipsoid already states,
and
[`gain_feasibility()`](https://FAkohoue.github.io/DesiredGainR/reference/gain_feasibility.md)
provides that answer exactly and without simulation. Over several cycles
the ranking can change, however, because truncation selection erodes the
genetic variance that response depends on, drift accumulates in a finite
population, and an antagonistic correlation can drive a secondary trait
past an unacceptable level. This function therefore searches for the
direction that performs best once those processes are represented.

## Usage

``` r
optimize_desired_gains(
  setup,
  n_cycles = 5L,
  mode = c("pareto", "economic", "target", "constrained"),
  budget = 60L,
  n_initial = NULL,
  n_replicates = 3L,
  economic_weights = NULL,
  target_gains = NULL,
  focal_trait = NULL,
  gain_floors = NULL,
  include_diversity = TRUE,
  non_negative = TRUE,
  stability_tolerance = 0.05,
  n_candidates = 2000L,
  checkpoint = NULL,
  seed = 42L,
  verbose = TRUE,
  ...
)
```

## Arguments

- setup:

  An object from
  [`founder_population()`](https://FAkohoue.github.io/DesiredGainR/reference/founder_population.md).

- n_cycles:

  Number of selection cycles per evaluation.

- mode:

  Ranking mode. See Details.

- budget:

  Total number of simulated directions, including the initial design.

- n_initial:

  Size of the quasi-random initial design. Defaults to ten times the
  number of free parameters.

- n_replicates:

  Simulation replicates per direction.

- economic_weights:

  Named weights required by `mode = "economic"`.

- target_gains:

  Named absolute targets required by `mode = "target"`.

- focal_trait, gain_floors:

  Objective and constraints required by `mode = "constrained"`.

- include_diversity:

  Whether to treat diversity, measured as the negated mean relationship
  among selected parents, as an additional Pareto objective. Family
  balancing and coancestry control can cost a large share of nominal
  gain, so the trade-off is reported rather than assumed away.

- non_negative:

  Whether to restrict the search to directions seeking improvement in
  every trait.

- stability_tolerance:

  Relative tolerance defining the stability region.

- n_candidates:

  Size of the quasi-random pool over which the acquisition function is
  maximised at each iteration.

- checkpoint:

  Optional file path. When supplied, the accumulated evaluations are
  written after every simulation and reloaded automatically on a
  subsequent call, so that an interrupted search is not lost.

- seed:

  Random seed. The caller's random number generator state is restored on
  exit.

- verbose:

  Whether to report progress.

- ...:

  Further arguments passed to
  [`simulate_selection_cycles()`](https://FAkohoue.github.io/DesiredGainR/reference/simulate_selection_cycles.md),
  such as `mating_system`, `n_parents` and `lower_is_better`.

## Value

An object of class `desiredgainr_optimisation` containing the evaluated
directions and outcomes, the posterior-mean Pareto set, the recommended
region, and the search diagnostics.

## The simulation is the objective, never the surrogate

Each candidate direction is evaluated by running
[`simulate_selection_cycles()`](https://FAkohoue.github.io/DesiredGainR/reference/simulate_selection_cycles.md).
A Gaussian process is fitted to the accumulated results, but it is used
only to decide where the next simulation should be spent. It never
filters, screens, or replaces an evaluation. Because the acquisition
function retains an exploration term across the whole search space, no
region can be permanently excluded by the surrogate, and a direction
that the surrogate expects to be poor can still be visited.

## Search space

Only the direction of a desired-gain vector affects the index, since the
attainable magnitude is fixed by selection intensity. The domain is
therefore the unit sphere, with `p - 1` free parameters for `p` traits,
which is small enough to cover densely for the trait numbers breeders
use. By default the search is restricted to the non-negative orthant,
meaning improvement is sought in every trait; set `non_negative = FALSE`
to admit directions that deliberately concede ground on a trait.

## Ranking modes

- `"pareto"`:

  Returns the non-dominated set of multi-cycle outcomes and the
  direction generating each. This is the default, because it is the only
  mode that does not require the economic weights that breeders find
  hardest to state; choosing a point on a frontier is an easier
  judgement, and it is weight elicitation by revealed preference. The
  search uses randomised augmented-Chebyshev scalarisation, which
  converges on the frontier while needing only single-objective
  improvement.

- `"economic"`:

  Maximises `sum(economic_weights * cumulative_gain)`.

- `"target"`:

  Minimises the distance between the cumulative gain and a stated
  absolute target, reporting the per-trait shortfall.

- `"constrained"`:

  Maximises the gain in `focal_trait` subject to `gain_floors` on the
  remaining traits, weighting improvement by the probability that every
  floor is met.

## Noise and the frontier

Simulation output is stochastic. Every direction is evaluated with the
same sequence of seeds, so that comparisons between directions share
their stochasticity, which removes far more comparison variance than
increasing the replicate count. The reported frontier is computed from
Gaussian-process posterior means rather than from raw replicate
averages, because a frontier built from raw draws is populated by
fortunate runs.

## Interpreting the result

No single best direction is returned. The optimum is conditional on the
supplied genetic covariance matrix, the founder germplasm and the
programme parameters, so the result reports a stability region: every
direction whose posterior outcome lies within `stability_tolerance` of
the best. Choose within that region on other grounds.

## See also

[`simulate_selection_cycles()`](https://FAkohoue.github.io/DesiredGainR/reference/simulate_selection_cycles.md),
[`gain_feasibility()`](https://FAkohoue.github.io/DesiredGainR/reference/gain_feasibility.md)
