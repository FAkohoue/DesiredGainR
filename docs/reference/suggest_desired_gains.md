# Suggest desired-gain directions from the current population

Gives a breeder who cannot state a complete desired-gain vector two
population-driven recommendations without inventing economic weights:

## Usage

``` r
suggest_desired_gains(
  setup,
  minimum_gains,
  lower_is_better = NULL,
  programme = list(),
  control = list()
)
```

## Arguments

- setup:

  A simulation setup from
  [`founder_population()`](https://FAkohoue.github.io/DesiredGainR/reference/founder_population.md).

- minimum_gains:

  Named non-negative vector of trait-specific minimum favourable gains,
  in estimated `G_target` genetic standard deviations.

- lower_is_better:

  Traits for which a reduction is favourable.

- programme:

  Named list of breeding-programme arguments accepted by
  [`simulate_selection_cycles()`](https://FAkohoue.github.io/DesiredGainR/reference/simulate_selection_cycles.md).
  Common entries are `mating_system`, `n_parents`, `n_crosses`, and
  `n_progeny_per_cross`.

- control:

  Named list of computational settings: `n_cycles` (5), `budget` (60),
  `n_initial` (`NULL`), `n_replicates` (50), `n_candidates` (2000),
  `screening_replicates` (50), `confirmation_replicates` (200),
  `confirmation_finalists` (5), `required_success_probability` (0.8),
  `search_starts` (3), `search_angle_tolerance` (10 degrees),
  `probability_level` (0.95), `uncertainty` (a named list), `model` (a
  named list), `seed` (42), `checkpoint` (`NULL`), and `verbose`
  (`TRUE`). The uncertainty list accepts `architecture_draws` (30),
  `covariance_draws` (0), `genetic_df`, and `residual_df`. A supported
  result remains explicitly conditional on the point covariance unless
  covariance draws and defensible degrees of freedom are supplied.
  Smaller simulation settings are useful for examples, not decisions.

## Value

An object of class `desiredgainr_gain_suggestion`. The two
recommendations contain a unit-norm desired-gain direction in favourable
genetic-SD space and a signed vector in original trait units. The object
also contains exact one-cycle feasibility, all simulated candidate
summaries, independent confirmation, and full search provenance.

## Details

1.  a **minimum-attainment direction**, chosen to maximise robust
    evidence that every trait reaches its own breeder-specified minimum;
    and

2.  a **maximum-balanced direction**, chosen to maximise the typical
    gain of the worst-responding trait in genetic standard deviations.

`minimum_gains` is always stated in favourable genetic standard
deviations calculated from `setup$G_target`, the biological
genetic-covariance estimate supplied by the breeder. `setup$G_realised`
is used only to diagnose whether the finite-QTL simulation reproduces
that target; using it for the units would make a breeder's threshold
depend on the simulator seed. Thus `c(yield = 1, disease = 0.5)` asks
for at least one estimated genetic SD of yield improvement and half an
SD of disease reduction when `lower_is_better = "disease"`. Values may
differ by trait and zero means non-decline, not omission. At least one
value must be positive.

The public call is deliberately short. `programme` contains arguments of
[`simulate_selection_cycles()`](https://FAkohoue.github.io/DesiredGainR/reference/simulate_selection_cycles.md)
that describe the breeding programme; `control` contains computational
settings. Unknown fields are rejected so a misspelling cannot silently
change a recommendation.

## Exact feasibility mathematics

For selection intensity \\i\\, genetic covariance \\G\\, phenotypic
covariance \\P\\, and index coefficients \\b\\, the response is \$\$R =
iGb / \sqrt{b^\mathsf{T}Pb}\$\$ and all attainable one-cycle responses
satisfy \$\$R^\mathsf{T}G^{-1}PG^{-1}R=i^2.\$\$ After orienting traits
and dividing each response by its estimated genetic SD, DesiredGainR
solves the strictly convex programme \$\$\min_z z^\mathsf{T}Bz \quad
\hbox{subject to}\quad z_j\ge m_j.\$\$ The exact active-set solver
checks the KKT conditions, which are necessary and sufficient because
\\B\\ is positive definite. It also solves the same problem with every
\\m_j=1\\; scaling its solution to the planned selection intensity gives
the mathematically greatest common lower bound attainable across traits
in one cycle.

## Multi-cycle recommendation

The analytical directions and a surrogate-assisted Pareto search are
evaluated through forward simulation. For each direction, joint success
is the event that every favourable cumulative response reaches its own
`minimum_gains` threshold at `control$n_cycles`. The table reports a
Jeffreys beta-binomial posterior probability and an exact, one-sided
Clopper-Pearson lower confidence bound with Bonferroni family-wise
coverage across all evaluated directions. No independence among
directions is required for that bound.

"Highest" has no unique meaning for several traits without preferences.
The least preferential definition used here is maximin: maximise the
worst-responding trait after all traits are expressed in genetic SD. The
recommendation is ranked by a simultaneous, distribution-free lower
confidence bound for the median worst-trait gain, obtained from binomial
order-statistic theory. It does not assume normal simulation outcomes.

Search and inference are separated in three stages. The adaptive search
discovers directions. Once that candidate set is locked, every direction
is rerun with new common random numbers; all confidence bounds and both
choices use this independent screening sample. Finally, the two locked
recommendations are rerun with a third seed stream. `confirmation` is
therefore an independent assessment of the reported vectors rather than
reuse of either optimisation or screening noise.

## References

Brown LD, Cai TT, DasGupta A (2001). Interval estimation for a binomial
proportion. *Statistical Science* 16:101–133.
[doi:10.1214/ss/1009213286](https://doi.org/10.1214/ss/1009213286)

Casella G, Berger RL (2002). *Statistical Inference*, 2nd ed. Duxbury.

Pesek J, Baker RJ (1969). Desired improvement in relation to selection
indices. *Canadian Journal of Plant Science* 49:803–804.
[doi:10.4141/cjps69-137](https://doi.org/10.4141/cjps69-137)

Yang W-N, Nelson BL (1991). Using common random numbers and control
variates in multiple-comparison procedures. *Operations Research*
39:583–591.
[doi:10.1287/opre.39.4.583](https://doi.org/10.1287/opre.39.4.583)

## See also

[`define_desired_gain_intervals()`](https://FAkohoue.github.io/DesiredGainR/reference/define_desired_gain_intervals.md),
[`optimize_desired_gains()`](https://FAkohoue.github.io/DesiredGainR/reference/optimize_desired_gains.md),
[`gain_feasibility()`](https://FAkohoue.github.io/DesiredGainR/reference/gain_feasibility.md)
