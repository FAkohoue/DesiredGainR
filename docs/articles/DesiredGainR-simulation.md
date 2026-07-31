# Comparing objectives over several cycles

## 1. What this layer is for, and what it is not for

The purpose is narrow and deliberately kept so: **to compare
desired-gain directions over several breeding cycles.** It is not a
crossing-plan tool, and parent, cross and mating allocation are outside
the package by design.

The heavy chunks below are marked `eval = FALSE` because a realistic run
takes minutes to hours. Output blocks shown as fixed text illustrate the
shape of each object; your numbers will depend on your founders and
parameters.

`AlphaSimR` is required and lives in `Suggests`. Every function fails
with an explicit installation message when it is absent.

------------------------------------------------------------------------

## 2. Why simulate at all

For any linear index the response is \\\mathbf{R} =
i\\\mathbf{G}\mathbf{b}/\sqrt{\mathbf{b}^\mathsf{T}\mathbf{P}\mathbf{b}}\\,
so attainable responses lie on the ellipsoid

\\\mathbf{R}^\mathsf{T}\mathbf{G}^{-1}\mathbf{P}\mathbf{G}^{-1}\mathbf{R}
= i^{2}.\\

[`gain_feasibility()`](https://FAkohoue.github.io/DesiredGainR/reference/gain_feasibility.md)
solves that exactly, without simulation. **A single-cycle calculation
therefore cannot distinguish desired-gain directions beyond what the
ellipsoid already states**, and any simulator that merely re-applies the
one-cycle formula \\T\\ times returns \\T\\ times the one-cycle answer.

Multi-cycle simulation earns its cost only by representing what makes
cycle five differ from cycle one:

1.  **Bulmer effect.** Truncation selection removes variance
    anisotropically, by \\k = i(i-x)\\ along the index direction. A
    direction that hammers one axis exhausts the variance that axis
    depends on.
2.  **Finite \\N_e\\.** Drift, and the fact that the most aggressive
    index usually erodes effective population size fastest.
3.  **Antagonistic correlations.** A secondary trait can cross an
    unacceptable threshold in cycle three, which no single-cycle
    response reveals.
4.  **Compounding estimation error**, when the index is rebuilt each
    cycle from that cycle’s own data.

If none of those matter for your question, use
[`gain_feasibility()`](https://FAkohoue.github.io/DesiredGainR/reference/gain_feasibility.md)
and stop.

------------------------------------------------------------------------

## 3. Founders come from your germplasm, not from a simulated genome

This is the design decision that determines whether the results mean
anything.

A simulator that generates founder genomes from a coalescent model
produces linkage disequilibrium, allele frequencies and population
structure belonging to an assumed demography. Optimising a desired-gain
direction against that is optimising for a fictional population.
DesiredGainR therefore requires the breeder’s own phased marker data.

### 3.1 The input contract

``` r
founders <- founder_haplotypes(dgr_hap1, dgr_hap2, dgr_map)
founders
#> <desiredgainr_founders>
#>   200 individuals, 450 variants, 3 chromosomes, diploid
#>   Map: converted from physical position at 1e-08 Morgan per base pair 
#>   Missing calls before resolution: 0.000% (policy: error)
```

`hap1` and `hap2` are marker-by-individual matrices coded 0, 1 or `NA`,
giving the allele carried on the first and second homologue. Row names
are variant identifiers, column names individual identifiers, and both
matrices must carry identical dimnames in identical order.

``` r
dim(dgr_hap1)
#> [1] 450 200
dgr_hap1[1:3, 1:4]
#>       CAND001 CAND002 CAND003 CAND004
#> M0001       0       0       1       0
#> M0002       0       0       1       1
#> M0003       0       1       1       1
head(dgr_map, 3)
#>   variant_id chromosome position_bp
#> 1      M0001          1       1e+06
#> 2      M0002          1       2e+06
#> 3      M0003          1       3e+06
```

### 3.2 Why dosage is refused

A haplotype is one chromosome copy, so at a biallelic locus it carries a
single allele rather than a count. Genotype dosage coded 0, 1, 2 is the
*sum* of the homologues:

\\\text{dosage} = \text{hap}\_1 + \text{hap}\_2.\\

It records how many alternative alleles an individual carries but not
which copy carries them, so it cannot say which alleles lie in *cis*
across loci. That co-occurrence is the linkage disequilibrium the
simulation exists to represent, and assigning heterozygous calls to
homologues at random would destroy it.

``` r
dosage <- dgr_hap1 + dgr_hap2
founder_haplotypes(dosage, dgr_hap2, dgr_map)
#> Error:
#> ! hap1 must contain only 0, 1 and NA, but also contains: 2. A haplotype is a single chromosome copy, so at a biallelic locus it carries one allele rather than an allele count. A 0/1/2 matrix is genotype dosage, which is the sum of the homologues and therefore discards phase; see haplotypes_from_inbred_dosage() when the material is highly inbred.
```

### 3.3 The exception: inbred material

In a diploid inbred line a dosage of 0 or 2 is unambiguous, so phase can
be derived without external phasing. Only heterozygous and missing calls
need a decision.

Measure the rates before converting, because a low overall level can
conceal a few badly affected individuals or markers:

``` r
diagnostics <- dosage_diagnostics(dosage)
diagnostics
#> <desiredgainr_dosage_diagnostics>
#>   450 variants x 200 individuals
#>   Heterozygous calls: 50.464%   Missing calls: 0.000%
#>   Per-individual heterozygosity, maximum 56.000% (CAND019)
#>   Per-marker heterozygosity, maximum 89.500%
```

``` r
converted <- haplotypes_from_inbred_dosage(
  dosage,
  heterozygous_policy = "drop_variant",
  missing_policy = "drop_variant"
)
founders <- founder_haplotypes(converted$hap1, converted$hap2, dgr_map)
```

**No threshold is imposed on residual heterozygosity**, because no
universally appropriate level exists: it depends on the generation,
mating history, crop, population type, genotyping error rate and
quality-control procedure. The rates are measured and reported instead,
and heterozygous calls are never resolved silently. The four policies
are `"error"`, `"drop_variant"`, `"drop_individual"` and `"mask"`. None
assigns a heterozygote to a homologue.

Use this route only for highly inbred diploid material. For outcrossing
or clonal germplasm, phase externally.

### 3.4 Diploid only

The contract is diploid and biallelic, and a polyploid dataset raises an
explanatory error rather than proceeding. Adding homologue matrices
would not make the simulation polyploid: autopolyploid and allopolyploid
meiosis require an explicit crop-specific pairing and recombination
model, including multivalent formation and double reduction. Applying a
diploid meiotic model to polyploid homologues would produce confidently
wrong results.

------------------------------------------------------------------------

## 4. Calibrating the trait architecture

Genome structure comes from your markers; trait variances and
correlations are calibrated to the covariance matrix you already
estimated.

``` r
setup <- founder_population(
  founders,
  G = dgr_G,
  h2 = stats::setNames(dgr_traits$heritability, traits),
  n_qtl_per_chromosome = 100L,
  seed = 42L
)
setup
```

    #> <desiredgainr_sim_setup>
    #>   Founders: 200 individuals, 3 chromosomes, 100 QTL per chromosome
    #>   Traits: GY, PHT, AD, ASI, EPP, GLS
    #>   Dominance simulated: no

Trait variances are taken from \\\operatorname{diag}(\mathbf{G})\\ and
correlations from the corresponding correlation matrix, so the
simulation reproduces the observed germplasm structure *and* the
observed trait covariances.

**For a clonal programme**, supply `dominance_degree`, and `G` should
then be the **genotypic** rather than the additive covariance, because
the unit of selection is the clone and dominance is inherited intact.

``` r
setup_clonal <- founder_population(
  founders,
  G = dgr_G,                       # genotypic covariance for a clonal crop
  h2 = stats::setNames(dgr_traits$heritability, traits),
  dominance_degree = stats::setNames(rep(0.3, length(traits)), traits),
  seed = 42L
)
```

------------------------------------------------------------------------

## 5. Simulating one objective

``` r
simulation <- simulate_selection_cycles(
  setup,
  desired_gains = c(GY = 1.0, PHT = 0.4, AD = 0.6,
                    ASI = 0.5, EPP = 0.4, GLS = 0.6),
  n_cycles = 5L,
  mating_system = "outcross",
  n_parents = 20L,
  n_crosses = 50L,
  n_progeny_per_cross = 10L,
  lower_is_better = lower_is_better,
  reestimate_index = TRUE,
  seed = 11L
)
simulation
```

    #> <desiredgainr_simulation>
    #>   outcross system, 5 cycles, 20 parents recycled
    #>   Index re-estimated each cycle: yes
    #>   Cumulative genetic gain:
    #>       GY     PHT      AD     ASI     EPP     GLS
    #>    1.842  -8.031  -1.226  -0.244  -0.019  -0.487
    #>   Final mean relationship among parents: 0.0913

### 5.1 The three mating systems

| System | Crops | Cycle structure |
|----|----|----|
| `"self"` | Wheat, rice, common bean, cowpea | Selected parents intercrossed, families advanced by selfing or doubled haploidy before evaluation. Recombination releases variance slowly. |
| `"outcross"` | Maize, sorghum, pearl millet | Random mating each cycle, so half the selection-induced disequilibrium decays per generation. |
| `"clonal"` | Cassava, sweetpotato, banana, potato | Selection acts on **total genetic value**, since the clone inherits dominance intact. Requires a setup built with `dominance_degree`. |

### 5.2 What each cycle records

``` r
simulation$cycle_table
```

    #>    cycle  trait genetic_mean cumulative_gain genetic_variance
    #> 1:     1     GY        0.402           0.402            0.548
    #> 2:     1    PHT       -1.907          -1.907          139.412
    #> ...
    #>    mean_relationship effective_size selection_intensity n_evaluated
    #>                0.021           23.9               1.755         500

Gain and the loss of diversity are recorded together deliberately. A
direction that maximises cumulative gain while collapsing effective
population size is not obviously preferable, and the trade-off should be
visible rather than inferred.

### 5.3 The re-estimation switch

`reestimate_index = TRUE`, the default, rebuilds the index from each
cycle’s own simulated data. That is what a breeding programme actually
does, and it propagates estimation error across cycles. Setting it
`FALSE` reuses the cycle one coefficients and isolates the effect of the
desired-gain direction alone.

**Run both.** A large divergence means the recommendation is sensitive
to estimation error rather than to the objective, and the breeder should
be told so.

### 5.4 The selection rule is deliberately simple

Truncation on the index, and nothing more. Optimal contribution
selection, family quotas and coancestry constraints are the province of
a mating-design tool. Duplicating them here would mean maintaining two
simulators with different genetic engines.

------------------------------------------------------------------------

## 6. Searching for the best objective

### 6.1 The search space

Only the *direction* of \\\mathbf{d}\\ affects the index, because the
attainable magnitude is fixed by selection intensity. The domain is
therefore the unit sphere,

\\\\\mathbf{d}\\\_2 = 1,\\

with \\p-1\\ free parameters for \\p\\ traits. For the three to eight
traits breeders use, that is small enough to cover densely. By default
the search is restricted to the non-negative orthant, meaning
improvement is sought in every trait; set `non_negative = FALSE` to
admit directions that deliberately concede ground somewhere.

### 6.2 The simulation is the objective function

``` r
optimisation <- optimize_desired_gains(
  setup,
  n_cycles = 5L,
  mode = "pareto",
  budget = 60L,
  n_initial = 30L,
  n_replicates = 3L,
  include_diversity = TRUE,
  mating_system = "outcross",
  n_parents = 20L, n_crosses = 50L, n_progeny_per_cross = 10L,
  lower_is_better = lower_is_better,
  checkpoint = "search_state.rds",
  seed = 42L
)
optimisation
```

    #> <desiredgainr_optimisation>
    #>   Mode: pareto   Cycles: 5   Evaluations: 60 (3 replicates each)
    #>   Objectives: GY, PHT, AD, ASI, EPP, GLS, diversity
    #>   Non-dominated directions: 9 of 60
    #>   No single direction is recommended; choose a point on the frontier.

A Gaussian process is fitted to the accumulated results, but **it
decides only where the next simulation should be spent.** It never
filters, screens, or substitutes for an evaluation, and because the
acquisition function retains an exploration term across the whole
domain, no region can be permanently excluded by the surrogate.

An earlier design used a deterministic infinitesimal-model recursion as
the surrogate’s prior mean. That was rejected: it would embed a second
genetic model whose assumptions could fail in precisely the situations
where simulation is most needed.

### 6.3 The four ranking modes

| Mode | Objective | Required arguments |
|----|----|----|
| `"pareto"` | Non-dominated set of multi-cycle outcomes | none |
| `"economic"` | Maximise \\\mathbf{w}^\mathsf{T}\mathbf{R}\_{\text{cum}}\\ | `economic_weights` |
| `"target"` | Minimise distance to stated absolute targets | `target_gains` |
| `"constrained"` | Maximise a focal trait subject to floors | `focal_trait`, `gain_floors` |

`"pareto"` is the default because it is the only mode that does not
require the economic weights breeders find hardest to state. Choosing a
point on a frontier is an easier judgement than stating weights, and it
is weight elicitation by revealed preference.

The Pareto search uses randomised augmented-Chebyshev scalarisation,
which converges on the frontier while requiring only single-objective
improvement. A fresh weight vector is drawn each iteration so the search
sweeps the frontier rather than one corner of it.

### 6.4 Noise, and why the frontier is smoothed

Simulation output is stochastic. Two properties matter for
interpretation.

**Common random numbers.** Every direction is evaluated with the same
sequence of seeds, so comparisons between directions share their
stochasticity. This removes far more comparison variance than increasing
the replicate count.

**The frontier comes from posterior means, not raw draws.** A frontier
built from raw replicate averages is populated by fortunate runs. The
reported non-dominated set is computed from Gaussian-process posterior
means.

``` r
optimisation$pareto_set[, c("d_GY", "d_ASI", "posterior_GY", "posterior_ASI")]
```

### 6.5 No single best direction is returned

For the scalar modes the result reports a **stability region**: every
direction whose posterior outcome lies within `stability_tolerance` of
the best.

``` r
optimisation$stability_proportion
optimisation$recommended_direction
```

The optimum is conditional on the supplied genetic covariance, the
founder germplasm and the programme parameters. Reporting one direction
to three decimal places would imply a precision the analysis does not
have.

### 6.6 Budget and resumption

A realistic budget takes hours. Supplying `checkpoint` writes the
accumulated evaluations after every simulation and reloads them
automatically on a later call, so an interruption does not discard the
work.

``` r
# Re-running the same call continues from where it stopped.
optimisation <- optimize_desired_gains(
  setup, n_cycles = 5L, mode = "pareto", budget = 120L,
  checkpoint = "search_state.rds", seed = 42L, ...
)
```

------------------------------------------------------------------------

## 7. What the results are conditional on

State these alongside any recommendation.

| Assumption | Consequence if wrong |
|----|----|
| Founder panel represents the breeding population | The optimised direction suits a different population |
| \\\mathbf{G}\\ is the covariance of the target population | Trait architecture is miscalibrated throughout |
| Quantitative trait loci are many and of small effect | The infinitesimal approximation degrades |
| The programme parameters match practice | Intensity and drift are misrepresented |
| Truncation is the operative selection rule | Coancestry-controlled programmes will differ |

The simulation is a decision-support tool, not a prediction. Its value
lies in ranking directions under stated assumptions, and it should be
reported with those assumptions attached.

------------------------------------------------------------------------

## 8. Function reference for this layer

| Function | Role |
|----|----|
| [`dosage_diagnostics()`](https://FAkohoue.github.io/DesiredGainR/reference/dosage_diagnostics.md) | Measure heterozygosity and missingness before converting |
| [`haplotypes_from_inbred_dosage()`](https://FAkohoue.github.io/DesiredGainR/reference/haplotypes_from_inbred_dosage.md) | Derive phase from inbred dosage, under an explicit policy |
| [`founder_haplotypes()`](https://FAkohoue.github.io/DesiredGainR/reference/founder_haplotypes.md) | Validate and package phased founders |
| [`founder_population()`](https://FAkohoue.github.io/DesiredGainR/reference/founder_population.md) | Build the AlphaSimR population, calibrated to \\\mathbf{G}\\ |
| [`simulate_selection_cycles()`](https://FAkohoue.github.io/DesiredGainR/reference/simulate_selection_cycles.md) | Run one desired-gain direction forward |
| [`optimize_desired_gains()`](https://FAkohoue.github.io/DesiredGainR/reference/optimize_desired_gains.md) | Search directions, surrogate-assisted |

------------------------------------------------------------------------

## 9. References

- Bulmer MG (1971). The effect of selection on genetic variability. *The
  American Naturalist* **105**:201-211.
- Gaynor RC, Gorjanc G, Hickey JM (2021). AlphaSimR: an R package for
  breeding program simulations. *G3 Genes\|Genomes\|Genetics*
  **11**:jkaa017. <https://doi.org/10.1093/g3journal/jkaa017>
- Knowles J (2006). ParEGO: a hybrid algorithm with on-line landscape
  approximation for expensive multiobjective optimization problems.
  *IEEE Transactions on Evolutionary Computation* **10**:50-66.
- Vieira RA, Nogueira APO, Fritsche-Neto R (2025). Optimizing the
  selection of quantitative traits in plant breeding using simulation.
  *Frontiers in Plant Science* **16**:1495662.
  <https://doi.org/10.3389/fpls.2025.1495662>
