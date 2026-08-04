# Working with other breeding software

## 1. Where DesiredGainR enters the analysis

DesiredGainR begins after genetic evaluation. Mixed-model or
genomic-prediction software supplies candidate breeding values and,
ideally, the fitted multi-trait covariance matrices. DesiredGainR then
defines and tests the objective, builds selection indices, ranks
candidates and reports uncertainty.

It does not analyse field designs, estimate marker effects, choose the
final parent set, perform optimal contribution selection (OCS), perform
optimum cross selection or allocate matings. Those downstream decisions
belong to [HapBlockR](https://github.com/FAkohoue/HapBlockR).

The division is intentional and one directional:

| Breeding question | Package responsible |
|----|----|
| What traits and response direction define the objective? | DesiredGainR |
| Is that objective attainable at the planned selection intensity? | DesiredGainR |
| What coefficients, scores and expected responses follow from the objective? | DesiredGainR |
| How stable is the ranking to uncertain weights or covariance matrices? | DesiredGainR |
| Which candidates should become parents under merit and relationship constraints? | HapBlockR |
| How much should each parent contribute under OCS? | HapBlockR |
| Which pairs should be crossed, and how should matings be allocated under optimum cross selection? | HapBlockR |

HapBlockR calls DesiredGainR. DesiredGainR does not call HapBlockR. This
avoids a circular dependency and keeps the scientific responsibilities
auditable.

------------------------------------------------------------------------

## 2. Minimum hand-off contract

Transfer the following objects explicitly:

1.  One stable candidate identifier used by every table and matrix.
2.  One GEBV or breeding-value column per trait.
3.  Trait names, units and favourable directions.
4.  Additive genetic covariance \\\mathbf{G}\\.
5.  Phenotypic or index-variable covariance \\\mathbf{P}\\.
6.  Prediction-error covariance or uncertainty draws when available.
7.  Population, environment and model provenance.

Never match candidates or traits by position alone.

``` r
stopifnot(!anyDuplicated(dgr_gebv$GenoID))
stopifnot(all(traits %in% names(dgr_gebv)))
stopifnot(identical(rownames(dgr_G), traits))
stopifnot(identical(colnames(dgr_G), traits))
```

------------------------------------------------------------------------

## 3. Import covariance from upstream software

[`import_covariance()`](https://FAkohoue.github.io/DesiredGainR/reference/import_covariance.md)
accepts matrices extracted from ASReml-R, `sommer`, `breedR`, BGLR or a
stage-wise workflow. The caller extracts the scientifically correct
model term. DesiredGainR validates its structure, orders the traits and
records its origin.

``` r
imported_G <- import_covariance(
  dgr_G[rev(traits), rev(traits)],
  trait_cols = traits,
  source = "matrix",
  estimand = "additive genetic covariance",
  P = dgr_P
)
imported_G
#> <desiredgainr_imported_covariance>
#>   Source: matrix 
#>   Estimand: additive genetic covariance 
#>   Traits: GY, PHT, AD, ASI, EPP, GLS 
#>   Condition number: 1.857e+04   Positive definite: yes
```

Because upstream interfaces change, passing the extracted covariance is
more stable than asking DesiredGainR to interpret an entire fitted model
object. The extraction script should remain with the analysis so the
chosen random term and covariance structure can be audited.

------------------------------------------------------------------------

## 4. Pass phased genomes to the simulation layer

[`haplotypes_from_inbred_dosage()`](https://FAkohoue.github.io/DesiredGainR/reference/haplotypes_from_inbred_dosage.md)
accepts inbred dosage data, while
[`founder_haplotypes()`](https://FAkohoue.github.io/DesiredGainR/reference/founder_haplotypes.md)
accepts phased homologues and a genetic map. Its input shape is
compatible with phased data returned by HapBlockR, but DesiredGainR does
not import HapBlockR and does not duplicate its haplotype analysis.

The resulting founder object can be passed to AlphaSimR through
[`founder_population()`](https://FAkohoue.github.io/DesiredGainR/reference/founder_population.md).
Read [Multi-cycle
simulation](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-simulation.md)
before interpreting a simulated difference between objectives.

------------------------------------------------------------------------

## 5. Hand the index to HapBlockR

HapBlockR’s `build_selection_index()` delegates `method = "dgsi"` and
`method = "qgsi"` to DesiredGainR. It returns a directionally aligned
`merit_score`, which HapBlockR can combine with genomic relationships,
haplotype complementarity, family representation and crossing
constraints. The intended flow is therefore:

``` text
genetic evaluation
  -> DesiredGainR objective, index, response and uncertainty
  -> HapBlockR parent selection, OCS, cross ranking and mate allocation
```

The main downstream choices are deliberately separate:

- `select_parents_ga_ts()` constructs a fixed-size shortlist while
  balancing genome-wide merit with complementary haplotype coverage.
- `select_parents_ocs()` performs OCS and mate allocation when a
  relationship matrix and a true OCS engine are supplied.
- `usefulness_criterion()` ranks candidate parent pairs by expected
  cross performance rather than parental merit alone.

OCS and optimum cross selection are related but not synonymous. OCS
chooses parental contributions while controlling population coancestry.
Optimum cross selection chooses the pairs and mating allocation under
cross-level merit, diversity and operational constraints. Both decisions
belong to HapBlockR. Neither changes the selection-index mathematics
owned by DesiredGainR.

``` r
# Run this stage in HapBlockR. For method = "dgsi" or "qgsi",
# build_selection_index() calls DesiredGainR internally.
index <- HapBlockR::build_selection_index(
  trait_values = candidate_gebv,
  genetic_cov = G,
  phenotypic_cov = P,
  desired_gains = desired_gains,
  directions = directions,
  method = "dgsi",
  n_select = analytical_n_select
)

# The returned merit score then feeds the appropriate HapBlockR decision.
ocs_plan <- HapBlockR::select_parents_ocs(
  merit = index$scores$merit_score,
  G = relationship_matrix,
  n_crosses = planned_crosses
)
```

Export identifiers, ranks, scores, per-trait GEBVs and the exact
objective used. Do not send only a sorted list: HapBlockR needs enough
information to apply diversity, relationship and crossing constraints
that are outside the index.

Also retain package versions, seeds, covariance provenance and the
selected-set stability analysis. A reproducible hand-off allows the next
breeding stage to distinguish a biological decision from a software
default.

### What `n_select` means in DesiredGainR

`n_select` is the number used for truncation-selection calculations. It
sets the selection proportion and hence the selection intensity. It also
defines a top-ranked set for response and stability diagnostics. It does
**not** prove that these candidates form the best parent set. The final
set can change when HapBlockR accounts for relatedness, contributions,
haplotype complementarity, cross usefulness and operational constraints.
