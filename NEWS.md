# DesiredGainR 0.5.0 (in development)

## Example data

Nine datasets replace the previous examples, which were uncorrelated standard
normal variates. That earlier material made the genetic covariance matrix
effectively an identity matrix, and an identity matrix cannot illustrate a
multi-trait selection index, because the purpose of an index is to resolve the
tension between correlated traits measured on different scales.

The new data describe one simulated tropical maize programme with six traits:
grain yield, plant height, anthesis date, the anthesis-silking interval, ears
per plant, and grey leaf spot severity. They carry three properties the earlier
examples lacked.

1. **Genuine antagonism.** Grain yield correlates positively with anthesis
   date, yet the objective requires yield to rise while the cycle shortens, so
   the two cannot both be pushed freely. Yield is also strongly and negatively
   correlated with the anthesis-silking interval. Consequently
   `gain_feasibility()` returns informative answers rather than declaring every
   target attainable.
2. **Realistic scale differences.** Genetic standard deviations span a factor
   of roughly one hundred and twenty, from ears per plant to plant height,
   which is what gives the standardisation and conditioning diagnostics
   something to detect.
3. **Internal consistency.** The markers explain the trait values. Genetic
   values are generated from ninety quantitative trait loci drawn from the
   marker panel, then whitened and recoloured so that their realised covariance
   equals `dgr_G` to numerical precision. Heritabilities are recoverable
   exactly as `diag(dgr_G) / diag(dgr_P)`.

The datasets are `dgr_traits`, `dgr_G`, `dgr_P`, `dgr_candidates`, `dgr_gebv`,
`dgr_history`, `dgr_hap1`, `dgr_hap2` and `dgr_map`. The haplotype matrices and
map satisfy the `founder_haplotypes()` contract directly, and `dgr_history`
records a previous cycle's decision produced by an undisclosed weight vector,
so that `retrospective_weights()` has something real to recover.

`tests/testthat/test-example-data.R` asserts each of these properties, since
data that quietly stopped demonstrating them would make every example
misleading.

The demonstration script `inst/examples/demo_pipeline.R` has been rewritten
around the new data and now runs the whole package in the order the decisions
actually arise: inspect the covariance matrices, translate between weights and
desired gains, test feasibility, recover the historical objective, compare
index families, check the effective weights, measure sensitivity to the stated
weights, fit the optimised desired-gain index, and finally score genomic
estimated breeding values with the quadratic index.

The former `inst/extdata/example_pheno.csv` and `inst/extdata/example_gebv.csv`
are superseded and may be deleted.

## Corrections found by running the demonstration

- **`selection_index()` returned the wrong candidates in `$selected`.** The
  ranking table is sorted by score after the selection flag is computed, and
  the returned subset was taken by filtering that sorted table with the
  unsorted flag. The rows were therefore chosen positionally rather than by
  merit. The `selected` column of `$ranking` was always correct, as were the
  observed selection differentials, so the fault was confined to the
  `$selected` table itself — which is nonetheless the element most users read
  first.

  Nothing in a suite of 356 tests detected it, because the object remained
  internally plausible: the right number of rows, the right columns, and
  consistent values. It surfaced only from noticing that two numbers printed by
  the demonstration script could not both be true, a Spearman correlation of
  0.938 between two indices alongside an overlap of one candidate in twenty
  between their selected sets. Tests now assert the relationship between those
  two quantities, and assert that `$selected` matches an independently
  recomputed top-n.

Four defects survived a passing test suite and were exposed only by reading
the output of `demo_pipeline.R` on the new data. Each is a case where the code
did what it was told and the result was nonetheless misleading.

- **The index coefficient of variation diverged on a centred index.**
  Standardising the traits places the index mean at zero, so `sd / mean`
  reported values of the order of 10^17. The quantity is undefined there, and
  it is now withheld with `CV_I_note` explaining why.

- **The objective-setting functions had no `lower_is_better`.** Every other
  entry point in the package accepts it, so `implied_economic_weights()`,
  `implied_desired_gains()` and `gain_feasibility()` silently interpreted their
  inputs in the raw trait direction. A positive desired gain for plant height
  was therefore a request to grow taller. The author of the demonstration
  script made exactly this mistake. All three now accept `lower_is_better` and
  orient the covariance matrices internally.

- **`implied_desired_gains()` could not return the units it was given.** A
  round trip through `implied_economic_weights(gain_units = "genetic_sd")` came
  back in trait units, which looks like a failed inversion. It now takes
  `gain_units`, and a test asserts the round trip is exact in all three.

- **`run_qgsi()` gave no warning when trait scales were mismatched.** Economic
  weights multiply genomic estimated breeding values directly, so on the
  original scale a trait with a large standard deviation dominates whatever
  weight it was given. In the demonstration, plant height overwhelmed the index
  and the expected gain in grain yield came out *negative* while selecting for
  higher yield. `run_qgsi()` now issues the same conditional scale warning that
  `run_dgsi()` already carried, naming the dominating trait and its share.

- **`selection_index()` gains `center_traits`, separate from `scale_traits`.**
  Centring was previously unconditional, which placed the index mean at zero
  in every call and therefore made the coefficient of variation permanently
  undefined. The two transformations are now controlled independently, as they
  already were in `run_qgsi()`, so a published result reporting `CV_I` on an
  index built from raw trait values can be reproduced. Centring never changes
  a ranking, since it shifts every score by the same constant.

- **Units are now reported in the units the caller used.** `gain_feasibility()`
  returned the attainable response in original trait units regardless of the
  `gain_units` the target was stated in, so a request of one genetic standard
  deviation of yield came back as 0.41 tonnes per hectare. Comparing the two
  invites misreading. `attainable_response_input_units` and
  `shortfall_input_units` are now returned alongside the trait-unit versions.

- **Coefficients are reported on a comparable footing.**
  `retrospective_weights()` returns `coefficients_per_sd`, and the print method
  shows that rather than the raw vector. A coefficient carries inverse trait
  units, so ears per plant, whose standard deviation is 0.10, received the
  largest raw coefficient in the example programme while representing a modest
  emphasis.

- `run_qgsi()` no longer repeats the citation for the expected-gain equation on
  every row of `expected_gain_per_trait`, where six identical sentences swamped
  the printed table. It is carried as a `basis` attribute on that table and as
  `expected_gain_basis` in the result.

- Minor: `gain_feasibility()` no longer reports a required selection of one
  candidate when the required proportion is smaller than a single candidate,
  since rounding up overstated what the population could deliver. A required
  intensity beyond the reach of normal truncation is reported as such rather
  than converted into a proportion.

## Corrections found by the test suite

- **The two desired-gain formulations are one index, not two.** For a square
  invertible `G`, the Yamada form reduces algebraically to the Pesek-Baker
  form, since `P⁻¹G(GP⁻¹G)⁻¹d = P⁻¹G(G⁻¹PG⁻¹)d = G⁻¹d`. Both `method` values
  are retained as documented routes to one estimand. They remain worth
  distinguishing, however, because only the Yamada route is available when `G`
  is rank deficient, and the two are not numerically equivalent under an
  ill-conditioned matrix. The documentation previously claimed they gave
  different coefficients and rankings; that claim was wrong.

- **Expected per-trait gains rest on two approximations, not one.** The
  identity `Cov(gamma, I) = Gamma w` is exact, because the third central
  moments of a zero-mean multivariate normal vanish. Two further steps are not.
  First, the achieved selection differential on the index is assumed to equal
  the normal-theory value `i * sd(I)`; the index is a quadratic form and is
  right skewed, so the differential actually achieved is larger. Second,
  `E[gamma | I]` is assumed linear in `I`, which also fails for a quadratic
  index.

  Monte Carlo testing shows these two errors partly cancel. Correcting only the
  first, by substituting the observed differential, makes the prediction worse
  rather than better. The total index standard deviation is therefore used
  because it is the correct linear-regression denominator, and not because it
  demonstrably fits simulation best: at modest curvature the residual error of
  five to nine per cent exceeds the two per cent difference between the two
  candidate denominators, so simulation cannot arbitrate between them. Both
  assumptions are now reported alongside the gains.

- `simulate_selection_cycles()` tracked genetic gain with `bv()`, whose values
  are expressed relative to the supplied population and therefore have a mean
  of zero in every cycle. Cumulative gain was consequently always zero. Gain is
  now tracked with `gv()`, which is anchored by the trait mean set when the
  trait was added. Covariance estimation for rebuilding the index continues to
  use breeding values for additive programmes, since centring does not affect a
  covariance.

- `effective_weights()` failed on a rank-deficient covariance matrix, where the
  shares are undefined, and that failure propagated out of `run_dgsi()`.
  Shares are now returned as `NA` in that case, and `run_dgsi()` treats the
  diagnostic as non-fatal.

- `effective_weights()` no longer warns by default, and its threshold adapts to
  the number of traits. A concentrated effective weight can arise from a
  deliberately asymmetric objective as readily as from mismatched trait scales,
  so warning on every asymmetry was noise. The scale mismatch that Crosbie et
  al. (1980) describe is still reported, by the separate `run_dgsi()` warning
  that fires only when trait standard deviations differ by more than fivefold.

This release begins the multi-cycle simulation layer. Its purpose is narrow and
deliberately kept so: it exists to compare desired-gain directions over several
breeding cycles, because a direction that maximises response in the first cycle
may exhaust the genetic variation that response depends on within a few more.
It is not a crossing-plan tool, and parent, cross and mating allocation remain
outside the scope of this package.

## Founders come from real phased marker data

DesiredGainR does not simulate founder genomes. Therefore the linkage
disequilibrium, allele frequency spectrum and population structure of every
simulation are those of the breeder's own germplasm rather than those of an
assumed demography, which removes the principal scientific weakness of
simulation-based recommendations.

- `founder_haplotypes()` validates and packages phased founder haplotypes. It
  accepts variant-by-individual matrices `hap1` and `hap2` of 0 and 1, together
  with a variant map giving chromosome and genetic position. This is the shape
  returned by phased variant-call-format readers, including
  `HapBlockR::read_phased_vcf()`.

  Genotype dosage in 0, 1, 2 coding is rejected, and the error explains why. A
  haplotype is one chromosome copy, so at a biallelic locus it carries a single
  allele rather than an allele count; dosage is the sum of the homologues and
  therefore records how many alternative alleles an individual carries but not
  which copy carries them. Assigning heterozygous calls at random would destroy
  the linkage disequilibrium that determines how favourable alleles segregate
  together. Hence phase must be supplied rather than inferred.

- Missing calls are permitted on input and resolved under an explicit
  `missing_policy` of `"error"`, `"drop_variant"` or `"drop_individual"`.
  Phased variant call format output legitimately carries missing values,
  whereas AlphaSimR requires complete haplotypes, so they must be resolved
  rather than imputed silently. The rate observed before resolution, and
  everything removed, are recorded in the returned object.

- The contract is deliberately **diploid and biallelic**, and a polyploid
  dataset raises an explanatory error. Extending the two homologue matrices to
  further copies would not make the simulation polyploid, because
  autopolyploid and allopolyploid meiosis require an explicit crop-specific
  pairing and recombination model, including multivalent formation and double
  reduction. Applying a diploid meiotic model to polyploid homologues would
  produce confidently wrong results, so polyploid support must wait for a
  dedicated homologue-level representation.

- `dosage_diagnostics()` reports overall, per-individual and per-marker
  heterozygosity and missingness. Run it before any conversion, because a low
  overall rate can conceal a small number of badly affected individuals or
  markers, and that distribution determines which conversion policy is
  cheaper.

- `haplotypes_from_inbred_dosage()` derives phase from dosage for highly inbred
  diploid material. A dosage of 0 or 2 is unambiguous, so doubled haploids,
  recombinant inbred lines and advanced selfed generations can be converted
  without external phasing.

  No threshold is imposed on residual heterozygosity, because no universally
  appropriate level exists: it depends on the generation, mating history, crop,
  population type, genotyping error rate and quality-control procedure, and
  must therefore be measured from each dataset. Heterozygous calls are never
  resolved silently and are never assigned to a homologue. The available
  policies are (i) `"error"`, (ii) `"drop_variant"`, (iii)
  `"drop_individual"`, and (iv) `"mask"`, which sets the call to missing for
  resolution by `founder_haplotypes()`. Every rate observed and every call
  removed or masked is recorded in the returned `conversion` attribute.

- `founder_population()` builds the AlphaSimR founder population and calibrates
  the trait architecture to the genetic covariance matrix the breeder has
  already estimated: variances from the diagonal of `G` and correlations from
  the corresponding correlation matrix.

## Recurrent selection across cycles

- `simulate_selection_cycles()` runs a recurrent-selection programme forward
  under a fixed desired-gain direction and records, for every cycle, the
  genetic mean and variance per trait, the mean relationship among selected
  parents, and the implied effective population size.

  Three mating systems are supported: (i) `"self"` for self-pollinated line
  development, advancing by selfing or doubled haploidy, (ii) `"outcross"` for
  recurrent selection in a random-mating population, and (iii) `"clonal"` for
  clonally propagated crops, where selection acts on total genetic value
  because the clone inherits dominance intact. The clonal system therefore
  requires a setup built with `dominance_degree`, and `G` should then be the
  genotypic rather than the additive covariance.

  `reestimate_index` controls whether the index is rebuilt from each cycle's
  own simulated data, which is what a breeding programme actually does and
  which propagates estimation error across cycles, or held fixed at the cycle
  one solution, which isolates the effect of the desired-gain direction itself.
  A large divergence between the two settings means the recommendation is
  sensitive to estimation error, and the breeder should be told so.

## Searching for the best desired-gain direction

- `optimize_desired_gains()` searches for the desired-gain direction that
  performs best once several cycles have been simulated. Only the direction is
  optimised, because the attainable magnitude is fixed by selection intensity,
  so the domain is the unit sphere with one fewer free parameter than there are
  traits. That is small enough to cover densely for the trait numbers breeders
  use.

  **The simulation is the objective function and is never replaced by a cheaper
  approximation.** A Gaussian process is fitted to the accumulated results, but
  it decides only where the next simulation should be spent; it never filters,
  screens, or substitutes for an evaluation, and because the acquisition
  function retains an exploration term across the whole domain, no region can
  be permanently excluded by the surrogate. An earlier design used a
  deterministic infinitesimal-model recursion as the surrogate's prior mean;
  that was rejected, because it would embed a second genetic model whose
  assumptions could fail in precisely the situations where simulation is most
  needed.

  Four ranking modes are available: (i) `"pareto"`, the default, returning the
  non-dominated set of multi-cycle outcomes without requiring economic weights,
  (ii) `"economic"`, maximising a weighted sum, (iii) `"target"`, minimising
  the distance to stated absolute targets, and (iv) `"constrained"`, maximising
  a focal trait subject to floors on the remainder.

  Two properties matter for interpretation. Every direction is evaluated with
  the same sequence of seeds, so comparisons between directions share their
  stochasticity, which removes more comparison variance than raising the
  replicate count. The reported frontier is computed from posterior means
  rather than raw replicate averages, because a frontier built from raw draws
  is populated by fortunate runs.

  No single best direction is returned for the scalar modes either. The optimum
  is conditional on the supplied genetic covariance, the founder germplasm and
  the programme parameters, so the result reports every direction whose
  posterior outcome lies within `stability_tolerance` of the best.

- Diversity, measured as the negated mean relationship among selected parents,
  is included as an additional Pareto objective by default. Family balancing
  and coancestry control can cost a large share of nominal gain, so the
  trade-off is reported rather than assumed away.

- Searches are resumable. Supplying `checkpoint` writes the accumulated
  evaluations after every simulation and reloads them automatically on a
  later call, because a realistic budget takes hours and an interruption
  should not discard the work.

## Dependencies

- `AlphaSimR` is added to `Suggests`. Every simulation function fails with an
  explicit installation message when it is absent, and the tests skip rather
  than fail. DesiredGainR does not depend on, import, or call HapBlockR, which
  keeps the dependency arrow one-directional.

# DesiredGainR 0.4.0

This release adds the objective-setting layer. Selection indices are
straightforward to compute; deciding what to select for is not, and the
literature identifies the specification of economic weights as the principal
obstacle to routine index use. Hence the new functions do not compute better
indices. They help a breeder state, test, and defend the objective an index is
built from.

Every addition is additive. No existing function, argument, result element, or
column name has been renamed or removed.

## Translating between the two ways of stating an objective

The Smith-Hazel economic index and the Pesek-Baker desired-gain index are two
parameterisations of one linear index. Therefore each desired-gain vector
corresponds to exactly one economic-weight vector.

- `implied_economic_weights()` returns the weights implied by desired gains,
  through `w = G^-1 P G^-1 d`.
- `implied_desired_gains()` performs the reverse translation, `d = G P^-1 G w`.
  This is the more useful direction in practice: a breeder who proposes weights
  can immediately see what response those weights actually request.

## Testing whether an objective is attainable

- `gain_feasibility()` reports the selection intensity a target response
  requires, `i = sqrt(d' G^-1 P G^-1 d)`, the selected proportion that delivers
  it, whether that proportion is reachable in the available population, and the
  per-trait shortfall when it is not. Covarrubias-Pazaran (2021) directs
  breeders to external software for this check; it is now available in R.

## Recovering an objective from past decisions

- `retrospective_weights()` recovers the linear rule implicit in a programme's
  historical selections, through `b = P^-1 s`. This is the method used to
  introduce a selection index at the International Maize and Wheat Improvement
  Center (CIMMYT) and the International Rice Research Institute (IRRI). The
  recovered weights are a starting point and reproduce past bias, so they
  should be inspected and adjusted before adoption.

## Quantifying how much the objective matters

- `weight_sensitivity()` perturbs the stated weights and reports how often the
  selected set survives unchanged. A high stability proportion means further
  effort refining the weights will not change the decision; a low one means the
  objective, rather than the index, is the binding uncertainty.
- `effective_weights()` reports each trait's contribution as `b_j * sigma_gj`
  rather than `b_j`, and warns when one trait dominates. Crosbie et al. (1980)
  showed that equal weights on unstandardised traits concentrate selection on
  whichever trait carries the largest genetic variance.
- `matrix_diagnostics()` reports the condition number of a covariance matrix.
  An index built by inverting an ill-conditioned matrix is numerically
  meaningless while remaining superficially plausible.

## Classical index families

- `selection_index()` builds the Smith-Hazel, base, Pesek-Baker, Yamada, and
  Mulamba-Mock indices, together with independent culling and tandem selection,
  from one interface, so that families and objectives can be compared on
  identical data. Guimaraes et al. (2021) found the Mulamba-Mock rank-sum
  index, which requires neither economic weights nor covariance matrices, to
  give the most balanced multi-trait response of the methods they compared.

  **Two indices share the name Pesek-Baker.** Pesek and Baker (1969) give
  `b = G^-1 d`, which is `method = "pesek_baker"` and what Rahimi and Debnath
  (2023) implement. Yamada et al. (1975) give `b = P^-1 G (G P^-1 G)^-1 d`,
  which is `method = "yamada"` and what Joukhadar et al. (2024) and `run_dgsi()`
  use. Both satisfy `G'b = d`, but they yield different coefficients and
  different rankings. Hence any comparison with published results must state
  which formulation was used.

## Evaluating an index

- `evaluate_index()` reports the four criteria used by Rahimi and Debnath
  (2023) and by the established selection-index software: (i) `R_HI`, the
  correlation between the index and net merit, (ii) `delta_H`, the expected
  gain in aggregate merit, (iii) `RE`, efficiency relative to direct selection
  on a main trait, and (iv) `CV_I`, the coefficient of variation of the index,
  together with the expected response for each trait.

  Relative efficiency below 1 is the intended trade-off of multi-trait
  selection, not a defect; every value reported by Rahimi and Debnath was below
  1.

## Other changes

- `estimate_genetic_covariance()` gains `method = "adjusted_means_surrogate"`,
  which returns the covariance of across-environment adjusted means labelled as
  a working surrogate for genetic covariance. Covarrubias-Pazaran (2021)
  recommends this for operational use when no fitted genetic covariance matrix
  is available. Earlier releases refused this workflow; it is now supported
  explicitly, with its assumptions recorded, rather than left to be performed
  silently outside the package.
- `run_dgsi()` now returns `non_iterated`, the classical Yamada solution using
  the supplied desired gains without iteration. Joukhadar et al. (2024)
  benchmark against exactly this comparator, and it costs nothing because it is
  the search starting point.
- `run_dgsi()` now returns `effective_weights` and warns when trait standard
  deviations differ by more than fivefold while `scale_traits = FALSE`. The
  warning is conditional on the scales actually differing, so analyses on
  comparable scales remain quiet.

# DesiredGainR 0.3.1

No function, argument, result-element, or column name was removed or renamed.
Existing code continues to run unchanged. Two corrections do change returned
*values*; both superseded quantities remain available under new names.

## Corrected results (numerical change)

- `run_qgsi()`: `expected_gain_per_trait$Expected_Genetic_Gain` now divides
  `Gamma %*% w` by the **total** index standard deviation,
  `sqrt(w'Gamma w + 2 tr(W Gamma W Gamma))`. Releases up to 0.3.0 divided by
  `sqrt(w'Gamma w)`, the purely linear (LGSI) index standard deviation, which
  overstated every per-trait gain whenever `W` was non-zero. The previous
  values are still returned in the new column
  `Expected_Genetic_Gain_LinearSD`. A Monte Carlo test against simulated
  truncation selection now pins the corrected formula
  (`tests/testthat/test-theory-simulation.R`).
- `run_qgsi()`: `theoretical_parameters$squared_index_merit_correlation` now
  uses the general definition `Cov(H, I)^2 / (Var(I) Var(H))`. Releases up to
  0.3.0 returned `Var(I) / Var(H)`, which equals the squared correlation only
  when the index is the MSPE-optimal predictor of net merit — not guaranteed,
  because `run_qgsi()` accepts arbitrary user weights. The former ratio is
  retained as `theoretical_parameters$variance_ratio_index_to_merit`.
  Unchanged when `true_G` is not supplied (both remain `NA`), and unchanged
  when `true_G == Gamma` (both equal 1).
- `run_dgsi()`: ties in the index score are now broken by ascending candidate
  identifier instead of by input row order, so the selected set no longer
  depends on how the rows were sorted. Only tied candidates are affected.

## Fixed

- `run_dgsi()` no longer leaves the caller's random number generator reseeded.
  The RNG state is saved on entry and restored on exit, so `seed` gives
  reproducible results without altering the user's stream.
- `run_qgsi()` now warns when `center_traits = FALSE` and `Gamma` is estimated
  internally. In that combination the scores use uncentred GEBVs while the
  estimated `Gamma` is a centred covariance, so the reported index mean,
  variances, and expected gains describe a different index from the one in
  `ranked_geno`. Supply `Gamma` explicitly, or use `center_traits = TRUE`.

## Documentation

- `run_dgsi()`: `dg` is now documented as being in **candidate
  standard-deviation units**, which is what the objective has always compared
  against, regardless of `scale_traits`. A new Details section gives the
  realised-response and objective formulas, explains the `0.25` and `0.05`
  floors that were previously undocumented constants, and states that
  best-of-`n_rep` selection is optimistically biased because it is made on the
  same candidates that are then selected.
- Added a worked `@examples` block to `run_dgsi()`.

## Internal

- `globalVariables()` is now declared once in `R/globals.R` instead of in two
  files with overlapping and partly obsolete entries.
- `.desiredgainr_z()` and `.validate_square_matrix()` are unused by any
  exported function. They are retained for now in case external code calls
  them via `:::`, and are marked for removal in a future release.
- `inst/CITATION` uses `bibentry()` instead of the deprecated `citEntry()`.
- `DESCRIPTION`: added `Language: en-GB` and `Depends: R (>= 4.1)`; dropped the
  hand-written `Author:`/`Maintainer:` fields, which duplicated `Authors@R` and
  could drift out of sync.
- `run_dgsi()` hoists the loop-invariant column means and standard deviations
  out of the search loop. Results are identical; the search is substantially
  faster on large candidate sets.
- Added `tests/testthat/test-theory-simulation.R`: Monte Carlo validation of
  the QGSI index mean, variance, per-trait gain and squared correlation, plus
  a Pesek-Baker consistency check, an RNG-restoration check, and a row-order
  invariance check.
- Added `tests/testthat/test-downstream-contract.R`, which pins the call
  surface and the result elements used by `hapblockr::build_selection_index()`
  (`coefficients`, `realised_response`, `ranked_geno$SelectionIndex`,
  `ranked_geno$QGSI`, `replicate_diagnostics$Chosen`, `component_summary`, the
  nine argument names that `dgsi_control` protects, and the
  `center_traits = TRUE` default). None of the 0.3.1 corrections touch these.

# DesiredGainR 0.3.0

- Removed the author-named DGSI wrapper; `run_dgsi()` is now the sole public
  desired-gain index entry point.
- Added `estimate_genetic_covariance()` with clearly distinguished prediction
  covariance, PEV-corrected, SE-diagonal-corrected, and
  relationship-adjusted estimators, together with assumptions and provenance.
- Implemented QGSI genomic covariance estimation from reference GEBVs using
  Supplementary Equations 19.1 and 19.2 of Cerón-Rojas et al. (2026).
- Added model-based QGSI mean, linear and quadratic variance, normal-selection
  intensity, expected net-merit response, and expected per-trait gains.
- Added optional simulation-only squared accuracy and mean squared prediction
  error calculations when the true genetic covariance matrix is supplied.
- Added exact top-ranked selection and clearly separated observed GEBV
  differentials from expected or realised genetic gain.
- Separated DGSI desired gains from QGSI linear and quadratic economic weights
  throughout the combined pipeline.
- Removed the unsupported constructor that fabricated quadratic economic
  weights from desired gains and empirical correlations.
- Added equation-level tests for scores, contributions, covariance estimators,
  response parameters, relationship-matrix alignment, and selection.

# DesiredGainR 0.2.0

- Renamed the package from DGQGSI to DesiredGainR.
- Added canonical `run_dgsi()` and `run_qgsi()` interfaces.
- Required an explicit genetic covariance matrix for DGSI and recorded the
  provenance of `P` and `G`.
- Changed threshold selection to eligibility followed by index ranking.
- Added automatic best-replicate selection, convergence, coefficient, ranking
  and selected-set stability diagnostics.
- Made missing-value handling explicit.
- Required an explicit symmetric quadratic-weight matrix for QGSI.
- Added candidate-specific linear, squared and cross-product contributions.

# DGQGSI 0.1.0

- Initial package structure.
