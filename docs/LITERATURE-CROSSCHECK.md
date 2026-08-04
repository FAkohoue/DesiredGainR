# DesiredGainR — cross-check against the Publications folder

## Architecture decision added 4 August 2026

Montesinos-López, Montesinos-López, Hernández-Suárez and Alemu (2026),
*A selection index with minimal genetic relatedness for multi-trait data
via binary quadratic programming* (Plant Methods 22:7,
<https://doi.org/10.1186/s13007-025-01484-4>), does **not** define
another DesiredGainR index family. Its QPMSI selects a fixed-size binary
subset by combining multi-trait merit with a genomic-relatedness
penalty. That is an operational parent-selection problem and belongs to
HapBlockR.

QPMSI is not true optimal contribution selection (OCS), because its
binary selection variables do not optimise continuous parental
contributions. It is also not a complete optimum-cross-selection method,
because it does not choose parent pairs or allocate matings. Those
distinctions remain HapBlockR’s responsibility. The intended integration
is that HapBlockR obtains the multi-trait merit score from DesiredGainR,
then applies parent selection, OCS or optimum cross selection as
required. DesiredGainR must not import HapBlockR or reproduce those
downstream optimisers.

Date 2026-07-30. Sources read in full or in relevant part:

1.  Covarrubias-Pazaran G. (2021) *Guideline: Practical implementation
    of selection indices.* CGIAR Excellence in Breeding, 18/06/2021.
2.  Covarrubias-Pazaran G. *Bringing a selection index into the
    CIMMYT-Maize programs.* EiB deck.
3.  Covarrubias-Pazaran G. *Bringing a selection index into the IRRI
    programs.* EiB deck.
4.  Joukhadar R. et al. (2024) *Optimising desired gain indices to
    maximise selection response.* Front. Plant Sci. 15:1337388.
5.  Guimarães P.H.R. et al. (2021) *Index selection can improve the
    selection efficiency in a rice recurrent selection population.*
    Euphytica 217:95.
6.  Rahimi M. & Debnath S. (2023) *Estimating optimum and base selection
    indices … new and simple SAS and R codes.* Sci. Rep. 13:18977.
7.  Vieira R.A., Nogueira A.P.O., Fritsche-Neto R. (2025) *Optimizing
    the selection of quantitative traits in plant breeding using
    simulation.* Front. Plant Sci. 16:1495662.

Headline: the package is **faithful to Joukhadar** on the mechanics it
implements, and **misaligned with the CGIAR guideline** on the workflow
its target users are told to follow. The largest single gap is that it
reports none of the index-evaluation criteria that every competing tool
reports.

------------------------------------------------------------------------

## MUST FIX

### F1. `scale_traits = FALSE` contradicts the guideline, and the current mix matches nothing

The guideline’s step 5 is *standardise the adjusted means*, and step 6
then defines `d = 1` as “select individuals that are 1 standard
deviation away from the population mean”. Joukhadar likewise used **the
correlation matrix of GEBVs as P** — i.e. a standardised space
throughout.

[`run_dgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_dgsi.md)
defaults to `scale_traits = FALSE`. It then computes `b` from *unscaled*
`P` and `G` but measures `realised_response` in *candidate SD units*.
That hybrid corresponds to neither reference. Coefficients live in raw
trait space while the objective lives in standardised space.

Given the `hapblockr` constraint, do not flip the default silently.
Instead: warn when `scale_traits = FALSE`, document that the
guideline-conformant workflow is `scale_traits = TRUE`, and flip the
default at the next major version with a deprecation cycle.

### F2. The objective function is not Joukhadar’s, but the package cites Joukhadar

Joukhadar’s goodness of fit is logarithmic, with `base = √(e^{dg})`, and
the accepted penalty is `q = Σᵢ (gofMAXᵢ − gofᵢ)/gofMAXᵢ`, minimised
across traits.

[`run_dgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_dgsi.md)
minimises `Σ vⱼ ((rⱼ − dⱼ)/max(|dⱼ|, 0.25))²` — a weighted relative
squared error with an undocumented floor. It is a defensible objective;
it is not the published one. Either add
`objective = c("relative_squared", "joukhadar")` with the published form
available, or stop implying fidelity in DESCRIPTION and README.

The rest of the algorithm *is* faithful: proposal recentring on the last
accepted sample, 1,000 iterations, 20 replicates, cross-replicate
correlation as a convergence check. Worth stating that explicitly in the
docs.

### F3. The Crosbie trap is exactly what the current defaults invite

Guimarães et al. quote Crosbie et al. (1980): when equal weights are
assigned across traits, *most weight is placed on the trait with the
largest genetic variance and least on the trait with the smallest*. The
guideline recommends starting from `d = 1` for every trait. With
`scale_traits = FALSE` that combination silently hands the index to
whichever trait has the biggest raw variance.

Add a diagnostic reporting each trait’s *effective* contribution to the
index (`bⱼ σ_gⱼ` as a share of the total), and warn when one trait
exceeds a threshold share. This is cheap and prevents the most common
failure mode in the cited literature.

### F4. The package blocks the guideline’s own recommended `G` surrogate

The guideline states plainly:

> “In practice G is not known and is enough to use across-environment
> adjusted means by a mixed model to scale, define the desired σ of
> progress and use the covariance of those adjusted means as a surrogate
> for G.”

[`estimate_genetic_covariance()`](https://FAkohoue.github.io/DesiredGainR/reference/estimate_genetic_covariance.md)
has no such method, and the package documentation states the opposite
position — that BLUEs and adjusted means cannot identify genetic
covariance. Both positions are defensible; the problem is that
DesiredGainR currently **refuses the workflow its target users have been
instructed to use**, with no supported on-ramp.

Add `method = "adjusted_means_surrogate"`, implementing exactly the
guideline’s recommendation, labelled with its estimand and its
limitations, and cited to Covarrubias-Pazaran (2021). This keeps the
package’s honesty while removing the wall. The existing provenance
machinery already carries the caveat.

### F5. No index-evaluation criteria — the biggest credibility gap

Rahimi & Debnath evaluate every index on four criteria, and name the
software breeders already use: **Genes, MIX, RIndSel, SelAction**.
DesiredGainR reports none of them:

| Criterion | Definition | Note |
|----|----|----|
| `RHI` | correlation between index and aggregate genotype, `r_HI` | the standard accuracy measure |
| `ΔH` | genetic gain of the aggregate |  |
| `RE` | efficiency relative to direct selection on the main trait | “index is better when RE \> 1” |
| `CV_I` | phenotypic coefficient of variation of the index |  |

For DGSI the package also reports no expected genetic response at all
(audit §3.2). Until these exist, a reviewer comparing DesiredGainR with
RIndSel will find it strictly less informative.

------------------------------------------------------------------------

## MUST HAVE

### H1. Retrospective weights, `b = P⁻¹s`

This appears in the guideline footnote *and* is the headline method in
both deployment decks:

> “If the selection differentials represent the breeder’s goal, then the
> weights determine the merit of individuals selected. Weights can be
> back calculated.”

At CIMMYT-Maize and IRRI this is how an index was actually bootstrapped:
take the selection differentials the programme has historically
achieved, recover the weights implied by them, then refine. It is
revealed-preference weight elicitation from data the breeder already
has, and it answers objective 1 better than anything I proposed in the
design document. It was absent from my plan.

`retrospective_weights(historical_selected, historical_population, P)` →
`b = P⁻¹s`, with the differentials also reported in σ units so they can
be used directly as desired gains for standardised traits.

### H2. Mulamba & Mock rank-sum index

In Guimarães et al. this index **beat Smith-Hazel and Tai on every
trait**, gave the best population shift, and requires **no economic
weights at all**. Crosbie et al.’s four cited advantages: unaffected by
unequal trait variances, requires no genetic parameters, simple, and
gives better selection differentials.

It is about fifteen lines of code and is the pragmatic baseline this
literature keeps returning to. Its absence is conspicuous.

### H3. Smith-Hazel and the base (Brim/Williams) index as first-class comparators

Rahimi & Debnath find the **base index** (`b = w` directly) sometimes
*outperforms* the optimum index on ΔH and RE. Joukhadar benchmarks
against the non-iterated Pesek-Baker. DesiredGainR can compute none of
these, so it cannot reproduce any comparison in any of these papers,
including the one it cites.

Note this is nearly free:
[`run_dgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_dgsi.md)
already computes the non-iterated solution as
`initial_b <- coefficient(dg)` before the search starts. Return it.

### H4. Feasibility check — replaces an external tool the guideline sends users to

The guideline instructs breeders to leave R and use **DESIRE**
(une.edu.au) to “understand if the progress of 1 or more σ is possible
given your current trait genetic correlations and selection intensity”.
Joukhadar’s INDEX3 targeted 4 σ and the paper notes this “was impossible
to achieve given that only 100 lines were selected from 3,005 (~3.3%) …
theoretically only 0.05% of the whole population could be above or below
4 standard deviations for a single trait.”

The exact test from the design document — `i_required = √(d'G⁻¹PG⁻¹d)`,
converted to a required selected proportion — closes this in-package.
Both sources independently establish demand for it. This is the highest
value-per-line feature in the whole plan.

### H5. Tai (1977) index

A generalisation of Pesek-Baker admitting secondary traits,
`b = G⁻¹Δg_d`, used as a main comparator in Guimarães et al. Directly
adjacent to what the package already does.

### H6. Restricted and predetermined-proportional-gain indices

Kempthorne & Nordskog (1959), Mallard (1972), Harville (1975), Tallis
(1985) — all cited by Joukhadar as the constrained family that
desired-gain indices are an alternative to. Already audit §6; the
literature confirms it is expected.

------------------------------------------------------------------------

## NICE TO HAVE

### N1. Gain as a percentage of the population mean

Guimarães et al. use the Vencovsky & Barriga convention throughout:
`SG% = ((X̄_selected − X̄_population) · h²) / X̄_population × 100`. This is
the reporting format in the Brazilian and CGIAR-adjacent literature.
DesiredGainR reports SD units only.

### N2. Family-structured index

The CIMMYT-Maize deck: “the index is calculated **among family and then
within family**.” The package’s own `example_pheno.csv` already carries
a `Family` column that nothing uses.

### N3. Categorical and marker-based gates

The IRRI deck combines quantitative BV thresholds with qualitative
R-gene status (`Xa21`, `Xa5`, `Pi54`, `Pita` = R/S) as hard filters
before indexing. `trait_min` is numeric-only, so this workflow cannot be
expressed.

### N4. Comparison against the programme’s existing culling algorithm

Both decks follow the same arc: encode the breeder’s current
sequential-culling rules as an algorithm, run the index, compare
selected sets and total merit (the maize deck reports merit 3.1327 vs
3.1345 and a genotype-by-genotype overlap table). A
`compare_to_culling()` helper would be the single most effective
adoption tool — it is how both CGIAR centres were persuaded.

Related: both decks show the index vs **independent culling** vs
**tandem selection** comparison. Those two baselines are trivial and
worth shipping.

### N5. Stage and pipeline metadata

The guideline’s step 1 is a table of traits × stages × surrogate ×
number of environments; the decks frame everything as *Program → Stage →
Product profile*. The package has no notion of a stage. Carrying this as
metadata would let results be filed against the scheme they came from.

------------------------------------------------------------------------

## Changes this forces on the 0.4–0.6 design document

### D1. The multi-cycle optimiser must compare directions under OCS, not only truncation selection

Vieira et al. report that **optimal contribution selection yields higher
index values than truncation selection while maintaining lower
inbreeding**, and that a single cycle of optimal haploid selection can
outperform both in the long term. If
[`optimize_desired_gains()`](https://FAkohoue.github.io/DesiredGainR/reference/optimize_desired_gains.md)
simulates truncation selection only, it optimises a strategy the
simulation literature already shows to be long-term suboptimal. OCS must
be a selectable selection rule in the cycle templates, not an
afterthought.

### D2. Number of parents recycled is a first-order lever

The same review: fewer parents favour short-term gain, more parents
favour long-term gain. This interacts directly with the desired-gain
direction and must be a simulator parameter exposed in the Pareto
output, not fixed.

### D3. Joukhadar’s stated limitation is precisely objective 2

> “more empirical testing is required using multi-breeding cycle data to
> ensure that the calculated indices are sufficiently powerful to
> maximise the correlation between the index and the net genetic merit.”

and

> “it is still not possible to determine if the method maximises the
> correlation between the selection index and net genetic merit …
> because the covariance between both is undefined given that it depends
> on the estimation of the economic weight.”

That second quotation is the same circularity flagged in the design
document — ranking multi-trait outcomes needs a merit definition, which
needs weights. It should be cited as motivation for the Pareto-frontier
default, which sidesteps it.

### D4. `P` convention

Joukhadar used the **correlation matrix of GEBVs** as `P`. The guideline
standardises before indexing. The design document should state the
standardised convention as the recommended default rather than leaving
`P` unconstrained.

------------------------------------------------------------------------

## Revised priority

| \# | Item | Effort | Source |
|----|----|----|----|
| 1 | Index-evaluation criteria (RHI, ΔH, RE, CV_I) + expected response | low | Rahimi; audit §3.2 |
| 2 | Feasibility check `i_required` | low | Guideline (DESIRE); Joukhadar |
| 3 | Retrospective weights `b = P⁻¹s` | low | Guideline; both decks |
| 4 | Mulamba & Mock, Smith-Hazel, base index, non-iterated Pesek-Baker | low | Guimarães; Rahimi |
| 5 | `adjusted_means_surrogate` for `G` | low | Guideline |
| 6 | Effective-weight diagnostic + `scale_traits` warning | low | Crosbie via Guimarães |
| 7 | Weight/desired-gain duality, sensitivity | medium | design doc |
| 8 | Joukhadar objective as an option | medium | Joukhadar |
| 9 | Tai, restricted, PPG indices | medium | Guimarães; Joukhadar |
| 10 | Simulation engine with OCS | high | Vieira |
| 11 | Surrogate-assisted optimiser | high | design doc |

Items 1–6 are all low effort, all directly demanded by the literature,
and together they would make DesiredGainR competitive with RIndSel and
consistent with the CGIAR guideline. None of them requires AlphaSimR.

------------------------------------------------------------------------

# Addendum: items arising from `00_integrated_review.md`

The integrated review corroborates every finding above. It also reports
numerical detail from the two EiB decks and from Rahimi & Debnath that
changes or adds to the list.

## F6 (CORRECTED). The two “Pesek-Baker” formulations are the same index

**This item was wrong as originally written, and the test suite caught
it.**

| Form | Coefficients | Used by |
|----|----|----|
| Pešek & Baker (1969), original | `b = G⁻¹d` | Rahimi & Debnath (2023); Tai index in Guimarães et al. |
| Yamada et al. (1975), phenotypic-criteria | `b = P⁻¹G(GP⁻¹G)⁻¹d` | Joukhadar et al. (2024); **DesiredGainR**; `hapblockr` |

For a square invertible `G` these are algebraically identical:

    P⁻¹G(GP⁻¹G)⁻¹d = P⁻¹G(G⁻¹PG⁻¹)d = G⁻¹d

The original claim that they “produce different coefficients and
different rankings” was false. `test-index-families.R` asserted the
difference, failed, and forced the correction.

Two consequences.

First, **the naming collision with `hapblockr` is far less serious than
reported**. Both packages compute the same index for an invertible `G`;
only the route differs. It remains worth aligning the names, but no user
will get contradictory numbers from it.

Second, the distinction is still real in two situations, and both are
worth documenting rather than hiding: (i) when `G` is singular or rank
deficient the direct inverse does not exist and only the Yamada route is
available, and (ii) the two routes are not numerically equivalent,
because one inverts `G` while the other inverts `P` and then `GP⁻¹G`.
Under an ill-conditioned matrix they diverge — which is the real
explanation for F7 below, and the reason both routes are retained with
conditioning reported.

## F7 (MUST FIX). Report matrix conditioning before inverting

Rahimi & Debnath report `R_HI = 0.0018` and `RE = 0.199` for their
Pesek-Baker index against `R_HI ≈ 0.989` for the optimum and base
indices — an almost total collapse. The integrated review attributes
this to a desired-gain vector “poorly aligned with the defined aggregate
objective.”

That is probably not the main cause. With `b = G⁻¹d` on **unscaled**
maize traits — plant height in cm, grain number as counts, 100-grain
weight, yield — `G` is severely ill-conditioned, so `d'G⁻¹PG⁻¹d`
inflates, `σ_I` explodes, and `R_HI = σ_HI/(σ_I σ_H)` collapses toward
zero. The failure is numerical, and it is the same scaling problem as F1
and F3 appearing in its most destructive form.

This is the strongest available argument for standardising traits, and
it justifies a hard requirement: report the condition number of `P` and
`G`, and warn before any inversion when it exceeds a threshold. Cheap,
and it would have flagged the published result.

## H7 (MUST HAVE). Index scope is program × stage × product profile, not universal

The IRRI deck reports retrospective coefficients per region, and the
**same resistance-gene indicator carries opposite signs in different
programmes**:

| Trait  | Bangladesh |       ESA |     India | Philippines |
|--------|-----------:|----------:|----------:|------------:|
| YLD_BV |      2.905 |     4.764 |     3.371 |       2.015 |
| ZNC_BV |      0.457 |     0.629 |     0.460 |       1.199 |
| Xa21   |     −3.148 | **1.003** |    −3.533 |      −0.844 |
| Xa5    |  **0.751** |    −0.436 | **0.915** |      −0.492 |

An index is therefore scoped to a programme, a stage, and a product
profile. This upgrades N5 from nice-to-have: results should carry that
scope as first-class metadata, and reusing an index outside its scope
should at minimum warn. It also reinforces the guideline’s and both
decks’ instruction that coefficients must not be read as biological
importance.

## H1 revised. Retrospective weights require an explicit fine-tuning step

The IRRI deck quantifies this: doubling the Philippines yield
coefficient moved final gain from 2.230 to 3.772. Retrospective weights
answer *“what linear rule best approximates past decisions?”*, not
*“what should we do next?”*

[`retrospective_weights()`](https://FAkohoue.github.io/DesiredGainR/reference/retrospective_weights.md)
must therefore return a **starting point paired with a
scenario-comparison workflow**, never a final answer. Pair it with the
weight-sensitivity tooling (§1.4 of the design document) so the breeder
can see what changing an emphasis does before committing.

## D1 reinforced. Diversity constraints are expensive and must be priced explicitly

The IRRI family-balanced comparison:

| Region      | Historical | Unrestricted index | Family-balanced index |
|-------------|-----------:|-------------------:|----------------------:|
| Bangladesh  |      2.477 |              5.457 |                 2.833 |
| ESA         |      3.481 |              7.018 |                 4.038 |
| India       |      3.178 |              5.617 |                 3.468 |
| Philippines |      1.919 |              3.361 |                 2.069 |

Family balancing costs roughly half the nominal index gain. Diversity
control cannot be an optional extra in the multi-cycle optimiser — the
trade-off is first-order, and the Pareto output must include a diversity
axis (N_e, family concentration, or coancestry) alongside the trait-gain
axes.

## Presentation caution: relative efficiency below 1 is not failure

Every RE value in Rahimi & Debnath was below 1 — no multi-trait index
beat direct selection on yield alone. That is the expected trade-off,
not a defect. When DesiredGainR reports RE it must frame it that way, or
users will read `RE < 1` as a broken index. Report RE alongside the
correlated responses that the direct-selection comparator sacrifices.

## Sobering baseline: the optimum index may buy nothing

In Rahimi & Debnath the base index (`b = a`) and the optimum index
(`b = P⁻¹Ga`) correlated at **0.99979**, and the base index scored
slightly better on ΔH and RE. This is dataset-specific, but it argues
strongly for shipping the base index as a routine comparator (H3): when
it matches the optimum index, the covariance-estimation effort bought
nothing, and that is worth reporting rather than hiding.

## Where I read the evidence differently from the integrated review

- **Rahimi’s Pesek-Baker collapse** — the review calls it a misaligned
  desired-gain vector; I read it as a conditioning artifact of unscaled
  traits (see F7). The distinction matters, because the review’s reading
  suggests desired-gain indices are fragile in principle, whereas mine
  says they are fragile *when applied to unstandardised traits* — which
  is fixable and is exactly what the guideline already prescribes.
- **Joukhadar’s penalty** — the review describes `q` as averaged across
  traits; the paper writes it as a sum over `n` traits. Immaterial to
  ranking, but the implementation should follow the paper.
- **Feasibility** is not a distinct step in the review’s ten-step
  workflow, even though the guideline explicitly sends breeders to
  external software for it. It deserves its own step (see H4).
