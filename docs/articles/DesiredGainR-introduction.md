# Introduction to DesiredGainR

## 1. Overview

Computing a multi-trait selection index is arithmetic. Deciding what the
index should select for is the difficult part, and it is the part that
determines whether the resulting selection is any good.

The literature is consistent on this point. Guimarães et al. (2021)
describe the assignment of meaningful economic weights as the principal
bottleneck to routine index use, and demonstrate the consequence
empirically: when they supplied arbitrary weights to a Smith-Hazel and
to a Tai index, the gains those indices had delivered under sensible
weights disappeared entirely. Covarrubias-Pazaran (2021), writing the
CGIAR implementation guideline, reaches the same conclusion from the
operational side and advises breeders not to interpret index
coefficients at all, because the desired response is the only decision
genuinely available to them.

DesiredGainR is therefore organised around the objective rather than
around the index. It provides the classical index families, but its
distinctive content is the layer that helps a breeder state an
objective, test whether it can be attained, and discover how much the
resulting decision depends on it.

This vignette orients the reader. The other articles cover each layer in
depth, and Section 6 provides a direct link to each one.

------------------------------------------------------------------------

## 2. What the package does, and what it does not

DesiredGainR begins after the genetic evaluation. It consumes genetic
predictions and covariance matrices produced by a fitted multi-trait
model elsewhere, and it never estimates them from raw phenotypes.

That boundary is deliberate. A covariance matrix of raw phenotypes
contains residual and environmental variation, and using it as a genetic
covariance silently mislabels the estimand. The package therefore
requires `G` explicitly, and where no fitted matrix exists,
[`estimate_genetic_covariance()`](https://FAkohoue.github.io/DesiredGainR/reference/estimate_genetic_covariance.md)
provides labelled approximations that record what they actually are.

| DesiredGainR does | DesiredGainR does not |
|----|----|
| Map economic weights and desired gains algebraically | Fit the multi-trait genetic model |
| Test whether a stated objective is attainable | Estimate breeding values or marker effects |
| Build and compare supported selection methods | Select the operational parent set |
| Report the standard index-evaluation criteria | Analyse field trials or adjust for design |
| Simulate desired-gain directions over declared scenarios | Perform OCS, rank crosses or allocate matings |

The last exclusion matters for anyone using both of the author’s
packages. Parent selection, optimal contribution selection (OCS),
optimum cross selection and mating allocation belong to
[HapBlockR](https://github.com/FAkohoue/HapBlockR). HapBlockR calls
DesiredGainR through `build_selection_index()` when DGSI or QGSI is
requested, then uses the resulting merit score for the operational
breeding decision. DesiredGainR does not call HapBlockR. See [Working
with other breeding
software](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-interoperation.md).

------------------------------------------------------------------------

## 3. The four layers

### 3.1 Defining the objective

The Smith-Hazel economic index and the Pesek-Baker desired-gain index
both produce linear scores. However, they start from different breeding
statements. Smith-Hazel starts from aggregate-value weights. Pesek-Baker
starts from a desired response direction.

For admissible covariance matrices, each desired-gain direction has an
implied weight vector that reproduces the same coefficient direction.
This is an algebraic translation. It does not create economic values for
free. The translation helps a breeder inspect the trade-offs implied by
an objective. It also shows the response implied by a proposed weight
vector.

Around it sit a feasibility test, a method for recovering the objective
a programme has been using implicitly, and a sensitivity analysis that
reports whether the decision depends on the weights at all.

### 3.2 Building an index

Multiple-trait selection includes multistage selection, tandem
selection, independent culling, and index selection. These strategies
allow different types of trade-off. The breeder must define those
trade-offs before choosing a method.

[`selection_index()`](https://FAkohoue.github.io/DesiredGainR/reference/selection_index.md)
provides Smith-Hazel, the base index, both desired-gain formulations,
the Mulamba-Mock rank-sum index, Elston’s multiplicative index,
independent culling, and a within-cohort sequential screen. One
interface lets these alternatives use identical data and settings.
[`comparison_objective()`](https://FAkohoue.github.io/DesiredGainR/reference/comparison_objective.md)
fixes one yardstick for all fitted methods.
[`compare_selection_methods()`](https://FAkohoue.github.io/DesiredGainR/reference/compare_selection_methods.md)
then checks fairness and compares responses, target attainment,
rankings, and selected sets across the classical, restricted, general,
iterative desired-gain and quadratic genomic results.
[`run_dgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_dgsi.md)
applies the iterative optimisation procedure of Joukhadar et al. (2024)
to the established desired-gain index.
[`run_qgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_qgsi.md)
fits the quadratic genomic selection index of Cerón-Rojas et al. (2026).

### 3.3 Evaluating an index

[`evaluate_index()`](https://FAkohoue.github.io/DesiredGainR/reference/evaluate_index.md)
reports the correlation between the index and one fixed definition of
aggregate genetic merit, \\R\_{HI}\\. It also gives expected change in
that merit, \\\Delta H\\. Expected response is shown for every trait as
\\\Delta_j\\. Relative efficiency (RE) refers to one declared main
trait. Index heritability describes the composite score. The coefficient
of variation of index scores, \\CV_I\\, is available when the score has
a meaningful non-zero mean.

These quantities are not interchangeable measures of index quality. In
particular, \\R\_{HI}\\ and \\\Delta H\\ require the same
aggregate-merit definition for every index, RE describes only the main
trait, and \\CV_I\\ is undefined for a centred index. Therefore,
DesiredGainR also directs the breeder to per-trait attainment,
response-direction alignment, selected-set stability, covariance
uncertainty and, for repeated cycles, diversity. The detailed
interpretation is given in [Multiple-trait
selection](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-index-families.md)
and the breeder guide.

### 3.4 Looking ahead

A single-cycle response calculation cannot distinguish desired-gain
directions beyond what the achievable-response ellipsoid already states.
Over several cycles the ranking can change, because truncation selection
erodes the variance that response depends on. The simulation layer
represents that, using founders built from the breeder’s own phased
marker data rather than a simulated genome.

------------------------------------------------------------------------

## 4. The example programme

The shipped data describe a simulated tropical maize programme with six
traits. The material is patterned on a real programme structure but
represents no real germplasm and no real trial.

``` r
dgr_traits
#>   trait               description      unit direction heritability   mean
#> 1    GY               Grain yield      t/ha  increase         0.35   5.20
#> 2   PHT              Plant height        cm  decrease         0.60 205.00
#> 3    AD             Anthesis date      days  decrease         0.70  66.00
#> 4   ASI Anthesis-silking interval      days  decrease         0.30   1.80
#> 5   EPP            Ears per plant     count  increase         0.40   0.98
#> 6   GLS   Grey leaf spot severity score 1-9  decrease         0.45   3.90
#>   genetic_sd
#> 1       0.75
#> 2      12.00
#> 3       2.50
#> 4       0.60
#> 5       0.10
#> 6       0.85
```

Three properties make these data useful for demonstration, and the
previous example data lacked all three.

**The objective is internally antagonistic.** Grain yield correlates
positively with anthesis date, yet the programme wants yield to rise
while the cycle shortens.

``` r
round(stats::cov2cor(dgr_G)["GY", ], 2)
#>    GY   PHT    AD   ASI   EPP   GLS 
#>  1.00  0.25  0.30 -0.55  0.45 -0.35
```

**The trait scales differ by a factor of about one hundred and twenty**,
from ears per plant to plant height. That is what makes standardisation
and matrix conditioning matter rather than being formalities.

``` r
round(range(dgr_traits$genetic_sd), 2)
#> [1]  0.1 12.0
```

**The population is internally consistent.** Genetic values were
generated from ninety quantitative trait loci drawn from the shipped
marker panel, then adjusted so that their realised covariance equals
`dgr_G` to numerical precision, and heritabilities are recoverable
exactly.

``` r
round(diag(dgr_G) / diag(dgr_P), 2)
#>   GY  PHT   AD  ASI  EPP  GLS 
#> 0.35 0.60 0.70 0.30 0.40 0.45
```

The nine datasets are documented together under
[`?"DesiredGainR-data"`](https://FAkohoue.github.io/DesiredGainR/reference/DesiredGainR-data.md).

------------------------------------------------------------------------

## 5. A short worked path

The full narrative is in the *Complete workflow* vignette. This is the
shape of it in four calls.

**Check the covariance matrices before building anything on them.**
Index coefficients come from inverting these matrices, and an
ill-conditioned matrix produces coefficients that are numerically
meaningless while looking plausible.

``` r
matrix_diagnostics(dgr_G, "G")$condition_number
#> [1] 18569.92
matrix_diagnostics(stats::cov2cor(dgr_G))$condition_number
#> [1] 8.75443
```

Almost all of that conditioning comes from the trait scales rather than
the correlations, which is why the package standardises by default.

**State the objective and see what it implies.**

``` r
desired_gains <- c(
  GY = 1.0, PHT = 0.4, AD = 0.6, ASI = 0.5,
  EPP = 0.4, GLS = 0.6
)
round(
  implied_economic_weights(
    desired_gains, dgr_G, dgr_P,
    lower_is_better = lower_is_better, gain_units = "genetic_sd"
  ),
  2
)
#>     GY    PHT     AD    ASI    EPP    GLS 
#>  14.52   0.03   2.32 -13.21 -15.15   0.68 
#> attr(,"provenance")
#> [1] "Implied by the supplied desired gains through w = G^-1 P G^-1 d; not an independently estimated economic value. Expressed in the favourable-direction space, so a positive weight favours movement in the breeder-defined direction."
```

**Ask whether it can be attained.**

``` r
gain_feasibility(
  desired_gains, dgr_G, dgr_P,
  n_candidates = nrow(dgr_candidates), n_select = 20L,
  lower_is_better = lower_is_better, gain_units = "genetic_sd"
)
#> <desiredgainr_feasibility>
#>   Required selection intensity: 3.2082
#>   Requires the top 0.1773%, which is fewer than one of 200 candidates
#>   Planned intensity: 1.7550 (top 10.0%)
#>   Feasible at planned intensity: no
#>   Feasible anywhere in this population: no
#>   Attainable fraction of the requested gain: 54.7%
```

It cannot, and not marginally: the target would require selecting fewer
than one candidate from two hundred. Roughly fifty-five per cent of the
requested gain is available at the intensity actually planned. Knowing
that before optimising is worth more than any subsequent refinement of
the index.

**Build the index.**

``` r
selection_index(
  dgr_candidates, traits,
  method = "smith_hazel",
  G = dgr_G, P = dgr_P,
  economic_weights = c(
    GY = 1.0, PHT = 0.2, AD = 0.5,
    ASI = 0.4, EPP = 0.3, GLS = 0.5
  ),
  lower_is_better = lower_is_better, n_select = 20L, main_trait = "GY"
)
#> <desiredgainr_index>
#>   Method: smith_hazel 
#>   Candidates: 200  Traits: 6 
#>   Traits standardised: yes 
#>   Selected: 20 (10.0%), intensity 1.755
#>   Coefficients:
#>     GY    PHT     AD    ASI    EPP    GLS 
#> 0.3550 0.0735 0.2695 0.2521 0.2013 0.2919 
#>   R_HI 0.7089  dH 1.2228  RE 0.7658
#>   CV_I undefined for a centred index
```

------------------------------------------------------------------------

## 6. Which vignette to read next

| If you want to | Read |
|----|----|
| Follow one programme end to end | [The complete workflow](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-workflow.md) |
| Decide what to select for | [Defining a breeding objective](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-objective.md) |
| Define ranges when exact gains are uncertain | [Defining desired gains and acceptable intervals](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-desired-gain-intervals.md) |
| Obtain or judge a covariance matrix | [Obtaining G and P](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-covariance.md) |
| Choose and compare multiple-trait methods | [Multiple-trait selection](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-index-families.md) |
| Use GEBVs, prediction errors, or distinct information traits | [Using predictions, prediction errors, and restricted responses](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-information.md) |
| Understand what iteration changes in a desired-gain index | [Iterative optimisation of the desired-gain index](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-dgsi.md) |
| Use the quadratic genomic index | [Quadratic genomic selection](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-qgsi.md) |
| Compare objectives over several cycles | [Multi-cycle simulation](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-simulation.md) |
| Inspect validation with real programmes | [Empirical validation](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-empirical-validation.md) |
| Combine this with other breeding software | [Working with other packages](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-interoperation.md) |

A plain-language companion for readers who will approve a selection
decision but do not run R is available through
[`open_desiredgain_guide()`](https://FAkohoue.github.io/DesiredGainR/reference/open_desiredgain_guide.md).

------------------------------------------------------------------------

## 7. Three conventions worth knowing before you start

**Direction is declared, not signed.** Every entry point accepts
`lower_is_better`, naming the traits that improve by falling. Gains and
weights are then stated as improvements throughout, and the package
orients the covariance matrices internally. Signing values by hand is
the commonest way to obtain a confidently wrong answer.

**Units are stated explicitly.** Desired gains can be given in original
trait units, genetic standard deviations, or phenotypic standard
deviations, through `gain_units`. Results are returned in the units they
were asked for.

**Coefficients are not comparable across traits.** A coefficient carries
inverse trait units, so a trait measured on a small scale receives a
large number for the same emphasis.
[`effective_weights()`](https://FAkohoue.github.io/DesiredGainR/reference/effective_weights.md)
reports the interpretable form, the coefficient multiplied by the
trait’s genetic standard deviation.

------------------------------------------------------------------------

## 8. References

- Cerón-Rojas JJ, Montesinos-López OA, Montesinos-López A, et
  al. (2026). Nonlinear genomic selection index accelerates multi-trait
  crop improvement. *Nature Communications* **17**:1991.
  <https://doi.org/10.1038/s41467-026-69890-3>
- Covarrubias-Pazaran G (2021). *Practical implementation of selection
  indices.* CGIAR Excellence in Breeding.
- Guimarães PHR, Melo PGS, Cordeiro ACC, Torga PP, Rangel PHN, de Castro
  AP (2021). Index selection can improve the selection efficiency in a
  rice recurrent selection population. *Euphytica* **217**:95.
  <https://doi.org/10.1007/s10681-021-02819-7>
- Joukhadar R, Li Y, Thistlethwaite R, Forrest KL, Tibbits JF, Trethowan
  R, Hayden MJ (2024). Optimising desired gain indices to maximise
  selection response. *Frontiers in Plant Science* **15**:1337388.
  <https://doi.org/10.3389/fpls.2024.1337388>
- Rahimi M, Debnath S (2023). Estimating optimum and base selection
  indices in plant and animal breeding programs by development new and
  simple SAS and R codes. *Scientific Reports* **13**:18977.
  <https://doi.org/10.1038/s41598-023-46368-6>
