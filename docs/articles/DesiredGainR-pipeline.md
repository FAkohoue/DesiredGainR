# Full pipeline: from trait values to a defended selection decision

### Purpose

This vignette assumes no prior knowledge of DesiredGainR. It runs the
entire package on the shipped example programme, one stage at a time,
showing every call, every output, and the decision each one feeds.
Nothing is skipped and nothing is left to be inferred from another
document.

Sixteen stages in four phases:

| Phase | Stages | Question answered |
|----|----|----|
| **Evidence** | 1-4 | What do we have, and can it be trusted? |
| **Objective** | 5-8 | What should we select for, and is it possible? |
| **Index** | 9-12 | Which index, and what does it actually do? |
| **Looking ahead** | 13-16 | Does this objective still hold over several cycles? |

Every stage names the vignette that carries its full statistical detail.
This document is the connective walkthrough. Those are the references.

The walkthrough stops before the operational crossing decision.
Throughout this vignette, `n_select` is an analytical truncation count
used to calculate selection intensity and response, while `n_parents`
and `n_crosses` in the simulation are scenario settings. None of them is
an optimised parent or mating plan. Use
[HapBlockR](https://github.com/FAkohoue/HapBlockR) after the index has
been built to choose parents, perform OCS, rank crosses and allocate
matings.

------------------------------------------------------------------------

## Phase 1: Evidence

### Stage 1 — What you need before you start

DesiredGainR begins *after* the genetic evaluation. It does not analyse
field trials, does not fit a mixed model, and does not estimate breeding
values. You must arrive with four things.

| Input | What it is | Where it comes from |
|----|----|----|
| Trait values | One adjusted mean, best linear unbiased prediction (BLUP) or genomic estimated breeding value (GEBV) per candidate per trait | Your multi-trait genetic evaluation |
| \\\mathbf{G}\\ | Genetic variance-covariance matrix | The same fitted model |
| \\\mathbf{P}\\ | Phenotypic variance-covariance matrix | The same fitted model |
| Trait directions | Which traits improve by rising, which by falling | The product profile |

If you have no fitted \\\mathbf{G}\\, Stage 3 shows the supported
approximations. If you have no trait directions written down, stop and
write them down: every later stage depends on them.

### Stage 2 — Load and inspect the data

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

Read this table before anything else. It carries the units, the
direction of improvement, the heritability and the genetic standard
deviation for each trait. Two columns decide most of what follows.

``` r
traits <- dgr_traits$trait
lower_is_better <- dgr_traits$trait[dgr_traits$direction == "decrease"]
lower_is_better
#> [1] "PHT" "AD"  "ASI" "GLS"
```

`lower_is_better` is passed to nearly every function in the package.
Supplying it means you state gains as *improvements* everywhere, and the
package handles the sign conventions internally. Signing values by hand
is the commonest source of a confidently wrong answer.

``` r
str(dgr_candidates, max.level = 1)
#> 'data.frame':    200 obs. of  8 variables:
#>  $ GenoID: chr  "CAND001" "CAND002" "CAND003" "CAND004" ...
#>  $ Family: chr  "F01" "F01" "F01" "F01" ...
#>  $ GY    : num  4.38 3.71 4.36 6.29 7.73 ...
#>  $ PHT   : num  200 189 182 190 216 ...
#>  $ AD    : num  67.3 67 65.4 59.7 63.6 ...
#>  $ ASI   : num  2.695 1.967 0.949 2.139 1.293 ...
#>  $ EPP   : num  1.154 0.752 0.891 1.061 1.091 ...
#>  $ GLS   : num  3.59 3.8 4.58 4.8 2.24 ...
head(dgr_candidates[, c("GenoID", "GY", "PHT", "GLS")], 3)
#>          GenoID     GY      PHT    GLS
#> CAND001 CAND001 4.3844 200.0641 3.5937
#> CAND002 CAND002 3.7096 188.6819 3.7963
#> CAND003 CAND003 4.3630 181.9463 4.5780
```

### Stage 3 — Check the covariance matrices

Covariance-based linear indices solve systems involving \\\mathbf{P}\\,
\\\mathbf{G}\\, or both. An ill-conditioned matrix can produce unstable
coefficients while the output still looks plausible. Therefore, check
the matrices before fitting a covariance-based index. Rank, threshold
and quadratic genomic methods use different inputs, which require their
own diagnostics.

``` r
diagnostics_G <- matrix_diagnostics(dgr_G, "G")
c(
  condition_number = diagnostics_G$condition_number,
  positive_definite = diagnostics_G$positive_definite,
  numerical_rank = diagnostics_G$numerical_rank
)
#>  condition_number positive_definite    numerical_rank 
#>          18569.92              1.00              6.00
```

A condition number near \\1.9 \times 10^{4}\\ looks alarming until it is
decomposed. Almost all of it comes from the trait scales rather than
from the correlations:

``` r
matrix_diagnostics(stats::cov2cor(dgr_G))$condition_number
#> [1] 8.75443
```

Standardising the traits removes three orders of magnitude of
ill-conditioning. This is the practical form of an observation by
Crosbie et al. (1980): equal weights applied to unstandardised traits
concentrate selection on whichever trait happens to carry the largest
variance. It is why DesiredGainR standardises by default.

**If you have no fitted \\\mathbf{G}\\.** Use
[`estimate_genetic_covariance()`](https://FAkohoue.github.io/DesiredGainR/reference/estimate_genetic_covariance.md),
which returns labelled approximations rather than relabelling a
covariance of predictions as a genetic covariance:

``` r
G_working <- estimate_genetic_covariance(
  values = dgr_candidates, trait_cols = traits,
  method = "adjusted_means_surrogate"
)
G_working$estimand
#> [1] "working surrogate for genetic covariance. covariance of across-environment adjusted means"
```

Four methods are available: `prediction_covariance`, `pev_corrected`,
`se_diagonal_corrected`, `relationship_adjusted`, and
`adjusted_means_surrogate`. The last implements what the CGIAR
implementation guideline recommends for operational use when no fitted
matrix exists.

### Stage 4 — Confirm units and heritabilities

The declared heritabilities must be recoverable from the two matrices.
If they are not, the metadata and the matrices are describing different
things.

``` r
round(diag(dgr_G) / diag(dgr_P), 2)
#>   GY  PHT   AD  ASI  EPP  GLS 
#> 0.35 0.60 0.70 0.30 0.40 0.45
dgr_traits$heritability
#> [1] 0.35 0.60 0.70 0.30 0.40 0.45
```

Check the spread of trait scales, because it determines whether
standardisation is optional or essential:

``` r
genetic_sd <- stats::setNames(dgr_traits$genetic_sd, dgr_traits$trait)
round(sort(genetic_sd, decreasing = TRUE), 2)
#>   PHT    AD   GLS    GY   ASI   EPP 
#> 12.00  2.50  0.85  0.75  0.60  0.10
round(max(genetic_sd) / min(genetic_sd))
#> [1] 120
```

A spread of 120-fold means standardisation is essential here. The
package warns automatically above fivefold.

------------------------------------------------------------------------

## Phase 2: Objective

### Stage 5 — State a first objective

There are two common ways to state a linear breeding objective. They are
algebraically related, but they have different biological meanings.

**Economic weights** \\\mathbf{w}\\ say what a unit of each trait is
worth. **Desired gains** \\\mathbf{d}\\ say how much of each trait you
want. The Smith-Hazel index uses the first, \\\mathbf{b} =
\mathbf{P}^{-1}\mathbf{G}\mathbf{w}\\. the Pesek-Baker desired-gain
index uses the second, \\\mathbf{b} =
\mathbf{P}^{-1}\mathbf{G}(\mathbf{G}\mathbf{P}^{-1}\mathbf{G})^{-1}\mathbf{d}\\.

Most breeders find desired gains easier, so start there. State them as
improvements, in genetic standard deviations, and let `lower_is_better`
handle direction.

``` r
desired_gains <- c(
  GY = 1.0, PHT = 0.4, AD = 0.6, ASI = 0.5,
  EPP = 0.4, GLS = 0.6
)
```

This says: one standard deviation of extra yield, four tenths of a
standard deviation shorter plants, six tenths earlier flowering, and so
on.

### Stage 6 — Calculate the algebraically implied weights

Setting the two coefficient formulae equal gives an exact algebraic
translation in both directions:

\\\mathbf{w} = \mathbf{G}^{-1}\mathbf{P}\mathbf{G}^{-1}\mathbf{d},
\qquad \mathbf{d} = \mathbf{G}\mathbf{P}^{-1}\mathbf{G}\mathbf{w}.\\

``` r
implied <- implied_economic_weights(
  desired_gains, dgr_G, dgr_P,
  lower_is_better = lower_is_better, gain_units = "genetic_sd"
)
round(implied, 2)
#>     GY    PHT     AD    ASI    EPP    GLS 
#>  14.52   0.03   2.32 -13.21 -15.15   0.68 
#> attr(,"provenance")
#> [1] "Implied by the supplied desired gains through w = G^-1 P G^-1 d; not an independently estimated economic value. Expressed in the favourable-direction space, so a positive weight favours movement in the breeder-defined direction."
```

These are implied weights. They reproduce the desired-gain coefficient
direction under the supplied covariance matrices. They are not observed
prices, profit coefficients, or breeder preferences. Use them to
diagnose the mathematical trade-offs. Use independently elicited
economic weights when the objective is aggregate economic merit.

**Two of these are negative, for traits you asked to improve.** That is
correct. Once oriented, grain yield correlates \\+0.55\\ with the
anthesis-silking interval and \\+0.45\\ with ears per plant, so
selecting hard for yield would carry both *past* the smaller gains you
requested. To deliver the ratio you asked for, the index must hold them
back. A negative weight never means the trait should get worse.

**The desired gains are already standardised. The returned weights use a
different scale.** Here `desired_gains` is expressed in genetic standard
deviations, but `dgr_G` and `dgr_P` are supplied in the original trait
units. Therefore,
[`implied_economic_weights()`](https://FAkohoue.github.io/DesiredGainR/reference/implied_economic_weights.md)
first converts the target to those original units and returns a weight
per original unit of each trait. Compare coefficients only after
conversion to a common trait scale.

The multiplication below does not standardise `desired_gains` a second
time. It converts each returned raw-unit weight into the contribution
associated with one genetic standard deviation of its trait. These
products are on a common genetic-standard-deviation scale and may be
compared as a diagnostic:

``` r
round(sort(implied * genetic_sd[names(implied)], decreasing = TRUE), 2)
#>    GY    AD   GLS   PHT   EPP   ASI 
#> 10.89  5.81  0.58  0.36 -1.51 -7.93
```

Ears per plant had the largest raw coefficient but one of the smaller
per-genetic-standard-deviation effects. If `G` and `P` had already been
standardised to genetic-standard-deviation units, the returned weights
would already be per genetic standard deviation and this multiplication
would be wrong. Standardise either the covariance space or the returned
weights, not both.

Confirm the translation inverts, which costs one line and catches unit
errors:

``` r
round(
  implied_desired_gains(
    implied, dgr_G, dgr_P,
    lower_is_better = lower_is_better, gain_units = "genetic_sd"
  ),
  3
)
#>  GY PHT  AD ASI EPP GLS 
#> 1.0 0.4 0.6 0.5 0.4 0.6 
#> attr(,"provenance")
#> [1] "Implied by the supplied economic weights through d = G P^-1 G w, expressed in genetic_sd units and up to the scalar set by selection intensity."
```

### Stage 7 — Test whether the objective is attainable

**This is the stage most often skipped, and the one that most often
changes the answer.**

Let \\\mathbf{b}\\ denote the linear-index coefficients, \\\mathbf{G}\\
the additive genetic covariance matrix, \\\mathbf{P}\\ the phenotypic
covariance matrix, and \\i\\ the standardised selection intensity. The
expected one-cycle response is

\\\mathbf{R} = i\\\frac{\mathbf{G}\mathbf{b}}
{\sqrt{\mathbf{b}^\mathsf{T}\mathbf{P}\mathbf{b}}}.\\

Consequently, every response available at intensity \\i\\ lies on the
achievable-response ellipsoid

\\\mathbf{R}^\mathsf{T}\mathbf{G}^{-1}\mathbf{P}
\mathbf{G}^{-1}\mathbf{R}=i^2.\\

A specified vector \\\mathbf{d}\\ lies on an ellipsoid whose intensity
is

\\i\_{\text{required}} =
\sqrt{\mathbf{d}^\mathsf{T}\mathbf{G}^{-1}\mathbf{P}\mathbf{G}^{-1}\mathbf{d}}.\\

``` r
feasibility <- gain_feasibility(
  desired_gains, dgr_G, dgr_P,
  n_candidates = nrow(dgr_candidates), n_select = 20L,
  lower_is_better = lower_is_better, gain_units = "genetic_sd"
)
feasibility
#> <desiredgainr_feasibility>
#>   Required selection intensity: 3.2082
#>   Requires the top 0.1773%, which is fewer than one of 200 candidates
#>   Planned intensity: 1.7550 (top 10.0%)
#>   Feasible at planned intensity: no
#>   Feasible anywhere in this population: no
#>   Attainable fraction of the requested gain: 54.7%
```

Read the output line by line. The required intensity of 3.2082
corresponds, under normal truncation selection, to retaining the best
0.1773% of an effectively continuous population. In a population of 200
candidates, this is 0.355 candidate. A breeder cannot select a fraction
of a candidate. Even selecting the single best candidate retains 0.5%,
which gives a lower intensity than 3.2082. Therefore, the exact vector
is unattainable in this candidate set.

The planned decision retains 20 of 200 candidates, or 10%. Its intensity
is 1.7550. The ratio \\1.7550/3.2082=0.547\\ means that 54.7% of the
requested vector can be delivered while preserving its exact trait
proportions:

``` r
round(feasibility$attainable_response_input_units, 3)
#>    GY   PHT    AD   ASI   EPP   GLS 
#> 0.547 0.219 0.328 0.274 0.219 0.328
```

The returned values are in the same genetic-standard-deviation units as
the input. Thus, the planned exact-ray response is 0.547 standard
deviations for grain yield, 0.219 for plant height, 0.328 for anthesis
date, 0.274 for the anthesis-silking interval, 0.219 for ears per plant,
and 0.328 for grey leaf spot resistance. They are not observed gains and
they are not six independent upper limits. They are the model-based
response at the planned intensity when the requested **ratio among
traits is held exactly**.

This interpretation has four consequences.

1.  The classical desired-gain index honours the direction, or ratio, of
    \\\mathbf{d}\\. Multiplying all elements by the same positive
    constant leaves the ranking unchanged. It changes only the intensity
    required to attain the stated magnitude.
2.  `feasible_at_planned_intensity = no` still allows the index to rank
    candidates. The complete stated vector lies beyond one-cycle
    response at the planned 10% selection proportion under the supplied
    \\\mathbf{G}\\ and \\\mathbf{P}\\.
3.  An exact desired-gain vector is not the same objective as a set of
    minimum gains. For example, another attainable response may exceed
    the grain-yield target and still meet every minimum without
    preserving the ratio \\1.0:0.4:0.6:0.5:0.4:0.6\\.
    [`gain_feasibility()`](https://FAkohoue.github.io/DesiredGainR/reference/gain_feasibility.md)
    tests the exact ray. Population-driven and interval searches assess
    trait-specific floors and must be interpreted through their
    worst-trait and joint-attainment results.
4.  The result is conditional on a linear index, one-cycle truncation
    selection, approximate normality, and the supplied covariance
    matrices. Therefore, covariance uncertainty and selected-set
    stability must be checked before a recommendation is treated as
    robust.

If the ratio itself is biologically required, scale the whole vector by
0.547 or relax the selection proportion. If each entry is a minimum
rather than a fixed ratio, use
[`suggest_desired_gains()`](https://FAkohoue.github.io/DesiredGainR/reference/suggest_desired_gains.md)
or interval optimisation and require the lower confidence bound for
joint attainment to support the recommendation.

### Stage 8 — Recover what the programme already selects for, then adjust

If your programme has been selecting without a formal index, the
differentials it achieved already encode its objective:

\\\mathbf{b} = \mathbf{P}^{-1}\mathbf{s}\\

where \\\mathbf{s}\\ holds the differences between the selected group
and the whole population.

``` r
recovered <- retrospective_weights(
  selected_values = dgr_candidates[dgr_history$selected, traits],
  population_values = dgr_candidates[, traits],
  trait_cols = traits
)
recovered
#> <desiredgainr_retrospective>
#>   40 selected from 200 candidates (20.0%)
#>   P: estimated from population_values 
#>   Recovered coefficients, per trait standard deviation:
#>      GY     PHT      AD     ASI     EPP     GLS 
#>  0.9324 -0.3000 -0.3986 -0.5244  0.1508 -0.3975 
#>   (Raw coefficients are in $coefficients; they carry inverse trait units and are not comparable across traits.)
```

``` r
round(recovered$selection_differential_sd, 2)
#>    GY   PHT    AD   ASI   EPP   GLS 
#>  1.11 -0.27 -0.42 -0.83  0.52 -0.42
```

This is how an index was introduced at both the International Maize and
Wheat Improvement Center (CIMMYT) and the International Rice Research
Institute (IRRI). Treat it as a starting point: it reproduces past
decisions including any bias in them, and must be adjusted deliberately
against the product profile before it becomes a forward-looking
objective.

------------------------------------------------------------------------

## Phase 3: Index

### Stage 9 — Build candidate indices

[`selection_index()`](https://FAkohoue.github.io/DesiredGainR/reference/selection_index.md)
provides six classical index families and two operational comparators
through one interface.

``` r
economic_weights <- c(
  GY = 1.0, PHT = 0.2, AD = 0.5,
  ASI = 0.4, EPP = 0.3, GLS = 0.5
)

smith_hazel <- selection_index(
  dgr_candidates, traits,
  method = "smith_hazel",
  G = dgr_G, P = dgr_P, economic_weights = economic_weights,
  lower_is_better = lower_is_better, n_select = 20L, main_trait = "GY"
)
rank_sum <- selection_index(
  dgr_candidates, traits,
  method = "mulamba_mock",
  lower_is_better = lower_is_better, n_select = 20L
)
base_index <- selection_index(
  dgr_candidates, traits,
  method = "base",
  G = dgr_G, P = dgr_P, economic_weights = economic_weights,
  lower_is_better = lower_is_better, n_select = 20L, main_trait = "GY"
)
smith_hazel
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

| Method | Needs weights | Needs \\\mathbf{G}\\, \\\mathbf{P}\\ | Use when |
|----|----|----|----|
| `smith_hazel` | yes | yes | Weights are defensible |
| `base` | yes | no | A quick comparator. It is often nearly as good. |
| `pesek_baker` | desired gains | \\\mathbf{G}\\ only | Gains easier than weights |
| `yamada` | desired gains | yes | Same index, better conditioned route |
| `mulamba_mock` | no | no | Weights unavailable or untrusted |
| `independent_culling` | thresholds | no | Non-compensatory requirements |
| `elston` | thresholds | no | Firm limits followed by balanced margins |
| `tandem` | trait order | no | Within-cohort sequential-screen comparator |
| [`restricted_index()`](https://FAkohoue.github.io/DesiredGainR/reference/restricted_index.md) | weights and restrictions | yes | One or more responses must be controlled |
| [`run_dgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_dgsi.md) | desired gains | yes | Iterative selected-differential search is required |
| [`run_qgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_qgsi.md) | linear and quadratic values | genomic covariance | Merit is curved or interactive |

### Stage 10 — Evaluate and compare

[`evaluate_index()`](https://FAkohoue.github.io/DesiredGainR/reference/evaluate_index.md)
reports the classical criteria used by Rahimi and Debnath (2023),
together with the response of every trait. Each quantity answers a
different question.

| Criterion | Full definition | Valid interpretation |
|----|----|----|
| \\R\_{HI}\\ | Correlation between the selection index \\I\\ and aggregate genetic merit \\H\\ | It ranges from -1 to 1. Negative values select against the stated merit, zero means no linear association, and values nearer 1 indicate closer agreement. It exists only when the same aggregate-weight vector defines \\H\\ for every index compared. |
| \\\Delta H\\ | Expected change in aggregate genetic merit | This is \\iR\_{HI}\sigma_H\\, where \\\sigma_H\\ is the genetic standard deviation of aggregate merit. Positive values improve the stated merit. Compare magnitudes only when all indices use the same merit definition, population and selection intensity. |
| \\\Delta_j\\ | Expected genetic response of trait \\j\\ | This is the primary biological result. Check its sign, magnitude and units for every trait. An aggregate statistic can conceal an unfavourable response. |
| RE | Relative efficiency for the declared main trait | RE = 1 matches direct response, 0.80 retains 80%, a value above 1 exceeds direct response through correlated information, and a negative value moves the main trait unfavourably. It says nothing about the other traits. |
| \\h_I^2\\ | Heritability of the index score treated as a composite trait | It ranges from 0 to 1 for compatible covariance matrices. Its square root is the accuracy with which the observed index predicts its own additive genetic component. It does not measure agreement with the breeding objective. |
| \\CV_I\\ | Coefficient of variation of index scores | This is \\100s_I/\|\bar I\|\\. It depends on the arbitrary zero of the score and is undefined for a centred index. Use it only for raw, uncentred scores. |

The software output name `delta_H` and the abbreviated print label `dH`
both mean \\\Delta H\\, the expected change in aggregate genetic merit.
They do not mean the response of an individual trait.

No classical scalar is sufficiently robust on its own. Start with
\\\Delta_j\\ for every trait. Then test exact-ray feasibility. For
minimum floors, inspect the worst-trait and joint-attainment
probabilities. Measure response alignment with
[`index_uncertainty()`](https://FAkohoue.github.io/DesiredGainR/reference/index_uncertainty.md).
Inspect rank correlation and selected-set overlap under perturbation.
For repeated cycles, finish with diversity or coancestry diagnostics.
This panel covers biological attainment, statistical stability, and the
retention of usable variation.

**Relative efficiency below 1 is not a failure.** It is the intended
trade of response in the main trait for response elsewhere. Every value
reported by Rahimi and Debnath was below 1. Conversely, a high relative
efficiency is not proof that the desired multi-trait objective was
attained.

Use
[`comparison_objective()`](https://FAkohoue.github.io/DesiredGainR/reference/comparison_objective.md)
to define one target and utility. Then use
[`compare_selection_methods()`](https://FAkohoue.github.io/DesiredGainR/reference/compare_selection_methods.md)
to verify the common comparison conditions and compare every supported
fitted family. Read response basis carefully. A QGSI with non-zero
curvature has approximate per-trait gains. The \\\mathbf{W}=0\\ special
case and the other linear-index responses are exact under their stated
covariance models. See [Multiple-trait
selection](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-index-families.md).

``` r
comparison <- compare_selection_methods(
  list(
    Smith_Hazel = smith_hazel,
    Base = base_index,
    Rank_sum = rank_sum
  )
)
comparison$fairness
#>                                             Condition Satisfied
#> 1                Same candidates and objective traits      TRUE
#> 2                    Same favourable trait directions      TRUE
#> 3                               Common response units      TRUE
#> 4                                Same number selected      TRUE
#> 5                            Same selection intensity      TRUE
#> 6 One common genetic covariance for response geometry     FALSE
#> 7             Common genetic covariance is invertible      TRUE
#> 8        Validation identifiers aligned when supplied      TRUE
#>                                                                        Interpretation
#> 1                                             Required and enforced by this function.
#> 2                                             Required and enforced by this function.
#> 3                                              Responses are reported in trait units.
#> 4                            Hard culling can retain fewer candidates than requested.
#> 5                        Model-based responses require one common selection pressure.
#> 6                Alignment uses fitted genetic covariance ; disagreement is reported.
#> 7 Mahalanobis alignment and Satoh projection require an invertible covariance matrix.
#> 8 No validation data were supplied. External validation evidence remains unavailable.
comparison$summary
#>         Method       Family Strategy N_selected Selected_fraction
#>         <char>       <char>   <char>      <int>             <num>
#> 1: Smith_Hazel  smith_hazel    index         20               0.1
#> 2:        Base         base    index         20               0.1
#> 3:    Rank_sum mulamba_mock rank_sum         20               0.1
#>    Selection_intensity      R_HI  Delta_H        RE Index_heritability
#>                  <num>     <num>    <num>     <num>              <num>
#> 1:            1.754983 0.7088841 1.222814 0.7658206          0.5124609
#> 2:            1.754983 0.6950918 1.199023 0.8431941          0.4831527
#> 3:            1.754983        NA       NA        NA                 NA
#>    Index_accuracy  MSPE Common_merit_response Common_merit_correlation
#>             <num> <num>                 <num>                    <num>
#> 1:      0.7158637    NA                    NA                       NA
#> 2:      0.6950918    NA                    NA                       NA
#> 3:             NA    NA                    NA                       NA
#>    Worst_expected_attainment Mean_expected_attainment Worst_observed_attainment
#>                        <num>                    <num>                     <num>
#> 1:                        NA                       NA                        NA
#> 2:                        NA                       NA                        NA
#> 3:                        NA                       NA                        NA
#>    Satoh_beta Mahalanobis_alignment Mahalanobis_residual
#>         <num>                 <num>                <num>
#> 1:         NA                    NA                   NA
#> 2:         NA                    NA                   NA
#> 3:         NA                    NA                   NA
#>    Validation_utility_response Validation_utility_rank_correlation
#>                          <num>                               <num>
#> 1:                          NA                                  NA
#> 2:                          NA                                  NA
#> 3:                          NA                                  NA
#>        Expected_response_basis
#>                         <char>
#> 1: exact linear-index response
#> 2: exact linear-index response
#> 3:                 unavailable
round(comparison$rank_correlation, 2)
#>             Smith_Hazel Base Rank_sum
#> Smith_Hazel        1.00 0.98     0.94
#> Base               0.98 1.00     0.91
#> Rank_sum           0.94 0.91     1.00
comparison$selected_overlap
#>             Smith_Hazel Base Rank_sum
#> Smith_Hazel          20   16       13
#> Base                 16   20       14
#> Rank_sum             13   14       20
```

Rank correlation uses every candidate. Selected-set overlap concerns the
final decision. Report both because close overall rankings can diverge
near the selection boundary.

Check which traits the index is really acting on:

``` r
smith_hazel$effective_weights[, c("Trait", "Genetic_share")]
#>     Trait Genetic_share
#>    <char>         <num>
#> 1:     GY    0.20205605
#> 2:    PHT    0.06339965
#> 3:     AD    0.23638957
#> 4:    ASI    0.15198589
#> 5:    EPP    0.12788741
#> 6:    GLS    0.21828143
```

### Stage 11 — Optimise the desired-gain vector

[`run_dgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_dgsi.md)
applies the iterative optimisation procedure of Joukhadar et al. (2024)
to the established desired-gain index. It samples desired-gain vectors,
builds the same index formula from each vector, and keeps the one whose
*realised* response in the selected set comes closest to your target.

The target above uses genetic-standard-deviation units.
[`run_dgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_dgsi.md)
uses candidate-standard-deviation units. Convert once before fitting.

``` r
candidate_sd <- vapply(
  dgr_candidates[, traits], stats::sd, numeric(1L)
)
target_original_units <- desired_gains * genetic_sd
dg_target <- target_original_units / candidate_sd

dgsi <- run_dgsi(
  init_data = dgr_candidates[, c("GenoID", "Family")],
  cand_data = dgr_candidates,
  trait_cols = traits,
  dg = dg_target,
  G = dgr_G, P = dgr_P,
  lower_is_better = lower_is_better,
  scale_traits = TRUE,
  n_select = 20L, n_iter = 200L, n_rep = 5L, seed = 42L
)
round(dgsi$realised_response, 3)
#>    GY   PHT    AD   ASI   EPP   GLS 
#> 0.961 0.549 1.357 0.283 0.277 0.287
```

Read the deviation from target against Stage 7:

``` r
round(dgsi$realised_response - dg_target, 3)
#>     GY    PHT     AD    ASI    EPP    GLS 
#>  0.422  0.221  0.858 -0.003  0.036 -0.139
```

`realised_response` is the selected-candidate differential in candidate
standard deviations. It is the quantity used by the iterative objective.
It is not transmitted genetic gain. Compare the model-based transmitted
response on the original genetic-standard-deviation scale as follows.

``` r
transmitted_genetic_sd <-
  dgsi$theoretical_response$original_units / genetic_sd
round(transmitted_genetic_sd, 3)
#>     GY    PHT     AD    ASI    EPP    GLS 
#>  0.259 -0.526 -1.014 -0.551  0.078  0.011
round(transmitted_genetic_sd / desired_gains, 3)
#>     GY    PHT     AD    ASI    EPP    GLS 
#>  0.259 -1.316 -1.690 -1.101  0.195  0.018
```

Feasibility constrains the requested response ray. The iterative search
is not forced to remain on that ray. Therefore, inspect every
trait-specific attainment value. Whether an uneven shortfall is
acceptable is a breeding decision.

Always inspect the replicate diagnostics. The search is stochastic, so
agreement across replicates is the evidence that the answer is stable:

``` r
dgsi$replicate_diagnostics
#>    Replicate Objective Iteration_of_best Selected Plateau
#>        <int>     <num>             <int>    <int>  <lgcl>
#> 1:         1  6.753815               170       20   FALSE
#> 2:         2  4.696348               196       20   FALSE
#> 3:         3  6.016223               172       20   FALSE
#> 4:         4  8.613511               162       20   FALSE
#> 5:         5  5.133730                93       20    TRUE
#>    Final_window_relative_improvement Chosen
#>                                <num> <lgcl>
#> 1:                        0.21956039  FALSE
#> 2:                        0.47600041  FALSE
#> 3:                        0.05548664  FALSE
#> 4:                        0.09781648  FALSE
#> 5:                        0.00000000   TRUE
```

And compare against the classical solution without iteration:

``` r
round(dgsi$non_iterated$realised_response, 3)
#>     GY    PHT     AD    ASI    EPP    GLS 
#>  1.468  0.425  0.536 -0.173  0.636  0.495
```

The anthesis-silking interval moves in the unfavourable direction in
this candidate set. Joukhadar et al. (2024) introduced their iterative
search to reduce this type of disagreement between a target and a
selected-candidate differential. This result concerns the current
candidates. It does not alter the exact proportional-response property
of the closed-form index under its covariance model.

### Stage 12 — Score genomic estimated breeding values

[`run_qgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_qgsi.md)
implements the quadratic genomic selection index of Cerón-Rojas et
al. (2026):

\\I\_{qg,i} = \mathbf{w}^{\mathsf T}\widehat{\mathbf{g}}\_i +
\widehat{\mathbf{g}}\_i^{\mathsf T}\mathbf{W}\widehat{\mathbf{g}}\_i\\

where \\\widehat{\mathbf{g}}\_i\\ is the candidate’s vector of genomic
estimated breeding values, \\\mathbf{w}\\ the linear economic weights,
and \\\mathbf{W}\\ a symmetric matrix of squared and cross-product
weights. Negative diagonal curvature in \\\mathbf{W}\\ favours
intermediate values. Positive curvature favours extremes.

``` r
quadratic_weights <- diag(c(
  GY = 0.05, PHT = -0.03, AD = -0.03, ASI = -0.04, EPP = 0.02, GLS = -0.03
))
dimnames(quadratic_weights) <- list(traits, traits)

qgsi <- run_qgsi(
  init_data = dgr_gebv["GenoID"], gebv_data = dgr_gebv,
  trait_cols = traits,
  linear_weights = economic_weights, W = quadratic_weights,
  lower_is_better = lower_is_better,
  scale_traits = TRUE, n_select = 20L
)
qgsi$expected_gain_per_trait
#>     Trait Expected_Genetic_Gain Expected_Genetic_Gain_LinearSD
#>    <char>                 <num>                          <num>
#> 1:     GY            1.38182165                     1.38527228
#> 2:    PHT            0.04857598                     0.04869728
#> 3:     AD            0.23800522                     0.23859955
#> 4:    ASI            1.21214644                     1.21517336
#> 5:    EPP            0.97186893                     0.97429584
#> 6:    GLS            0.93598981                     0.93832713
```

Compare QGSI with a linear genomic selection index by fitting the same
call with `W` set to a zero matrix. Keep the GEBVs, linear weights,
genomic covariance, transformation, and selected count unchanged.
Compare per-trait response, rank correlation, selected overlap, and
validated nonlinear utility. The detailed code is given in
[Multiple-trait
selection](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-index-families.html#compare-a-linear-and-quadratic-genomic-index).

**Standardisation matters more here than anywhere else in the package**,
because economic weights multiply the breeding values directly. Left
unstandardised on these data, plant height carries half the linear index
despite holding nearly the smallest weight, and the expected gain in
grain yield comes out *negative* while the index selects for more of it.
[`run_qgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_qgsi.md)
warns when the trait scales make this likely.

------------------------------------------------------------------------

## Phase 4: Looking ahead

A single-cycle calculation cannot distinguish desired-gain directions
beyond what the feasibility ellipsoid already tells you. Over several
cycles the ranking can change, because truncation selection erodes the
genetic variance that the response depends on. Phase 4 represents that.

The simulation is calibrated to your own germplasm: founders are built
from real phased marker data, never from a simulated genome.

### Stage 13 — Build founders from phased marker data

``` r
founders <- founder_haplotypes(dgr_hap1, dgr_hap2, dgr_map)
founders
#> <desiredgainr_founders>
#>   200 individuals, 450 variants, 3 chromosomes, diploid
#>   Map: converted from physical position at 1e-08 Morgan per base pair 
#>   Missing calls before resolution: 0.000% (policy: error)
```

`hap1` and `hap2` are marker-by-individual matrices coded 0 or 1, giving
the allele on each chromosome copy. This is what a phase-preserving
variant call format reader returns.

**Dosage coded 0, 1, 2 is not accepted**, because it records how many
alternative alleles an individual carries but not which copy carries
them, and that co-occurrence is the linkage disequilibrium the
simulation exists to represent.

For highly inbred diploid material, dosage 0 and 2 are unambiguous, so
phase can be derived without external phasing. Measure the
heterozygosity first:

``` r
dosage <- dgr_hap1 + dgr_hap2
diagnostics <- dosage_diagnostics(dosage)
c(
  heterozygosity = diagnostics$overall_heterozygosity,
  missing = diagnostics$overall_missing
)
#> heterozygosity        missing 
#>      0.5046444      0.0000000
```

``` r
converted <- haplotypes_from_inbred_dosage(
  dosage,
  heterozygous_policy = "drop_variant"
)
```

No threshold is imposed on residual heterozygosity, because no
universally appropriate level exists. It is measured and reported
instead. Heterozygous calls are never resolved silently.

### Stage 14 — Calibrate the simulation

Genome structure comes from your markers. Trait architecture is
calibrated to the covariance matrix you already estimated, so the
simulation reproduces both. Quantitative trait loci (QTL) are sampled on
that structure for the simulation model.

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

For a **clonal** programme add `dominance_degree`, and supply the
*genotypic* rather than the additive covariance as `G`, because the
clone inherits dominance intact.

### Stage 15 — Simulate one objective across cycles

``` r
simulation <- simulate_selection_cycles(
  setup,
  desired_gains = desired_gains,
  n_cycles = 5L,
  mating_system = "outcross",
  n_parents = 20L, n_crosses = 50L, n_progeny_per_cross = 10L,
  lower_is_better = lower_is_better,
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

Output above is illustrative of the object’s shape. Your numbers depend
on your founders and parameters.

Three mating systems are supported, and they differ in more than naming:

| System | Crops | What differs |
|----|----|----|
| `"self"` | Wheat, rice, bean | Advanced by selfing or doubled haploidy. Recombination releases variance slowly. |
| `"outcross"` | Maize, sorghum, millet | Random mating occurs each cycle. Half the selection-induced disequilibrium decays per generation. |
| `"clonal"` | Cassava, sweetpotato, banana | Selection on *total* genetic value, because dominance is inherited intact |

The per-cycle table records the genetic mean and variance for each
trait, the mean relationship among selected parents, and the implied
effective population size, so that gain and the loss of diversity can be
read together.

This is a scenario comparison, not OCS or optimum cross selection. The
simulation asks whether one desired-gain objective remains attractive
under a declared reproduction rule. It does not claim that its recycled
parents or crosses are the operational optimum. Pass the index merit
score to HapBlockR when the real parent, contribution and mating
decisions are made.

`reestimate_index = TRUE`, the default, rebuilds the index from each
cycle’s own data, which is what a programme actually does and which
propagates estimation error across cycles. Setting it `FALSE` isolates
the effect of the desired-gain direction alone. A large divergence
between the two means the recommendation is sensitive to estimation
error.

### Stage 16 — Search for the best objective over cycles

Because only the *direction* of \\\mathbf{d}\\ matters, the search space
is the unit sphere with \\p-1\\ free parameters for \\p\\ traits — small
enough to cover densely.

``` r
optimisation <- optimize_desired_gains(
  setup,
  n_cycles = 5L,
  mode = "pareto",
  budget = 60L,
  n_replicates = 3L,
  mating_system = "outcross",
  n_parents = 20L, n_crosses = 50L, n_progeny_per_cross = 10L,
  lower_is_better = lower_is_better,
  checkpoint = "search.rds",
  seed = 42L
)
optimisation
```

    #> <desiredgainr_optimisation>
    #>   Mode: pareto   Cycles: 5   Evaluations: 60 (3 replicates each)
    #>   Objectives: GY, PHT, AD, ASI, EPP, GLS, diversity
    #>   Non-dominated directions: 9 of 60
    #>   No single direction is recommended. Choose a point on the frontier.

**The simulation is the objective function and is never replaced by a
cheaper approximation.** A Gaussian process decides only where to spend
the next simulation. It never filters or substitutes for one. Because
the acquisition function keeps exploring, no region can be permanently
excluded.

Four ranking modes are available:

| Mode | Objective | Needs |
|----|----|----|
| `"pareto"` | Non-dominated outcomes | nothing beyond the setup |
| `"economic"` | Maximise \\\mathbf{w}^\mathsf{T}\mathbf{R}\\ | `economic_weights` |
| `"target"` | Minimise distance to absolute targets | `target_gains` |
| `"constrained"` | Maximise one trait subject to floors | `focal_trait`, `gain_floors` |

`"pareto"` is the default because it is the only mode that does not
require the weights breeders find hardest to state. Choosing a point on
a frontier is an easier judgement, and it is weight elicitation by
revealed preference.

A realistic budget takes hours, so supply `checkpoint` to make the
search resumable.

------------------------------------------------------------------------

### Sign-off checklist

Complete this before promoting any output to a breeding recommendation.

| Check | Where it came from |
|----|----|
| Covariance matrices positive definite, conditioning inspected | Stage 3 |
| Heritabilities consistent with \\\mathbf{G}\\ and \\\mathbf{P}\\ | Stage 4 |
| Trait directions declared, not signed by hand | Stage 2 |
| Objective feasible at the planned intensity, or knowingly not | Stage 7 |
| Index family chosen and justified against a comparator | Stages 9-10 |
| Effective weights inspected. No unintended trait dominance | Stage 10 |
| Replicate diagnostics agree | Stage 11 |
| Traits standardised, or the warning read and accepted | Stages 3, 12 |
| Decision robust to plausible weight perturbation | *Objective vignette* |
| Final parents, contributions and crosses decided in HapBlockR, not inferred from `n_select` | [Working with other breeding software](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-interoperation.md) |
| Random seed, package version and session recorded | below |

``` r
packageVersion("DesiredGainR")
#> [1] '0.5.0'
sessionInfo()
#> R version 4.5.0 (2025-04-11 ucrt)
#> Platform: x86_64-w64-mingw32/x64
#> Running under: Windows 11 x64 (build 26200)
#> 
#> Matrix products: default
#>   LAPACK version 3.12.1
#> 
#> locale:
#> [1] LC_COLLATE=English_United States.utf8 
#> [2] LC_CTYPE=English_United States.utf8   
#> [3] LC_MONETARY=English_United States.utf8
#> [4] LC_NUMERIC=C                          
#> [5] LC_TIME=English_United States.utf8    
#> 
#> time zone: America/Bogota
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] DesiredGainR_0.5.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] digest_0.6.39     desc_1.4.3        R6_2.6.1          fastmap_1.2.0    
#>  [5] xfun_0.57         cachem_1.1.0      knitr_1.51        htmltools_0.5.9  
#>  [9] rmarkdown_2.31    lifecycle_1.0.5   cli_3.6.6         sass_0.4.10      
#> [13] pkgdown_2.2.0     data.table_1.18.4 textshaping_1.0.5 jquerylib_0.1.4  
#> [17] systemfonts_1.3.2 compiler_4.5.0    rstudioapi_0.18.0 tools_4.5.0      
#> [21] ragg_1.5.2        bslib_0.11.0      evaluate_1.0.5    yaml_2.3.12      
#> [25] otel_0.2.0        jsonlite_2.0.0    htmlwidgets_1.6.4 rlang_1.2.0      
#> [29] fs_2.1.0
```

------------------------------------------------------------------------

### Where each stage is covered in depth

| Stage | Vignette |
|----|----|
| 1-4 | [Obtaining G and P](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-covariance.md) |
| 5-8 | [Defining a breeding objective](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-objective.md) |
| 9-10 | [Multiple-trait selection](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-index-families.md) |
| 11 | [Iterative optimisation of the desired-gain index](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-dgsi.md) |
| 12 | [Quadratic genomic selection](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-qgsi.md) |
| 13-16 | [Comparing objectives over several cycles](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-simulation.md) |
| Downstream decision | [Working with other breeding software](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-interoperation.md) |

------------------------------------------------------------------------

### References

- Cerón-Rojas JJ, Montesinos-López OA, Montesinos-López A, et
  al. (2026). Nonlinear genomic selection index accelerates multi-trait
  crop improvement. *Nature Communications* **17**:1991.
- Covarrubias-Pazaran G (2021). *Practical implementation of selection
  indices.* CGIAR Excellence in Breeding.
- Guimarães PHR et al. (2021). Index selection can improve the selection
  efficiency in a rice recurrent selection population. *Euphytica*
  **217**:95.
- Joukhadar R et al. (2024). Optimising desired gain indices to maximise
  selection response. *Frontiers in Plant Science* **15**:1337388.
- Rahimi M, Debnath S (2023). Estimating optimum and base selection
  indices. *Scientific Reports* **13**:18977.
