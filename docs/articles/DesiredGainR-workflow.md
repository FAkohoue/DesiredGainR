# The complete workflow

## 1. Purpose

Every other vignette is a deep dive into one layer. This one runs the
whole package, start to finish, on the shipped example programme, so
that the relationships between the stages are visible in one place
rather than left to be inferred across eight documents.

The order below is the order in which the decisions actually arise. It
is not the order of the function reference, and it deliberately places
the objective before the index, because an index built on an objective
that was never tested is a well-computed answer to an unexamined
question.

Each stage gives the call and the decision it feeds. The vignette named
at the end of each stage carries the full argument reference and the
statistical detail.

| Phase      | Stages | Question                                         |
|------------|--------|--------------------------------------------------|
| Evidence   | 1-2    | What do we have, and can it be trusted?          |
| Objective  | 3-5    | What should we select for, and is it possible?   |
| Index      | 6-8    | Which index, and what does it actually do?       |
| Confidence | 9-10   | Would a different objective change the decision? |

------------------------------------------------------------------------

## 2. Stage 1 — Inspect the covariance matrices

Index coefficients are obtained by inverting `P` and `G`. An
ill-conditioned matrix yields coefficients that are numerically
meaningless while remaining superficially plausible, so conditioning is
checked before anything is built on top of it.

``` r
matrix_diagnostics(dgr_G, "G")$condition_number
#> [1] 18569.92
matrix_diagnostics(dgr_P, "P")$condition_number
#> [1] 10611.84
```

Those numbers look alarming until they are decomposed. Almost all of the
ill-conditioning comes from the trait scales rather than from the
correlations:

``` r
matrix_diagnostics(stats::cov2cor(dgr_G))$condition_number
#> [1] 8.75443
```

Standardising removes three orders of magnitude. This is the practical
form of the observation by Crosbie et al. (1980) that equal weights on
unstandardised traits concentrate selection on whichever trait carries
the largest variance, and it is why the package standardises by default.

See [Obtaining G and
P](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-covariance.md).

------------------------------------------------------------------------

## 3. Stage 2 — Confirm what the matrices represent

`G` should be a genetic covariance matrix from a fitted multi-trait
model. Where none exists,
[`estimate_genetic_covariance()`](https://FAkohoue.github.io/DesiredGainR/reference/estimate_genetic_covariance.md)
provides labelled approximations, including the covariance of
across-environment adjusted means that the CGIAR guideline recommends
for operational use. Each method records its estimand rather than
relabelling a covariance of predictions as a genetic covariance.

``` r
G_working <- estimate_genetic_covariance(
  values = dgr_candidates, trait_cols = traits,
  method = "adjusted_means_surrogate"
)
G_working$estimand
#> [1] "working surrogate for genetic covariance; covariance of
#>      across-environment adjusted means"
```

The example data ship a known `G`, so the rest of this vignette uses it
directly.

------------------------------------------------------------------------

## 4. Stage 3 — State the objective

Desired gains are stated as improvements, and `lower_is_better` declares
which traits improve by falling. The covariance matrices are oriented
internally, so nothing needs signing by hand.

``` r
desired_gains <- c(
  GY = 1.0, PHT = 0.4, AD = 0.6, ASI = 0.5,
  EPP = 0.4, GLS = 0.6
)

implied <- implied_economic_weights(
  desired_gains, dgr_G, dgr_P,
  lower_is_better = lower_is_better, gain_units = "genetic_sd"
)
round(implied, 2)
#>     GY    PHT     AD    ASI    EPP    GLS 
#>  14.52   0.03   2.32 -13.21 -15.15   0.68 
#> attr(,"provenance")
#> [1] "Implied by the supplied desired gains through w = G^-1 P G^-1 d; not an independently estimated economic value. Expressed in the favourable-direction space, so a positive weight always means the trait matters."
```

Two of those weights are negative for traits the objective asks to
improve, which is correct rather than an error. In addition, the target
values are already genetic-standard-deviation gains, whereas the
returned weights are per original trait unit because `dgr_G` and `dgr_P`
are in original units. The following multiplication expresses each
weight per genetic standard deviation; it does not standardise the
target again:

``` r
genetic_sd <- stats::setNames(dgr_traits$genetic_sd, dgr_traits$trait)
round(sort(implied * genetic_sd[names(implied)], decreasing = TRUE), 2)
#>    GY    AD   GLS   PHT   EPP   ASI 
#> 10.89  5.81  0.58  0.36 -1.51 -7.93
```

If the covariance matrices had already been transformed to
genetic-standard- deviation units, this multiplication would be
incorrect because the weights would already be on that scale.

See [Defining a breeding
objective](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-objective.md).

------------------------------------------------------------------------

## 5. Stage 4 — Test whether the objective is attainable

This stage exists because an optimisation can return a numerical answer
even when the stated one-cycle target is outside the achievable-response
ellipsoid.
[`gain_feasibility()`](https://FAkohoue.github.io/DesiredGainR/reference/gain_feasibility.md)
tests the exact desired-gain ray before any shortfall is interpreted.

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

The required intensity is 3.2082, equivalent under normal truncation to
the best 0.1773%. Among 200 candidates, that proportion represents 0.355
candidate. It cannot be rounded to one because selecting one retains
0.5% and gives a lower intensity. Therefore, the exact vector is not
attainable in this finite candidate set.

The planned selection of 20 candidates has intensity 1.7550. Its ratio
to the required intensity is 0.547, so 54.7% of every target component
is available while the requested trait proportions are preserved.

``` r
round(feasibility$attainable_response_input_units, 3)
#>    GY   PHT    AD   ASI   EPP   GLS 
#> 0.547 0.219 0.328 0.274 0.219 0.328
```

These values are in the same genetic-standard-deviation units as the
input. They describe the point on the requested ray, not independent
upper limits for the six traits. If the entries are minimum floors
rather than a fixed ratio, use population-driven or interval
optimisation and inspect worst-trait and joint attainment. The
conclusion remains conditional on the linear model, normal truncation
selection and the estimated covariance matrices.

Carry this number forward to Stage 8, where the optimiser will be seen
to distribute the shortfall unevenly rather than uniformly.

------------------------------------------------------------------------

## 6. Stage 5 — Recover what the programme already selects for

Where an index has never been used, the differentials a programme has
achieved already encode its objective.

``` r
recovered <- retrospective_weights(
  selected_values = dgr_candidates[dgr_history$selected, traits],
  population_values = dgr_candidates[, traits],
  trait_cols = traits
)
round(recovered$selection_differential_sd, 2)
#>    GY   PHT    AD   ASI   EPP   GLS 
#>  1.11 -0.27 -0.42 -0.83  0.52 -0.42
```

The recovered rule reproduces past decisions, including any bias in
them, so it is a starting point rather than a forward-looking objective.
The objective vignette verifies the recovery against the weights that
generated the historical decision.

------------------------------------------------------------------------

## 7. Stage 6 — Build and compare index families

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

Relative efficiency below one is the intended trade of response in the
main trait for response elsewhere in the objective, not a defect. Every
value reported by Rahimi and Debnath (2023) was below one.

The Mulamba-Mock rank-sum index needs neither weights nor covariance
matrices, which is why Guimarães et al. (2021) found it competitive.
Compare the two by rank correlation *and* by set overlap, because the
two must be consistent:

``` r
paired <- merge(
  smith_hazel$ranking[, c("id", "score")],
  rank_sum$ranking[, c("id", "score")],
  by = "id", suffixes = c("_sh", "_mm")
)
round(stats::cor(paired$score_sh, paired$score_mm, method = "spearman"), 3)
#> [1] 0.938
length(intersect(smith_hazel$selected$id, rank_sum$selected$id))
#> [1] 13
```

A high correlation paired with a nearly empty intersection would
indicate a fault rather than a finding. That check found a genuine
row-alignment bug during development.

See [Choosing an
index](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-index-families.md).

------------------------------------------------------------------------

## 8. Stage 7 — Check which traits the index acts on

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

------------------------------------------------------------------------

## 9. Stage 8 — Optimise the desired-gain vector

[`run_dgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_dgsi.md)
implements the iterative method of Joukhadar et al. (2024), searching
for a desired-gain vector whose realised response in the selected set
approaches the stated target.

``` r
dgsi <- run_dgsi(
  init_data = dgr_candidates[, c("GenoID", "Family")],
  cand_data = dgr_candidates,
  trait_cols = traits,
  dg = desired_gains,
  G = dgr_G, P = dgr_P,
  lower_is_better = lower_is_better,
  scale_traits = TRUE,
  n_select = 20L, n_iter = 200L, n_rep = 5L, seed = 42L
)
round(dgsi$realised_response - desired_gains, 3)
#>     GY    PHT     AD    ASI    EPP    GLS 
#> -1.455 -0.117 -0.738  0.566  0.007 -0.038
```

Read that against Stage 4. Feasibility constrains the requested
*proportion*; the optimiser is not bound to it, and distributes the
shortfall unevenly, overshooting grain yield while undershooting
anthesis date. Whether that trade is acceptable is a breeding decision
rather than a numerical one.

The classical solution without iteration is returned as a comparator:

``` r
round(dgsi$non_iterated$realised_response, 3)
#>     GY    PHT     AD    ASI    EPP    GLS 
#>  1.526  0.274  0.441 -0.256  0.562  0.228
```

The anthesis-silking interval moves in the wrong direction there.
Joukhadar et al. reported exactly this failure of the unoptimised
method, and it reproduces here on independent data.

See [Optimising desired
gains](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-dgsi.md).

------------------------------------------------------------------------

## 10. Stage 9 — Score genomic estimated breeding values

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

Standardisation matters more here than anywhere else in the package. On
the original scale, plant height carries half the linear index despite
holding nearly the smallest weight, and the expected gain in grain yield
turns *negative* while the index selects for more of it.
[`run_qgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_qgsi.md)
warns when the trait scales make that likely.

See [Quadratic genomic
selection](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-qgsi.md).

------------------------------------------------------------------------

## 11. Stage 10 — Ask how much the objective matters

``` r
weight_sensitivity(
  economic_weights, dgr_candidates, dgr_G, dgr_P,
  n_select = 20L, trait_cols = traits, n_draws = 100L
)
#> <desiredgainr_sensitivity>
#>   100 draws, log-scale perturbation SD 0.25, selecting 20
#>   Median selected-set agreement: 0.905
#>   Decisions reproducing the original set (agreement >= 0.90): 56.0%
#>   Median rank correlation with the stated objective: 0.994
#>   Weights the decision is most sensitive to:
#>   PHT    GY   ASI 
#> 0.516 0.475 0.046 
#>   (Sensitivity of the decision, not share of the index; see effective_weights().)
```

Note that the weight the decision is most sensitive to is not the trait
carrying most of the index. Influence on the decision and contribution
to the index value are different quantities, and reading them together
is the point.

------------------------------------------------------------------------

## 12. What comes after

Two questions lie beyond a single cycle and beyond this package.

**Does this objective still look right in five cycles?** Truncation
selection erodes the genetic variance that response depends on, so a
direction that maximises first-cycle response can exhaust itself. The
simulation layer compares directions over successive cycles, with
founders built from the breeder’s own phased marker data. See [Comparing
objectives over several
cycles](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-simulation.md).

**Which candidates should actually be crossed?** Parent, cross and
mating allocation, coancestry control and operational feasibility are
outside this package’s scope by design. See [Working with other breeding
software](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-interoperation.md).

------------------------------------------------------------------------

## 13. Where each stage is covered in depth

| Stage | Vignette |
|----|----|
| 1-2 | [Obtaining G and P](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-covariance.md) |
| 3-5 | [Defining a breeding objective](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-objective.md) |
| 6-7 | [Choosing an index](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-index-families.md) |
| 8 | [Optimising desired gains](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-dgsi.md) |
| 9 | [Quadratic genomic selection](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-qgsi.md) |
| 10 | [Defining a breeding objective](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-objective.md) |
| Beyond | [Comparing objectives over several cycles](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-simulation.md); [Working with other breeding software](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-interoperation.md) |

``` r
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
