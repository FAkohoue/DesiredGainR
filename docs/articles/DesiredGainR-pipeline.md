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
This document is the connective walkthrough; those are the references.

------------------------------------------------------------------------

## Phase 1: Evidence

### Stage 1 — What you need before you start

DesiredGainR begins *after* the genetic evaluation. It does not analyse
field trials, does not fit a mixed model, and does not estimate breeding
values. You must arrive with four things.

| Input | What it is | Where it comes from |
|----|----|----|
| Trait values | One adjusted mean, BLUP or GEBV per candidate per trait | Your multi-trait genetic evaluation |
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

Every index in this package is computed by inverting \\\mathbf{P}\\,
\\\mathbf{G}\\, or both. An ill-conditioned matrix produces coefficients
that are numerically meaningless while looking entirely plausible, so
this check comes before anything is built on top of them.

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
#> [1] "working surrogate for genetic covariance; covariance of across-environment adjusted means"
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

There are two ways to state a breeding objective, and they are
equivalent.

**Economic weights** \\\mathbf{w}\\ say what a unit of each trait is
worth. **Desired gains** \\\mathbf{d}\\ say how much of each trait you
want. The Smith-Hazel index uses the first, \\\mathbf{b} =
\mathbf{P}^{-1}\mathbf{G}\mathbf{w}\\; the Pesek-Baker desired-gain
index uses the second, \\\mathbf{b} =
\mathbf{P}^{-1}\mathbf{G}(\mathbf{G}\mathbf{P}^{-1}\mathbf{G})^{-1}\mathbf{d}\\.

Most breeders find desired gains easier, so start there. State them as
improvements, in genetic standard deviations, and let `lower_is_better`
handle direction.

``` r
desired_gains <- c(GY = 1.0, PHT = 0.4, AD = 0.6, ASI = 0.5,
                   EPP = 0.4, GLS = 0.6)
```

This says: one standard deviation of extra yield, four tenths of a
standard deviation shorter plants, six tenths earlier flowering, and so
on.

### Stage 6 — Translate to economic weights

Setting the two index formulae equal gives an exact translation in both
directions:

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
#> [1] "Implied by the supplied desired gains through w = G^-1 P G^-1 d; not an independently estimated economic value. Expressed in the favourable-direction space, so a positive weight always means the trait matters."
```

**Two of these are negative, for traits you asked to improve.** That is
correct. Once oriented, grain yield correlates \\+0.55\\ with the
anthesis-silking interval and \\+0.45\\ with ears per plant, so
selecting hard for yield would carry both *past* the smaller gains you
requested. To deliver the ratio you asked for, the index must hold them
back. A negative weight never means the trait should get worse.

**The magnitudes are not comparable.** A weight carries inverse trait
units, so a trait on a small scale gets a big number for the same
emphasis. Always rescale before comparing:

``` r
round(sort(implied * genetic_sd[names(implied)], decreasing = TRUE), 2)
#>    GY    AD   GLS   PHT   EPP   ASI 
#> 10.89  5.81  0.58  0.36 -1.51 -7.93
```

Ears per plant had the largest raw number and is one of the smaller
effects.

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

For any linear index the response is \\\mathbf{R} =
i\\\mathbf{G}\mathbf{b}/\sqrt{\mathbf{b}^\mathsf{T}\mathbf{P}\mathbf{b}}\\,
so the attainable responses lie on an ellipsoid and a stated target
requires exactly

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

The request looked unremarkable and requires selecting fewer than one
candidate from two hundred. About fifty-five per cent of it is available
at the intensity actually planned:

``` r
round(feasibility$attainable_response_input_units, 3)
#>    GY   PHT    AD   ASI   EPP   GLS 
#> 0.547 0.219 0.328 0.274 0.219 0.328
```

**Two consequences worth understanding now.** First, the classical
desired-gain index honours only the *direction* of \\\mathbf{d}\\;
multiplying every element by a constant changes nothing, because the
magnitude is fixed by selection intensity. Second, an antagonistic
correlation structure can make a modest-looking target unreachable, and
no amount of optimisation will fix that.

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

This is how an index was introduced at both CIMMYT and IRRI. Treat it as
a starting point: it reproduces past decisions including any bias in
them, and must be adjusted deliberately against the product profile
before it becomes a forward-looking objective.

------------------------------------------------------------------------

## Phase 3: Index

### Stage 9 — Build candidate indices

[`selection_index()`](https://FAkohoue.github.io/DesiredGainR/reference/selection_index.md)
provides every classical family through one interface.

``` r
economic_weights <- c(GY = 1.0, PHT = 0.2, AD = 0.5,
                      ASI = 0.4, EPP = 0.3, GLS = 0.5)

smith_hazel <- selection_index(
  dgr_candidates, traits, method = "smith_hazel",
  G = dgr_G, P = dgr_P, economic_weights = economic_weights,
  lower_is_better = lower_is_better, n_select = 20L, main_trait = "GY"
)
rank_sum <- selection_index(
  dgr_candidates, traits, method = "mulamba_mock",
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

| Method | Needs weights | Needs \\\mathbf{G}\\, \\\mathbf{P}\\ | Use when |
|----|----|----|----|
| `smith_hazel` | yes | yes | Weights are defensible |
| `base` | yes | no | A quick comparator; often nearly as good |
| `pesek_baker` | desired gains | \\\mathbf{G}\\ only | Gains easier than weights |
| `yamada` | desired gains | yes | Same index, better conditioned route |
| `mulamba_mock` | no | no | Weights unavailable or untrusted |
| `independent_culling` | thresholds | no | Non-compensatory requirements |
| `tandem` | trait order | no | Comparator for sequential selection |

### Stage 10 — Evaluate and compare

The four criteria reported are those of Rahimi and Debnath (2023), which
the established selection-index software also reports.

| Criterion    | Meaning                                                   |
|--------------|-----------------------------------------------------------|
| \\R\_{HI}\\  | Correlation between index and net merit                   |
| \\\Delta H\\ | Expected gain in aggregate merit                          |
| RE           | Efficiency relative to direct selection on the main trait |
| \\CV_I\\     | Coefficient of variation of the index                     |

**Relative efficiency below 1 is not a failure.** It is the intended
trade of response in the main trait for response elsewhere. Every value
reported by Rahimi and Debnath was below 1.

Compare families by rank correlation *and* by overlap of the selected
sets, because the two must be consistent with each other:

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

A high correlation with a nearly empty intersection indicates a fault,
not a finding. That exact check caught a row-alignment bug during
development.

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
implements the iterative method of Joukhadar et al. (2024). It samples
desired-gain vectors, builds the index from each, and keeps the one
whose *realised* response in the selected set comes closest to your
target.

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
round(dgsi$realised_response, 3)
#>    GY   PHT    AD   ASI   EPP   GLS 
#> 1.572 0.303 0.258 0.528 0.449 0.532
```

Read the deviation from target against Stage 7:

``` r
round(dgsi$realised_response - desired_gains, 3)
#>     GY    PHT     AD    ASI    EPP    GLS 
#>  0.572 -0.097 -0.342  0.028  0.049 -0.068
```

Feasibility constrained the requested *proportion*. The optimiser is not
bound to that proportion, so it distributes the shortfall unevenly,
overshooting grain yield while undershooting anthesis date. Whether that
trade is acceptable is a breeding decision, not a numerical one.

Always inspect the replicate diagnostics. The search is stochastic, so
agreement across replicates is the evidence that the answer is stable:

``` r
dgsi$replicate_diagnostics
#>    Replicate Objective Iteration_of_best Selected Plateau
#>        <int>     <num>             <int>    <int>  <lgcl>
#> 1:         1 0.7427606               166       20   FALSE
#> 2:         2 1.3406629               190       20   FALSE
#> 3:         3 1.0778294                35       20    TRUE
#> 4:         4 0.8216413               120       20   FALSE
#> 5:         5 2.0393576               142       20   FALSE
#>    Final_window_relative_improvement Chosen
#>                                <num> <lgcl>
#> 1:                        0.71931914   TRUE
#> 2:                        0.30952774  FALSE
#> 3:                        0.00000000  FALSE
#> 4:                        0.51939104  FALSE
#> 5:                        0.01409831  FALSE
```

And compare against the classical solution without iteration:

``` r
round(dgsi$non_iterated$realised_response, 3)
#>     GY    PHT     AD    ASI    EPP    GLS 
#>  1.526  0.274  0.441 -0.256  0.562  0.228
```

The anthesis-silking interval moves the *wrong way* there. Joukhadar et
al. reported precisely this failure of the unoptimised method.

### Stage 12 — Score genomic estimated breeding values

[`run_qgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_qgsi.md)
implements the quadratic genomic selection index of Cerón-Rojas et
al. (2026):

\\I\_{qg,i} = \mathbf{w}^\mathsf{T}\hat{\boldsymbol{\gamma}}\_i +
\hat{\boldsymbol{\gamma}}\_i^\mathsf{T}\mathbf{W}\hat{\boldsymbol{\gamma}}\_i\\

where \\\hat{\boldsymbol{\gamma}}\_i\\ is the candidate’s vector of
genomic estimated breeding values, \\\mathbf{w}\\ the linear economic
weights, and \\\mathbf{W}\\ a symmetric matrix of squared and
cross-product weights. Negative diagonal curvature in \\\mathbf{W}\\
favours intermediate values; positive curvature favours extremes.

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
  dosage, heterozygous_policy = "drop_variant"
)
```

No threshold is imposed on residual heterozygosity, because no
universally appropriate level exists; it is measured and reported
instead, and heterozygous calls are never resolved silently.

### Stage 14 — Calibrate the simulation

Genome structure comes from your markers. Trait architecture is
calibrated to the covariance matrix you already estimated, so the
simulation reproduces both.

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

Output above is illustrative of the object’s shape; your numbers depend
on your founders and parameters.

Three mating systems are supported, and they differ in more than naming:

| System | Crops | What differs |
|----|----|----|
| `"self"` | Wheat, rice, bean | Advanced by selfing or doubled haploidy; recombination releases variance slowly |
| `"outcross"` | Maize, sorghum, millet | Random mating each cycle; half the selection-induced disequilibrium decays per generation |
| `"clonal"` | Cassava, sweetpotato, banana | Selection on *total* genetic value, because dominance is inherited intact |

The per-cycle table records the genetic mean and variance for each
trait, the mean relationship among selected parents, and the implied
effective population size, so that gain and the loss of diversity can be
read together.

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
    #>   No single direction is recommended; choose a point on the frontier.

**The simulation is the objective function and is never replaced by a
cheaper approximation.** A Gaussian process decides only where to spend
the next simulation; it never filters or substitutes for one, and
because the acquisition function keeps exploring, no region can be
permanently excluded.

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
| Effective weights inspected; no unintended trait dominance | Stage 10 |
| Replicate diagnostics agree | Stage 11 |
| Traits standardised, or the warning read and accepted | Stages 3, 12 |
| Decision robust to plausible weight perturbation | *Objective vignette* |
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

| Stage | Vignette                                                   |
|-------|------------------------------------------------------------|
| 1-4   | *Obtaining the genetic and phenotypic covariance matrices* |
| 5-8   | *Defining a breeding objective*                            |
| 9-10  | *Choosing an index*                                        |
| 11    | *Optimising desired gains*                                 |
| 12    | *Quadratic genomic selection*                              |
| 13-16 | *Comparing objectives over several cycles*                 |

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
