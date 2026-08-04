# Multiple-trait selection: methods, choice, and comparison

## 1. Learning objectives

Plant breeding programmes usually improve several traits at the same
time. Yield, quality, adaptation, maturity, and disease resistance may
all influence the release decision. Progress in one trait can also
reduce progress in another. Therefore, the method must represent both
the breeding objective and the genetic relationships among traits.

This vignette has four objectives. They are to: (i) define
multiple-trait selection, (ii) explain its main methods, (iii) show when
each method is appropriate, and (iv) provide a comprehensive comparison
analysis.

At the end, the breeder should be able to justify the selected method.
The justification should refer to the objective, the available
information, the expected responses, and the stability of the final
ranking.

## 2. What is multiple-trait selection?

Multiple-trait selection is the choice of candidates using two or more
traits. The traits can enter the decision at the same stage or at
different stages. The breeder must first decide how trade-offs are
allowed.

Four broad strategies are recognised by Beavis, Lamkey, Mahama, and Suza
(2023).

| Strategy | Definition | Suitable situation | Main limitation |
|----|----|----|----|
| Multistage selection | Different traits are evaluated at different stages. | Expensive traits are measured after inexpensive screening. | Early rejection can remove candidates that would perform well later. |
| Tandem selection | One trait is improved first. Another trait becomes the focus later. | One urgent trait dominates a specific breeding stage. | Other traits can deteriorate while attention is placed elsewhere. |
| Independent culling | Every candidate must satisfy a limit for every required trait. | Release standards or biological limits are firm. | Excellence elsewhere gives no compensation for a failed limit. |
| Index selection | Traits are combined into one score of overall merit. | Trade-offs are allowed and simultaneous improvement is required. | The score is useful only when its objective and inputs are defensible. |

These strategies answer different questions. Hence, a comparison should
never treat them as interchangeable algorithms. The breeder’s decision
rule comes first. The statistical method follows.

## 3. The quantities used by an index

Let \\\mathbf{x}\\ contain the information recorded for a candidate. Let
\\\mathbf{g}\\ contain the additive genetic values of the objective
traits. The two vectors can contain the same traits. They can also
differ. For example, records from relatives, testcrosses, environments,
or correlated traits can help predict a smaller set of objective traits.

Three covariance matrices describe this problem.

| Symbol | Definition | Breeding meaning |
|----|----|----|
| \\\mathbf{P}=\operatorname{Var}(\mathbf{x})\\ | Covariance among information variables | Describes the variation and redundancy in the evidence used for selection. |
| \\\mathbf{C}=\operatorname{Cov}(\mathbf{x},\mathbf{g})\\ | Covariance between information and objective traits | Describes how strongly each information source predicts genetic value. |
| \\\mathbf{G}=\operatorname{Var}(\mathbf{g})\\ | Genetic covariance among objective traits | Describes available genetic variation and correlated response. |

The aggregate genetic merit is

\\ H=\mathbf{a}^{\mathsf T}\mathbf{g}, \\

where \\\mathbf{a}\\ contains economic weights. A linear selection index
is

\\ I=\mathbf{b}^{\mathsf T}\mathbf{x}, \\

where \\\mathbf{b}\\ contains index coefficients. Under truncation
selection, the expected response vector is

\\ \Delta =i\frac{\mathbf{C}^{\mathsf T}\mathbf{b}}
{\sqrt{\mathbf{b}^{\mathsf T}\mathbf{P}\mathbf{b}}}. \\

The term \\i\\ is the selection intensity. It increases when a smaller
proportion of candidates is selected. When information and objective
traits are identical, \\\mathbf{C}=\mathbf{G}\\. This gives the familiar
expression used throughout classical selection-index theory.

Use
[`selection_information()`](https://FAkohoue.github.io/DesiredGainR/reference/selection_information.md)
and
[`generalized_index()`](https://FAkohoue.github.io/DesiredGainR/reference/generalized_index.md)
when the information and objective traits differ. Read [Using
predictions, prediction errors, and restricted
responses](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-information.md)
for the complete formulation.

## 4. The methods available in DesiredGainR

### 4.1 Economic-value methods

The Smith-Hazel index maximises expected response in aggregate genetic
merit. Its coefficients are

\\ \mathbf{b}=\mathbf{P}^{-1}\mathbf{G}\mathbf{a}. \\

Choose this method when the relative economic values are credible. It is
also the correct reference when every method is compared against one
fixed definition of merit.

The base index uses \\\mathbf{b}=\mathbf{a}\\. It applies the economic
weights directly. It provides a useful benchmark because it avoids
covariance inversion. A small advantage for Smith-Hazel over the base
index suggests that the estimated covariance structure contributes
little to the decision.

### 4.2 Desired-gain methods

The Pesek-Baker and Yamada formulations use a desired response direction
\\\mathbf{d}\\. Their coefficient expressions are

\\ \mathbf{b}=\mathbf{G}^{-1}\mathbf{d} \\

and

\\ \mathbf{b}=\mathbf{P}^{-1}\mathbf{G}
(\mathbf{G}\mathbf{P}^{-1}\mathbf{G})^{-1}\mathbf{d}. \\

They give the same coefficient direction when \\\mathbf{G}\\ is square
and invertible. Therefore, they are two computational routes to the same
linear index. A material difference between their answers indicates
numerical instability or a rank-deficient matrix.

Choose this family when the breeder can state relative desired gains
more confidently than economic weights. The vector fixes response
proportions. The selection intensity fixes the attainable magnitude.
Multiplying every desired gain by the same positive constant leaves the
ranking unchanged.

[`run_dgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_dgsi.md)
implements the iterative desired-gain selection index (DGSI). It
searches for a candidate ranking whose selected differential is close to
the requested vector. The closed-form desired-gain index remains its
mathematical reference. Read [Optimising desired
gains](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-dgsi.md)
before using the iterative method.

### 4.3 Rank and threshold methods

The Mulamba-Mock rank-sum method ranks candidates within traits and
combines the ranks. It requires neither economic weights nor covariance
matrices when equal rank weights are used. It is valuable as a simple
comparator. Its score does not provide a closed-form predicted genetic
response in DesiredGainR.

Independent culling applies firm acceptance limits. Supply each limit in
the original trait units. Use a minimum for traits where larger values
are favourable. Use a maximum for traits listed in `lower_is_better`.

The Elston multiplicative index also applies firm limits. It then ranks
the eligible candidates using their joint margins above those limits. It
is useful when balanced superiority among eligible candidates matters.

Tandem selection focuses on one trait at each stage. DesiredGainR
includes it as a historical and operational comparator.

### 4.4 Restricted methods

A restricted index maximises merit while controlling response in
specified traits. DesiredGainR provides five forms.

| Method | Breeding question |
|----|----|
| Kempthorne-Nordskog | Which index holds specified traits at zero expected change? |
| Restricted Smith-Hazel | Which projection form gives the same zero-response restriction? |
| Tallis | Which index preserves stated response proportions? |
| Mallard | Which index targets stated response amounts? |
| Harville | Which index balances merit against a finite restriction penalty? |

Restrictions usually reduce response in unconstrained net merit. That
loss is the cost of respecting the biological or commercial condition.
Report both the restricted response and the loss relative to the
unrestricted index.

Satoh (2024) gives a geometric interpretation of proportional
restriction. The scalar \\\beta\\ measures progress along the desired
direction. The Mahalanobis residual measures departure from that
direction after accounting for genetic covariance. DesiredGainR reports
both quantities through
[`evaluate_restricted_response()`](https://FAkohoue.github.io/DesiredGainR/reference/evaluate_restricted_response.md).

### 4.5 Non-linear economic value

The quadratic genomic selection index (QGSI) represents linear, squared,
and cross-product economic values. Choose it when biological or economic
merit is genuinely curved. A desired gain is a response target. It is a
different quantity from a quadratic weight. Read [Quadratic genomic
selection](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-qgsi.md)
for the full demonstration.

## 5. How to choose a method

Begin with the breeder’s statement. Then choose the simplest method that
represents it accurately.

| Breeder’s statement | Starting method | Required evidence |
|----|----|----|
| “I can defend relative economic values.” | Smith-Hazel | \\\mathbf{G}\\, \\\mathbf{P}\\, and economic weights |
| “I can state the response direction.” | Pesek-Baker or Yamada | \\\mathbf{G}\\, and \\\mathbf{P}\\ for Yamada |
| “I can state acceptable gain intervals.” | [`suggest_desired_gains()`](https://FAkohoue.github.io/DesiredGainR/reference/suggest_desired_gains.md) | \\\mathbf{G}\\, \\\mathbf{P}\\, intervals, and candidates |
| “Every candidate must satisfy these limits.” | Independent culling | Candidate values and limits |
| “Eligible candidates should have balanced margins.” | Elston | Candidate values and limits |
| “Some traits should remain unchanged.” | Restricted index | \\\mathbf{G}\\, \\\mathbf{P}\\, merit weights, and restrictions |
| “One trait is screened at each stage.” | Tandem or multistage selection | Stage order and stage-specific evidence |
| “I need a simple comparator with few assumptions.” | Base or Mulamba-Mock | Weights for the base index, or candidate ranks |
| “Merit contains curvature or trait interactions.” | QGSI | Genomic estimated breeding values and quadratic economic values |

Use the following decision sequence.

1.  Define the objective traits and their favourable directions.
2.  Identify firm limits that allow no trade-off.
3.  Decide whether the remaining objective is economic, desired-gain, or
    non-linear.
4.  Confirm the covariance matrices and units required by that method.
5.  Fit one simple benchmark and one objective-matched method.
6.  Compare biological response before comparing scalar criteria.
7.  Examine ranking stability and uncertainty.
8.  Use simulation when repeated cycles or non-linear decisions matter.

## 6. Conditions for a fair comparison

A comprehensive comparison requires common conditions. Otherwise, a
numerical difference can reflect the setup rather than the method.

Use the same:

1.  candidates and trait records;
2.  trait directions and units;
3.  genetic and phenotypic covariance estimates;
4.  selected proportion or selection intensity;
5.  aggregate-merit definition when comparing \\R\_{HI}\\ or \\\Delta
    H\\;
6.  training, validation, and test partitions;
7.  simulation scenarios, replicate seeds, and computational budget.

The comparison also needs several outcomes. No single statistic is
sufficient.

| Outcome | Question answered |
|----|----|
| Expected response for every trait | What genetic change is predicted? |
| Target attainment | How much of each requested gain is reached? |
| Feasibility | Can the complete target be reached at the planned intensity? |
| \\R\_{HI}\\ | How closely does the index rank one fixed aggregate merit? |
| \\\Delta H\\ | How much change in that fixed merit is expected? |
| Relative efficiency (RE) | How much response in the declared main trait remains relative to direct selection? |
| Index heritability | How repeatable is the composite score under the fitted covariance model? |
| Mahalanobis alignment | Does the response follow the requested multivariate direction? |
| Selected-set overlap | Do two methods choose the same candidates? |
| Rank correlation | Do two methods order all candidates similarly? |
| Selection probability | Does candidate choice remain stable after prediction uncertainty is propagated? |
| Information-deletion efficiency | Which records contribute useful selection information? |
| Multi-cycle response and diversity | Does the apparent advantage persist across breeding cycles? |

The coefficient of variation of the index, \\CV_I\\, is retained for
reproduction of published analyses. It depends on the score origin and
becomes undefined after centring. Therefore, it should never decide the
preferred method.

## 7. A comprehensive comparison analysis

### 7.1 State the objective

The example seeks gains in grain yield (GY), ears per plant (EPP), and
lower values for plant height (PHT), anthesis date (AD),
anthesis-silking interval (ASI), and grey leaf spot score (GLS). The
target is expressed in genetic standard deviations.

``` r
desired_gains_sd <- c(
  GY = 1.0,
  PHT = 0.4,
  AD = 0.6,
  ASI = 0.5,
  EPP = 0.4,
  GLS = 0.6
)
```

First test feasibility. This separates an unrealistic target from a poor
method.

``` r
feasibility <- gain_feasibility(
  desired_gains_sd,
  dgr_G,
  dgr_P,
  n_candidates = nrow(dgr_candidates),
  n_select = 20L,
  lower_is_better = lower_is_better,
  gain_units = "genetic_sd"
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

### 7.2 Put every method on one explicit scale

The following transformation expresses candidates and covariance
matrices in genetic-standard-deviation units.
[`selection_index()`](https://FAkohoue.github.io/DesiredGainR/reference/selection_index.md)
then receives `scale_traits = FALSE`. Hence, the desired-gain vector and
every reported response use the same units.

``` r
genetic_sd <- sqrt(diag(dgr_G))[traits]
inverse_sd <- diag(1 / genetic_sd, nrow = length(traits))
dimnames(inverse_sd) <- list(traits, traits)

G_sd <- inverse_sd %*% dgr_G %*% inverse_sd
P_sd <- inverse_sd %*% dgr_P %*% inverse_sd

candidates_sd <- dgr_candidates
centred <- sweep(
  as.matrix(dgr_candidates[, traits]),
  2L,
  colMeans(dgr_candidates[, traits]),
  "-"
)
candidates_sd[, traits] <- sweep(centred, 2L, genetic_sd, "/")
```

The economic weights below are implied by the desired-gain direction.
This creates one common merit definition for the mathematical
comparison.

``` r
common_weights <- implied_economic_weights(
  desired_gains_sd,
  G_sd,
  P_sd,
  lower_is_better = lower_is_better,
  gain_units = "trait"
)
round(common_weights, 3)
#>     GY    PHT     AD    ASI    EPP    GLS 
#> 10.887  0.362  5.808 -7.928 -1.515  0.576 
#> attr(,"provenance")
#> [1] "Implied by the supplied desired gains through w = G^-1 P G^-1 d; not an independently estimated economic value. Expressed in the favourable-direction space, so a positive weight always means the trait matters."
```

Negative implied weights can occur. They are mathematically valid under
correlated response. The desired response remains the biological
statement to interpret.

### 7.3 Fit objective-matched methods and benchmarks

``` r
common_arguments <- list(
  values = candidates_sd,
  trait_cols = traits,
  id_col = "GenoID",
  G = G_sd,
  P = P_sd,
  lower_is_better = lower_is_better,
  center_traits = FALSE,
  scale_traits = FALSE,
  n_select = 20L,
  main_trait = "GY"
)

smith_hazel <- do.call(selection_index, c(
  common_arguments,
  list(method = "smith_hazel", economic_weights = common_weights)
))

pesek_baker <- do.call(selection_index, c(
  common_arguments,
  list(
    method = "pesek_baker",
    desired_gains = desired_gains_sd,
    aggregate_weights = common_weights
  )
))

base_index <- do.call(selection_index, c(
  common_arguments,
  list(method = "base", economic_weights = common_weights)
))

rank_sum <- do.call(selection_index, c(
  common_arguments,
  list(method = "mulamba_mock")
))
```

The first two methods should agree because the economic weights were
derived from the desired gains. This is a mathematical verification. The
base and rank-sum methods provide simpler benchmarks.

Threshold methods require limits rather than a response vector. These
limits use the units of `candidates_sd`. A maximum is supplied for each
trait listed in `lower_is_better`.

``` r
limits_sd <- c(
  GY = -0.50,
  PHT = 0.75,
  AD = 0.75,
  ASI = 0.75,
  EPP = -0.75,
  GLS = 0.75
)

culling <- do.call(selection_index, c(
  common_arguments,
  list(
    method = "independent_culling",
    culling_thresholds = limits_sd
  )
))

elston <- do.call(selection_index, c(
  common_arguments,
  list(method = "elston", culling_thresholds = limits_sd)
))

tandem <- do.call(selection_index, c(
  common_arguments,
  list(method = "tandem", tandem_order = c("GY", "GLS", "AD"))
))
```

### 7.4 Assemble the comparison

``` r
comparison <- compare_selection_methods(
  list(
    Smith_Hazel = smith_hazel,
    Pesek_Baker = pesek_baker,
    Base = base_index,
    Rank_sum = rank_sum,
    Culling = culling,
    Elston = elston,
    Tandem = tandem
  ),
  target_gains = desired_gains_sd
)
comparison
#> <desiredgainr_method_comparison>
#>   Methods: Smith_Hazel, Pesek_Baker, Base, Rank_sum, Culling, Elston, Tandem 
#>   All recorded comparison conditions are satisfied.
#>       Method N_selected      R_HI  Delta_H        RE Worst_expected_attainment
#>  Smith_Hazel         20 0.4117023 5.630275 0.5268772                 0.5470366
#>  Pesek_Baker         20 0.4117023 5.630275 0.5268772                 0.5470366
#>         Base         20 0.3779319 5.168444 0.3454239                -0.1915515
#>     Rank_sum         20        NA       NA        NA                        NA
#>      Culling         20        NA       NA        NA                        NA
#>       Elston         20        NA       NA        NA                        NA
#>       Tandem         20        NA       NA        NA                        NA
#>  Mahalanobis_alignment
#>                1.00000
#>                1.00000
#>                0.83361
#>                     NA
#>                     NA
#>                     NA
#>                     NA
```

Read the fairness table first.

``` r
comparison$fairness
#>                                              Condition Satisfied
#> 1                      Same candidates and trait order      TRUE
#> 2                Same direction, centring, and scaling      TRUE
#> 3                                 Same number selected      TRUE
#> 4                             Same selection intensity      TRUE
#> 5 Common aggregate merit among merit-based comparisons      TRUE
#>                                                             Interpretation
#> 1                                  Required and enforced by this function.
#> 2                                  Required and enforced by this function.
#> 3  A culling method may select fewer candidates when few pass every limit.
#> 4                Expected responses require one common selection pressure.
#> 5 R_HI and Delta_H answer one common question only under one merit vector.
```

A culling method may select fewer than 20 candidates. That result is
useful. It shows that the limits and the requested selection quota
conflict.

### 7.5 Compare biological response

``` r
comparison$responses[
  , c(
    "Method", "Trait", "Expected_response",
    "Observed_differential", "Expected_attainment"
  )
]
#>          Method  Trait Expected_response Observed_differential
#>          <char> <char>             <num>                 <num>
#>  1: Smith_Hazel     GY        0.54703663            2.71839267
#>  2: Smith_Hazel    PHT        0.21881465            0.51969887
#>  3: Smith_Hazel     AD        0.32822198            0.64456380
#>  4: Smith_Hazel    ASI        0.27351832           -0.30315333
#>  5: Smith_Hazel    EPP        0.21881465            1.05538000
#>  6: Smith_Hazel    GLS        0.32822198            0.69842471
#>  7: Pesek_Baker     GY        0.54703663            2.71839267
#>  8: Pesek_Baker    PHT        0.21881465            0.51969887
#>  9: Pesek_Baker     AD        0.32822198            0.64456380
#> 10: Pesek_Baker    ASI        0.27351832           -0.30315333
#> 11: Pesek_Baker    EPP        0.21881465            1.05538000
#> 12: Pesek_Baker    GLS        0.32822198            0.69842471
#> 13:        Base     GY        0.35864060            2.68933267
#> 14:        Base    PHT        0.04922368            0.30407054
#> 15:        Base     AD        0.09127695            0.62602380
#> 16:        Base    ASI       -0.09577576           -0.58177000
#> 17:        Base    EPP        0.06607095            0.77743000
#> 18:        Base    GLS        0.09854367            0.13686588
#> 19:    Rank_sum     GY                NA            1.91381933
#> 20:    Rank_sum    PHT                NA            0.66151429
#> 21:    Rank_sum     AD                NA            0.80376780
#> 22:    Rank_sum    ASI                NA            1.14847167
#> 23:    Rank_sum    EPP                NA            1.24483000
#> 24:    Rank_sum    GLS                NA            1.10650706
#> 25:     Culling     GY                NA            1.48661267
#> 26:     Culling    PHT                NA            0.59683721
#> 27:     Culling     AD                NA            0.60975580
#> 28:     Culling    ASI                NA            1.07222167
#> 29:     Culling    EPP                NA            0.60343000
#> 30:     Culling    GLS                NA            0.62050118
#> 31:      Elston     GY                NA            1.37939933
#> 32:      Elston    PHT                NA            0.69501637
#> 33:      Elston     AD                NA            0.50001980
#> 34:      Elston    ASI                NA            1.07734667
#> 35:      Elston    EPP                NA            0.67248000
#> 36:      Elston    GLS                NA            0.84341294
#> 37:      Tandem     GY                NA            1.66673933
#> 38:      Tandem    PHT                NA           -0.03648821
#> 39:      Tandem     AD                NA            0.91060980
#> 40:      Tandem    ASI                NA            0.13530500
#> 41:      Tandem    EPP                NA            0.55378000
#> 42:      Tandem    GLS                NA            1.55563059
#>          Method  Trait Expected_response Observed_differential
#>          <char> <char>             <num>                 <num>
#>     Expected_attainment
#>                   <num>
#>  1:           0.5470366
#>  2:           0.5470366
#>  3:           0.5470366
#>  4:           0.5470366
#>  5:           0.5470366
#>  6:           0.5470366
#>  7:           0.5470366
#>  8:           0.5470366
#>  9:           0.5470366
#> 10:           0.5470366
#> 11:           0.5470366
#> 12:           0.5470366
#> 13:           0.3586406
#> 14:           0.1230592
#> 15:           0.1521282
#> 16:          -0.1915515
#> 17:           0.1651774
#> 18:           0.1642395
#> 19:                  NA
#> 20:                  NA
#> 21:                  NA
#> 22:                  NA
#> 23:                  NA
#> 24:                  NA
#> 25:                  NA
#> 26:                  NA
#> 27:                  NA
#> 28:                  NA
#> 29:                  NA
#> 30:                  NA
#> 31:                  NA
#> 32:                  NA
#> 33:                  NA
#> 34:                  NA
#> 35:                  NA
#> 36:                  NA
#> 37:                  NA
#> 38:                  NA
#> 39:                  NA
#> 40:                  NA
#> 41:                  NA
#> 42:                  NA
#>     Expected_attainment
#>                   <num>
```

Expected response and observed differential answer different questions.
Expected response predicts transmitted change under the covariance
model. Observed differential describes the candidates selected from this
data set. Rank, threshold, and tandem methods retain the second quantity
because a closed-form expected response is unavailable here.

Target attainment equals expected response divided by the requested
gain. A value of 1 reaches the target for that trait. A value of 0.70
reaches 70%. Negative values move in an unfavourable direction. Read the
worst trait before the mean because one severe failure can be hidden by
strong gains elsewhere.

### 7.6 Compare summary criteria

``` r
comparison$summary[
  , c(
    "Method", "N_selected", "R_HI", "Delta_H", "RE",
    "Index_heritability", "Worst_expected_attainment",
    "Mahalanobis_alignment", "Mahalanobis_residual"
  )
]
#>         Method N_selected      R_HI  Delta_H        RE Index_heritability
#>         <char>      <int>     <num>    <num>     <num>              <num>
#> 1: Smith_Hazel         20 0.4117023 5.630275 0.5268772          0.2439162
#> 2: Pesek_Baker         20 0.4117023 5.630275 0.5268772          0.2439162
#> 3:        Base         20 0.3779319 5.168444 0.3454239          0.1428325
#> 4:    Rank_sum         20        NA       NA        NA                 NA
#> 5:     Culling         20        NA       NA        NA                 NA
#> 6:      Elston         20        NA       NA        NA                 NA
#> 7:      Tandem         20        NA       NA        NA                 NA
#>    Worst_expected_attainment Mahalanobis_alignment Mahalanobis_residual
#>                        <num>                 <num>                <num>
#> 1:                 0.5470366               1.00000         2.087841e-15
#> 2:                 0.5470366               1.00000         1.064591e-15
#> 3:                -0.1915515               0.83361         3.663562e-01
#> 4:                        NA                    NA                   NA
#> 5:                        NA                    NA                   NA
#> 6:                        NA                    NA                   NA
#> 7:                        NA                    NA                   NA
```

Use each quantity for its defined purpose.

1.  Read `R_HI` and `Delta_H` only where the same aggregate merit was
    supplied.
2.  Read `RE` as efficiency for GY alone.
3.  Read index heritability as repeatability of the composite score.
4.  Read worst attainment as protection against hidden trait failure.
5.  Read Mahalanobis alignment as agreement with the complete desired
    direction.
6.  Read the residual as departure from the requested proportions.

### 7.7 Compare rankings and selected candidates

``` r
round(comparison$rank_correlation, 2)
#>             Smith_Hazel Pesek_Baker Base Rank_sum Culling Elston Tandem
#> Smith_Hazel        1.00        1.00 0.92     0.63    0.31   0.31   0.41
#> Pesek_Baker        1.00        1.00 0.92     0.63    0.31   0.31   0.41
#> Base               0.92        0.92 1.00     0.33    0.17   0.17   0.35
#> Rank_sum           0.63        0.63 0.33     1.00    0.43   0.44   0.40
#> Culling            0.31        0.31 0.17     0.43    1.00   1.00   0.24
#> Elston             0.31        0.31 0.17     0.44    1.00   1.00   0.24
#> Tandem             0.41        0.41 0.35     0.40    0.24   0.24   1.00
round(comparison$selected_jaccard, 2)
#>             Smith_Hazel Pesek_Baker Base Rank_sum Culling Elston Tandem
#> Smith_Hazel        1.00        1.00 0.74     0.33    0.18   0.18   0.25
#> Pesek_Baker        1.00        1.00 0.74     0.33    0.18   0.18   0.25
#> Base               0.74        0.74 1.00     0.21    0.11   0.11   0.21
#> Rank_sum           0.33        0.33 0.21     1.00    0.21   0.25   0.33
#> Culling            0.18        0.18 0.11     0.21    1.00   0.74   0.21
#> Elston             0.18        0.18 0.11     0.25    0.74   1.00   0.21
#> Tandem             0.25        0.25 0.21     0.33    0.21   0.21   1.00
comparison$selected_overlap
#>             Smith_Hazel Pesek_Baker Base Rank_sum Culling Elston Tandem
#> Smith_Hazel          20          20   17       10       6      6      8
#> Pesek_Baker          20          20   17       10       6      6      8
#> Base                 17          17   20        7       4      4      7
#> Rank_sum             10          10    7       20       7      8     10
#> Culling               6           6    4        7      20     17      7
#> Elston                6           6    4        8      17     20      7
#> Tandem                8           8    7       10       7      7     20
```

Rank correlation uses every candidate. Jaccard similarity uses the
selected sets. A high rank correlation can coexist with a different
parent list near the selection boundary. Therefore, report both.

## 8. Add uncertainty and validation

The point comparison above is the beginning of the analysis. A breeding
recommendation also needs stability evidence.

### 8.1 Covariance and weight uncertainty

Use
[`index_uncertainty()`](https://FAkohoue.github.io/DesiredGainR/reference/index_uncertainty.md)
to resample the genetic and phenotypic covariance matrices. Use
[`weight_sensitivity()`](https://FAkohoue.github.io/DesiredGainR/reference/weight_sensitivity.md)
when economic weights are uncertain. Report coefficient intervals,
response intervals, rank stability, and selected-set stability.

### 8.2 Candidate prediction uncertainty

Genomic estimated breeding values (GEBVs) have prediction error. When
candidate-specific prediction error covariance (PEC) is available, use
[`candidate_score_uncertainty()`](https://FAkohoue.github.io/DesiredGainR/reference/candidate_score_uncertainty.md).
Report each candidate’s score interval and selection probability. Give
particular attention to candidates near the selection boundary.

### 8.3 Value of each information source

Use
[`index_information_efficiency()`](https://FAkohoue.github.io/DesiredGainR/reference/index_information_efficiency.md)
for a Smith-Hazel or general economic index. It applies Cunningham’s
deletion-efficiency principle. A trait with little information value can
then be removed from future phenotyping plans if the saving is
operationally important.

### 8.4 Out-of-sample validation

Fit the index with training data. Fix its coefficients and
transformations. Apply
[`predict()`](https://rdrr.io/r/stats/predict.html) to a later cohort or
a held-out set. Compare predicted responses with realised breeding
values, progeny means, or later-stage trial performance. This separates
genuine prediction from reuse of the fitting data.

### 8.5 Multi-cycle validation

Use
[`simulate_selection_cycles()`](https://FAkohoue.github.io/DesiredGainR/reference/simulate_selection_cycles.md)
when the decision concerns several cycles. Compare cumulative response,
joint target attainment, genetic variance, and diversity. Use matched
scenarios and replicate seeds. Read [Multi-cycle
simulation](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-simulation.md)
for calibration and stress testing.

## 9. Choosing the preferred method

A method is preferred when five conditions are satisfied.

1.  It represents the breeder’s actual objective.
2.  Its required inputs have adequate quality.
3.  Its expected response respects every critical trait.
4.  Its ranking remains stable under credible uncertainty.
5.  Its advantage persists in validation or relevant multi-cycle
    scenarios.

The method with the largest single statistic can fail this standard. For
example, a high aggregate-merit response can hide failure of a disease
limit. A stable rank can also follow the wrong biological direction.
Therefore, the final recommendation should report the complete decision
panel.

Close the analysis with a practical statement. Name the selected method.
State why it matches the objective. Report the expected response for
every trait. State the main uncertainty. Then describe how the decision
will be validated in the breeding programme.

## 10. References

- Beavis W, Lamkey K, Mahama AA, Suza W (2023). Multiple Trait
  Selection. In Suza WP and Lamkey KR, editors, *Quantitative Genetics
  for Plant Breeding*. Iowa State University Digital Press.
- Cunningham EP (1969). The relative efficiencies of selection indexes.
  *Acta Agriculturae Scandinavica* **19**:45-48.
- Elston RC (1963). A weight-free index for the purpose of ranking or
  selection with respect to several traits at a time. *Biometrics*
  **19**:85-97.
- Kempthorne O, Nordskog AW (1959). Restricted selection indices.
  *Biometrics* **15**:10-19.
- Mulamba NN, Mock JJ (1978). Improvement of yield potential of the Eto
  Blanco maize population by breeding for plant traits. *Egyptian
  Journal of Genetics and Cytology* **7**:40-51.
- Pesek J, Baker RJ (1969). Desired improvement in relation to selection
  indices. *Canadian Journal of Plant Science* **49**:803-804.
- Satoh M (2024). Characteristics of restricted selection indices and
  geometrical interpretation of restricted breeding values. *Journal of
  Animal Breeding and Genetics* **141**:353-363.
- Smith HF (1936). A discriminant function for plant selection. *Annals
  of Eugenics* **7**:240-250.
- Yamada Y, Yokouchi K, Nishida A (1975). Selection index when genetic
  gains of individual traits are of primary concern. *Japanese Journal
  of Genetics* **50**:33-41.
