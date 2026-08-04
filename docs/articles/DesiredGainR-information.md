# Using predictions, prediction errors, and restricted responses

## 1. Separate the evidence from the breeding objective

A selection index uses two vectors. The information vector contains the
measurements used to rank candidates. The objective vector contains the
genetic quantities that the breeder seeks to improve.

The two vectors often contain the same traits. In that case,
\\\mathbf{C}=\mathbf{G}\\. However, mixed-model and genomic analyses can
provide several sources of information for a smaller set of objective
traits.

Define:

1.  \\\mathbf{P}=\operatorname{Var}(\mathbf{x})\\, the covariance of the
    information vector.
2.  \\\mathbf{C}=\operatorname{Cov}(\mathbf{x},\mathbf{g})\\, the
    covariance between information and objective traits.
3.  \\\mathbf{G}=\operatorname{Var}(\mathbf{g})\\, the genetic
    covariance of the objective traits.

The three matrices must form a valid joint covariance:

\\ \begin{bmatrix} \mathbf{P} & \mathbf{C}\\ \mathbf{C}^{\mathsf T} &
\mathbf{G} \end{bmatrix}. \\

[`selection_information()`](https://FAkohoue.github.io/DesiredGainR/reference/selection_information.md)
checks this condition.

``` r
information <- selection_information(
  values = dgr_candidates,
  P = dgr_P,
  C = dgr_G,
  G = dgr_G
)
information
#> <desiredgainr_information>
#>   Candidates: 200 
#>   Information variables: 6 
#>   Objective traits: 6 
#>   Joint minimum eigenvalue: 0.004384
```

This example uses the same six traits in both vectors. Therefore it
reproduces the classical case.

## 2. Fit economic and desired-gain objectives

For economic weights \\\mathbf{a}\\, the general solution is

\\ \mathbf{b}=\mathbf{P}^{-1}\mathbf{C}\mathbf{a}. \\

For desired gains \\\mathbf{d}\\, the solution is

\\ \mathbf{b}= \mathbf{P}^{-1}\mathbf{C} \left( \mathbf{C}^{\mathsf
T}\mathbf{P}^{-1}\mathbf{C} \right)^{-1}\mathbf{d}. \\

The expected response in the objective traits is

\\ \Delta\mathbf{g} = \frac{i\mathbf{C}^{\mathsf T}\mathbf{b}}
{\sqrt{\mathbf{b}^{\mathsf T}\mathbf{P}\mathbf{b}}}, \\

where \\i\\ is the standardised selection intensity.

``` r
desired_gains_sd <- c(
  GY = 1.0, PHT = 0.4, AD = 0.6,
  ASI = 0.5, EPP = 0.4, GLS = 0.6
)
desired_gains <- desired_gains_sd * sqrt(diag(dgr_G))

general_dg <- generalized_index(
  information,
  objective = desired_gains,
  method = "desired_gain",
  n_select = 20L,
  lower_is_better = lower_is_better
)
general_dg
#> <desiredgainr_generalized_index>
#>   Method: desired_gain 
#>   Information coefficients:
#>       GY      PHT       AD      ASI      EPP      GLS 
#>  2.37568 -0.03506 -0.47920  1.24645 -1.45683 -0.53687 
#>   Expected objective response:
#>       GY      PHT       AD      ASI      EPP      GLS 
#>  0.41028 -2.62578 -0.82055 -0.16411  0.02188 -0.27899
```

The coefficients refer to the information variables. The response refers
to the objective traits. These names can differ. The requested gains
began in genetic standard deviations. Multiplication by
`sqrt(diag(dgr_G))` converted them to the original trait units used by
\\\mathbf{P}\\, \\\mathbf{C}\\, and \\\mathbf{G}\\. The
`lower_is_better` argument then oriented the reductions.

## 3. Use genomic predictions correctly

A genomic estimated breeding value (GEBV) is a prediction. Its
covariance describes the spread of predictions. The prediction error
covariance (PEC) describes the remaining uncertainty.

For calibrated best linear unbiased prediction (BLUP):

\\ \operatorname{Cov}(\widehat{\mathbf{g}},\mathbf{g}) =
\operatorname{Var}(\widehat{\mathbf{g}}). \\

Hence a GEBV index can use:

1.  \\\mathbf{P}\_{x}=\operatorname{Var}(\widehat{\mathbf{g}})\\.
2.  \\\mathbf{C}=\operatorname{Var}(\widehat{\mathbf{g}})\\.
3.  \\\mathbf{G}=\operatorname{Var}(\widehat{\mathbf{g}})
    +\overline{\mathbf{PEC}}\\.

This identity requires calibrated predictions. Use cross-validation or a
later cohort to inspect calibration.

The full PEC remains valuable after fitting the index.
Candidate-specific PEC matrices provide score standard errors:

\\ \operatorname{SE}(I_i) = \sqrt{\mathbf{b}^{\mathsf
T}\mathbf{PEC}\_i\mathbf{b}}. \\

``` r
example_pec <- diag(stats::setNames(rep(0.10, length(traits)), traits))

score_uncertainty <- candidate_score_uncertainty(
  general_dg,
  prediction_error_covariance = example_pec,
  n_draws = 1000L,
  seed = 91L
)
head(score_uncertainty$summary)
#>         id     score score_se     lower     upper selection_probability
#>     <char>     <num>    <num>     <num>     <num>                 <num>
#> 1: CAND051 -17.00498 0.991915 -19.06810 -15.12936                 1.000
#> 2: CAND092 -19.30859 0.991915 -21.34802 -17.39222                 1.000
#> 3: CAND135 -19.85767 0.991915 -21.70744 -17.93782                 1.000
#> 4: CAND130 -19.94099 0.991915 -21.90581 -17.86567                 1.000
#> 5: CAND087 -20.67578 0.991915 -22.63881 -18.72685                 0.984
#> 6: CAND005 -20.89164 0.991915 -22.72547 -19.10280                 0.978
```

The selection probability shows how often a candidate enters the
selected set after prediction error is propagated. A candidate near the
cut-off can carry a wide score interval. That result helps the breeder
decide where further phenotyping has the greatest value.

## 4. Measure the value of each information source

Cunningham’s deletion efficiency measures the contribution of each
information variable. It compares the full optimum index with the index
formed after deleting one variable.

``` r
economic_weights <- c(
  GY = 1.0, PHT = 0.2, AD = 0.5,
  ASI = 0.4, EPP = 0.3, GLS = 0.5
)

smith_hazel <- selection_index(
  dgr_candidates,
  traits,
  method = "smith_hazel",
  G = dgr_G,
  P = dgr_P,
  economic_weights = economic_weights,
  scale_traits = TRUE,
  scale_by = "phenotypic",
  n_select = 20L
)

index_information_efficiency(smith_hazel)
#>    Information Full_index_efficiency Efficiency_after_deletion
#>         <char>                 <num>                     <num>
#> 1:          AD                     1                 0.7169894
#> 2:          GY                     1                 0.9440188
#> 3:         EPP                     1                 0.9756465
#> 4:         GLS                     1                 0.9788647
#> 5:         PHT                     1                 0.9823905
#> 6:         ASI                     1                 0.9994423
#>    Proportional_loss Response_variance_lost
#>                <num>                  <num>
#> 1:      0.2830106424            0.485926261
#> 2:      0.0559811928            0.108828492
#> 3:      0.0243535198            0.048113946
#> 4:      0.0211353254            0.041823949
#> 5:      0.0176095022            0.034908910
#> 6:      0.0005576911            0.001115071
```

An efficiency of 0.99 retains 99% of the full index standard deviation.
Such a variable may add little selection value. Its field cost can still
justify retention for another purpose. Therefore, combine this statistic
with cost, heritability, and programme needs.

## 5. Interpret restricted breeding values

Satoh (2024) showed that a restricted breeding value (RBV) is a
projection of an ordinary breeding value (BV). For a desired-gain
direction \\\mathbf{d}\\:

\\ \mathbf{g}\_{R} = \mathbf{d} \left( \mathbf{d}^{\mathsf
T}\mathbf{G}^{-1}\mathbf{d} \right)^{-1} \mathbf{d}^{\mathsf
T}\mathbf{G}^{-1}\mathbf{g} = \beta\mathbf{d}. \\

The scalar \\\beta\\ has a direct interpretation. It measures progress
along the requested direction. Larger values indicate greater progress.

``` r
rbv <- restricted_breeding_values(
  dgr_gebv,
  G = dgr_G,
  direction = desired_gains,
  lower_is_better = lower_is_better
)
rbv
#> <desiredgainr_restricted_bv>
#>   Restriction: proportional 
#>   Candidates: 200 
#>   Largest numerical violation: 2.665e-15 
#>   Beta range: -1.306 to 1.165
head(rbv$restricted)
#>               GY        PHT         AD         ASI          EPP         GLS
#> row_1 -0.3205120  2.0512769  0.6410240  0.12820481 -0.017093974  0.21794817
#> row_2 -0.1481996  0.9484776  0.2963992  0.05927985 -0.007903980  0.10077574
#> row_3 -0.3342387  2.1391276  0.6684774  0.13369547 -0.017826063  0.22728230
#> row_4  0.2730533 -1.7475409 -0.5461065 -0.10922131  0.014562841 -0.18567622
#> row_5 -0.0847987  0.5427117  0.1695974  0.03391948 -0.004522598  0.05766312
#> row_6  0.3438428 -2.2005937 -0.6876855 -0.13753711  0.018338281 -0.23381308
head(rbv$beta)
#>      row_1      row_2      row_3      row_4      row_5      row_6 
#> -0.4273494 -0.1975995 -0.4456516  0.3640710 -0.1130649  0.4584570
```

Every projected vector follows the desired proportions. The candidate
ranking can therefore use \\\beta\\ when the objective is a complete
proportional direction.

An achieved population response can be evaluated in the same space.

``` r
response_check <- evaluate_restricted_response(
  response = general_dg$expected_response,
  direction = desired_gains,
  G = dgr_G,
  lower_is_better = lower_is_better
)
response_check
#> <desiredgainr_restricted_response>
#>   Satoh beta: 0.54704 
#>   Mahalanobis alignment: 1 
#>   Mahalanobis residual: 8.1766e-15
```

The output provides three complementary quantities:

1.  `beta` measures progress along the desired-gain direction.
2.  `mahalanobis_alignment` measures directional agreement.
3.  `mahalanobis_residual` measures departure from the requested
    proportions.

## 6. Connect Yamada and restricted selection

Satoh (2024) proved that the Yamada desired-gain index is a special case
of the Kempthorne and Nordskog restricted index. A complete proportional
direction leaves a one-dimensional response space. Economic weights then
have no effect on the ranking.

``` r
proportional <- restricted_index(
  dgr_candidates,
  traits,
  method = "harville",
  G = dgr_G,
  P = dgr_P,
  target_gains = desired_gains,
  lower_is_better = lower_is_better,
  scale_traits = FALSE,
  n_select = 20L
)

proportional$constraint$satoh_response
#> <desiredgainr_restricted_response>
#>   Satoh beta: 0.54704 
#>   Mahalanobis alignment: 1 
#>   Mahalanobis residual: 2.8438e-14
head(proportional$restricted_breeding_values$beta)
#>     CAND001     CAND002     CAND003     CAND004     CAND005     CAND006 
#> -0.17752587 -0.27390211 -0.02502757  0.68023578  0.80895456 -0.58831815
```

This result gives a direct bridge between desired-gain selection,
restricted selection, and simulation.

## 7. Use the Elston index as a careful comparator

The Elston multiplicative index ranks candidates by the product of their
margins above breeder-defined floors. DesiredGainR first applies every
floor. It then ranks eligible candidates. This ordering prevents a large
value in two traits from hiding a failed floor in another trait.

``` r
elston <- selection_index(
  dgr_candidates,
  trait_cols = c("GY", "EPP"),
  method = "elston",
  culling_thresholds = c(GY = 0, EPP = 0),
  center_traits = TRUE,
  scale_traits = TRUE,
  n_select = 20L
)

head(elston$ranking)
#>         id    score  rank selected
#>     <char>    <num> <int>   <lgcl>
#> 1: CAND092 3.843129     1     TRUE
#> 2: CAND103 3.809303     2     TRUE
#> 3: CAND051 3.798083     3     TRUE
#> 4: CAND193 3.779833     4     TRUE
#> 5: CAND137 3.655901     5     TRUE
#> 6: CAND066 3.652627     6     TRUE
```

The Elston index provides a useful weight-light comparator. Its response
theory is limited. Therefore, use the observed selected differential and
simulation when it changes the operational decision.

## 8. Practical sequence

The recommended sequence is:

1.  Fit the multi-trait or multi-environment mixed model.
2.  Export GEBVs, \\\mathbf{G}\\, and full PEC information.
3.  Define \\\mathbf{P}\\, \\\mathbf{C}\\, and \\\mathbf{G}\\
    explicitly.
4.  Fit the economic or desired-gain index.
5.  Propagate candidate-specific prediction error.
6.  Inspect information deletion efficiency.
7.  Evaluate restricted response with Satoh’s criterion.
8.  Validate the decision in a later cohort.

This sequence preserves the distinction between evidence, objective,
expected response, and uncertainty.

## References

- Beavis WD, Lamkey K, Mahama AA, Suza W (2023). Multiple Trait
  Selection. In *Quantitative Genetics for Plant Breeding*. Iowa State
  University Digital Press.
- Beavis WD, Mahama AA, Suza W (2023). Multi Environment Trials: Linear
  Mixed Models. In *Quantitative Genetics for Plant Breeding*. Iowa
  State University Digital Press.
- Cunningham EP (1969). The relative efficiencies of selection indexes.
  *Acta Agriculturae Scandinavica* **19**:45-48.
- Elston RC (1963). A weight-free index for ranking or selection with
  respect to several traits at a time. *Biometrics* **19**:85-97.
- Henderson CR, Quaas RL (1976). Multiple trait evaluation using
  relatives’ records. *Journal of Animal Science* **43**:1188-1197.
- Satoh M (2024). Characteristics of restricted selection indices and
  geometrical interpretation of restricted breeding values. *Journal of
  Animal Breeding and Genetics* **141**:353-363.
  <https://doi.org/10.1111/jbg.12845>
