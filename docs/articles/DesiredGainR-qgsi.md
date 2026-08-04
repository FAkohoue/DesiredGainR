# Quadratic genomic selection

## 1. The decision that QGSI represents

The quadratic genomic selection index (QGSI) is useful when a straight
line gives an inadequate representation of a trait combination’s value.
Three examples are an intermediate optimum, diminishing benefit beyond a
useful level, and complementarity between two traits. If every
additional unit has the same value, a linear index is clearer and should
be preferred.

For candidate \\i\\, DesiredGainR calculates

\\ I\_{qg,i}=\mathbf{w}^{\mathsf T}\widehat{\mathbf{g}}\_i+
\widehat{\mathbf{g}}\_i^{\mathsf T}\mathbf{W} \widehat{\mathbf{g}}\_i,
\\

where \\\widehat{\mathbf{g}}\_i\\ contains genomic estimated breeding
values (GEBVs), \\\mathbf{w}\\ contains linear economic weights, and
\\\mathbf{W}\\ is a symmetric matrix of squared and cross-product
economic weights. The package ranks larger QGSI scores first.

**Desired gains are not QGSI weights.** A desired gain states the
response a breeder wants. A QGSI weight states the value assigned to a
candidate’s linear, squared or cross-product genetic merit. DesiredGainR
never substitutes one for the other.

------------------------------------------------------------------------

## 2. Translate the biology into weights

For one standardised, favourable-direction trait \\x_j\\, its
contribution is

\\u_j(x_j)=w_jx_j+W\_{jj}x_j^2.\\

The signs have precise meanings.

| Specification | Interpretation |
|----|----|
| \\W\_{jj}=0\\ | The trait has a linear value. |
| \\W\_{jj}\<0\\ | The value is concave: improvement has diminishing benefit and eventually reaches an optimum. |
| \\W\_{jj}\>0\\ | Extreme values in either direction receive increasing value. Use this only when both extremes are genuinely valuable. |
| \\W\_{jk}\>0\\ | Same-sign deviations are rewarded and opposite-sign deviations are penalised. This includes both the jointly favourable and jointly unfavourable quadrants. |
| \\W\_{jk}\<0\\ | Opposite-sign deviations are rewarded and same-sign deviations are penalised. |

When \\W\_{jj}\<0\\ and \\w_j\>0\\, the single-trait optimum is
\\-w_j/(2W\_{jj})\\ standard deviations above the reference mean. This
identity provides a defensible way to set curvature. For example,
\\w_j=0.40\\ and \\W\_{jj}=-0.10\\ place the optimum at two standard
deviations.

For an off-diagonal entry, the coefficient of \\x_jx_k\\ in the score is
\\2W\_{jk}\\ because a symmetric quadratic form contains both
\\W\_{jk}\\ and \\W\_{kj}\\. Enter half of the intended cross-product
coefficient in each of the two symmetric cells.

------------------------------------------------------------------------

## 3. Prepare the GEBVs

The example contains 200 candidates and six traits. Grain yield and ears
per plant improve by increasing. Plant height, flowering time,
anthesis-silking interval and grey leaf spot severity improve by
decreasing.

``` r
dgr_traits[, c("trait", "description", "direction", "unit")]
#>   trait               description direction      unit
#> 1    GY               Grain yield  increase      t/ha
#> 2   PHT              Plant height  decrease        cm
#> 3    AD             Anthesis date  decrease      days
#> 4   ASI Anthesis-silking interval  decrease      days
#> 5   EPP            Ears per plant  increase     count
#> 6   GLS   Grey leaf spot severity  decrease score 1-9
head(dgr_gebv[, c("GenoID", traits)])
#>    GenoID      GY      PHT      AD     ASI     EPP     GLS
#> 1 CAND001  0.1040   5.6096  1.3217 -0.3749 -0.0256  0.1115
#> 2 CAND002  0.4459   1.5528  3.7261 -0.0895  0.0301 -0.8196
#> 3 CAND003  0.2930  16.4992  2.9357  0.0730  0.0255 -0.2172
#> 4 CAND004 -0.3515 -12.4867 -2.9110 -0.0533 -0.0232  0.0952
#> 5 CAND005 -0.1015  -7.5274  1.8098  0.1244 -0.0493 -0.6216
#> 6 CAND006  0.2835  -4.4497 -0.3108  0.1581 -0.0988  0.3139
```

`lower_is_better` declares the unfavourable measurement directions. The
function then changes their signs internally, so all weights below are
stated in a favourable-direction space. A positive weight for grey leaf
spot therefore rewards greater resistance, not greater disease severity.

We use `scale_traits = TRUE`. Consequently, both \\\mathbf{w}\\ and
\\\mathbf{W}\\ refer to GEBVs measured in standard-deviation units.
Without scaling, their units would depend on the original measurement
scale: \\w_j\\ has inverse trait units, \\W\_{jj}\\ has inverse squared
units, and \\W\_{jk}\\ has inverse cross-trait units.

------------------------------------------------------------------------

## 4. Define one explicit quadratic objective

The following objective rewards improvement in all six traits,
introduces diminishing returns, and adds complementarity between grain
yield and disease resistance. These numbers demonstrate the mechanics.
They are not universal recommendations.

``` r
linear_weights <- c(
  GY = 1.00, PHT = 0.15, AD = 0.10,
  ASI = 0.20, EPP = 0.40, GLS = 0.40
)

W <- diag(c(
  GY = -0.08, PHT = -0.04, AD = -0.04,
  ASI = -0.06, EPP = -0.08, GLS = -0.08
))
dimnames(W) <- list(traits, traits)

# The complete cross-product coefficient is 2 * 0.03 = 0.06.
W["GY", "GLS"] <- W["GLS", "GY"] <- 0.03
W
#>        GY   PHT    AD   ASI   EPP   GLS
#> GY  -0.08  0.00  0.00  0.00  0.00  0.03
#> PHT  0.00 -0.04  0.00  0.00  0.00  0.00
#> AD   0.00  0.00 -0.04  0.00  0.00  0.00
#> ASI  0.00  0.00  0.00 -0.06  0.00  0.00
#> EPP  0.00  0.00  0.00  0.00 -0.08  0.00
#> GLS  0.03  0.00  0.00  0.00  0.00 -0.08
```

The diagonal entries make the objective concave along each trait axis.
The positive `GY` by `GLS` term gives additional value to candidates
that are above the reference mean for both yield and resistance. It also
raises the score of candidates below the mean for both traits, because
the product of two negative deviations is positive. The linear terms and
negative diagonal curvature are expected to oppose that second effect
here, but this must be verified from the candidate-specific
contributions. Before using a cross-product term, the breeding programme
should state why the interaction has value beyond the sum of the
separate effects and inspect all four sign combinations.

------------------------------------------------------------------------

## 5. Fit the index and select candidates

``` r
qgsi <- run_qgsi(
  init_data = dgr_gebv["GenoID"],
  gebv_data = dgr_gebv,
  trait_cols = traits,
  linear_weights = linear_weights,
  W = W,
  Gamma = dgr_G,
  lower_is_better = lower_is_better,
  center_traits = TRUE,
  scale_traits = TRUE,
  n_select = 20L
)
qgsi
#> <quadratic_genomic_index>
#>   Candidates: 200 
#>   Traits: 6 
#>   Gamma: user supplied 
#>   Selected: 20 (10.0%)
#>   Model index SD: 2.7387
```

`Gamma = dgr_G` supplies the genomic covariance used by the theoretical
calculations. When `Gamma` is absent,
[`run_qgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_qgsi.md)
estimates the covariance of the supplied GEBVs and records that
estimand. It does not relabel a covariance of predictions as the total
additive genetic covariance.

Inspect the ranking together with its linear and quadratic parts.

``` r
head(qgsi$ranked_geno[, c(
  "GenoID", "LinearPart", "QuadraticPart", "QGSI", "Rank", "Selected"
)])
#>     GenoID LinearPart QuadraticPart     QGSI  Rank Selected
#>     <char>      <num>         <num>    <num> <num>   <lgcl>
#> 1: CAND148   3.912298    -1.0326361 2.879662     1     TRUE
#> 2: CAND095   3.055209    -0.6356583 2.419550     2     TRUE
#> 3: CAND150   2.710777    -0.5128079 2.197969     3     TRUE
#> 4: CAND190   3.054647    -0.8718337 2.182813     4     TRUE
#> 5: CAND102   2.616351    -0.7059104 1.910440     5     TRUE
#> 6: CAND142   2.311149    -0.4366633 1.874486     6     TRUE
```

The `LinearPart` and `QuadraticPart` are components of the same score.
Their signs do not say that a candidate is biologically good or bad by
themselves. the selection decision follows their sum, `QGSI`.

------------------------------------------------------------------------

## 6. Interpret the predicted response

``` r
qgsi$expected_gain_per_trait
#>     Trait Expected_Genetic_Gain Expected_Genetic_Gain_LinearSD
#>    <char>                 <num>                          <num>
#> 1:     GY             2.9024196                      3.0085687
#> 2:    PHT            -0.3908794                     -0.4051749
#> 3:     AD            -0.3029565                     -0.3140364
#> 4:    ASI             2.2565622                      2.3390906
#> 5:    EPP             1.7877263                      1.8531082
#> 6:    GLS             1.2331510                      1.2782506
qgsi$observed_selection_differential
#>     Trait      Mean_all Mean_selected Observed_GEBV_differential
#>    <char>         <num>         <num>                      <num>
#> 1:     GY -1.854853e-17     1.4174281                  1.4174281
#> 2:    PHT -1.275022e-17    -0.1595742                 -0.1595742
#> 3:     AD -1.405126e-18    -0.2105603                 -0.2105603
#> 4:    ASI  5.655199e-18     0.8812971                  0.8812971
#> 5:    EPP  2.389582e-18     0.9695206                  0.9695206
#> 6:    GLS  4.397524e-18     0.7940917                  0.7940917
```

These two tables answer different questions.

- `expected_gain_per_trait` is the model-based expected genetic
  response. It uses the genomic covariance and divides by the standard
  deviation of the complete index, including its quadratic variance.
- `observed_selection_differential` is the difference between the
  selected candidates’ mean GEBV and the full candidate mean. It
  describes the selected group. It is not automatically the response
  transmitted to progeny.

The expected-gain calculation is a linear-regression approximation. A
quadratic score is not normally distributed, so the normal-theory
selection differential becomes less accurate as curvature becomes
stronger. Report the assumptions stored with the result rather than
presenting the approximation as an observed outcome.

Under the centred multivariate-normal model,
\\\operatorname{Cov}(\widehat{\mathbf{g}}, \widehat{\mathbf{g}}^{\mathsf
T}\mathbf{W} \widehat{\mathbf{g}})=0\\ because the relevant third
central moments are zero. Therefore, the numerator of the
linear-regression gain is \\\mathbf{\Gamma}\mathbf{w}\\. \\\mathbf{W}\\
enters through the total index variance. This is a property of that
model, not proof that curvature has no biological effect. When curvature
materially changes the selected set, inspect the observed GEBV
differential and use simulation or later-cohort validation to assess the
non-linear decision.

``` r
qgsi$theoretical_parameters[c(
  "linear_index_variance",
  "quadratic_index_variance",
  "total_index_variance",
  "index_standard_deviation",
  "selection_intensity"
)]
#> $linear_index_variance
#> [1] 6.980418
#> 
#> $quadratic_index_variance
#> [1] 0.5199212
#> 
#> $total_index_variance
#> [1] 7.50034
#> 
#> $index_standard_deviation
#> [1] 2.738675
#> 
#> $selection_intensity
#> [1] 1.754983
```

Accuracy and mean squared prediction error require a true genetic
covariance, which is normally available only in simulation. They are
deliberately left unavailable in empirical applications rather than
being estimated by treating the GEBV covariance as truth.

------------------------------------------------------------------------

## 7. Identify what changed each candidate’s score

Candidate-specific contributions make the non-linearity auditable.

``` r
head(qgsi$linear_contributions)
#>     GenoID  Linear_GY  Linear_PHT   Linear_AD  Linear_ASI Linear_EPP
#>     <char>      <num>       <num>       <num>       <num>      <num>
#> 1: CAND001  0.2406880 -0.08109141 -0.06309848  0.23986008 -0.1954300
#> 2: CAND002  1.0776148 -0.01630549 -0.17594006  0.04885063  0.1732933
#> 3: CAND003  0.7033355 -0.25499517 -0.13884557 -0.05990562  0.1428422
#> 4: CAND004 -0.8743168  0.20790126  0.13554757  0.02462308 -0.1795424
#> 5: CAND005 -0.2623492  0.12870268 -0.08600564 -0.09430606 -0.3523194
#> 6: CAND006  0.6800807  0.07955270  0.01351684 -0.11686043 -0.6799999
#>     Linear_GLS
#>          <num>
#> 1: -0.05457917
#> 2:  0.55730321
#> 3:  0.16142957
#> 4: -0.04386745
#> 5:  0.42718538
#> 6: -0.18758851
head(qgsi$quadratic_contributions[, seq_len(
  min(8L, ncol(qgsi$quadratic_contributions))
)])
#> [1] 1 2 3 4 5 6
```

Use these tables to answer concrete questions: was a candidate selected
because of yield, because it approached an intermediate plant-height
optimum, or because it combined yield and resistance? A large
contribution should be checked against the biological justification for
the corresponding weight.

------------------------------------------------------------------------

## 8. Determine whether the quadratic term matters

Fit the linear special case by setting \\\mathbf{W}=0\\, then compare
rankings and selected sets. If the decision barely changes, the
quadratic complexity is not buying a materially different decision.

``` r
linear_only <- run_qgsi(
  init_data = dgr_gebv["GenoID"],
  gebv_data = dgr_gebv,
  trait_cols = traits,
  linear_weights = linear_weights,
  W = matrix(0, length(traits), length(traits),
    dimnames = list(traits, traits)
  ),
  Gamma = dgr_G,
  lower_is_better = lower_is_better,
  center_traits = TRUE,
  scale_traits = TRUE,
  n_select = 20L
)

rank_comparison <- merge(
  as.data.frame(qgsi$ranked_geno[, c("GenoID", "Rank", "Selected")]),
  as.data.frame(linear_only$ranked_geno[, c("GenoID", "Rank", "Selected")]),
  by = "GenoID", suffixes = c("_quadratic", "_linear")
)

c(
  spearman_rank_correlation = stats::cor(
    rank_comparison$Rank_quadratic,
    rank_comparison$Rank_linear,
    method = "spearman"
  ),
  selected_in_both = sum(
    rank_comparison$Selected_quadratic & rank_comparison$Selected_linear
  )
)
#> spearman_rank_correlation          selected_in_both 
#>                 0.9909548                15.0000000
```

This comparison is more informative than asking whether the quadratic
term is statistically non-zero. The breeder’s question is whether it
changes which candidates are retained, whether that change is stable,
and whether the non-linear value judgement is defensible.

------------------------------------------------------------------------

## 9. Minimum evidence before operational use

Do not adopt a QGSI because it is more flexible than a linear index.
Adopt it only when the programme can defend the curvature.

1.  State the biological or economic reason for every non-zero element
    of \\\mathbf{W}\\.
2.  Record whether GEBVs were centred and scaled, and define the
    reference population used for both operations.
3.  Confirm that \\\mathbf{W}\\ is symmetric and that trait names and
    directions match the GEBV columns.
4.  Compare the quadratic ranking and selected set with the linear
    special case.
5.  Inspect candidate-specific contributions for implausible dominance
    by one trait or interaction.
6.  Examine the response for every trait, not only the total QGSI score.
7.  Test rank and selected-set stability under plausible weight and
    covariance perturbations.
8.  Validate the decision in later cohorts before assigning realised
    breeding value to the selected set.

The practical output is not “QGSI is better”. It is a documented
statement that an explicit non-linear objective changed the decision in
a useful, stable and biologically defensible way.

------------------------------------------------------------------------

## 10. Related documentation

- [Multiple-trait
  selection](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-index-families.md)
  places QGSI beside the linear and desired-gain families.
- [Optimising desired
  gains](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-dgsi.md)
  explains the iterative desired-gain selection index, which solves a
  different breeding problem.
- [`run_qgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_qgsi.md)
  gives the complete argument and return-value contract.
- [`compare_dg_and_qgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/compare_dg_and_qgsi.md)
  compares rankings descriptively without implying that DGSI and QGSI
  optimise the same objective.

## References

Ceron-Rojas JJ, Montesinos-Lopez OA, Montesinos-Lopez A, et al. (2026).
Nonlinear genomic selection index accelerates multi-trait crop
improvement. *Nature Communications*, 17, 1991.
<https://doi.org/10.1038/s41467-026-69890-3>
