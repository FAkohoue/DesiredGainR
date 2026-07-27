# DesiredGainR: DGSI and QGSI Workflows

## Scope

DesiredGainR provides two distinct multi-trait selection indices. DGSI
calibrates a linear score towards breeder-defined desired responses.
QGSI predicts a quadratic economic merit from genomic estimated breeding
values (GEBVs). The package does not infer genetic values, fit the
upstream genomic model, or treat desired gains as economic weights.

``` r
traits <- c("yield", "protein", "disease")
n <- 60
candidates <- data.frame(
  GenoID = paste0("G", seq_len(n)),
  Family = rep(paste0("F", 1:6), each = 10),
  yield = rnorm(n),
  protein = rnorm(n),
  disease = rnorm(n)
)
G <- cov(candidates[traits])
dimnames(G) <- list(traits, traits)
```

## Optimised desired-gain index

[`run_dgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_dgsi.md)
requires the genetic covariance matrix `G`. If `P` is omitted, the
function estimates an empirical working covariance from the reference
values and reports that provenance explicitly.

``` r
dgsi <- run_dgsi(
  init_data = candidates[c("GenoID", "Family")],
  cand_data = candidates,
  trait_cols = traits,
  dg = c(yield = 0.6, protein = 0.3, disease = 0.4),
  G = G,
  lower_is_better = "disease",
  n_select = 10,
  n_iter = 100,
  n_rep = 5,
  seed = 42
)

dgsi$best_replicate
#> [1] 1
dgsi$replicate_diagnostics
#>    Replicate Objective Iteration_of_best Selected Plateau
#>        <int>     <num>             <int>    <int>  <lgcl>
#> 1:         1  1.098514                66       10   FALSE
#> 2:         2  1.098514                11       10   FALSE
#> 3:         3  1.289913                22       10   FALSE
#> 4:         4  1.289913                51       10   FALSE
#> 5:         5  1.098514                35       10   FALSE
#>    Final_window_relative_improvement Chosen
#>                                <num> <lgcl>
#> 1:                         0.8017427   TRUE
#> 2:                         0.6996890  FALSE
#> 3:                         0.7671994  FALSE
#> 4:                         0.7671994  FALSE
#> 5:                         0.8017427  FALSE
head(dgsi$ranked_geno)
#>    GenoID Family SelectionIndex Eligible Selected  Rank
#>    <char> <char>          <num>   <lgcl>   <lgcl> <num>
#> 1:     G9     F1       5.024595     TRUE     TRUE     1
#> 2:    G12     F2       4.547937     TRUE     TRUE     2
#> 3:     G7     F1       3.401556     TRUE     TRUE     3
#> 4:    G53     F6       2.997152     TRUE     TRUE     4
#> 5:    G25     F3       2.787374     TRUE     TRUE     5
#> 6:    G48     F5       2.686831     TRUE     TRUE     6
```

All stochastic replicates run internally. DesiredGainR returns the
replicate with the smallest response objective automatically. The
breeder reviews the objective, plateau, coefficient, rank-correlation,
and selected-set agreement diagnostics but does not choose among runs
manually.

In `select_mode = "eligible_top_n"`, thresholds define the eligible pool
after traits have been oriented so that larger is favourable. The index
then ranks eligible candidates. Missing data stop the analysis unless an
explicit complete-case or reference-mean-imputation policy is requested.

## Approximating genetic covariance

The covariance of BLUPs or GEBVs is a covariance of predictions. It is
generally smaller than total genetic covariance because predictions are
shrunken. DesiredGainR does not silently use it as `G`.

[`estimate_genetic_covariance()`](https://FAkohoue.github.io/DesiredGainR/reference/estimate_genetic_covariance.md)
makes an approximation explicit and returns its estimand, provenance,
assumptions, component matrices, and numerical diagnostics.

``` r
G_working <- estimate_genetic_covariance(
  values = candidates,
  trait_cols = traits,
  method = "prediction_covariance"
)
G_working
#> DesiredGainR covariance estimate
#>   Method: prediction_covariance 
#>   Estimand: covariance of supplied genetic predictions 
#>   Genotypes: 60  Traits: 3
G_working$provenance
#> [1] "sample covariance of supplied genetic predictions"
```

When full cross-trait prediction-error covariance (PEV) is available for
compatible multivariate BLUPs, the package can apply

``` math

\widehat G =
\operatorname{Cov}(\widehat u) + \overline{\operatorname{PEV}}.
```

``` r
pev <- diag(c(yield = 0.10, protein = 0.08, disease = 0.12))
dimnames(pev) <- list(traits, traits)

G_approx <- estimate_genetic_covariance(
  values = candidates,
  trait_cols = traits,
  method = "pev_corrected",
  prediction_error_covariance = pev
)
G_approx$G
#>               yield     protein     disease
#> yield    1.42539038  0.01834906 -0.07791071
#> protein  0.01834906  0.91525829 -0.03908052
#> disease -0.07791071 -0.03908052  0.89497981
G_approx$assumptions
#> [1] "Rows represent the target population of candidate genotypes."                                       
#> [2] "Trait columns are genetic predictions on a common, interpretable scale."                            
#> [3] "Values are compatible multivariate BLUPs and PEV matrices are from the same fitted model and scale."
#> [4] "The BLUP covariance identity is applicable to the target candidate sample."
```

Pass the complete `G_approx` object as `G` in
[`run_dgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_dgsi.md)
to retain its estimation method, estimand, assumptions, and diagnostics
in `dgsi$covariance_provenance$G`. A genuine covariance matrix from an
upstream model can still be passed directly.

A common trait-by-trait PEV matrix or a trait-by-trait-by-genotype PEV
array can be supplied. If only prediction standard errors are available,
`method = "se_diagonal_corrected"` corrects the diagonal but cannot
recover cross-trait prediction-error covariance.
`method = "relationship_adjusted"` accounts for a named genotype
relationship matrix but does not undo prediction shrinkage.

BLUEs, adjusted phenotypic means, and raw phenotypes contain residual or
environmental covariance. Their covariance is not a genetic covariance
estimate unless an upstream genetic model or suitable genetic and
replication structure separates those components.

## Quadratic genomic selection index

For each candidate, QGSI evaluates

``` math

I_{qg} = \mathbf w^\mathsf{T}\hat{\boldsymbol\gamma}
+ \hat{\boldsymbol\gamma}^\mathsf{T}
\mathbf W\hat{\boldsymbol\gamma}.
```

The breeder supplies linear economic weights $`\mathbf w`$ and a
symmetric matrix $`\mathbf W`$ of squared and cross-product economic
weights. The off-diagonal coefficient multiplying $`\gamma_i\gamma_j`$
is $`2W_{ij}`$. Negative diagonal curvature represents stabilising
selection around the chosen origin, whereas positive curvature favours
extremes. DesiredGainR does not infer this economic objective from
desired gains or correlations.

``` r
W <- matrix(
  c(
     0.10, 0.02, -0.01,
     0.02, 0.05,  0.00,
    -0.01, 0.00, -0.08
  ),
  3,
  dimnames = list(traits, traits)
)

qgsi <- run_qgsi(
  init_data = candidates[c("GenoID", "Family")],
  gebv_data = candidates,
  trait_cols = traits,
  linear_weights = c(yield = 1, protein = 0.5, disease = 0.7),
  W = W,
  lower_is_better = "disease",
  n_select = 10
)

qgsi
#> <quadratic_genomic_index>
#>   Candidates: 60 
#>   Traits: 3 
#>   Gamma: estimated from reference GEBVs 
#>   Selected: 10 (16.7%)
#>   Model index SD: 1.4425
head(qgsi$ranked_geno)
#>    GenoID Family LinearPart QuadraticPart     QGSI  Rank Selected
#>    <char> <char>      <num>         <num>    <num> <num>   <lgcl>
#> 1:     G9     F1   3.222095   0.382769476 3.604864     1     TRUE
#> 2:    G12     F2   2.472549   0.494057578 2.966607     2     TRUE
#> 3:    G53     F6   2.440532  -0.021194995 2.419337     3     TRUE
#> 4:     G7     F1   1.994999   0.222529931 2.217529     4     TRUE
#> 5:     G1     F1   2.118862  -0.004260127 2.114602     5     TRUE
#> 6:     G4     F1   1.916362   0.094384172 2.010746     6     TRUE
```

GEBVs are centred by default because QGSI theory assumes zero means.
`lower_is_better` changes the trait direction before scoring. Economic
weights must describe the resulting favourable-direction trait space.
When scaling is requested, the weights must also refer to
standard-deviation units.

### Estimating genomic covariance

Without a supplied `Gamma`, the function estimates the genomic
covariance matrix from reference GEBVs with Supplementary Equation 19.2
of Cerón-Rojas et al. (2026):

``` math

\hat{\Gamma} = g^{-1}\hat{\gamma}^{\mathsf T}\hat{\gamma}.
```

Candidate GEBVs serve as the reference when `reference_gebv_data` is
absent. A separate reference population may be supplied. With a named
genomic relationship matrix $`\Phi`$, the function instead evaluates
Supplementary Equation 19.1:

``` math

\hat{\Gamma} =
g^{-1}\hat{\gamma}^{\mathsf T}\Phi^{-1}\hat{\gamma}.
```

A spectral Moore–Penrose inverse is used and reported for a
rank-deficient positive-semidefinite relationship matrix.

``` r
qgsi$Gamma
#>              yield    protein    disease
#> yield   1.30330054 0.01804324 0.07661220
#> protein 0.01804324 0.82133732 0.03842917
#> disease 0.07661220 0.03842917 0.76206348
qgsi$covariance_provenance
#> $source
#> [1] "estimated from reference GEBVs"
#> 
#> $equation
#> [1] "Ceron-Rojas et al. (2026), Supplementary Equation 19.2"
#> 
#> $reference_n
#> [1] 60
#> 
#> $divisor
#> [1] 60
#> 
#> $covariance_centring
#> [1] "reference GEBV column means removed"
#> 
#> $relationship_rank
#> [1] NA
#> 
#> $relationship_dimension
#> [1] NA
```

### Model-based parameters

DesiredGainR reports the model expected index, its linear and quadratic
variance components, normal-selection intensity, expected quadratic
net-merit response, and expected per-trait gains from Supplementary
Equation 16.

``` r
qgsi$theoretical_parameters
#> $model_expected_index
#> [1] 0.1096213
#> 
#> $linear_index_variance
#> [1] 2.034247
#> 
#> $quadratic_index_variance
#> [1] 0.0466331
#> 
#> $total_index_variance
#> [1] 2.08088
#> 
#> $index_standard_deviation
#> [1] 1.442525
#> 
#> $selection_intensity
#> [1] 1.499106
#> 
#> $expected_net_merit_response
#> [1] 2.162498
#> 
#> $true_merit_variance
#> [1] NA
#> 
#> $merit_index_covariance
#> [1] NA
#> 
#> $squared_index_merit_correlation
#> [1] NA
#> 
#> $mean_squared_prediction_error
#> [1] NA
#> 
#> $accuracy_and_mspe_available
#> [1] FALSE
#> 
#> $assumptions
#> [1] "GEBVs represent the stated transformed trait space."                                                                    
#> [2] "Gamma is the genomic covariance of those GEBVs."                                                                        
#> [3] "Normal-selection response uses the reported selection intensity; it is a model-based expectation."                      
#> [4] "Accuracy and MSPE are not reported because true_G was not supplied; Gamma is not relabelled as true genetic covariance."
qgsi$expected_gain_per_trait
#>      Trait Expected_Genetic_Gain
#>     <char>                 <num>
#> 1:   yield             1.4357058
#> 2: protein             0.4788791
#> 3: disease             0.6614063
#>                                                                                                     Basis
#>                                                                                                    <char>
#> 1: Ceron-Rojas et al. (2026), Supplementary Equation 16; linear regression of breeding value on the index
#> 2: Ceron-Rojas et al. (2026), Supplementary Equation 16; linear regression of breeding value on the index
#> 3: Ceron-Rojas et al. (2026), Supplementary Equation 16; linear regression of breeding value on the index
```

Squared index–merit correlation and mean squared prediction error
require `true_G`. That matrix is generally known only in simulation.
These quantities remain unavailable for empirical data instead of being
calculated under the unverifiable assumption that estimated `Gamma` is
the true covariance.

The selected candidates’ mean GEBV shift is returned separately:

``` r
qgsi$observed_selection_differential
#>      Trait      Mean_all Mean_selected Observed_GEBV_differential
#>     <char>         <num>         <num>                      <num>
#> 1:   yield -2.289835e-17    1.28040702                 1.28040702
#> 2: protein  5.551115e-18    0.08028095                 0.08028095
#> 3: disease  7.864080e-18    1.02581750                 1.02581750
```

This is an observed differential among supplied predictions, not
realised genetic gain. Realised gain requires later-cycle observations.

Candidate-specific contribution tables reproduce the linear and
quadratic score:

``` r
head(qgsi$linear_contributions)
#>    GenoID Linear_yield Linear_protein Linear_disease
#>    <char>        <num>          <num>          <num>
#> 1:     G1   1.39761928     -0.2265997     0.94784224
#> 2:     G2  -0.53803734      0.0496329     0.93160971
#> 3:     G3   0.38978924      0.2479295    -0.18498698
#> 4:     G4   0.65952344      0.6568860     0.59995208
#> 5:     G5   0.43092916     -0.4066284    -0.09641948
#> 6:     G6  -0.07946368      0.6082889     0.20208591
head(qgsi$quadratic_contributions)
#>    GenoID Quadratic_yield_x_yield Quadratic_yield_x_protein
#>    <char>                   <num>                     <num>
#> 1:     G1            0.1953339652              -0.025336009
#> 2:     G2            0.0289484178              -0.002136348
#> 3:     G3            0.0151935655               0.007731220
#> 4:     G4            0.0434971165               0.034658539
#> 5:     G5            0.0185699938              -0.014018243
#> 6:     G6            0.0006314477              -0.003866950
#>    Quadratic_protein_x_protein Quadratic_yield_x_disease
#>                          <num>                     <num>
#> 1:                 0.010269485              -0.037849217
#> 2:                 0.000492685               0.014321166
#> 3:                 0.012293806               0.002060170
#> 4:                 0.086299853              -0.011305213
#> 5:                 0.033069332               0.001187142
#> 6:                 0.074003086               0.000458814
#>    Quadratic_protein_x_disease Quadratic_disease_x_disease
#>                          <num>                       <num>
#> 1:                           0                -0.146678352
#> 2:                           0                -0.141697412
#> 3:                           0                -0.005586969
#> 4:                           0                -0.058766123
#> 5:                           0                -0.001517831
#> 6:                           0                -0.006667545
```

They explain each candidate’s score. A nonlinear QGSI contribution is
not one global additive marker-effect vector.

## Comparing the methods

The combined pipeline requires separate objectives:

``` r
both <- run_dgsi_qgsi_pipeline(
  mode = "both",
  init_data = candidates[c("GenoID", "Family")],
  trait_cols = traits,
  dg = c(yield = 0.6, protein = 0.3, disease = 0.4),
  cand_data = candidates,
  G = G,
  n_select = 10,
  n_iter = 50,
  n_rep = 3,
  gebv_data = candidates,
  qgsi_linear_weights = c(yield = 1, protein = 0.5, disease = 0.7),
  W = W,
  qgsi_n_select = 10,
  lower_is_better = "disease",
  debug = FALSE
)

head(both$comparison_result$comparison_table)
#>    GenoID Family DG_SelectionIndex Eligible DG_Selected DG_Rank QGSI_LinearPart
#>    <char> <char>             <num>   <lgcl>      <lgcl>   <num>           <num>
#> 1:     G9     F1         1.7358426     TRUE        TRUE       1       3.2220947
#> 2:    G12     F2         1.5119273     TRUE        TRUE       2       2.4725495
#> 3:     G7     F1         1.1559742     TRUE        TRUE       3       1.9949993
#> 4:    G53     F6         1.0204345     TRUE        TRUE       4       2.4405321
#> 5:     G4     F1         0.9673805     TRUE        TRUE       5       1.9163616
#> 6:    G58     F6         0.9059897     TRUE        TRUE       6       0.8766554
#>    QGSI_QuadraticPart     QGSI QGSI_Rank QGSI_Selected RankDiff_DG_minus_QGSI
#>                 <num>    <num>     <num>        <lgcl>                  <num>
#> 1:         0.38276948 3.604864         1          TRUE                      0
#> 2:         0.49405758 2.966607         2          TRUE                      0
#> 3:         0.22252993 2.217529         4          TRUE                     -1
#> 4:        -0.02119499 2.419337         3          TRUE                      1
#> 5:         0.09438417 2.010746         6          TRUE                     -1
#> 6:         0.30852975 1.185185        14         FALSE                     -8
#>    AbsRankDiff_DG_vs_QGSI
#>                     <num>
#> 1:                      0
#> 2:                      0
#> 3:                      1
#> 4:                      1
#> 5:                      1
#> 6:                      8
```

Rank agreement is descriptive. DGSI targets desired responses; QGSI
predicts a quadratic economic merit. Agreement therefore does not
validate either method.

## Interpretation and reporting

Report the upstream genetic model, target population and environments,
covariance provenance, trait directions, transformations, selection
intensity, DGSI desired responses, QGSI economic weights, random seed,
optimisation diagnostics, and independent later-cycle evaluation.

The algorithms calculate results from the supplied genetic values and
parameters. Their breeding relevance is driven by whether those inputs
represent the target programme and whether the breeder’s objectives are
scientifically and economically defensible.

``` r
packageVersion("DesiredGainR")
#> [1] '0.3.0'
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
#> [1] DesiredGainR_0.3.0
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
