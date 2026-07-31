# Fit and evaluate a quadratic genomic selection index

`run_qgsi()` implements the quadratic genomic selection index (QGSI)
\$\$I\_{qg,i} = w^\mathsf{T}\hat{\gamma}\_i +
\hat{\gamma}\_i^\mathsf{T}W\hat{\gamma}\_i,\$\$ where
\\\hat{\gamma}\_i\\ is the vector of trait genomic estimated breeding
values (GEBVs), \\w\\ contains the linear economic weights, and \\W\\ is
the symmetric matrix of economic weights for squared and cross-product
genetic merit. Desired gains are not economic weights and are not used
to construct either \\w\\ or \\W\\.

## Usage

``` r
run_qgsi(
  init_data,
  gebv_data,
  trait_cols,
  linear_weights,
  W,
  id_col = "GenoID",
  reference_gebv_data = NULL,
  Gamma = NULL,
  relationship_matrix = NULL,
  true_G = NULL,
  lower_is_better = NULL,
  center_traits = TRUE,
  scale_traits = FALSE,
  missing_policy = c("error", "complete_cases", "mean_impute"),
  n_select = NULL,
  selection_proportion = NULL,
  return_contributions = TRUE,
  relationship_tolerance = 1e-08,
  symmetry_tolerance = sqrt(.Machine$double.eps),
  debug = FALSE
)

run_qgsi_desired_gain(
  init_data,
  gebv_data,
  trait_cols,
  id_col = "GenoID",
  dg,
  W_d = NULL,
  quadratic_diag_weights = NULL,
  lower_is_better = NULL,
  center_traits = FALSE,
  scale_traits = FALSE,
  impute_missing = FALSE,
  return_components = TRUE,
  debug = TRUE
)
```

## Arguments

- init_data:

  Candidate identifiers and metadata.

- gebv_data:

  Candidate identifiers and one GEBV column per trait.

- trait_cols:

  Unique names of the trait GEBV columns.

- linear_weights:

  Named linear economic-weight vector \\w\\.

- W:

  Symmetric matrix of economic weights for squared and cross-product
  terms. For traits `i` and `j`, the coefficient multiplying
  \\\gamma_i\gamma_j\\ is `2 * W[i, j]`.

- id_col:

  Candidate identifier column.

- reference_gebv_data:

  Optional reference-population GEBVs used to estimate \\\Gamma\\. The
  candidate GEBVs are used when this is `NULL`.

- Gamma:

  Optional genomic covariance matrix in the original trait coordinates.
  When supplied, it is transformed consistently with the GEBVs.

- relationship_matrix:

  Optional named genomic relationship matrix \\\Phi\\ for the reference
  genotypes. It is used only when `Gamma` is `NULL`; row and column
  names must match reference genotype identifiers.

- true_G:

  Optional true genetic covariance matrix. This is normally available
  only in simulation. Supplying it enables MSPE and squared index-merit
  correlation calculations; these quantities are not fabricated from
  empirical data when `true_G` is absent.

- lower_is_better:

  Traits for which smaller original values are favourable. Their GEBVs
  are multiplied by -1. `linear_weights` and `W` must describe this
  favourable-direction trait space.

- center_traits:

  Whether to centre GEBVs using reference means. The QGSI theory assumes
  zero-mean GEBVs, so `TRUE` is the default. When this is `FALSE` and
  `Gamma` is estimated internally, the scores use uncentred GEBVs while
  the estimated `Gamma` is a centred covariance; the reported
  theoretical parameters then describe a different index from the one
  that produced `ranked_geno`, and a warning is issued.

- scale_traits:

  Whether to divide GEBVs by reference standard deviations. Supplied
  weights must refer to the resulting scale.

- missing_policy:

  Explicit missing-value policy.

- n_select:

  Optional exact number of top-ranked candidates to select.

- selection_proportion:

  Optional selected proportion in `(0, 1]`. Specify at most one of
  `n_select` and `selection_proportion`.

- return_contributions:

  Whether to return candidate-specific linear and quadratic contribution
  tables.

- relationship_tolerance:

  Relative eigenvalue tolerance used for the Moore-Penrose inverse of
  `relationship_matrix`.

- symmetry_tolerance:

  Relative tolerance used when validating symmetric matrices.

- debug:

  Whether to print progress messages.

- dg:

  Deprecated compatibility name for `linear_weights`. It is interpreted
  as linear economic weights, not desired genetic gains.

- W_d:

  Deprecated compatibility name for `W`.

- quadratic_diag_weights:

  Optional explicit diagonal economic weights used to create `W_d` only
  when a full matrix is absent.

- impute_missing:

  Deprecated logical missing-value switch.

- return_components:

  Whether to return candidate-specific contributions.

## Value

An object of class `quadratic_genomic_index`. It contains rankings,
optional selection decisions, score contributions, the genomic
covariance estimate and provenance, and QGSI theoretical parameters.
Expected per-trait gains are the linear-regression gains in
Supplementary Equation 16 of Ceron-Rojas et al. (2026);
`observed_selection_differential` instead reports the selected
candidates' GEBV shift and is not labelled realised genetic gain.

## Details

The genomic covariance matrix \\\Gamma\\ may be supplied, estimated by
\$\$\hat{\Gamma} = g^{-1}\hat{\gamma}^{\mathsf{T}}\hat{\gamma},\$\$ or
estimated with a supplied genomic relationship matrix \\\Phi\\ as
\$\$\hat{\Gamma} =
g^{-1}\hat{\gamma}^{\mathsf{T}}\Phi^{-1}\hat{\gamma}.\$\$ A
Moore-Penrose inverse is used, and reported, when \\\Phi\\ is
positive-semidefinite but rank deficient.

## Changes in 0.3.1

Two reported quantities were corrected. Argument names, result element
names, and column names are unchanged, and the superseded quantities are
still returned under new names.

- `expected_gain_per_trait$Expected_Genetic_Gain` now divides \\\Gamma
  w\\ by the total index standard deviation \\\sqrt{w^\mathsf{T}\Gamma
  w + 2\operatorname{tr}(W\Gamma W\Gamma)}\\. Earlier releases divided
  by \\\sqrt{w^\mathsf{T}\Gamma w}\\, the purely linear index standard
  deviation, which inflated every gain when `W` was non-zero. The
  previous values are returned as `Expected_Genetic_Gain_LinearSD`.

- `theoretical_parameters$squared_index_merit_correlation` now uses the
  general \\\operatorname{Cov}(H, I)^2 /
  (\operatorname{Var}(I)\operatorname{Var}(H))\\. Earlier releases
  returned \\\operatorname{Var}(I)/\operatorname{Var}(H)\\, which equals
  the squared correlation only when the index is the MSPE-optimal
  predictor of net merit. That ratio is retained as
  `variance_ratio_index_to_merit`.

## References

Ceron-Rojas JJ, Montesinos-Lopez OA, Montesinos-Lopez A, et al. (2026).
Nonlinear genomic selection index accelerates multi-trait crop
improvement. *Nature Communications*, 17, 1991.
[doi:10.1038/s41467-026-69890-3](https://doi.org/10.1038/s41467-026-69890-3)

## Examples

``` r
traits <- c("yield", "height")
gebv <- data.frame(
  GenoID = paste0("G", 1:8),
  yield = c(-1.2, -0.7, -0.2, 0.1, 0.3, 0.7, 1.0, 1.4),
  height = c(0.8, 0.3, -0.1, -0.4, 0.5, -0.7, -0.2, -0.6)
)
W <- matrix(
  c(0.10, -0.02, -0.02, -0.05),
  2, dimnames = list(traits, traits)
)
fit <- run_qgsi(
  init_data = gebv["GenoID"],
  gebv_data = gebv,
  trait_cols = traits,
  linear_weights = c(yield = 1, height = 0.2),
  W = W,
  n_select = 2
)
fit$theoretical_parameters
#> $model_expected_index
#> [1] 0.0658125
#> 
#> $linear_index_variance
#> [1] 0.544475
#> 
#> $quadratic_index_variance
#> [1] 0.01014943
#> 
#> $total_index_variance
#> [1] 0.5546244
#> 
#> $index_standard_deviation
#> [1] 0.7447311
#> 
#> $selection_intensity
#> [1] 1.271106
#> 
#> $expected_net_merit_response
#> [1] 0.9466324
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
#> $variance_ratio_index_to_merit
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
#> [4] "Expected per-trait gains divide Gamma %*% w by the total index standard deviation, which includes the quadratic variance."                                                                                                                
#> [5] "Per-trait gains are a linear-regression approximation. The index is a quadratic form and therefore not normally distributed, so the normal-theory selection differential is inexact and the residual error grows with the curvature in W."
#> [6] "Accuracy and MSPE are not reported because true_G was not supplied; Gamma is not relabelled as true genetic covariance."                                                                                                                  
#> 
fit$expected_gain_per_trait
#>     Trait Expected_Genetic_Gain Expected_Genetic_Gain_LinearSD
#>    <char>                 <num>                          <num>
#> 1:  yield             1.0187456                       1.028197
#> 2: height            -0.4471813                      -0.451330
fit$observed_selection_differential
#>     Trait      Mean_all Mean_selected Observed_GEBV_differential
#>    <char>         <num>         <num>                      <num>
#> 1:  yield -3.122502e-17         1.025                      1.025
#> 2: height  1.040834e-17        -0.350                     -0.350
```
