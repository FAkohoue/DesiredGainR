# Obtaining G and P

## 1. What the two matrices mean

The genetic covariance matrix \\\mathbf{G}\\ describes additive genetic
variation and covariation in the target population. The phenotypic or
index-variable covariance matrix \\\mathbf{P}\\ describes the variables
on which the index is calculated. They answer different questions and
are not interchangeable.

The strongest source is a fitted multi-trait genetic model that reports
the trait covariance components and their uncertainty. DesiredGainR
begins after that model has been fitted. It does not estimate breeding
values from raw field records.

Never use the covariance of raw phenotypes as \\\mathbf{G}\\. It
contains environmental, residual and design variation. Likewise, the
covariance of best linear unbiased predictions or genomic estimated
breeding values (GEBVs) is a covariance of predictions and is usually
shrunken relative to total additive genetic covariance.

### 1.1 When the information and objective traits differ

The classical notation assumes that the same traits appear in the
information vector and the breeding objective. Chapter 12 supports a
wider formulation.

Define:

1.  \\\mathbf{P}\_{x}=\operatorname{Var}(\mathbf{x})\\ for the selection
    information
2.  \\\mathbf{C}=\operatorname{Cov}(\mathbf{x},\mathbf{g})\\ between
    information and objective traits
3.  \\\mathbf{G}\_{H}=\operatorname{Var}(\mathbf{g})\\ for the objective
    traits

The economic index is

\\ \mathbf{b}=\mathbf{P}\_{x}^{-1}\mathbf{C}\mathbf{a}. \\

The desired-gain index is

\\ \mathbf{b} = \mathbf{P}\_{x}^{-1}\mathbf{C} \left(
\mathbf{C}^{\mathsf T}\mathbf{P}\_{x}^{-1}\mathbf{C}
\right)^{-1}\mathbf{d}. \\

Use
[`selection_information()`](https://FAkohoue.github.io/DesiredGainR/reference/selection_information.md)
and
[`generalized_index()`](https://FAkohoue.github.io/DesiredGainR/reference/generalized_index.md)
for this case. This formulation can combine field records, family means,
genomic estimated breeding values, and environment-specific predictions.

### 1.2 Retain prediction error covariance

The prediction error variance (PEV) measures uncertainty for one trait.
The prediction error covariance (PEC) extends it across traits. Full PEC
matrices preserve the uncertainty in trait combinations.

Candidate-specific PEC matrices support:

1.  index score standard errors
2.  uncertainty intervals for candidate scores
3.  probabilities of entering the selected set
4.  targeted collection of further phenotypic evidence

Use
[`candidate_score_uncertainty()`](https://FAkohoue.github.io/DesiredGainR/reference/candidate_score_uncertainty.md)
after fitting the index. A diagonal PEV approximation remains useful
when full PEC is unavailable. State that approximation in the analysis
report.

------------------------------------------------------------------------

## 2. Import a fitted covariance safely

[`import_covariance()`](https://FAkohoue.github.io/DesiredGainR/reference/import_covariance.md)
matches trait names, restores the requested order, distinguishes
correlations from covariances and records the estimand.

``` r
shuffled_G <- dgr_G[rev(traits), rev(traits)]

imported <- import_covariance(
  shuffled_G,
  trait_cols = traits,
  source = "matrix",
  estimand = "additive genetic covariance"
)
imported
#> <desiredgainr_imported_covariance>
#>   Source: matrix 
#>   Estimand: additive genetic covariance 
#>   Traits: GY, PHT, AD, ASI, EPP, GLS 
#>   Condition number: 1.857e+04   Positive definite: yes
```

If software reports a correlation matrix and separate variances, declare
that explicitly.

``` r
imported_from_correlation <- import_covariance(
  stats::cov2cor(dgr_G),
  trait_cols = traits,
  source = "matrix",
  estimand = "additive genetic covariance",
  is_correlation = TRUE,
  variances = diag(dgr_G)
)
max(abs(imported_from_correlation$covariance - dgr_G))
#> [1] 1.776357e-15
```

------------------------------------------------------------------------

## 3. When the fitted G matrix is unavailable

[`estimate_genetic_covariance()`](https://FAkohoue.github.io/DesiredGainR/reference/estimate_genetic_covariance.md)
provides labelled working alternatives. The method name is part of the
scientific conclusion.

| Method | What it estimates | Principal limitation |
|----|----|----|
| `prediction_covariance` | Covariance of supplied genetic predictions | Does not undo prediction shrinkage |
| `pev_corrected` | Prediction covariance plus mean prediction-error covariance | Requires full, compatible cross-trait prediction-error covariance |
| `se_diagonal_corrected` | Prediction covariance with diagonal error correction | Cannot recover cross-trait prediction-error covariance |
| `relationship_adjusted` | Prediction covariance adjusted for relationships | Still does not undo prediction shrinkage |
| `adjusted_means_surrogate` | Covariance of adjusted means | Working surrogate, not additive genetic covariance |

For example:

``` r
working_G <- estimate_genetic_covariance(
  dgr_gebv,
  trait_cols = traits,
  method = "prediction_covariance"
)
working_G
#> DesiredGainR covariance estimate
#>   Method: prediction_covariance 
#>   Estimand: covariance of supplied genetic predictions 
#>   Genotypes: 200  Traits: 6
```

The printed estimand should appear in every analysis report. Do not
shorten it to “genetic covariance” unless that is what the upstream
model estimated.

------------------------------------------------------------------------

## 4. Diagnose both matrices

``` r
matrix_diagnostics(dgr_G, "G")
#> $name
#> [1] "G"
#> 
#> $dimension
#> [1] 6
#> 
#> $eigenvalues
#> [1] 1.448651e+02 5.512711e+00 9.091787e-01 4.697793e-01 1.404278e-01
#> [6] 7.801062e-03
#> 
#> $minimum_eigenvalue
#> [1] 0.007801062
#> 
#> $maximum_eigenvalue
#> [1] 144.8651
#> 
#> $reciprocal_condition
#> [1] 5.385053e-05
#> 
#> $condition_number
#> [1] 18569.92
#> 
#> $numerical_rank
#> [1] 6
#> 
#> $positive_definite
#> [1] TRUE
matrix_diagnostics(dgr_P, "P")
#> $name
#> [1] "P"
#> 
#> $dimension
#> [1] 6
#> 
#> $eigenvalues
#> [1] 240.78332142   8.27015401   1.86206326   1.40179746   1.02624362
#> [6]   0.02269007
#> 
#> $minimum_eigenvalue
#> [1] 0.02269007
#> 
#> $maximum_eigenvalue
#> [1] 240.7833
#> 
#> $reciprocal_condition
#> [1] 9.423438e-05
#> 
#> $condition_number
#> [1] 10611.84
#> 
#> $numerical_rank
#> [1] 6
#> 
#> $positive_definite
#> [1] TRUE
```

Check the following before fitting an index:

1.  Row and column names contain exactly the same traits in the analysis
    order.
2.  Units agree with the candidate values and objective.
3.  Both matrices are symmetric and finite.
4.  Required matrices are positive definite, or any rank deficiency is
    understood.
5.  \\\mathbf{P}-\mathbf{G}\\ is positive semidefinite when the matrices
    describe compatible phenotypic and additive genetic variation.
6.  The condition number is not so large that small input errors
    dominate the coefficients.

[`bend_covariance()`](https://FAkohoue.github.io/DesiredGainR/reference/bend_covariance.md)
can repair a small numerical definiteness failure. It cannot turn an
inappropriate estimand into an appropriate one. Always report the
original minimum eigenvalue and the size of the adjustment.

------------------------------------------------------------------------

## 5. What to report

For each matrix, record its source, population, traits, units, model
term, estimation method, software version, uncertainty information and
any repair. If only a working surrogate is available, repeat the
analysis once a fitted multi-trait covariance becomes available and
examine whether the ranking and selected set change.

The next step is [Multiple-trait
selection](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-index-families.md).
