# Estimate a working genetic covariance matrix

Estimates a covariance matrix for use in a desired-gain selection index
when a covariance matrix from a fitted multivariate genetic model is not
available. The estimand depends on `method`; the function never relabels
the covariance of predictions as total genetic covariance.

## Usage

``` r
estimate_genetic_covariance(
  values,
  trait_cols = NULL,
  method = c("prediction_covariance", "pev_corrected", "se_diagonal_corrected",
    "relationship_adjusted", "adjusted_means_surrogate"),
  prediction_error_covariance = NULL,
  prediction_se = NULL,
  relationship_matrix = NULL,
  ids = NULL,
  eigen_tolerance = 1e-08,
  symmetry_tolerance = 1e-08
)
```

## Arguments

- values:

  A data frame or numeric matrix containing one row per genotype and one
  column per trait.

- trait_cols:

  Character vector naming and ordering the trait columns. For a matrix,
  `colnames(values)` are used when `trait_cols` is omitted.

- method:

  Estimation method. See Details.

- prediction_error_covariance:

  Prediction-error covariance (PEV) information for
  `method = "pev_corrected"`. Supply either one common trait-by-trait
  matrix or a trait-by-trait-by-genotype array.

- prediction_se:

  Prediction standard errors for `method = "se_diagonal_corrected"`.
  Supply a named vector with one value per trait or a genotype-by-trait
  matrix or data frame.

- relationship_matrix:

  Genotype relationship matrix for `method = "relationship_adjusted"`.

- ids:

  Genotype identifiers in the row order of `values`. Required when
  `relationship_matrix` is supplied; the matrix is matched by its
  dimnames.

- eigen_tolerance:

  Relative eigenvalue tolerance used for positive semidefiniteness and
  numerical-rank checks.

- symmetry_tolerance:

  Relative tolerance used to check matrix symmetry.

## Value

An object of class `desiredgainr_covariance_estimate` containing `G`,
the estimated matrix; `estimand`; `provenance`; `assumptions`; component
matrices; and numerical diagnostics.

## Details

`method = "prediction_covariance"` returns the sample covariance of the
supplied predictions. For BLUPs or GEBVs this is a covariance of
predicted genetic values and is generally shrunken relative to total
genetic covariance.

`method = "pev_corrected"` uses \$\$\widehat{G} =
\operatorname{Cov}(\widehat{u}) + \overline{\operatorname{PEV}}.\$\$
This follows the BLUP identity \\\operatorname{PEV} =
\operatorname{Var}(u) - \operatorname{Var}(\widehat{u})\\. It
approximates total genetic covariance when the supplied values are
compatible multivariate BLUPs, their PEV matrices are on the same scale,
and the candidate sample and prediction model support that identity.
Full cross-trait PEV matrices are required.

`method = "se_diagonal_corrected"` adds mean squared prediction standard
errors to the diagonal of the prediction covariance. Because standard
errors contain no cross-trait prediction-error covariance, off-diagonal
entries remain covariances of predictions. The result is therefore a
partial working approximation, not a complete estimate of total genetic
covariance.

`method = "adjusted_means_surrogate"` returns the covariance of
across-environment adjusted means and labels it as a working surrogate
for genetic covariance. Covarrubias-Pazaran (2021) recommends this for
operational use when no genetic covariance matrix is available, on the
grounds that a mixed-model adjustment removes most environmental and
design variation. However, the adjusted means still carry estimation
error, so the result is not an estimate of additive genetic covariance
and will differ from the covariance governing response in the next
generation. Hence it should be reported as a surrogate and revisited
once a fitted multi-trait model becomes available.

`method = "relationship_adjusted"` estimates the covariance of
predictions after accounting for a supplied relationship matrix:
\$\$\widehat{G}\_{pred} = X_c^\mathsf{T} K^+ X_c /
\operatorname{rank}(K).\$\$ It does not undo BLUP or GEBV shrinkage.
Estimating total genetic covariance for related candidates requires
prediction-error covariance from the fitted multivariate model,
including the relevant relationship structure.

Raw phenotypes contain residual and environmental covariance. Their
covariance cannot identify genetic covariance without a genetic model,
replication, or relationship information, and should not be passed to
this function as genetic predictions. Best linear unbiased estimates
(BLUEs) and adjusted means occupy an intermediate position: they are not
genetic covariance estimates either, but
`method = "adjusted_means_surrogate"` makes their use explicit and
auditable rather than leaving breeders to substitute them silently.

## References

Covarrubias-Pazaran G (2021). *Practical implementation of selection
indices.* CGIAR Excellence in Breeding.
