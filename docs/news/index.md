# Changelog

## DesiredGainR 0.3.0

- Removed the author-named DGSI wrapper;
  [`run_dgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_dgsi.md)
  is now the sole public desired-gain index entry point.
- Added
  [`estimate_genetic_covariance()`](https://FAkohoue.github.io/DesiredGainR/reference/estimate_genetic_covariance.md)
  with clearly distinguished prediction covariance, PEV-corrected,
  SE-diagonal-corrected, and relationship-adjusted estimators, together
  with assumptions and provenance.
- Implemented QGSI genomic covariance estimation from reference GEBVs
  using Supplementary Equations 19.1 and 19.2 of Cerón-Rojas et
  al. (2026).
- Added model-based QGSI mean, linear and quadratic variance,
  normal-selection intensity, expected net-merit response, and expected
  per-trait gains.
- Added optional simulation-only squared accuracy and mean squared
  prediction error calculations when the true genetic covariance matrix
  is supplied.
- Added exact top-ranked selection and clearly separated observed GEBV
  differentials from expected or realised genetic gain.
- Separated DGSI desired gains from QGSI linear and quadratic economic
  weights throughout the combined pipeline.
- Removed the unsupported constructor that fabricated quadratic economic
  weights from desired gains and empirical correlations.
- Added equation-level tests for scores, contributions, covariance
  estimators, response parameters, relationship-matrix alignment, and
  selection.

## DesiredGainR 0.2.0

- Renamed the package from DGQGSI to DesiredGainR.
- Added canonical
  [`run_dgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_dgsi.md)
  and
  [`run_qgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_qgsi.md)
  interfaces.
- Required an explicit genetic covariance matrix for DGSI and recorded
  the provenance of `P` and `G`.
- Changed threshold selection to eligibility followed by index ranking.
- Added automatic best-replicate selection, convergence, coefficient,
  ranking and selected-set stability diagnostics.
- Made missing-value handling explicit.
- Required an explicit symmetric quadratic-weight matrix for QGSI.
- Added candidate-specific linear, squared and cross-product
  contributions.
