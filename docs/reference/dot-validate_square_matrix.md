# Validate a square matrix

Superseded internally by `.dgr_covariance()` and
`.dgr_qg_symmetric_matrix()`, which also check names, finiteness, and
symmetry. Retained only so that external code calling
`DesiredGainR:::.validate_square_matrix()` does not break. Scheduled for
removal once no reverse dependency relies on it.

## Usage

``` r
.validate_square_matrix(M, p, nm = "matrix")
```

## Arguments

- M:

  Matrix.

- p:

  Expected dimension.

- nm:

  Object name for messages.
