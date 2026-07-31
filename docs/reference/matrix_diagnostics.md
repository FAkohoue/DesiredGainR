# Report the numerical conditioning of a covariance matrix

Selection-index coefficients are obtained by inverting covariance
matrices, and an ill-conditioned matrix produces coefficients that are
numerically meaningless while remaining superficially plausible. Rahimi
and Debnath (2023) report \\R\_{HI} = 0.0018\\ for a desired-gain index
on unstandardised maize traits spanning plant height in centimetres,
grain counts, and yield; that collapse is consistent with an inversion
performed on a severely ill-conditioned genetic covariance matrix.

## Usage

``` r
matrix_diagnostics(M, name = "matrix")
```

## Arguments

- M:

  Symmetric covariance matrix.

- name:

  Object name used in messages.

## Value

A list giving the eigenvalues, the reciprocal condition number, the
numerical rank, and a logical flag for positive definiteness.

## Details

Therefore inspect conditioning before trusting any index built on
original trait scales, and standardise the traits when the reciprocal
condition number is small.
