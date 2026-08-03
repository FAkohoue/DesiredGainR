# Propagate sampling error in the covariance matrices into an index

Every quantity a selection index reports treats \\\mathbf{G}\\ and
\\\mathbf{P}\\ as known constants. They are not. They are estimates,
often from a trial with far fewer independent families than the
covariance matrix has free parameters, and the index coefficients are a
non-linear function of them. An index whose coefficients cannot be
distinguished from zero is not a different index from one whose
coefficients are sharply determined, yet the point estimate presents
them identically.

## Usage

``` r
index_uncertainty(
  index,
  genetic_df,
  residual_df = NULL,
  n_draws = 500L,
  level = 0.95,
  seed = NULL
)
```

## Arguments

- index:

  A fitted object from
  [`selection_index()`](https://FAkohoue.github.io/DesiredGainR/reference/selection_index.md)
  using one of `"smith_hazel"`, `"base"`, `"pesek_baker"` or `"yamada"`.

- genetic_df:

  Degrees of freedom for the genetic covariance matrix. See the guidance
  above; this is the number of independent genetic units.

- residual_df:

  Degrees of freedom for the residual covariance matrix \\\mathbf{P} -
  \mathbf{G}\\. When `NULL`, `P` is held fixed and the reported
  uncertainty is attributable to `G` alone.

- n_draws:

  Number of resampling draws.

- level:

  Width of the reported intervals.

- seed:

  Optional seed. The random number generator state is restored on exit.

## Value

An object of class `desiredgainr_uncertainty` holding the coefficient
summary, the estimation-cost summary, the rank stability, and the draws
themselves.

## Details

This function resamples the covariance matrices from their approximate
sampling distribution, refits the index on each draw, and reports how
much the answer moves.

## What is resampled

A covariance matrix estimated from \\\nu\\ independent multivariate
normal vectors has a Wishart sampling distribution, so draws are taken
as

\$\$\mathbf{G}^{\*} \sim \mathcal{W}\_p(\nu_g, \mathbf{G}/\nu_g),\$\$

which has expectation \\\mathbf{G}\\. When `residual_df` is supplied the
residual matrix \\\mathbf{E} = \mathbf{P} - \mathbf{G}\\ is drawn the
same way and the phenotypic matrix reassembled as \\\mathbf{P}^{\*} =
\mathbf{G}^{\*} + \mathbf{E}^{\*}\\. That structure matters:
\\\mathbf{G}\\ and \\\mathbf{P}\\ are not independently estimated
quantities, because both are built from the same records, and drawing
them independently would break a correlation that damps the movement in
\\\mathbf{P}^{-1}\mathbf{G}\\.

## Choosing the degrees of freedom

This is the input that governs the width of every interval below, and it
is the one users get wrong.

`genetic_df` is the number of **independent genetic units** that
informed the genetic covariance, not the number of plots and not the
number of observations. In a half-sib trial it is on the order of the
number of families; in a diallel, the number of parents; in a clonally
replicated trial, the number of distinct clones. Replication within a
genotype sharpens the residual matrix, not the genetic one.

Supplying the number of plots where the number of families was meant
will understate every interval by roughly the square root of the
replication. When the design does not yield a clean answer, run the
function twice at the plausible extremes and report both.

## What the Wishart assumption buys and what it costs

It is an approximation. A restricted maximum likelihood estimate from an
unbalanced trial with a relationship matrix is not a sample covariance
from independent normal vectors, and its true sampling distribution has
no closed form. The Wishart matches the mean exactly, matches the
variance of a balanced design, and enforces positive semidefiniteness on
every draw, which a naive normal perturbation of the elements does not.

Where the fitting software reports an asymptotic variance-covariance
matrix of the variance components, that is a better basis than this.
Where it does not, which is the common case, this is a defensible
substitute provided the result is described as what it is.

## What is reported

The coefficient intervals answer how well determined the index is. The
second block answers the question that actually matters, which is what
the estimation error costs.

For `smith_hazel` and `base`, that cost is the classical relative
efficiency: under each drawn truth, the accuracy of the fitted index
divided by the accuracy of the index that would have been optimal for
that truth.

For `pesek_baker` and `yamada` it is sharper. The Pesek-Baker
coefficients \\\mathbf{b} = \mathbf{G}^{-1}\mathbf{d}\\ produce a
correlated response proportional to \\\mathbf{G}\mathbf{b} =
\mathbf{d}\\ exactly, so an index fitted on the true \\\mathbf{G}\\
delivers the desired gains in the requested proportions by construction.
Any departure under a drawn truth is therefore attributable to
estimation error alone, and `alignment`, the cosine between the achieved
response \\\mathbf{G}^{\*}\hat{\mathbf{b}}\\ and the desired gains
\\\mathbf{d}\\, measures it directly. An alignment interval whose lower
bound sits well below one means the desired gains were specified more
precisely than the data can deliver them.

## References

Hayes, J.F. and Hill, W.G. (1981) Modification of estimates of
parameters in the construction of genetic selection indices.
*Biometrics* 37, 483-493.

Sales, J. and Hill, W.G. (1976) Effect of sampling errors on efficiency
of selection indices. *Animal Production* 22, 1-17.

## See also

[`selection_index()`](https://FAkohoue.github.io/DesiredGainR/reference/selection_index.md),
[`bend_covariance()`](https://FAkohoue.github.io/DesiredGainR/reference/bend_covariance.md),
[`matrix_diagnostics()`](https://FAkohoue.github.io/DesiredGainR/reference/matrix_diagnostics.md),
[`predict.desiredgainr_index()`](https://FAkohoue.github.io/DesiredGainR/reference/predict.desiredgainr_index.md)

## Examples

``` r
set.seed(1)
traits <- c("yield", "protein")
G <- matrix(c(1.0, 0.2, 0.2, 0.5), 2, dimnames = list(traits, traits))
P <- matrix(c(2.5, 0.4, 0.4, 1.2), 2, dimnames = list(traits, traits))
values <- as.data.frame(matrix(
  stats::rnorm(80),
  ncol = 2, dimnames = list(paste0("G", 1:40), traits)
))
fit <- selection_index(
  values, traits,
  method = "pesek_baker", G = G, P = P,
  desired_gains = c(yield = 1, protein = 0.5), n_select = 8
)
index_uncertainty(fit, genetic_df = 30, residual_df = 120, n_draws = 200)
#> <desiredgainr_uncertainty>
#>   Method: pesek_baker 
#>   Draws: 200 usable of 200  Genetic df: 30  Residual df: 120
#>   95% intervals
#>   Mode: joint resampling of G and the residual 
#> 
#>   Coefficients:
#>      Trait  Estimate        SD     Lower    Upper Sign_stability Relative_SD
#>     <char>     <num>     <num>     <num>    <num>          <num>       <num>
#> 1:   yield 0.6760846 0.2602587 0.4152548 1.401397          1.000   0.3849499
#> 2: protein 0.5747822 0.3367344 0.1063244 1.308021          0.995   0.5858471
#> 
#>   alignment: 0.992 [0.959, 1.000]
#>   Rank correlation with the point estimate: 0.964 (lower 0.877)
#>   Selected set retained: 93.1% of 8 (lower 74.7%)
```
