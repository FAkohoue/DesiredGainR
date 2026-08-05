# Reproducing published results

## Why this vignette exists

Every other test in this package checks it against itself, or against an
equation it also implements. That establishes internal consistency and
nothing else. A formula can be implemented consistently and still be the
wrong formula.

This vignette checks the package against numbers produced by other
people, with other software, on data this package has never seen: Rahimi
and Debnath (2023), *Scientific Reports* **13**, 18977, computed in SAS
PROC IML and published in tables.

The reference values are frozen in
[`rahimi_debnath_2023()`](https://FAkohoue.github.io/DesiredGainR/reference/rahimi_debnath_2023.md).
Preserve them exactly. The moment they are derived from this package
they stop being an external check.

``` r
reference <- rahimi_debnath_2023()
reference$traits
#> [1] "plant_height"         "number_of_grain"      "number_of_row"       
#> [4] "row_length"           "leaf_length"          "hundred_grain_weight"
#> [7] "yield"
```

## Scope of the reproduction

The article’s genetic and phenotypic covariance matrices are in
supplementary data files that are absent from the publication. The index
**coefficients** therefore remain unavailable for direct recomputation.
This limits the quantities that can be checked.

It is not, because the Pesek–Baker index has an exact algebraic property
that the published table either satisfies or does not.

### The identity

Pesek–Baker sets \\\mathbf{b} = \mathbf{G}^{-1}\mathbf{d}\\. The
expected correlated response to selection is

\\\Delta = i\\\frac{\mathbf{G}\mathbf{b}}
{\sqrt{\mathbf{b}^\mathsf{T}\mathbf{P}\mathbf{b}}} =
i\\\frac{\mathbf{G}\mathbf{G}^{-1}\mathbf{d}}
{\sqrt{\mathbf{b}^\mathsf{T}\mathbf{P}\mathbf{b}}} =
\frac{i}{\sqrt{\mathbf{b}^\mathsf{T}\mathbf{P}\mathbf{b}}}\\\mathbf{d}.\\

The response is **exactly proportional to the desired gains**, with a
single constant of proportionality shared by every trait. No
approximation is involved.

So: divide the article’s published gains by the article’s published
\\d\\ vector. If the authors’ implementation and this package’s
understanding of the method agree, the seven ratios from the unrounded
calculations must be the same number. Ratios calculated from the printed
table can differ because both the numerator and denominator were
rounded.

``` r
ratio <- reference$pesek_baker_gain / reference$desired_gains
round(ratio, 5)
#>         plant_height      number_of_grain        number_of_row 
#>              0.33585              0.33586              0.33585 
#>           row_length          leaf_length hundred_grain_weight 
#>              0.33527              0.33583              0.33585 
#>                yield 
#>              0.33587
```

``` r
c(
  mean = mean(ratio),
  sd = stats::sd(ratio),
  relative_sd = stats::sd(ratio) / mean(ratio)
)
#>         mean           sd  relative_sd 
#> 0.3357687812 0.0002213001 0.0006590848
```

The displayed ratios do not agree to three or four significant figures.
Six round to 0.336 at three significant figures. `row_length` rounds to
0.335. All seven round to 0.34 at two significant figures. The relative
standard deviation is about 0.00066. This small spread is descriptive
evidence, but it does not define the reproduction tolerance.

The appropriate check propagates the printed decimal precision. The
desired gains are printed to three decimal places. Their rounding
half-unit is 0.0005. The responses are printed to four decimal places.
Their rounding half-unit is 0.00005.

``` r
desired_half_unit <- 0.0005
response_half_unit <- 0.00005

ratio_lower <-
  (reference$pesek_baker_gain - response_half_unit) /
    (reference$desired_gains + desired_half_unit)
ratio_upper <-
  (reference$pesek_baker_gain + response_half_unit) /
    (reference$desired_gains - desired_half_unit)

common_ratio_interval <- c(
  lower = max(ratio_lower),
  upper = min(ratio_upper)
)
round(common_ratio_interval, 7)
#>     lower     upper 
#> 0.3358452 0.3358549
```

The seven intervals share the common range 0.3358452 to 0.3358549. The
value 0.33585 lies inside this range. One exact proportionality constant
can therefore produce every published value after rounding.

### Rounding rules used by the test

| Published quantity | Printed precision | Rounding interval used |
|----|---:|---:|
| Desired gain | Three decimal places | Printed value \\\pm 0.0005\\ |
| Expected response | Four decimal places | Printed value \\\pm 0.00005\\ |
| Common proportionality constant | Intersection across all seven ratio intervals | 0.3358452 to 0.3358549 |

These bounds come directly from the printed decimal places. They were
not estimated from the observed ratio spread. The joint interval check
is stronger than a loose comparison tolerance. It verifies the defining
Pesek-Baker response identity against seven published traits. Direct
coefficient reproduction still requires the unpublished covariance
matrices.

## Why the reported \\R\_{HI}=0.0018\\ does not invalidate Pesek-Baker

The article reports \\R\_{HI}=0.0018\\ for Pesek-Baker and \\0.9887\\
for the Smith-Hazel optimum index. These values answer a net-merit
question. They do not measure whether Pesek-Baker achieved its desired
response direction.

Pesek-Baker starts from a desired-gain vector \\\mathbf{d}\\. This
vector states the relative genetic changes sought across traits. By
contrast, \\R\_{HI}\\ is the correlation between the index \\I\\ and an
aggregate merit \\H=\mathbf{a}^{\mathsf T}\mathbf{g}\\. It requires a
separate vector of merit weights \\\mathbf{a}\\.

The article used \\\mathbf{a}=\mathbf{d}\\. With
\\\mathbf{b}=\mathbf{G}^{-1}\mathbf{d}\\, this gives

\\R\_{HI} = \frac{\mathbf{d}^\mathsf{T}\mathbf{d}}
{\sqrt{\mathbf{d}^\mathsf{T}\mathbf{G}^{-1}\mathbf{P}\mathbf{G}^{-1}\mathbf{d}}
\\\sqrt{\mathbf{d}^\mathsf{T}\mathbf{G}\mathbf{d}}}.\\

This expression is mathematically the correlation between the
Pesek-Baker index and the constructed merit \\H\_{\mathbf
d}=\mathbf{d}^{\mathsf T}\mathbf{g}\\. The difficulty is biological
interpretation. Desired gains carry response units. Merit weights carry
value per trait unit. Substituting one for the other creates a merit
definition that changes when trait measurement units change.

The Pesek-Baker decision itself has the correct unit invariance. If
traits are rescaled by a diagonal matrix \\\mathbf{S}\\, then
\\\mathbf{d}^{\*}=\mathbf{S}\mathbf{d}\\ and
\\\mathbf{b}^{\*}=\mathbf{S}^{-1}\mathbf{b}\\. Therefore,
\\{\mathbf{b}^{\*}}^{\mathsf T}\mathbf{x}^{\*}=\mathbf{b}^{\mathsf
T}\mathbf{x}\\. The candidate ranking is unchanged. However, treating
desired gains as merit weights gives \\H\_{\mathbf
d}^{\*}=\mathbf{d}^{\mathsf T}\mathbf{S}^{2}\mathbf{g}\\. That merit
definition changes with the units.

``` r
range(reference$desired_gains)
#> [1]  0.224 45.046
max(reference$desired_gains) / min(reference$desired_gains)
#> [1] 201.0982
```

The numerical desired gains span a factor of about 200. This scale
disparity makes the substitution especially fragile. The value 0.0018
therefore gives little information about attainment of the desired-gain
objective.

### How to evaluate a desired-gain index

Desired gains define a valid multivariate response objective. Economic
weights are unnecessary for fitting that objective. The primary
evaluation should report:

1.  expected genetic response for every trait;
2.  proportionality or Mahalanobis alignment with \\\mathbf{d}\\;
3.  the attainable fraction at the planned selection intensity;
4.  worst-trait and joint attainment when targets are minimum floors;
5.  observed selection differentials for the current candidates; and
6.  stability under covariance and prediction uncertainty.

\\R\_{HI}\\ and \\\Delta H\\ remain useful when the breeding programme
also defines one common aggregate merit vector \\\mathbf{a}\\. They then
provide an additional economic comparison under the same merit
definition for every method. They are not prerequisites for a valid
desired-gain analysis.

### How DesiredGainR handles the distinction

When `aggregate_weights` is absent, DesiredGainR still fits the
desired-gain index and reports its per-trait response and attainment
diagnostics. It returns `NA` for \\R\_{HI}\\ and \\\Delta H\\ because
the separate aggregate merit \\H\\ has not been defined.

When the breeder supplies one common `aggregate_weights` vector,
DesiredGainR reports \\R\_{HI}\\ and \\\Delta H\\ against that same
merit for every fitted method. The comparison then has a common economic
yardstick.

The implied weights

\\\mathbf{w}=\mathbf{G}^{-1}\mathbf{P}\mathbf{G}^{-1}\mathbf{d}\\

identify the Smith-Hazel objective that reproduces a given Pesek-Baker
coefficient direction.
[`implied_economic_weights()`](https://FAkohoue.github.io/DesiredGainR/reference/implied_economic_weights.md)
provides this translation as a diagnostic. Using a different implied
vector for each desired-gain direction would also change the comparison
yardstick. DesiredGainR therefore keeps that choice explicit.

## The quadratic index: a reference implementation

The linear families are anchored above.
[`run_qgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_qgsi.md)
implements a different paper — Cerón-Rojas et al. (2026) — and needs its
own anchor.

That paper’s replication deposit contains the authors’ own R scripts.
The intermediate data files those scripts read were absent from the
deposit. Their published numbers therefore remain unavailable for direct
recomputation. The **algorithm** can still be reproduced. A reference
implementation provides a stronger comparison than a table because it
can be run on any input rather than only on theirs.

`tests/testthat/test-reference-implementation.R` transcribes their code
verbatim and runs both implementations on identical inputs. The
quantities compared are:

| Reference | This package |
|----|----|
| `b = P^-1 G w` | `selection_index(method = "smith_hazel")` coefficients |
| `corrHI = sqrt(b'Pb) / sqrt(w'Gw)` | `evaluate_index()$R_HI` |
| `Rs = ks * sqrt(b'Pb)` | `evaluate_index()$delta_H` |
| `Vpq = b'Pb + 2 tr(BB P BB P)` | [`run_qgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_qgsi.md) index variance |
| `VHq = w'Gw + 2 tr(A G A G)` | QGSI merit variance |

The second row is worth noting because the two expressions do not look
alike. The reference computes accuracy as a ratio of standard
deviations. This package computes a correlation. For the optimum index
\\\mathbf{b} = \mathbf{P}^{-1}\mathbf{G}\mathbf{w}\\ we have
\\\mathbf{b}'\mathbf{G}\mathbf{w} = \mathbf{b}'\mathbf{P}\mathbf{b}\\,
so they are the same number by an identity that is itself asserted as a
separate test. Agreement establishes that both parties mean the same
thing by accuracy, which a coefficient comparison alone would not.

## What this reproduction does not establish

Being explicit, because a reproduction vignette invites more confidence
than it earns.

- **The Rahimi and Debnath coefficients are not yet compared.** Their
  article points to BioStudies accession `S-BSST853`. The current
  vignette freezes the printed tables but does not yet redistribute or
  reanalyse that archive, so it verifies the defining property rather
  than claiming a complete rerun.
- **The Cerón-Rojas comparison is against code, not against published
  numbers.** It establishes that two independent implementations of the
  same algorithm agree. It does not establish that either reproduces the
  figures in the paper, because the inputs behind those figures were not
  deposited.
- **Neither anchor covers the simulation layer.** Forward simulation
  over breeding cycles is checked only against its own internal
  consistency.
- **The Monte Carlo tests remain the primary evidence** for the rest of
  the package. This vignette adds two external anchors. It does not
  replace them.

The qualitative orderings that *are* checkable, and hold:

``` r
c(
  optimum_RE = reference$optimum_method1_criteria[["RE"]],
  pesek_baker_RE = reference$pesek_baker_criteria[["RE"]]
)
#>     optimum_RE pesek_baker_RE 
#>         0.5504         0.1986
c(
  optimum_yield_gain = reference$optimum_method1_gain[["yield"]],
  pesek_baker_yield_gain = reference$pesek_baker_gain[["yield"]]
)
#>     optimum_yield_gain pesek_baker_yield_gain 
#>                 2.6488                 0.9559
```

An index using economic weights achieves more response in the main trait
than one constrained to a fixed proportional direction. That is the
expected trade-off, not a defect of either method.

## References

Rahimi, M. and Debnath, S. (2023) Estimating optimum and base selection
indices in plant and animal breeding programs by development new and
simple SAS and R codes. *Scientific Reports* **13**, 18977.

Cerón-Rojas, J.J., Montesinos-López, O.A., Montesinos-López, A., Vitale,
P., Pérez-Rodríguez, P., Fernandes, S.B., Ortiz, R. and Crossa, J.
(2026) Replication Data for: Nonlinear Genomic Selection Index
Accelerates Multi-Trait Crop Improvement. CIMMYT Research Data &
Software Repository Network, V1.
