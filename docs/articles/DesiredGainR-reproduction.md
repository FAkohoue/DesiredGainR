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
understanding of the method agree, all seven ratios must be the same
number.

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

All seven agree to four significant figures. The relative spread is
around one part in a thousand, which is what rounding to four figures in
a printed table produces.

This is a stronger check than comparing coefficients would have been.
Seven independent ratios cannot coincide by accident, and the test would
fail if either party had the method wrong.

### Tolerances

| Quantity | Tolerance | Why |
|----|----|----|
| Relative spread of the seven ratios | \\2\times10^{-3}\\ | Table values are printed to four significant figures. The smallest, `row_length` = 0.224, carries the most rounding. |
| Common ratio against 0.33585 | \\10^{-3}\\ | As above |
| Each trait’s ratio individually | \\5\times10^{-3}\\ | Dominated by `row_length`, where a half-unit in the fourth figure is 0.2 per cent |

These are set by the precision of the published table, not fitted to
make the test pass.

## The \\R\_{HI} = 0.0018\\ anomaly

The article reports \\R\_{HI} = 0.0018\\ for Pesek–Baker against
\\0.9887\\ for the optimum index on the same data. Taken at face value
that says the desired-gain index is uncorrelated with merit, which would
make the entire family useless.

It is not an arithmetic error. It follows from using the desired-gain
vector \\\mathbf{d}\\ as the aggregate weights \\\mathbf{a}\\, which is
what their code does, and what this package did until version 0.5.0. The
quantity computed is then

\\R\_{HI} = \frac{\mathbf{d}^\mathsf{T}\mathbf{d}}
{\sqrt{\mathbf{d}^\mathsf{T}\mathbf{G}^{-1}\mathbf{P}\mathbf{G}^{-1}\mathbf{d}}
\\\sqrt{\mathbf{d}^\mathsf{T}\mathbf{G}\mathbf{d}}},\\

which is not the correlation between the index and net merit under any
definition of merit, and collapses when the trait scales differ by
orders of magnitude.

``` r
range(reference$desired_gains)
#> [1]  0.224 45.046
max(reference$desired_gains) / min(reference$desired_gains)
#> [1] 201.0982
```

The traits span a factor of 200. That is the whole explanation.

### What this package does instead

There is no theory under which desired gains define net merit.
DesiredGainR therefore withholds \\R\_{HI}\\ and \\\Delta H\\ unless the
analysis supplies one **common aggregate objective**. Implied economic
weights

\\\mathbf{w} = \mathbf{G}^{-1}\mathbf{P}\mathbf{G}^{-1}\mathbf{d},\\

remain available through
[`implied_economic_weights()`](https://FAkohoue.github.io/DesiredGainR/reference/implied_economic_weights.md)
as a diagnostic, but are not substituted silently because every
desired-gain direction would then be evaluated against a different
definition of merit.

``` r
set.seed(23)
traits <- reference$traits
p <- length(traits)
sd_vector <- reference$desired_gains

correlation <- diag(p)
correlation[lower.tri(correlation)] <- correlation[upper.tri(correlation)] <- 0.25
G <- outer(sd_vector, sd_vector) * correlation
dimnames(G) <- list(traits, traits)
P <- G + diag(sd_vector^2 * 1.5)
dimnames(P) <- list(traits, traits)

values <- as.data.frame(matrix(
  stats::rnorm(p * 60, sd = rep(sd_vector, each = 60)),
  ncol = p,
  dimnames = list(paste0("g", 1:60), traits)
))

fit <- selection_index(
  values, traits,
  method = "pesek_baker", G = G, P = P,
  desired_gains = sd_vector,
  aggregate_weights = reference$economic_weights_method1,
  scale_traits = FALSE, n_select = 12
)
round(fit$evaluation$R_HI, 4)
#> [1] 0.6897
```

On trait scales as disparate as this dataset’s, \\R\_{HI}\\ is an
ordinary correlation rather than a number near zero.

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
