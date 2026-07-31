# DesiredGainR

**Computing a selection index is arithmetic. Deciding what it should
select for is the difficult part.**

DesiredGainR is organised around the breeding objective rather than
around the index. The literature is consistent that specifying economic
weights, not computing coefficients, is what obstructs routine index
use: Guimarães et al. (2021) showed that supplying arbitrary weights
eliminated the gains the same indices delivered under sensible ones, and
Covarrubias-Pazaran (2021) advises breeders not to interpret index
coefficients at all, because the desired response is the only decision
genuinely available to them.

The package therefore translates exactly between economic weights and
desired gains, tests whether a stated objective can be attained at the
planned selection intensity, recovers the objective a programme has been
applying implicitly, and reports how far a decision depends on weights
that are themselves uncertain. It provides the classical index families,
the iterative desired-gain index of Joukhadar et al. (2024), and the
quadratic genomic selection index of Cerón-Rojas et al. (2026), together
with a multi-cycle simulation layer whose founders are built from the
breeder’s own phased marker data.

DesiredGainR begins **after** the genetic evaluation. It does not
analyse field trials, fit mixed models, or estimate breeding values, and
it does not design crossing plans.

------------------------------------------------------------------------

## Installation

DesiredGainR requires R 4.1.0 or later.

``` r
install.packages("remotes")
remotes::install_github("FAkohoue/DesiredGainR",
  build_vignettes = TRUE,
  dependencies = TRUE
)
```

Set `build_vignettes = FALSE` to skip building the vignettes locally;
they remain available on the package website.

The multi-cycle simulation layer additionally requires **AlphaSimR**,
which is not installed automatically:

``` r
install.packages("AlphaSimR")
```

------------------------------------------------------------------------

## Documentation

Full documentation is at **<https://FAkohoue.github.io/DesiredGainR/>**.

| Vignette | Covers |
|----|----|
| [Introduction](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-introduction.html) | Orientation, input formats, definitions, conventions |
| [Full pipeline](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-pipeline.html) | Sixteen-stage walkthrough assuming no prior knowledge |
| [Complete workflow](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-workflow.html) | The same path condensed, with cross-references |
| [Defining a breeding objective](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-objective.html) | How to state, test and defend an objective |
| [Multi-cycle simulation](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-simulation.html) | Comparing objectives over several breeding cycles |

The [function
reference](https://FAkohoue.github.io/DesiredGainR/reference/) documents
every exported function.
[`open_desiredgain_guide()`](https://FAkohoue.github.io/DesiredGainR/reference/open_desiredgain_guide.md)
opens a plain-language Breeder’s Guide for readers who approve a
selection decision without running R.

------------------------------------------------------------------------

## Community

| Document | Purpose |
|----|----|
| [Support](https://FAkohoue.github.io/DesiredGainR/SUPPORT.md) | Where to ask, and what to record first |
| [Contributing](https://FAkohoue.github.io/DesiredGainR/CONTRIBUTING.md) | How to report a bug or propose a change |
| [Code of conduct](https://FAkohoue.github.io/DesiredGainR/CODE_OF_CONDUCT.md) | Expected behaviour, including scientific conduct |
| [Security](https://FAkohoue.github.io/DesiredGainR/SECURITY.md) | Private disclosure, including defects that return wrong numbers |
| [Citation](https://FAkohoue.github.io/DesiredGainR/CITATION.cff) | Machine-readable citation metadata |

Report reproducible bugs through the [issue
tracker](https://github.com/FAkohoue/DesiredGainR/issues). Do not place
confidential germplasm or phenotype records in public issues.

------------------------------------------------------------------------

## Citation

``` r
citation("DesiredGainR")
```

Machine-readable metadata are in
[CITATION.cff](https://FAkohoue.github.io/DesiredGainR/CITATION.cff). A
digital object identifier has not yet been assigned; do not cite a
provisional one.

Cite the method as well as the software when a specific index is used.
The relevant primary references are listed below.

------------------------------------------------------------------------

## References

- Cerón-Rojas JJ, Montesinos-López OA, Montesinos-López A, et
  al. (2026). Nonlinear genomic selection index accelerates multi-trait
  crop improvement. *Nature Communications* **17**:1991.
  <https://doi.org/10.1038/s41467-026-69890-3>
- Covarrubias-Pazaran G (2021). *Practical implementation of selection
  indices.* CGIAR Excellence in Breeding.
- Gaynor RC, Gorjanc G, Hickey JM (2021). AlphaSimR: an R package for
  breeding program simulations. *G3 Genes\|Genomes\|Genetics*
  **11**:jkaa017. <https://doi.org/10.1093/g3journal/jkaa017>
- Guimarães PHR, Melo PGS, Cordeiro ACC, Torga PP, Rangel PHN, de Castro
  AP (2021). Index selection can improve the selection efficiency in a
  rice recurrent selection population. *Euphytica* **217**:95.
  <https://doi.org/10.1007/s10681-021-02819-7>
- Joukhadar R, Li Y, Thistlethwaite R, Forrest KL, Tibbits JF, Trethowan
  R, Hayden MJ (2024). Optimising desired gain indices to maximise
  selection response. *Frontiers in Plant Science* **15**:1337388.
  <https://doi.org/10.3389/fpls.2024.1337388>
- Pesek J, Baker RJ (1969). Desired improvement in relation to selection
  indices. *Canadian Journal of Plant Science* **49**:803-804.
  <https://doi.org/10.4141/cjps69-137>
- Rahimi M, Debnath S (2023). Estimating optimum and base selection
  indices in plant and animal breeding programs by development new and
  simple SAS and R codes. *Scientific Reports* **13**:18977.
  <https://doi.org/10.1038/s41598-023-46368-6>
- Smith HF (1936). A discriminant function for plant selection. *Annals
  of Eugenics* **7**:240-250.
- Yamada Y, Yokouchi K, Nishida A (1975). Selection index when genetic
  gains of individual traits are of primary concern. *Japanese Journal
  of Genetics* **50**:33-41. <https://doi.org/10.1266/jjg.50.33>

------------------------------------------------------------------------

## Licence

MIT © Félicien Akohoue
