# Contributing to DesiredGainR

DesiredGainR welcomes focused bug reports, documentation corrections,
tests, and scientifically justified method improvements.

Before opening an issue, please search the existing issues and reduce
the problem to a reproducible example. Include the DesiredGainR version,
R version, operating system, relevant optional-package versions, random
seed, and non-sensitive input dimensions. Do not attach confidential
germplasm, genotype, phenotype, pedigree, or trial data.

For a code contribution:

1.  open an issue for substantial behavioural or interface changes;
2.  create a branch from the current development branch;
3.  preserve existing public interfaces unless a deprecation has been
    agreed, because at least one package depends on this one;
4.  add tests covering the corrected behaviour and its failure modes;
5.  run
    [`devtools::document()`](https://devtools.r-lib.org/reference/document.html),
    [`devtools::test()`](https://devtools.r-lib.org/reference/test.html),
    and `devtools::check(args = "--as-cran")`;
6.  use British English in prose while preserving established function
    and argument names; and
7.  describe scientific assumptions, data exclusions, fallbacks, and any
    change to a breeder-facing recommendation.

## Contributions that need extra care

This package reports quantities that a breeder will act on, so three
classes of change carry more risk than their size suggests.

**Anything that alters a returned number.** State in the pull request
which outputs change and by how much, and retain the superseded quantity
under a clearly named element where a user might reasonably have
depended on it.

**Anything involving units, orientation, or row alignment.** Most
defects found during development were of this kind rather than
mathematical: a coefficient reported in inverse trait units, a
covariance matrix not oriented by `lower_is_better`, or a table filtered
by a vector in a different row order. The formula was right in every
case. Add a test that compares two derived quantities which must agree
with one another, since a single quantity in isolation rarely looks
wrong.

**Anything that adds a threshold or default constant.** Prefer measuring
and reporting to imposing a cut-off. Where a constant is unavoidable,
document what it is for and why that value.

## Tests

Tests belong in the most closely related file under `tests/testthat/`;
create a clearly named new file when no suitable file exists. Prefer
tests that check a formula against an independent construction, or
against Monte Carlo simulation, over tests that restate the
implementation.

Never commit external executables, credentials, confidential data,
generated check directories, or compiled objects.

By contributing, you agree that your contribution is distributed under
the package licence and that you will follow
[CODE_OF_CONDUCT.md](https://FAkohoue.github.io/DesiredGainR/CODE_OF_CONDUCT.md).
