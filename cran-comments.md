# cran-comments.md

## Submission

This file describes the planned first submission of DesiredGainR 0.5.0.

## Test environments

- Local: Windows 11, R 4.5.0
- Configured GitHub Actions matrix: ubuntu-latest (R-devel, release,
  oldrel-1, 4.1), windows-latest (release), macos-latest (release)
- Source tarball checked with `R CMD check --as-cran`

The GitHub matrix must be green on the release commit before submission; this
file does not treat a configured job as evidence that an uncommitted working
tree has passed remotely.

## R CMD check results

Current local source-tarball check (`--as-cran --no-manual`):

0 errors | 0 warnings | 2 notes

One incoming-feasibility note contains the expected `New submission` message
and reports the three newly generated pkgdown article URLs as unavailable until
the rebuilt site is committed and deployed. The second note states that the
offline check could not verify the system clock. Package installation, code,
examples, tests, vignettes and vignette rebuilding all pass.

The local TinyTeX installation is missing the Courier font metric required by
R's no-index PDF-manual fallback. The indexed manual builds successfully. The
full Linux, Windows and macOS CI checks—including the PDF manual—must be green
on the release commit before this section is replaced by the submission result.

The expected release-commit note is:

> New submission

## Notes for the reviewer

**Suggested rather than imported dependencies.** `AlphaSimR` and `ASRgenomics`
are in `Suggests`. The simulation layer and `bend_covariance()` are the only
consumers, and neither should be able to make the package uninstallable.
Every example and test that needs them is guarded, and continuous integration
includes a job with both installed that rejects unexpected skips. Tests of
mutually exclusive dependency-absence branches are explicitly allowlisted.

**Long-running examples.** Examples that call `simulate_selection_cycles()` or
`propagate_covariance_uncertainty()` are wrapped in `\donttest{}`. They are
genuinely slow rather than merely cautious: the second rebuilds a founder
population for every covariance draw.

**The Breeder's Guide.** `inst/extdata` contains a rendered HTML guide built
from `inst/guide/DesiredGainR_Breeder_Guide.Rmd` by
`data-raw/build_breeder_guide.R`. It is a standalone document rather than a
vignette because it is intended for readers who do not run R, so it carries no
vignette engine and is located through `system.file()`.

**Reproduction of published results.** `vignette("DesiredGainR-reproduction")`
checks the package against tables published in Rahimi and Debnath (2023),
*Scientific Reports* 13:18977. The reference values are transcribed and frozen
in `rahimi_debnath_2023()`; no data from that article is redistributed.

## Downstream dependencies

`HapBlockR` calls `run_dgsi()` through `build_selection_index()`. The argument
names and result elements it reads are unchanged in this release. Two
behaviours it does not depend on did change, and are documented in `NEWS.md`:
the desired-gain unit conversion when `scale_traits = FALSE`, and the meaning
of a simulation cycle.
