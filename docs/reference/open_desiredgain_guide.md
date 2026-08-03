# Locate or Open the DesiredGainR Breeder's Guide

The PDF and HTML editions of the DesiredGainR Breeder's Guide are
standalone companions to the *Defining a breeding objective* vignette.
They cover the same material in plain language: what a breeding
objective is, why economic weights are difficult to state, how to read a
feasibility result, how to choose among the index families, and how to
judge whether a recommendation is robust. No R code is required to read
the guide.

## Usage

``` r
open_desiredgain_guide(open = TRUE, format = c("html", "pdf"))
```

## Arguments

- open:

  Logical, default `TRUE`. When `TRUE` and the session is interactive,
  the file is opened with the operating system's default application
  through [`browseURL`](https://rdrr.io/r/utils/browseURL.html). When
  `FALSE`, or in a non-interactive session, only the path is returned.

- format:

  Character, one of `"html"` (default) or `"pdf"`. HTML is the default
  because it needs only pandoc, whereas the PDF edition also requires a
  LaTeX installation and is therefore the more likely of the two to be
  absent.

## Value

The file path, invisibly, to the selected edition in the installed
package. An explicit error is raised when the requested edition is
unavailable, rather than a silent fallback to the other format.

## Details

It is intended for breeders, programme managers and reviewers who will
use or approve a selection decision but do not necessarily run R
themselves.

The PDF provides a fixed-layout reference edition and the HTML edition
renders the same content for on-screen reading. These standalone files
have no vignette engine, so they are not indexed by
[`vignette()`](https://rdrr.io/r/utils/vignette.html)/[`browseVignettes()`](https://rdrr.io/r/utils/browseVignettes.html)
like the package's `.Rmd` vignettes. Their version-controlled source is
`inst/guide/DesiredGainR_Breeder_Guide.Rmd`; both generated files ship
in `inst/extdata/` and are located via
[`system.file`](https://rdrr.io/r/base/system.file.html).

## See also

The *Defining a breeding objective* vignette
([`vignette("DesiredGainR-objective", package = "DesiredGainR")`](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-objective.md))
for the same material with full statistical detail and runnable code;
[`gain_feasibility`](https://FAkohoue.github.io/DesiredGainR/reference/gain_feasibility.md),
[`implied_economic_weights`](https://FAkohoue.github.io/DesiredGainR/reference/implied_economic_weights.md),
[`weight_sensitivity`](https://FAkohoue.github.io/DesiredGainR/reference/weight_sensitivity.md).

## Examples

``` r
# Locate without opening, which is what a script should do.
path <- try(open_desiredgain_guide(open = FALSE), silent = TRUE)
#> [open_desiredgain_guide] Guide located at: C:/Users/fakohoue/AppData/Local/R/win-library/4.5/DesiredGainR/extdata/DesiredGainR_Breeder_Guide.html
if (!inherits(path, "try-error")) file.exists(path)
#> [1] TRUE
```
