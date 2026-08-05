# Published selection-index results from Rahimi and Debnath (2023)

Reference values transcribed from Tables 2 and 3 of Rahimi and Debnath
(2023), for seven traits measured on 28 maize inbred lines evaluated in
a randomised complete block design with three replications.

## Usage

``` r
rahimi_debnath_2023()
```

## Format

A list with components:

- `traits`:

  The seven trait names, in the article's order.

- `desired_gains`:

  The \\d\\ vector used for the Pesek-Baker index, which the article
  states is \\\sqrt{\mathrm{diag}(\mathbf{G})}\\, the genetic standard
  deviations (Table 2, final column).

- `genetic_variances`:

  The implied diagonal of \\\mathbf{G}\\, obtained as `desired_gains^2`.

- `economic_weights_method1`:

  Unit economic weights (Table 2).

- `pesek_baker_gain`:

  Expected genetic advance per trait for the Pesek-Baker index (Table
  3).

- `pesek_baker_criteria`:

  The article's reported \\R\_{HI}\\, \\\Delta H\\, RE and \\CV_I\\ for
  the Pesek-Baker index.

- `optimum_method1_gain`:

  Expected genetic advance for the optimum (Smith-Hazel) index under
  Method 1 weights (Table 3).

- `optimum_method1_criteria`:

  Its reported criteria.

## Source

Rahimi, M. and Debnath, S. (2023) Estimating optimum and base selection
indices in plant and animal breeding programs by development new and
simple SAS and R codes. *Scientific Reports* 13, 18977.
[doi:10.1038/s41598-023-46368-6](https://doi.org/10.1038/s41598-023-46368-6)

## Details

This object freezes the article's printed table values. It is a
transcription check, not a full data reanalysis; see the reproduction
vignette for that distinction and for the archived-data accession.

See
[`vignette("DesiredGainR-reproduction")`](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-reproduction.md)
for the reproduction and its stated tolerances.

## Examples

``` r
reference <- rahimi_debnath_2023()
# The Pesek-Baker property: expected response is exactly proportional to the
# desired gains. Rounded table values show small ratio differences. The
# reproduction vignette propagates their printed decimal precision.
round(reference$pesek_baker_gain / reference$desired_gains, 5)
#>         plant_height      number_of_grain        number_of_row 
#>              0.33585              0.33586              0.33585 
#>           row_length          leaf_length hundred_grain_weight 
#>              0.33527              0.33583              0.33585 
#>                yield 
#>              0.33587 
```
