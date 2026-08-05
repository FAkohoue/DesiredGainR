# Compare an iteratively optimised DGSI with a QGSI

Merges results from
[`run_dgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_dgsi.md)
and
[`run_qgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_qgsi.md)
by candidate identifier. It reports scores, ranks, and the selection
flag from each fit. DGSI means desired-gain selection index. The
[`run_dgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_dgsi.md)
object contains the established desired-gain index after the iterative
search applied by Joukhadar et al. (2024). QGSI means quadratic genomic
selection index.

## Usage

``` r
compare_dg_and_qgsi(
  dg_result,
  qgsi_result,
  id_col = "GenoID",
  include_metadata = TRUE,
  compute_rank_differences = TRUE,
  add_correlation_summary = TRUE,
  sort_by = c("DG_rank", "QGSI_rank", "DG", "QGSI", "none"),
  debug = FALSE
)
```

## Arguments

- dg_result:

  A result returned by
  [`run_dgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_dgsi.md).

- qgsi_result:

  A result returned by
  [`run_qgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_qgsi.md).

- id_col:

  Candidate identifier column.

- include_metadata:

  Whether to retain DGSI metadata columns.

- compute_rank_differences:

  Whether to add signed and absolute rank differences.

- add_correlation_summary:

  Whether to calculate Pearson score and Spearman rank correlations.

- sort_by:

  Output sorting rule.

- debug:

  Whether to print progress messages.

## Value

A list with `comparison_table`, `correlation_summary`, and
`decision_summary`. The last component reports the common candidate
count, selected counts, selected overlap, and Jaccard similarity.

## Details

The comparison is descriptive. DGSI targets a desired response
direction. QGSI represents linear, squared, and cross-product economic
value. Agreement does not validate either objective and disagreement
does not identify a winner. The breeder must first decide which
objective represents the programme.
