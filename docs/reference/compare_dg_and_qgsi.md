# Compare DGSI and QGSI candidate rankings

Merges canonical
[`run_dgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_dgsi.md)
and
[`run_qgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_qgsi.md)
results by candidate identifier and reports score and rank agreement.
The comparison is descriptive: DGSI and QGSI optimise different breeding
objectives, so agreement is not interpreted as validation of either
method.

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

A list with `comparison_table` and `correlation_summary`.
