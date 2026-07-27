# Run DGSI, QGSI, or both workflows

This convenience wrapper keeps the two breeding objectives separate.
`dg` contains desired responses for
[`run_dgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_dgsi.md),
whereas `qgsi_linear_weights` and `W` contain the linear and quadratic
economic weights for
[`run_qgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_qgsi.md).
A desired-gain vector is never reused as a QGSI economic-weight vector.

## Usage

``` r
run_dgsi_qgsi_pipeline(
  mode = c("both", "dg", "qgsi"),
  init_data,
  trait_cols,
  id_col = "GenoID",
  dg = NULL,
  cand_data = NULL,
  ref_data = NULL,
  P = NULL,
  G = NULL,
  select_mode = c("top_n", "trait_thresholds"),
  n_select = 100L,
  trait_min_sd = NULL,
  fallback_to_top_n = TRUE,
  n_iter = 1000L,
  n_rep = 20L,
  sd_scale = 1,
  seed = 42L,
  ridge_P = 1e-06,
  ridge_M = 1e-06,
  dg_scale_traits = FALSE,
  gebv_data = NULL,
  qgsi_linear_weights = NULL,
  W = NULL,
  qgsi_reference_gebv_data = NULL,
  Gamma = NULL,
  relationship_matrix = NULL,
  true_G = NULL,
  qgsi_center_traits = TRUE,
  qgsi_scale_traits = FALSE,
  qgsi_missing_policy = c("error", "complete_cases", "mean_impute"),
  qgsi_n_select = NULL,
  qgsi_selection_proportion = NULL,
  lower_is_better = NULL,
  merge_outputs = TRUE,
  compare_sort_by = c("DG_rank", "QGSI_rank", "DG", "QGSI", "none"),
  debug = TRUE
)
```

## Arguments

- mode:

  One of `"both"`, `"dg"`, or `"qgsi"`.

- init_data:

  Candidate identifiers and metadata.

- trait_cols:

  Trait columns shared by the requested workflows.

- id_col:

  Candidate identifier column.

- dg:

  Named desired-gain vector required for DGSI modes.

- cand_data:

  Candidate genetic values required for DGSI modes.

- ref_data:

  Optional DGSI reference values.

- P:

  Optional DGSI phenotypic or index-variable covariance matrix.

- G:

  Genetic covariance matrix required for DGSI modes.

- select_mode:

  DGSI selection mode.

- n_select:

  DGSI selection count.

- trait_min_sd:

  Optional favourable-direction DGSI eligibility thresholds.

- fallback_to_top_n:

  Whether DGSI falls back when no candidate is eligible.

- n_iter, n_rep, sd_scale, seed, ridge_P, ridge_M:

  DGSI optimisation controls passed to
  [`run_dgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_dgsi.md).

- dg_scale_traits:

  Whether DGSI values are centred and scaled.

- gebv_data:

  Candidate GEBVs required for QGSI modes.

- qgsi_linear_weights:

  Named QGSI linear economic weights.

- W:

  Symmetric QGSI squared/cross-product economic-weight matrix.

- qgsi_reference_gebv_data:

  Optional reference GEBVs for estimating `Gamma`.

- Gamma:

  Optional genomic covariance matrix.

- relationship_matrix:

  Optional named genomic relationship matrix.

- true_G:

  Optional true genetic covariance matrix for simulation-based QGSI
  accuracy and MSPE calculations.

- qgsi_center_traits, qgsi_scale_traits:

  QGSI trait transformations.

- qgsi_missing_policy:

  QGSI missing-value policy.

- qgsi_n_select:

  Optional QGSI selection count.

- qgsi_selection_proportion:

  Optional QGSI selected proportion.

- lower_is_better:

  Traits for which smaller original values are favourable.

- merge_outputs:

  Whether to merge rankings when `mode = "both"`.

- compare_sort_by:

  Sorting rule passed to
  [`compare_dg_and_qgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/compare_dg_and_qgsi.md).

- debug:

  Whether to print progress messages.

## Value

A list containing the requested `dg_result`, `qgsi_result`, and, when
requested, `comparison_result`.

## Examples

``` r
set.seed(7)
traits <- c("yield", "height")
values <- data.frame(
  GenoID = paste0("G", 1:20),
  yield = rnorm(20),
  height = rnorm(20)
)
G <- stats::cov(values[traits])
dimnames(G) <- list(traits, traits)
W <- matrix(
  c(0.05, 0, 0, -0.03), 2,
  dimnames = list(traits, traits)
)
result <- run_dgsi_qgsi_pipeline(
  mode = "both",
  init_data = values["GenoID"],
  trait_cols = traits,
  dg = c(yield = 0.5, height = 0.2),
  cand_data = values,
  G = G,
  n_select = 4,
  n_iter = 20,
  n_rep = 2,
  gebv_data = values,
  qgsi_linear_weights = c(yield = 1, height = 0.2),
  W = W,
  qgsi_n_select = 4,
  debug = FALSE
)
```
