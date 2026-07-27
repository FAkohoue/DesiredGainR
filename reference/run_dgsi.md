# Fit an optimised desired-gain selection index

`run_dgsi()` fits a desired-gain selection index (DGSI) from candidate
genetic values. All traits are first oriented so that larger values are
favourable. Independent stochastic replicates are run internally, and
the replicate with the smallest objective is returned automatically.

## Usage

``` r
run_dgsi(
  init_data,
  cand_data,
  trait_cols,
  dg,
  P = NULL,
  G,
  ref_data = NULL,
  id_col = "GenoID",
  scale_traits = FALSE,
  lower_is_better = NULL,
  select_mode = c("top_n", "eligible_top_n"),
  n_select = 100L,
  trait_min = NULL,
  empty_eligibility = c("error", "fallback_top_n"),
  missing_policy = c("error", "complete_cases", "mean_impute"),
  n_iter = 1000L,
  n_rep = 20L,
  sd_scale = 1,
  seed = 42L,
  ridge_P = 1e-06,
  ridge_M = 1e-06,
  objective_weights = NULL,
  plateau_window = 100L,
  plateau_tolerance = 1e-04,
  validation_data = NULL,
  return_all_reps = TRUE,
  debug = FALSE
)
```

## Arguments

- init_data:

  Data frame containing candidate identifiers and metadata.

- cand_data:

  Data frame containing candidate identifiers and one column per trait.

- trait_cols:

  Character vector naming the trait columns.

- dg:

  Named numeric vector of desired gains in the analysis scale.

- P:

  Phenotypic or index-variable covariance matrix. If `NULL`, an
  empirical working covariance matrix is estimated from `ref_data`; its
  provenance is reported and it is not described as a genetic
  covariance.

- G:

  Genetic covariance matrix or a complete object returned by
  [`estimate_genetic_covariance()`](https://FAkohoue.github.io/DesiredGainR/reference/estimate_genetic_covariance.md).
  This input is required. Passing the complete estimate object preserves
  its method, estimand, assumptions, and diagnostics in the DGSI result.

- ref_data:

  Optional reference data used for scaling and, when needed, estimation
  of the working `P` matrix.

- id_col:

  Candidate identifier column.

- scale_traits:

  Whether to centre and scale traits using `ref_data`.

- lower_is_better:

  Traits for which smaller original values are favourable.

- select_mode:

  Either `"top_n"` or `"eligible_top_n"`. In `"eligible_top_n"` mode,
  thresholds define eligibility and the index ranks eligible candidates;
  thresholds never replace index ranking.

- n_select:

  Maximum number of candidates retained by the optimisation.

- trait_min:

  Named favourable-direction eligibility thresholds in the transformed
  analysis scale.

- empty_eligibility:

  Action when no candidate is eligible.

- missing_policy:

  Missing-value policy. The default rejects missing values. Imputation
  is never performed silently.

- n_iter:

  Search iterations per replicate.

- n_rep:

  Number of independent stochastic replicates.

- sd_scale:

  Perturbation scale for sampled desired-gain vectors.

- seed:

  Random seed.

- ridge_P, ridge_M:

  Non-negative ridge constants used in matrix solves.

- objective_weights:

  Optional named non-negative weights for the trait-specific squared
  response deviations.

- plateau_window:

  Number of final iterations used to diagnose a plateau.

- plateau_tolerance:

  Maximum relative improvement in the plateau window.

- validation_data:

  Optional independent candidate values used only to evaluate the
  winning coefficients after optimisation.

- return_all_reps:

  Whether to retain full replicate results.

- debug:

  Whether to print progress messages.

## Value

An object of class `desired_gain_index`. The `best_replicate` component
identifies the replicate selected automatically. The
`replicate_diagnostics`, `rank_correlation`, `coefficient_stability`,
and `selected_set_agreement` components describe optimisation stability.
