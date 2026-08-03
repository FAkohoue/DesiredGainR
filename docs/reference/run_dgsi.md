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
  allow_incompatible_estimated_P = FALSE,
  replicate_selection = c("holdout", "training"),
  holdout_fraction = 0.3,
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

  Named numeric vector of desired gains, expressed in **candidate
  standard-deviation units of the favourable-direction trait space**.
  This is true regardless of `scale_traits`: the realised response that
  `dg` is compared against is always divided by the candidate column
  standard deviations, so `dg = c(yield = 0.5)` requests a
  half-standard-deviation shift in the selected mean, never half a tonne
  per hectare. To express a target in original trait units, divide it by
  that trait's candidate standard deviation before passing it here. See
  Details.

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

  Optional independent candidate values. With
  `replicate_selection = "holdout"`, these data choose the stochastic
  replicate and are never used to fit its coefficients. When absent, an
  internal candidate holdout is split off before optimisation.

- return_all_reps:

  Whether to retain full replicate results.

- allow_incompatible_estimated_P:

  When `P` is estimated from `ref_data` and falls below `G` by more than
  sampling error explains, the two matrices cannot describe the same
  population and the call stops. Set this to `TRUE` to continue anyway;
  the override and the offending eigenvalue are recorded in
  `covariance_provenance$compatibility`.

- replicate_selection:

  How the winning stochastic search is chosen. `"holdout"`, the default,
  fits every replicate without the held-out candidates and chooses among
  them on that untouched set. `"training"` reproduces the earlier
  behaviour of taking the lowest training objective, which is a minimum
  over `n_rep` draws and is biased downward.

- holdout_fraction:

  Fraction of candidates reserved for choosing among replicates when
  `replicate_selection = "holdout"`.

- debug:

  Whether to print progress messages.

## Value

An object of class `desired_gain_index`. The `best_replicate` component
identifies the replicate selected automatically. The
`replicate_diagnostics`, `rank_correlation`, `coefficient_stability`,
and `selected_set_agreement` components describe optimisation stability.

## Units of `dg` and of `realised_response`

For each trait the realised response is \$\$r_j =
\frac{\bar{x}\_{j,\mathrm{selected}} -
\bar{x}\_{j,\mathrm{all}}}{s_j},\$\$ where \\s_j\\ is the standard
deviation of trait \\j\\ across all candidates in the
favourable-direction analysis space. Both `dg` and `realised_response`
are therefore in candidate standard-deviation units.

This holds whatever `scale_traits` is set to. The coefficient solve
happens in the units of `G` and `P`, which are raw trait units when
`scale_traits = FALSE`, so `dg` is converted from standard-deviation
units to that space before the solve and the coefficients are invariant
to the units the traits happen to be recorded in. Releases up to 0.5.0
passed the standard-deviation vector into the solve directly, which
asked for \\d_j\\ *raw* units per trait; with candidate standard
deviations of 10 and 1 and equal requested gains, that delivered
standardised response in the ratio 1:10 rather than 1:1.

## Empirical differential against transmitted response

`realised_response` is a property of the candidates that were selected:
the standardised differential between the selected group and the whole
set. It is not what the next generation inherits.

`theoretical_response` reports the model-based expectation \$\$\Delta G
= i \frac{\mathbf{G}\mathbf{b}}
{\sqrt{\mathbf{b}^\mathsf{T}\mathbf{P}\mathbf{b}}}\$\$ in the analysis
space, in the original trait units, and standardised. This is the
criterion on which a desired-gain index can be compared with the
classical families, which report the same quantity through
[`evaluate_index()`](https://FAkohoue.github.io/DesiredGainR/reference/evaluate_index.md).
Keep the two apart when reporting: a large differential among candidates
with a small transmitted response means selection is acting mostly on
non-heritable variation.

## Objective function

The search minimises \$\$\sum_j v_j \left(\frac{r_j -
d_j}{\max(\|d_j\|,\\0.25)}\right)^2,\$\$ with `objective_weights`
\\v_j\\. The floor of `0.25` on the denominator prevents traits with a
near-zero desired gain from dominating the objective through division by
a vanishing scale factor. Proposal standard deviations are
`sd_scale * pmax(abs(centre), 0.05)`, where the `0.05` floor keeps the
search from stalling when a component of the desired-gain vector is at
or near zero. Both floors are fixed constants in this release.

## Optimisation and reproducibility

The search is a stochastic hill climb over perturbed desired-gain
vectors: each proposal is mapped to coefficients by \\b =
P^{-1}G(G^\mathsf{T}P^{-1}G)^{-1}d\\ and accepted only if it lowers the
objective. The objective is a step function of \\b\\, because it depends
on \\b\\ only through the identity of the selected set, so a
derivative-free search is used. With the default holdout rule, the split
is made before optimisation and the replicate is chosen on candidates
that did not contribute to fitting. `validation_data`, when supplied, is
preferred to an internal split and can therefore select the replicate as
well as evaluate it. The final coefficients are applied to all
candidates only after that choice; they are not refitted. `run_dgsi()`
seeds the RNG from `seed` and restores the caller's RNG state before
returning.

Ties in the index score are broken by ascending candidate identifier, so
the selected set does not depend on input row order.

## Examples

``` r
set.seed(3)
traits <- c("yield", "disease")
candidates <- data.frame(
  GenoID = paste0("G", seq_len(40)),
  yield = rnorm(40),
  disease = rnorm(40)
)
# In practice G comes from a fitted multi-trait genetic model, or from
# estimate_genetic_covariance(); it is not the covariance of raw phenotypes.
G <- matrix(c(1.0, -0.3, -0.3, 0.8), 2, dimnames = list(traits, traits))
P <- G + diag(c(0.8, 0.6))
dimnames(P) <- list(traits, traits)

fit <- run_dgsi(
  init_data = candidates["GenoID"],
  cand_data = candidates,
  trait_cols = traits,
  dg = c(yield = 0.6, disease = 0.4),
  G = G, P = P,
  lower_is_better = "disease",
  n_select = 8,
  n_iter = 50,
  n_rep = 3,
  seed = 3
)
fit$coefficients
#>      yield    disease 
#> 0.09362514 0.04617554 
fit$realised_response
#>     yield   disease 
#> 0.9357724 0.8574219 
fit$replicate_diagnostics
#>    Replicate Objective Iteration_of_best Selected Plateau
#>        <int>     <num>             <int>    <int>  <lgcl>
#> 1:         1 0.3393454                21        8   FALSE
#> 2:         2 0.3393454                 3        8   FALSE
#> 3:         3 0.3393454                35        8   FALSE
#>    Final_window_relative_improvement Chosen
#>                                <num> <lgcl>
#> 1:                         0.7996234   TRUE
#> 2:                         0.6075916  FALSE
#> 3:                         0.7996234  FALSE
```
