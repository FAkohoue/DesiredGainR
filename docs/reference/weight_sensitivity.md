# Assess how sensitive a selection decision is to the stated objective

Economic weights and desired gains are estimates, and Guimaraes et al.
(2021) demonstrated that a poorly chosen weight vector can eliminate the
gain an index would otherwise deliver. However, a decision is often far
less sensitive to the objective than breeders fear. This function
quantifies which of the two situations applies.

## Usage

``` r
weight_sensitivity(
  economic_weights,
  values,
  G,
  P,
  n_select,
  trait_cols = names(economic_weights),
  relative_sd = 0.25,
  n_draws = 200L,
  agreement_threshold = 0.9,
  seed = 42L
)
```

## Arguments

- economic_weights:

  Named numeric vector of economic weights.

- values:

  Numeric matrix or data frame of candidate trait values.

- G:

  Genetic variance-covariance matrix, named by trait.

- P:

  Phenotypic variance-covariance matrix, named by trait.

- n_select:

  Number of candidates selected.

- trait_cols:

  Character vector naming and ordering the trait columns. Defaults to
  the names of `economic_weights`.

- relative_sd:

  Standard deviation of the log-scale perturbation applied to each
  weight.

- n_draws:

  Number of perturbed weight vectors evaluated.

- agreement_threshold:

  Selected-set overlap above which a draw is counted as reproducing the
  original decision.

- seed:

  Random seed. The caller's random number generator state is restored on
  exit.

## Value

An object of class `desiredgainr_sensitivity` giving the distribution of
selected-set agreement, the rank correlation with the original index,
the stability proportion, and the per-trait weight ratios that most
strongly drive disagreement.

## Details

Weight vectors are perturbed multiplicatively on the log scale, the
index is rebuilt for each draw, and the resulting selected sets are
compared with the set obtained from the stated objective.

## Influence is not the same as contribution

`weight_influence` measures how strongly perturbing each weight disturbs
the selected set. It is not the share of the index that the trait
contributes, which is what
[`effective_weights()`](https://FAkohoue.github.io/DesiredGainR/reference/effective_weights.md)
reports, and the two routinely disagree.

A trait carrying only a small part of the index can still be the largest
lever on the decision, because it moves the index in a direction the
remaining traits do not already cover, and the candidates near the
selection threshold are therefore reordered by it. A trait that
dominates the index may conversely be a weak lever, because the ranking
already follows it and scaling it further changes little. Read the two
diagnostics together: the first says which weights are worth arguing
about, the second says which traits the index is acting on.
