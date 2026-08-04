# Define the information used by a selection index

Mixed-model and genomic analyses can use records, family means, genomic
estimated breeding values, and environment-specific predictions as
selection information. The breeding objective can contain a different
set of genetic quantities. Let \\P=Var(x)\\, \\C=Cov(x,g)\\, and
\\G=Var(g)\\. This function stores that model and checks its joint
covariance.

## Usage

``` r
selection_information(values, P, C, G, id_col = NULL)
```

## Arguments

- values:

  Candidate-by-information matrix or data frame.

- P:

  Covariance matrix of the information variables.

- C:

  Covariance between information variables and objective traits.

- G:

  Genetic covariance matrix of the objective traits.

- id_col:

  Optional candidate identifier column in values.

## Value

An object of class desiredgainr_information.

## References

Beavis WD, Lamkey K, Mahama AA, Suza W (2023). Multiple Trait Selection.
In *Quantitative Genetics for Plant Breeding*. Iowa State University
Digital Press.

Henderson CR, Quaas RL (1976). Multiple trait evaluation using
relatives' records. *Journal of Animal Science* 43:1188-1197.

## See also

[`generalized_index()`](https://FAkohoue.github.io/DesiredGainR/reference/generalized_index.md),
[`candidate_score_uncertainty()`](https://FAkohoue.github.io/DesiredGainR/reference/candidate_score_uncertainty.md)
