# Optimising desired gains

## 1. What the iterative method changes

The classical desired-gain coefficient solve fixes a response direction
under the covariance model. The iterative desired-gain selection index
(DGSI) of Joukhadar et al. searches nearby desired-gain vectors and
evaluates the selected candidates directly. Its purpose is to reduce
disagreement between the requested standardised gain vector and the
differential of the selected group.

[`run_dgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_dgsi.md)
does not discover the breeding objective. The breeder must still state
the direction, minimum requirements or acceptable intervals. Test
attainability before interpreting an optimisation result.

### DGSI is not fused into `selection_index()`

The functions have separate contracts.

| Function | What it fits |
|----|----|
| `selection_index(method = "pesek_baker")` | The closed-form Pesek–Baker desired-gain index. |
| `selection_index(method = "yamada")` | The closed-form Yamada desired-gain index. For a square, invertible \\\mathbf{G}\\ it has the same coefficient direction as Pesek–Baker. |
| [`run_dgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_dgsi.md) | An iterative search over proposed desired-gain directions, evaluated through the candidates selected by each proposal. |

[`selection_index()`](https://FAkohoue.github.io/DesiredGainR/reference/selection_index.md)
never starts the iterative search. Conversely,
[`run_dgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_dgsi.md)
retains its own non-iterated Yamada solution in `dgsi$non_iterated`, so
the breeder can determine what iteration changed.

For a proposed direction \\\mathbf{d}\\, the coefficient map is

\\ \mathbf{b}(\mathbf{d})= \mathbf{P}^{-1}\mathbf{G}
(\mathbf{G}^{\mathsf T}\mathbf{P}^{-1}\mathbf{G})^{-1}\mathbf{d}. \\

After ranking candidates with \\\mathbf{b}(\mathbf{d})\\, the search
compares the standardised differential \\r_j\\ with the requested gain
\\d_j\\ through

\\ L(\mathbf{d})=\sum_j v_j
\left\\\frac{r_j-d_j}{\max(\|d_j\|,0.25)}\right\\^{2}, \\

where \\v_j\\ is the optional `objective_weights` value. The denominator
floor prevents a near-zero target from dominating the search. Because
the selected set changes in steps as coefficients move, the objective is
not smooth and the implementation uses a reproducible derivative-free
search.

``` r
desired_gains <- c(
  GY = 1.0, PHT = 0.4, AD = 0.6,
  ASI = 0.5, EPP = 0.4, GLS = 0.6
)

gain_feasibility(
  desired_gains, dgr_G, dgr_P,
  n_candidates = nrow(dgr_candidates), n_select = 20L,
  lower_is_better = lower_is_better,
  gain_units = "genetic_sd"
)
#> <desiredgainr_feasibility>
#>   Required selection intensity: 3.2082
#>   Requires the top 0.1773%, which is fewer than one of 200 candidates
#>   Planned intensity: 1.7550 (top 10.0%)
#>   Feasible at planned intensity: no
#>   Feasible anywhere in this population: no
#>   Attainable fraction of the requested gain: 54.7%
```

If the exact ray is unattainable, optimisation can find a closer
selected set; it cannot create unavailable genetic variation. If the
entries are minimum floors rather than an exact ratio, use
[`suggest_desired_gains()`](https://FAkohoue.github.io/DesiredGainR/reference/suggest_desired_gains.md)
or interval optimisation and assess joint attainment instead.

------------------------------------------------------------------------

## 2. Fit independent searches

The example below is deliberately small enough for a vignette.
Operational analyses should use more iterations and independent
replicates, followed by holdout or external validation.

``` r
dgsi <- run_dgsi(
  init_data = dgr_candidates["GenoID"],
  cand_data = dgr_candidates,
  trait_cols = traits,
  dg = desired_gains,
  G = dgr_G,
  P = dgr_P,
  lower_is_better = lower_is_better,
  n_select = 20L,
  n_iter = 60L,
  n_rep = 3L,
  seed = 742L,
  replicate_selection = "holdout"
)
dgsi
#> <desired_gain_index>
#>   Traits: 6   Candidates: 200   Selected: 20
#>   Selection rule: top_n 
#>   Coefficients:
#>      GY     PHT      AD     ASI     EPP     GLS 
#>  3.6396  0.0659  0.7469 -0.6898  4.2241 -0.0277 
#>   Realised response against the desired gains:
#>     GY    PHT     AD    ASI    EPP    GLS 
#> 1.5623 0.4525 0.5510 0.1162 0.7475 0.2092 
#>   Objective 2.78141, chosen from 3 replicates.
#>   Replicate selection used: internal pre-fit holdout 
#>   The comparison observations were excluded before fitting.
#>   P provenance: user supplied
```

`dg` and `realised_response` use candidate standard-deviation units. The
coefficient solve still respects the original units of \\\mathbf{G}\\
and \\\mathbf{P}\\; the package performs the conversion internally.

The operational output is the ranked candidate table, not the
coefficient vector by itself.

``` r
head(dgsi$ranked_geno[, c(
  "GenoID", "SelectionIndex", "Rank", "Eligible", "Selected"
)])
#>     GenoID SelectionIndex  Rank Eligible Selected
#>     <char>          <num> <num>   <lgcl>   <lgcl>
#> 1: CAND051      -23.77901     1     TRUE     TRUE
#> 2: CAND092      -25.63016     2     TRUE     TRUE
#> 3: CAND193      -26.44297     3     TRUE     TRUE
#> 4: CAND115      -26.49760     4     TRUE     TRUE
#> 5: CAND135      -26.62914     5     TRUE     TRUE
#> 6: CAND130      -26.91844     6     TRUE     TRUE
dgsi$selected_geno[, c("GenoID", "SelectionIndex", "Rank")]
#>      GenoID SelectionIndex  Rank
#>      <char>          <num> <num>
#>  1: CAND051      -23.77901     1
#>  2: CAND092      -25.63016     2
#>  3: CAND193      -26.44297     3
#>  4: CAND115      -26.49760     4
#>  5: CAND135      -26.62914     5
#>  6: CAND130      -26.91844     6
#>  7: CAND160      -27.36982     7
#>  8: CAND005      -28.09595     8
#>  9: CAND004      -28.15358     9
#> 10: CAND158      -28.48562    10
#> 11: CAND146      -28.56305    11
#> 12: CAND087      -28.97873    12
#> 13: CAND137      -29.24266    13
#> 14: CAND103      -29.26149    14
#> 15: CAND007      -29.42558    15
#> 16: CAND075      -29.44933    16
#> 17: CAND065      -29.72210    17
#> 18: CAND199      -30.12385    18
#> 19: CAND036      -30.73888    19
#> 20: CAND071      -30.76035    20
#>      GenoID SelectionIndex  Rank
#>      <char>          <num> <num>
```

`Eligible` matters only when `select_mode = "eligible_top_n"`. In the
default `top_n` mode every complete candidate is eligible and the index
selects the requested number with deterministic tie-breaking.

------------------------------------------------------------------------

## 3. Separate the selected differential from transmitted response

``` r
round(rbind(
  requested = dgsi$desired_gain,
  optimised_direction = dgsi$optimised_d,
  selected_differential = dgsi$realised_response,
  expected_genetic_response_sd = dgsi$theoretical_response$standardised
), 3)
#>                                 GY   PHT    AD   ASI   EPP   GLS
#> requested                    1.000 0.400 0.600 0.500 0.400 0.600
#> optimised_direction          1.040 0.610 1.019 0.859 0.789 0.337
#> selected_differential        1.562 0.452 0.551 0.116 0.747 0.209
#> expected_genetic_response_sd 0.366 0.215 0.359 0.302 0.278 0.119
```

The selected differential describes the candidates retained. The
theoretical response is the model-based expectation transmitted to the
next generation. They are not synonyms. A large selected differential
with a small theoretical response indicates that much of the observed
separation is not expected to be inherited.

------------------------------------------------------------------------

## 4. Inspect optimisation stability

``` r
dgsi$replicate_diagnostics
#>    Replicate Objective Iteration_of_best Selected Plateau
#>        <int>     <num>             <int>    <int>  <lgcl>
#> 1:         1 1.4841447                60       20   FALSE
#> 2:         2 0.9569090                23       20   FALSE
#> 3:         3 0.3441077                33       20   FALSE
#>    Final_window_relative_improvement Chosen
#>                                <num> <lgcl>
#> 1:                         0.3077807   TRUE
#> 2:                         0.5536885  FALSE
#> 3:                         0.8395049  FALSE
dgsi$rank_correlation
#>           [,1]      [,2]      [,3]
#> [1,] 1.0000000 0.8608820 0.8692597
#> [2,] 0.8608820 1.0000000 0.9736353
#> [3,] 0.8692597 0.9736353 1.0000000
dgsi$selected_set_agreement
#>             replicate_1 replicate_2 replicate_3
#> replicate_1   1.0000000   0.5384615   0.4814815
#> replicate_2   0.5384615   1.0000000   0.5384615
#> replicate_3   0.4814815   0.5384615   1.0000000
```

The winning replicate is selected on held-out candidates by default.
Report:

1.  the training and holdout objectives;
2.  whether the search reached a plateau;
3.  coefficient and rank stability across replicates;
4.  selected-set agreement;
5.  the response for every trait; and
6.  the covariance source and its uncertainty.

If independent replicates produce different selected sets with similar
objectives, the data do not identify one operational recommendation.
Increase the evidence or present the stable candidate core rather than
hiding the instability.

------------------------------------------------------------------------

## 5. Compare with the non-iterated solution

``` r
comparison <- data.frame(
  Trait = traits,
  Requested = unname(dgsi$desired_gain[traits]),
  Non_iterated = unname(dgsi$non_iterated$realised_response[traits]),
  Iterative = unname(dgsi$realised_response[traits])
)
comparison$Absolute_error_non_iterated <- abs(
  comparison$Non_iterated - comparison$Requested
)
comparison$Absolute_error_iterative <- abs(
  comparison$Iterative - comparison$Requested
)
comparison[-1] <- round(comparison[-1], 3)
comparison
#>   Trait Requested Non_iterated Iterative Absolute_error_non_iterated
#> 1    GY       1.0        1.526     1.562                       0.526
#> 2   PHT       0.4        0.274     0.452                       0.126
#> 3    AD       0.6        0.441     0.551                       0.159
#> 4   ASI       0.5       -0.256     0.116                       0.756
#> 5   EPP       0.4        0.562     0.747                       0.162
#> 6   GLS       0.6        0.228     0.209                       0.372
#>   Absolute_error_iterative
#> 1                    0.562
#> 2                    0.052
#> 3                    0.049
#> 4                    0.384
#> 5                    0.347
#> 6                    0.391
```

Compute both objective values from the full candidate set so that the
comparison uses the same observations. The stored `dgsi$objective` can
instead be a holdout objective and must not be compared directly with a
full-sample quantity.

``` r
dgsi_loss <- function(response, target) {
  sum(((response - target) / pmax(abs(target), 0.25))^2)
}

c(
  non_iterated_full_sample = dgsi_loss(
    dgsi$non_iterated$realised_response, dgsi$desired_gain
  ),
  iterative_full_sample = dgsi_loss(
    dgsi$realised_response, dgsi$desired_gain
  )
)
#> non_iterated_full_sample    iterative_full_sample 
#>                 3.279795                 2.108015
```

In this reproducible run, the full-sample loss falls from about 3.280 to
2.108, a reduction of approximately 36%. That improvement is not uniform
across traits. The anthesis-silking interval changes from an
unfavourable differential to a favourable one, and plant height and
anthesis date move much closer to their requests. Grain yield, ears per
plant and grey leaf spot do not all become closer individually. DGSI
minimises the declared multi-trait loss; it does not guarantee a smaller
error for every trait. The trait table must therefore accompany the
total objective.

The iterative method has value only if the improvement is meaningful on
held-out candidates and stable across searches. A smaller training
objective alone is not sufficient evidence.

The internal holdout is used to choose among stochastic replicates. It
is not an independent breeding cycle and does not establish
transportability. When a later cohort is available, pass it through
`validation_data` and retain it outside coefficient fitting.

------------------------------------------------------------------------

## 6. What constitutes a sufficient DGSI demonstration

A defensible analysis reports all of the following:

1.  the requested gain vector and whether it represents an exact ratio
    or minimum floors;
2.  the exact-ray feasibility result at the planned selection intensity;
3.  the non-iterated Yamada response;
4.  the optimised direction, selected differential and expected
    transmitted response;
5.  the selected candidates and their index ranks;
6.  training and holdout objectives;
7.  rank, coefficient and selected-set stability across replicates; and
8.  validation in a later cohort when one is available.

The demonstration above now covers items 1–7 with the shipped
population. Item 8 necessarily requires data from a genuinely later
breeding cohort.

For objectives stated as intervals rather than one exact vector, read
[Defining desired gains and acceptable
intervals](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-desired-gain-intervals.md).

## References

Joukhadar R, Daetwyler HD, Bansal U, Gendall AR, Hayden MJ (2024). An
iterative desired-gain selection index for simultaneous improvement of
multiple traits.
