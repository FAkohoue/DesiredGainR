# Iterative optimisation of the desired-gain index

## 1. What the iterative method changes

The classical desired-gain selection index (DGSI) fixes a response
direction under the covariance model. Joukhadar et al. (2024) retained
that established index formula and added an iterative search over nearby
desired-gain vectors. The search evaluates the selected candidates
directly. Its purpose is to reduce disagreement between the requested
standardised gain vector and the differential of the selected group.

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

The feasibility target uses genetic-standard-deviation units.
[`run_dgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_dgsi.md)
uses candidate-standard-deviation units. Convert the same biological
target before fitting.

``` r
genetic_sd <- sqrt(diag(dgr_G))[traits]
candidate_sd <- vapply(
  dgr_candidates[, traits], stats::sd, numeric(1L)
)
target_original_units <- desired_gains * genetic_sd
dg_target <- target_original_units / candidate_sd
round(dg_target, 3)
#>    GY   PHT    AD   ASI   EPP   GLS 
#> 0.540 0.327 0.499 0.286 0.241 0.426
```

If the exact ray is unattainable, optimisation can find a closer
selected set. It cannot create unavailable genetic variation. If the
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
  dg = dg_target,
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
#>       GY      PHT       AD      ASI      EPP      GLS 
#>  13.5687  -0.1837   4.0918  -3.2739 -38.2414   0.1771 
#>   Realised response against the desired gains:
#>      GY     PHT      AD     ASI     EPP     GLS 
#>  1.2522  0.0620  0.9686 -0.0270  0.2308  0.0663 
#>   Objective 8.69789, chosen from 3 replicates.
#>   Replicate selection used: internal pre-fit holdout 
#>   The comparison observations were excluded before fitting.
#>   P provenance: user supplied
```

`dg` and `realised_response` use candidate standard-deviation units. The
coefficient solve still respects the original units of \\\mathbf{G}\\
and \\\mathbf{P}\\. The package performs the conversion internally.

The operational output is the ranked candidate table, not the
coefficient vector by itself.

``` r
head(dgsi$ranked_geno[, c(
  "GenoID", "SelectionIndex", "Rank", "Eligible", "Selected"
)])
#>     GenoID SelectionIndex  Rank Eligible Selected
#>     <char>          <num> <num>   <lgcl>   <lgcl>
#> 1: CAND130      -133.2926     1     TRUE     TRUE
#> 2: CAND135      -151.7086     2     TRUE     TRUE
#> 3: CAND092      -152.4793     3     TRUE     TRUE
#> 4: CAND005      -153.6863     4     TRUE     TRUE
#> 5: CAND051      -155.1229     5     TRUE     TRUE
#> 6: CAND087      -156.8744     6     TRUE     TRUE
dgsi$selected_geno[, c("GenoID", "SelectionIndex", "Rank")]
#>      GenoID SelectionIndex  Rank
#>      <char>          <num> <num>
#>  1: CAND130      -133.2926     1
#>  2: CAND135      -151.7086     2
#>  3: CAND092      -152.4793     3
#>  4: CAND005      -153.6863     4
#>  5: CAND051      -155.1229     5
#>  6: CAND087      -156.8744     6
#>  7: CAND004      -158.6038     7
#>  8: CAND150      -159.4416     8
#>  9: CAND147      -159.9039     9
#> 10: CAND158      -160.0289    10
#> 11: CAND160      -161.5619    11
#> 12: CAND065      -161.6259    12
#> 13: CAND146      -162.4855    13
#> 14: CAND115      -162.8150    14
#> 15: CAND075      -166.1193    15
#> 16: CAND153      -166.5070    16
#> 17: CAND074      -166.7513    17
#> 18: CAND175      -166.9187    18
#> 19: CAND054      -167.3839    19
#> 20: CAND186      -167.4782    20
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
  requested_candidate_sd = dgsi$desired_gain,
  optimised_direction = dgsi$optimised_d,
  selected_differential = dgsi$realised_response,
  expected_response_candidate_sd = dgsi$theoretical_response$standardised
), 3)
#>                                   GY    PHT    AD    ASI   EPP   GLS
#> requested_candidate_sd         0.540  0.327 0.499  0.286 0.241 0.426
#> optimised_direction            2.651 -0.351 5.142  2.738 0.008 0.967
#> selected_differential          1.252  0.062 0.969 -0.027 0.231 0.066
#> expected_response_candidate_sd 0.250 -0.033 0.484  0.258 0.001 0.091
```

The selected differential describes the candidates retained. The
theoretical response is the model-based expectation transmitted to the
next generation. They are not synonyms. A large selected differential
with a small theoretical response indicates that much of the observed
separation is not expected to be inherited.

To compare transmitted response with the original feasibility target,
return to genetic-standard-deviation units.

``` r
transmitted_genetic_sd <-
  dgsi$theoretical_response$original_units / genetic_sd
round(rbind(
  requested_genetic_sd = desired_gains,
  transmitted_genetic_sd = transmitted_genetic_sd,
  attainment = transmitted_genetic_sd / desired_gains
), 3)
#>                           GY   PHT     AD    ASI   EPP    GLS
#> requested_genetic_sd   1.000 0.400  0.600  0.500 0.400  0.600
#> transmitted_genetic_sd 0.462 0.040 -0.582 -0.451 0.001 -0.128
#> attainment             0.462 0.101 -0.970 -0.902 0.004 -0.214
```

------------------------------------------------------------------------

## 4. Inspect optimisation stability

``` r
dgsi$replicate_diagnostics
#>    Replicate Objective Iteration_of_best Selected Plateau
#>        <int>     <num>             <int>    <int>  <lgcl>
#> 1:         1  4.763623                48       20   FALSE
#> 2:         2  3.481674                30       20   FALSE
#> 3:         3  3.370560                23       20   FALSE
#>    Final_window_relative_improvement Chosen
#>                                <num> <lgcl>
#> 1:                       0.005523655  FALSE
#> 2:                       0.273149398  FALSE
#> 3:                       0.296346133   TRUE
dgsi$rank_correlation
#>           [,1]      [,2]      [,3]
#> [1,] 1.0000000 0.9535618 0.8320378
#> [2,] 0.9535618 1.0000000 0.8824196
#> [3,] 0.8320378 0.8824196 1.0000000
dgsi$selected_set_agreement
#>             replicate_1 replicate_2 replicate_3
#> replicate_1   1.0000000   0.6000000   0.3333333
#> replicate_2   0.6000000   1.0000000   0.4814815
#> replicate_3   0.3333333   0.4814815   1.0000000
```

The winning replicate is selected on held-out candidates by default.
Report:

1.  The training and holdout objectives.
2.  Whether the search reached a plateau.
3.  Coefficient and rank stability across replicates.
4.  Selected-set agreement.
5.  The response for every trait.
6.  The covariance source and its uncertainty.

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
#> 1    GY     0.540        1.468     1.252                       0.928
#> 2   PHT     0.327        0.425     0.062                       0.098
#> 3    AD     0.499        0.536     0.969                       0.037
#> 4   ASI     0.286       -0.173    -0.027                       0.459
#> 5   EPP     0.241        0.636     0.231                       0.395
#> 6   GLS     0.426        0.495     0.066                       0.070
#>   Absolute_error_iterative
#> 1                    0.712
#> 2                    0.265
#> 3                    0.469
#> 4                    0.313
#> 5                    0.010
#> 6                    0.359
```

Compute both objective values from the full candidate set so that the
comparison uses the same observations. The stored `dgsi$objective` can
instead be a holdout objective. Compare it only with other holdout
quantities.

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
#>                 8.150162                 5.194453
```

The two values show whether iteration reduced the declared full-sample
loss. The improvement can differ among traits. Therefore, always inspect
the trait-specific errors as well as the total. The iterative procedure
minimises the declared multi-trait loss. It does not guarantee a smaller
error for every trait.

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

1.  The requested gain vector and whether it represents an exact ratio
    or minimum floors.
2.  The exact-ray feasibility result at the planned selection intensity.
3.  The non-iterated Yamada response.
4.  The optimised direction, selected differential and expected
    transmitted response.
5.  The selected candidates and their index ranks.
6.  Training and holdout objectives.
7.  Rank, coefficient and selected-set stability across replicates.
8.  Validation in a later cohort when one is available.

The demonstration above now covers items 1–7 with the shipped
population. Item 8 necessarily requires data from a genuinely later
breeding cohort.

For objectives stated as intervals rather than one exact vector, read
[Defining desired gains and acceptable
intervals](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-desired-gain-intervals.md).

## References

Joukhadar R, Li Y, Thistlethwaite R, Forrest KL, Tibbits JF, Trethowan
R, Hayden MJ (2024). Optimising desired gain indices to maximise
selection response. *Frontiers in Plant Science* **15**:1337388.
<https://doi.org/10.3389/fpls.2024.1337388>
