# Defining a breeding objective

## 1. Why this layer exists

An index answers the question “given this objective, which candidates?”.
Nothing in index theory answers “what should the objective be?”, and
that is the question breeders find hard.

Guimarães et al. (2021) quantified the cost of getting it wrong. Working
in a lowland rice recurrent-selection population, they applied
Smith-Hazel and Tai indices under three weight schemes: genetic standard
deviations, a variation-based index, and an arbitrary vector. Under the
arbitrary weights the Smith-Hazel index produced *no* gain for grain
yield, leaf scald or grain discoloration, and unfavourable change in
three of the six traits. The index machinery was identical in all three
runs. Only the objective differed.

Covarrubias-Pazaran (2021) reaches the same place from the operational
side. The CGIAR guideline instructs breeders to state a desired
differential and then adds, in the same paragraph, that the resulting
coefficients should not be interpreted, because they lack meaning under
strong genetic correlations. The desired response is the only decision
genuinely available.

This vignette covers the six tools DesiredGainR provides for that
decision.

------------------------------------------------------------------------

## 2. The duality: two ways of saying one thing

### 2.1 The result

The Smith-Hazel index sets `b = P⁻¹Gw` from economic weights `w`. The
desired-gain index sets `b = P⁻¹G(GP⁻¹G)⁻¹d` from desired gains `d`.
Setting the two equal gives

\\w = G^{-1}PG^{-1}d, \qquad d = GP^{-1}Gw.\\

So they are not competing methods. They are two coordinate systems on
one index, and every objective stated in either can be read in the
other.

### 2.2 Why it matters practically

A breeder who cannot put a currency value on grey leaf spot resistance
can usually say how much of it they want relative to yield. The duality
turns that into weights.

The reverse direction is the more useful one in practice, and is
under-exploited. When a programme proposes weights, showing them the
response those weights actually request often reveals that the objective
is not what anyone intended.

``` r
desired_gains <- c(
  GY = 1.0, PHT = 0.4, AD = 0.6, ASI = 0.5,
  EPP = 0.4, GLS = 0.6
)

implied <- implied_economic_weights(
  desired_gains, dgr_G, dgr_P,
  lower_is_better = lower_is_better, gain_units = "genetic_sd"
)
round(implied, 3)
#>      GY     PHT      AD     ASI     EPP     GLS 
#>  14.516   0.030   2.323 -13.213 -15.150   0.678 
#> attr(,"provenance")
#> [1] "Implied by the supplied desired gains through w = G^-1 P G^-1 d; not an independently estimated economic value. Expressed in the favourable-direction space, so a positive weight always means the trait matters."
```

The translation is exactly invertible when the units match, and checking
that is worth the one line it costs:

``` r
round(
  implied_desired_gains(
    implied, dgr_G, dgr_P,
    lower_is_better = lower_is_better, gain_units = "genetic_sd"
  ),
  3
)
#>  GY PHT  AD ASI EPP GLS 
#> 1.0 0.4 0.6 0.5 0.4 0.6 
#> attr(,"provenance")
#> [1] "Implied by the supplied economic weights through d = G P^-1 G w, expressed in genetic_sd units and up to the scalar set by selection intensity."
```

### 2.3 Two ways to misread the output

Both of the following mislead readers who have not met them before, and
both appear in the vector above.

**Negative weights arise for traits the objective asks to improve.**
Here the anthesis-silking interval and ears per plant both receive
negative weights, although the objective requests gains in both. This is
correct. Once oriented, grain yield correlates +0.55 with the
anthesis-silking interval and +0.45 with ears per plant, so selecting
for a full standard deviation of yield would carry both *past* the
smaller gains requested for them. Delivering the stated ratio therefore
requires the index to hold them back.

A negative weight never means the trait should decrease.
`lower_is_better` has already applied direction, so larger is better
throughout.

**Magnitudes are not comparable across traits.** A weight carries
inverse trait units, so a trait on a small scale receives a large number
for the same emphasis. Rescale before comparing:

``` r
genetic_sd <- stats::setNames(dgr_traits$genetic_sd, dgr_traits$trait)
round(sort(implied * genetic_sd[names(implied)], decreasing = TRUE), 2)
#>    GY    AD   GLS   PHT   EPP   ASI 
#> 10.89  5.81  0.58  0.36 -1.51 -7.93
```

Ears per plant carried the largest raw number in the vector and is one
of the smaller effects. Reading the unscaled vector would send a breeder
to argue about the wrong trait.

------------------------------------------------------------------------

## 3. Feasibility: can the objective be attained at all?

### 3.1 The achievable-response ellipsoid

For any linear index the response vector is `R = i·Gb/√(bᵀPb)`. Writing
`u = b/√(bᵀPb)` so that `uᵀPu = 1`, the attainable responses are
`{i·Gu}`, which is the ellipsoid

\\R^\mathsf{T}G^{-1}PG^{-1}R = i^{2}.\\

A stated target `d` therefore requires exactly

\\i\_{\mathrm{required}} = \sqrt{d^\mathsf{T}G^{-1}PG^{-1}d},\\

and the selected proportion delivering that intensity follows from the
normal-truncation relationship.

### 3.2 Two consequences

**The classical desired-gain index honours only the direction of `d`.**
Multiplying every desired gain by a constant leaves the index unchanged,
because the attainable magnitude is fixed by selection intensity.
Joukhadar et al. (2024) note this in passing when their INDEX1 and
INDEX3, targeting 0.5 and 4.0 standard deviations respectively, produced
identical rankings without iteration. Breeders who set absolute targets
and then read back a realised gain are being misled unless the point is
made explicit.

**An antagonistic correlation structure can make a modest-looking target
unreachable.** This is not a marginal effect.

``` r
feasibility <- gain_feasibility(
  desired_gains, dgr_G, dgr_P,
  n_candidates = nrow(dgr_candidates), n_select = 20L,
  lower_is_better = lower_is_better, gain_units = "genetic_sd"
)
feasibility
#> <desiredgainr_feasibility>
#>   Required selection intensity: 3.2082
#>   Requires the top 0.1773%, which is fewer than one of 200 candidates
#>   Planned intensity: 1.7550 (top 10.0%)
#>   Feasible at planned intensity: no
#>   Feasible anywhere in this population: no
#>   Attainable fraction of the requested gain: 54.7%
```

The request looks unremarkable — between 0.4 and 1.0 genetic standard
deviations across six traits — and requires selecting fewer than one
candidate from two hundred.

``` r
round(feasibility$attainable_response_input_units, 3)
#>    GY   PHT    AD   ASI   EPP   GLS 
#> 0.547 0.219 0.328 0.274 0.219 0.328
```

Fifty-five per cent of the request, in the requested proportion, is what
the planned intensity delivers.

### 3.3 Where this check belongs

The CGIAR guideline directs breeders to separate software to determine
whether “the progress of 1 or more σ is possible given your current
trait genetic correlations and selection intensity”. That check belongs
in the same place as the index, and it belongs *before* the optimisation
rather than after it disappoints.

------------------------------------------------------------------------

## 4. Recovering an objective from past decisions

### 4.1 The method

Where a programme has selected for years without a formal index, the
differentials it achieved already encode its objective.
Covarrubias-Pazaran recovers the corresponding linear rule as

\\b = P^{-1}s,\\

where `s` holds the differentials between the selected group and the
full population. This is how an index was introduced at both the
International Maize and Wheat Improvement Center (CIMMYT) and the
International Rice Research Institute (IRRI).

``` r
recovered <- retrospective_weights(
  selected_values = dgr_candidates[dgr_history$selected, traits],
  population_values = dgr_candidates[, traits],
  trait_cols = traits
)
recovered
#> <desiredgainr_retrospective>
#>   40 selected from 200 candidates (20.0%)
#>   P: estimated from population_values 
#>   Recovered coefficients, per trait standard deviation:
#>      GY     PHT      AD     ASI     EPP     GLS 
#>  0.9324 -0.3000 -0.3986 -0.5244  0.1508 -0.3975 
#>   (Raw coefficients are in $coefficients; they carry inverse trait units and are not comparable across traits.)
```

### 4.2 Checking the recovery

The example data record the weights that generated the historical
decision, so the recovery can be verified rather than asserted.
Comparing on a common scale, oriented and normalised to grain yield:

``` r
direction <- ifelse(dgr_traits$direction == "increase", 1, -1)
names(direction) <- traits
oriented <- recovered$coefficients_per_sd * direction
generating <- attr(dgr_history, "generating_weights")

data.frame(
  trait = traits,
  recovered = round(oriented / oriented[["GY"]], 2),
  generating = round(generating / generating[["GY"]], 2)
)
#>     trait recovered generating
#> GY     GY      1.00       1.00
#> PHT   PHT      0.32       0.25
#> AD     AD      0.43       0.60
#> ASI   ASI      0.56       0.45
#> EPP   EPP      0.16       0.30
#> GLS   GLS      0.43       0.55
```

The recovery captures the shape of the objective without reproducing it
exactly, which is the honest result. `b = P⁻¹s` returns the linear rule
that best reproduces the *observed differentials*, and those
differentials are themselves one realisation of a selection decision in
a finite population.

### 4.3 What it is for, and what it is not for

A retrospective weight vector answers “what linear rule best
approximates the decisions already made?”. It does not answer “what
should this programme select for next?”.

The distinction is operational, not pedantic. The IRRI deployment
reports retrospective coefficients for four regional programmes in which
the *same* resistance-gene indicator carries opposite signs, and the
fine-tuning step that followed — doubling the Philippines yield
coefficient — raised the reported final gain from 2.230 to 3.772.
Recovered weights reproduce past behaviour, including any bias in it.
Treat them as a starting point, inspect them, and adjust deliberately
against the target product profile.

------------------------------------------------------------------------

## 5. How much does the objective matter?

### 5.1 The question

Weights are estimates, and refining them costs meetings. Before spending
that effort it is worth knowing whether the decision would change.

``` r
economic_weights <- c(
  GY = 1.0, PHT = 0.2, AD = 0.5,
  ASI = 0.4, EPP = 0.3, GLS = 0.5
)
sensitivity <- weight_sensitivity(
  economic_weights, dgr_candidates, dgr_G, dgr_P,
  n_select = 20L, trait_cols = traits, n_draws = 100L
)
sensitivity
#> <desiredgainr_sensitivity>
#>   100 draws, log-scale perturbation SD 0.25, selecting 20
#>   Median selected-set agreement: 0.905
#>   Decisions reproducing the original set (agreement >= 0.90): 56.0%
#>   Median rank correlation with the stated objective: 0.994
#>   Weights the decision is most sensitive to:
#>   PHT    GY   ASI 
#> 0.516 0.475 0.046 
#>   (Sensitivity of the decision, not share of the index; see effective_weights().)
```

Weights are perturbed multiplicatively on the log scale, the index is
rebuilt for each draw, and the resulting selected sets are compared with
the set obtained from the stated objective.

### 5.2 Reading the result

| Observation | Interpretation |
|----|----|
| High stability, high rank correlation | The decision is robust; further argument about weights will not change it |
| Low stability, high rank correlation | The ranking is stable but the threshold is contested; consider selecting more candidates |
| Low stability, low rank correlation | The objective, not the index, is the binding uncertainty |

### 5.3 Influence is not contribution

`weight_influence` measures how strongly perturbing each weight disturbs
the selected set. It is a different quantity from the share of the index
a trait contributes, and the two routinely disagree.

In the output above, plant height is the weight the decision is most
sensitive to, while contributing about six per cent of the index:

``` r
smith_hazel <- selection_index(
  dgr_candidates, traits,
  method = "smith_hazel",
  G = dgr_G, P = dgr_P, economic_weights = economic_weights,
  lower_is_better = lower_is_better, n_select = 20L, main_trait = "GY"
)
smith_hazel$effective_weights[, c("Trait", "Genetic_share")]
#>     Trait Genetic_share
#>    <char>         <num>
#> 1:     GY    0.20205605
#> 2:    PHT    0.06339965
#> 3:     AD    0.23638957
#> 4:    ASI    0.15198589
#> 5:    EPP    0.12788741
#> 6:    GLS    0.21828143
```

There is no contradiction. A trait carrying little of the index can
still be the largest lever on the decision, because it moves the index
in a direction the remaining traits do not cover, and therefore reorders
the candidates near the selection threshold. A trait that dominates the
index may conversely be a weak lever, because the ranking already
follows it.

Read the two together. The first says which weights are worth arguing
about; the second says which traits the index is acting on.

------------------------------------------------------------------------

## 6. The Crosbie trap

Crosbie et al. (1980), as reported by Guimarães et al. (2021), observed
that assigning equal weights across unstandardised traits places most of
the selection pressure on whichever trait carries the largest genetic
variance. This is the single most common way for an objective to be
quietly subverted, and the CGIAR guideline’s recommendation to start
from a desired differential of one for every trait makes it easy to walk
into.

The example programme is built to show it. Trait standard deviations
span a factor of about one hundred and twenty:

``` r
round(sort(genetic_sd, decreasing = TRUE), 2)
#>   PHT    AD   GLS    GY   ASI   EPP 
#> 12.00  2.50  0.85  0.75  0.60  0.10
```

The consequence is visible in the conditioning of the covariance
matrices. Most of it comes from the scales rather than from the
correlations:

``` r
c(
  covariance = matrix_diagnostics(dgr_G)$condition_number,
  correlation = matrix_diagnostics(stats::cov2cor(dgr_G))$condition_number
)
#>  covariance correlation 
#> 18569.91997     8.75443
```

Standardising removes three orders of magnitude of ill-conditioning.
[`run_dgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_dgsi.md)
and
[`run_qgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_qgsi.md)
therefore warn when trait standard deviations differ by more than
fivefold and standardisation has been switched off, naming the trait
that dominates and its share.

------------------------------------------------------------------------

## 7. Summary of the tools

| Function | Question it answers |
|----|----|
| [`implied_economic_weights()`](https://FAkohoue.github.io/DesiredGainR/reference/implied_economic_weights.md) | What weights do these desired gains imply? |
| [`implied_desired_gains()`](https://FAkohoue.github.io/DesiredGainR/reference/implied_desired_gains.md) | What response do these weights actually request? |
| [`gain_feasibility()`](https://FAkohoue.github.io/DesiredGainR/reference/gain_feasibility.md) | Is this objective attainable, and at what intensity? |
| [`retrospective_weights()`](https://FAkohoue.github.io/DesiredGainR/reference/retrospective_weights.md) | What has this programme been selecting for? |
| [`weight_sensitivity()`](https://FAkohoue.github.io/DesiredGainR/reference/weight_sensitivity.md) | Would a different objective change the decision? |
| [`effective_weights()`](https://FAkohoue.github.io/DesiredGainR/reference/effective_weights.md) | Which traits is the index really acting on? |
| [`matrix_diagnostics()`](https://FAkohoue.github.io/DesiredGainR/reference/matrix_diagnostics.md) | Can these covariance matrices be inverted safely? |

------------------------------------------------------------------------

## 8. What is not here

Two routes to an objective are outside the package’s current scope, and
both are legitimate.

**A profit function.** The only principled definition of an economic
weight is the derivative of profit with respect to the trait mean. Where
a programme can write down prices, discounts and costs, differentiating
that function is better than any elicitation. DesiredGainR does not yet
provide a helper for it.

**Structured elicitation.** Marginal rates of substitution are weight
ratios directly, so a short sequence of trade-off questions pins the
vector up to scale, with a consistency check. This is standard decision
analysis and would fit the package well.

Both are recorded as intended additions rather than omissions of
principle.

------------------------------------------------------------------------

## 9. References

- Covarrubias-Pazaran G (2021). *Practical implementation of selection
  indices.* CGIAR Excellence in Breeding.
- Crosbie TM, Mock JJ, Smith OS (1980), as discussed in Guimarães et al.
  (2021).
- Guimarães PHR, Melo PGS, Cordeiro ACC, Torga PP, Rangel PHN, de Castro
  AP (2021). Index selection can improve the selection efficiency in a
  rice recurrent selection population. *Euphytica* **217**:95.
  <https://doi.org/10.1007/s10681-021-02819-7>
- Joukhadar R, Li Y, Thistlethwaite R, Forrest KL, Tibbits JF, Trethowan
  R, Hayden MJ (2024). Optimising desired gain indices to maximise
  selection response. *Frontiers in Plant Science* **15**:1337388.
  <https://doi.org/10.3389/fpls.2024.1337388>
- Pesek J, Baker RJ (1969). Desired improvement in relation to selection
  indices. *Canadian Journal of Plant Science* **49**:803-804.
  <https://doi.org/10.4141/cjps69-137>
- Smith HF (1936). A discriminant function for plant selection. *Annals
  of Eugenics* **7**:240-250.
