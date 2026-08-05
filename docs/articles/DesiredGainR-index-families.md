# Multiple-trait selection: methods, choice, and comparison

## 1. Learning objectives

Plant breeding programmes usually improve several traits at the same
time. Yield, quality, adaptation, maturity, and disease resistance may
all influence the release decision. Progress in one trait can also
reduce progress in another. Therefore, the method must represent both
the breeding objective and the genetic relationships among traits.

This vignette has four objectives. They are to: (i) define
multiple-trait selection, (ii) explain its main methods, (iii) show when
each method is appropriate, and (iv) provide a comprehensive comparison
analysis.

At the end, the breeder should be able to justify the selected method.
The justification should refer to the objective, the available
information, the expected responses, and the stability of the final
ranking.

## 2. What is multiple-trait selection?

Multiple-trait selection is the choice of candidates using two or more
traits. The traits can enter the decision at the same stage or at
different stages. The breeder must first decide how trade-offs are
allowed.

Four broad strategies are recognised by Beavis, Lamkey, Mahama, and Suza
(2023).

| Strategy | Definition | Suitable situation | Main limitation |
|----|----|----|----|
| Multistage selection | Different traits are evaluated at different stages. | Expensive traits are measured after inexpensive screening. | Early rejection can remove candidates that would perform well later. |
| Tandem selection | One trait is improved across one or more cycles before another trait becomes the focus. | One trait has temporary programme priority. | Other traits can deteriorate before the focus changes. |
| Independent culling | Every candidate must satisfy a limit for every required trait. | Release standards or biological limits are firm. | Excellence elsewhere gives no compensation for a failed limit. |
| Index selection | Traits are combined into one score of overall merit. | Trade-offs are allowed and simultaneous improvement is required. | The score is useful only when its objective and inputs are defensible. |

These strategies answer different questions. Hence, a comparison should
never treat them as interchangeable algorithms. The breeder’s decision
rule comes first. The statistical method follows.

### 2.1 Choose the strategy before choosing the index

No strategy is universally best. The suitable strategy follows from four
programme decisions.

1.  **Are any limits non-negotiable?** Use independent culling for those
    limits. An index can rank the candidates that pass.
2.  **Are traits measured at different stages or costs?** Use a
    multistage design. Place inexpensive or highly heritable screens
    early. Reserve costly traits for the reduced set.
3.  **Must one trait dominate several breeding cycles?** Use a true
    tandem schedule. Monitor correlated deterioration in the other
    traits before the focal trait changes.
4.  **Can progress in one trait compensate for less progress in
    another?** Use index selection. Then choose the index family from
    the biological objective.

These choices can be combined. For example, disease susceptibility can
define a firm culling limit. A desired-gain index can then rank the
eligible lines for yield, quality, and maturity. This hybrid rule
represents the programme more accurately than forcing every trait into
one score.

### 2.2 What DesiredGainR implements

The four strategies have different software coverage. The distinction
below is important.

| Strategy | Current package support | What can be compared directly |
|----|----|----|
| Multistage selection | No dedicated stage-allocation optimiser | Analyse the candidate set retained after each stage. Record attrition, cost, and final response separately. |
| Classical tandem selection across cycles | No dedicated alternating-cycle simulator | Define separate cycle scenarios outside the one-call comparison. Monitor response in every trait. |
| Within-cohort sequential screening | `selection_index(method = "tandem")` | Rankings, selected candidates, and observed differentials within one cohort. |
| Independent culling | `selection_index(method = "independent_culling")` | Passing rate, failed limits, selected candidates, and observed differentials. |
| Index selection | [`selection_index()`](https://FAkohoue.github.io/DesiredGainR/reference/selection_index.md), [`restricted_index()`](https://FAkohoue.github.io/DesiredGainR/reference/restricted_index.md), [`generalized_index()`](https://FAkohoue.github.io/DesiredGainR/reference/generalized_index.md), [`run_dgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_dgsi.md), and [`run_qgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_qgsi.md) | The available criteria depend on the fitted index family. They are defined in Section 6. |

The option named `method = "tandem"` is retained for compatibility. It
applies successive trait screens to one candidate cohort. Classical
tandem selection changes the focal trait across cycles or generations.
These procedures should not be interpreted as the same strategy.

### 2.3 Compare strategies on programme outcomes

Strategy comparison requires more than a common candidate score. Use the
same initial population, selection budget, time horizon, and release
constraints. Then compare the following outcomes.

| Outcome | Why it matters |
|----|----|
| Candidates retained after each decision | Shows where useful material is lost. |
| Cost per selected candidate | Gives multistage selection its operational value. |
| Expected transmitted response by trait | Measures the genetic consequence of selection. |
| Frequency of failed mandatory limits | Protects non-negotiable requirements. |
| Cycle time | Distinguishes rapid simultaneous selection from a tandem schedule. |
| Rank and selected-set agreement | Shows whether strategies change the immediate decision. |
| Multi-cycle response and genetic variance | Shows whether short-term gain persists. |
| Validation in a later cohort | Tests whether the rule transfers beyond the fitted candidates. |

Therefore, the preferred strategy is the simplest rule that respects
every firm constraint and gives the strongest validated response within
the programme’s cost and time limits. A larger index statistic alone
cannot choose among the four strategies.

## 3. The quantities used by an index

Let \\\mathbf{x}\\ contain the information recorded for a candidate. Let
\\\mathbf{g}\\ contain the additive genetic values of the objective
traits. The two vectors can contain the same traits. They can also
differ. For example, records from relatives, testcrosses, environments,
or correlated traits can help predict a smaller set of objective traits.

Three covariance matrices describe this problem.

| Symbol | Definition | Breeding meaning |
|----|----|----|
| \\\mathbf{P}=\operatorname{Var}(\mathbf{x})\\ | Covariance among information variables | Describes the variation and redundancy in the evidence used for selection. |
| \\\mathbf{C}=\operatorname{Cov}(\mathbf{x},\mathbf{g})\\ | Covariance between information and objective traits | Describes how strongly each information source predicts genetic value. |
| \\\mathbf{G}=\operatorname{Var}(\mathbf{g})\\ | Genetic covariance among objective traits | Describes available genetic variation and correlated response. |

The aggregate genetic merit is

\\ H=\mathbf{a}^{\mathsf T}\mathbf{g}, \\

where \\\mathbf{a}\\ contains economic weights. A linear selection index
is

\\ I=\mathbf{b}^{\mathsf T}\mathbf{x}, \\

where \\\mathbf{b}\\ contains index coefficients. Under truncation
selection, the expected response vector is

\\ \Delta =i\frac{\mathbf{C}^{\mathsf T}\mathbf{b}}
{\sqrt{\mathbf{b}^{\mathsf T}\mathbf{P}\mathbf{b}}}. \\

The term \\i\\ is the selection intensity. It increases when a smaller
proportion of candidates is selected. When information and objective
traits are identical, \\\mathbf{C}=\mathbf{G}\\. This gives the familiar
expression used throughout classical selection-index theory.

Use
[`selection_information()`](https://FAkohoue.github.io/DesiredGainR/reference/selection_information.md)
and
[`generalized_index()`](https://FAkohoue.github.io/DesiredGainR/reference/generalized_index.md)
when the information and objective traits differ. Read [Using
predictions, prediction errors, and restricted
responses](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-information.md)
for the complete formulation.

## 4. The methods available in DesiredGainR

The table below links each statistical method to its exact package call.
The value in the `method` column is supplied to the stated function.
This makes the mathematical name and the R interface explicit.

| Index family | Function and method value | Main input that defines the objective |
|----|----|----|
| Smith-Hazel | `selection_index(method = "smith_hazel")` | Economic weights |
| Base index | `selection_index(method = "base")` | Economic weights |
| Pesek-Baker | `selection_index(method = "pesek_baker")` | Desired gains |
| Yamada | `selection_index(method = "yamada")` | Desired gains |
| Mulamba-Mock rank sum | `selection_index(method = "mulamba_mock")` | Rank weights, or equal weights |
| Elston multiplicative index | `selection_index(method = "elston")` | Trait-specific culling limits |
| Independent culling | `selection_index(method = "independent_culling")` | Trait-specific culling limits |
| Within-cohort sequential screen | `selection_index(method = "tandem")` | Ordered trait screens in one cohort |
| Kempthorne-Nordskog | `restricted_index(method = "kempthorne_nordskog")` | Economic weights and zero-response traits |
| Restricted Smith-Hazel | `restricted_index(method = "restricted_smith_hazel")` | Economic weights and zero-response traits |
| Tallis | `restricted_index(method = "tallis")` | Economic weights and proportional target gains |
| Mallard | `restricted_index(method = "mallard")` | Economic weights and target gain amounts |
| Harville | `restricted_index(method = "harville")` | Economic weights, target gains, and a finite penalty |
| General economic index | `generalized_index(method = "economic")` | Economic weights and distinct information sources |
| General desired-gain index | `generalized_index(method = "desired_gain")` | Desired gains and distinct information sources |
| Desired-gain index with iterative search | [`run_dgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_dgsi.md) | Desired gains and candidate values |
| QGSI | [`run_qgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_qgsi.md) | Linear, squared, and cross-product economic values |

DGSI means desired-gain selection index. QGSI means quadratic genomic
selection index. The general index applies when the information
variables and the objective traits differ. For example, testcross
records may predict the additive genetic values of parental lines.

### 4.1 Economic-value methods

The Smith-Hazel index maximises expected response in aggregate genetic
merit. Its coefficients are

\\ \mathbf{b}=\mathbf{P}^{-1}\mathbf{G}\mathbf{a}. \\

Choose this method when the relative economic values are credible. It is
also the correct reference when every method is compared against one
fixed definition of merit.

The base index uses \\\mathbf{b}=\mathbf{a}\\. It applies the economic
weights directly. It provides a useful benchmark because it avoids
covariance inversion. A small advantage for Smith-Hazel over the base
index suggests that the estimated covariance structure contributes
little to the decision.

### 4.2 Desired-gain methods

The Pesek-Baker and Yamada formulations use a desired response direction
\\\mathbf{d}\\. Their coefficient expressions are

\\ \mathbf{b}=\mathbf{G}^{-1}\mathbf{d} \\

and

\\ \mathbf{b}=\mathbf{P}^{-1}\mathbf{G}
(\mathbf{G}\mathbf{P}^{-1}\mathbf{G})^{-1}\mathbf{d}. \\

They give the same coefficient direction when \\\mathbf{G}\\ is square
and invertible. Therefore, they are two computational routes to the same
linear index. A material difference between their answers indicates
numerical instability or a rank-deficient matrix.

Choose this family when the breeder can state relative desired gains
more confidently than economic weights. The vector fixes response
proportions. The selection intensity fixes the attainable magnitude.
Multiplying every desired gain by the same positive constant leaves the
ranking unchanged.

[`run_dgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_dgsi.md)
uses the established desired-gain selection index (DGSI) formula. It
adds the iterative search applied by Joukhadar et al. (2024). The search
changes the input desired-gain vector. It seeks a candidate ranking
whose selected differential is close to the breeder’s target. Joukhadar
et al. did not introduce a separate index family. Read [Optimising
desired
gains](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-dgsi.md)
before using the iterative search.

The closed-form Pesek-Baker and Yamada methods define a response
direction from covariance matrices. DGSI adds candidate-level iteration.
It searches for a ranking whose realised selected differential
approaches the requested direction. Fit both when candidate-level
attainment matters. Agreement gives a useful implementation check. A
difference directs attention to finite candidate supply, selection
intensity, and the shape of the candidate cloud.

### 4.3 Rank and threshold methods

The Mulamba-Mock rank-sum method ranks candidates within traits and
combines the ranks. It requires neither economic weights nor covariance
matrices when equal rank weights are used. It is valuable as a simple
comparator. Its score does not provide a closed-form predicted genetic
response in DesiredGainR.

Independent culling applies firm acceptance limits. Supply each limit in
the original trait units. Use a minimum for traits where larger values
are favourable. Use a maximum for traits listed in `lower_is_better`.

The Elston multiplicative index also applies firm limits. It then ranks
the eligible candidates using their joint margins above those limits. It
is useful when balanced superiority among eligible candidates matters.

`selection_index(method = "tandem")` applies ordered trait screens
within one candidate cohort. It is useful as an operational comparator.
It does not represent classical tandem selection across breeding cycles.

### 4.4 Restricted methods

A restricted index maximises merit while controlling response in
specified traits. DesiredGainR provides five forms.

| Method | Breeding question |
|----|----|
| Kempthorne-Nordskog | Which index holds specified traits at zero expected change? |
| Restricted Smith-Hazel | Which projection form gives the same zero-response restriction? |
| Tallis | Which index preserves stated response proportions? |
| Mallard | Which index targets stated response amounts? |
| Harville | Which index balances merit against a finite restriction penalty? |

Restrictions usually reduce response in unconstrained net merit. That
loss is the cost of respecting the biological or commercial condition.
Report both the restricted response and the loss relative to the
unrestricted index.

Satoh (2024) gives a geometric interpretation of proportional
restriction. The scalar \\\beta\\ measures progress along the desired
direction. The Mahalanobis residual measures departure from that
direction after accounting for genetic covariance. DesiredGainR reports
both quantities through
[`evaluate_restricted_response()`](https://FAkohoue.github.io/DesiredGainR/reference/evaluate_restricted_response.md).

All five restricted forms use the same interface. The restriction
determines which form has a clear interpretation. Use
Kempthorne-Nordskog or restricted Smith-Hazel for zero expected change.
Use Tallis for a response direction. Use Mallard for stated response
amounts. Use Harville when the programme accepts a measured departure
from the target in exchange for greater aggregate merit.

### 4.5 Non-linear economic value

The quadratic genomic selection index (QGSI) represents linear, squared,
and cross-product economic values. Choose it when biological or economic
merit is genuinely curved. A desired gain is a response target. It is a
different quantity from a quadratic weight. Read [Quadratic genomic
selection](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-qgsi.md)
for the full demonstration.

[`run_qgsi_desired_gain()`](https://FAkohoue.github.io/DesiredGainR/reference/run_qgsi.md)
remains as a compatibility entry point. It raises an informative error
because desired gains and quadratic economic values define different
objectives. Use
[`run_dgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_dgsi.md)
for desired gains. Use
[`run_qgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_qgsi.md)
for curved economic merit.

## 5. How to choose a strategy and an index family

Begin with the programme structure. Then choose the simplest index
family that represents the allowed trade-offs.

| Breeder’s statement | Starting method | Required evidence |
|----|----|----|
| “I can defend relative economic values.” | Smith-Hazel | \\\mathbf{G}\\, \\\mathbf{P}\\, and economic weights |
| “I can state the response direction.” | Pesek-Baker or Yamada | \\\mathbf{G}\\, and \\\mathbf{P}\\ for Yamada |
| “I can state acceptable gain intervals.” | [`suggest_desired_gains()`](https://FAkohoue.github.io/DesiredGainR/reference/suggest_desired_gains.md) | \\\mathbf{G}\\, \\\mathbf{P}\\, intervals, and candidates |
| “Every candidate must satisfy these limits.” | Independent culling | Candidate values and limits |
| “Eligible candidates should have balanced margins.” | Elston | Candidate values and limits |
| “Some traits should remain unchanged.” | Restricted index | \\\mathbf{G}\\, \\\mathbf{P}\\, merit weights, and restrictions |
| “Different traits become available at defined stages.” | Multistage selection | Measurement cost, reliability, stage capacity, and attrition |
| “One trait must dominate successive cycles.” | Classical tandem schedule | Cycle schedule and correlated response in every trait |
| “One cohort must pass ordered trait screens.” | `selection_index(method = "tandem")` | Screen order and candidate values |
| “I need a simple comparator with few assumptions.” | Base or Mulamba-Mock | Weights for the base index, or candidate ranks |
| “Merit contains curvature or trait interactions.” | QGSI | Genomic estimated breeding values and quadratic economic values |

Use the following decision sequence.

1.  Define the objective traits and their favourable directions.
2.  Identify firm limits that allow no compensation.
3.  Record when each trait becomes available and what it costs to
    measure.
4.  Decide whether selection is simultaneous, multistage, or spread
    across cycles.
5.  For simultaneous selection, define an economic, desired-gain,
    restricted, rank-based, threshold, or quadratic objective.
6.  Confirm the covariance matrices and units required by that
    objective.
7.  Fit one simple benchmark and one objective-matched method.
8.  Compare biological response before comparing scalar criteria.
9.  Examine ranking stability and uncertainty.
10. Use validation and simulation when the decision extends beyond one
    cohort.

## 6. What can be compared and how

### 6.1 Exact scope of the comparison functions

[`compare_selection_methods()`](https://FAkohoue.github.io/DesiredGainR/reference/compare_selection_methods.md)
accepts fitted objects from
[`selection_index()`](https://FAkohoue.github.io/DesiredGainR/reference/selection_index.md),
[`restricted_index()`](https://FAkohoue.github.io/DesiredGainR/reference/restricted_index.md),
[`generalized_index()`](https://FAkohoue.github.io/DesiredGainR/reference/generalized_index.md),
[`run_dgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_dgsi.md),
and
[`run_qgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_qgsi.md).
Internal adapters harmonise candidates, favourable trait directions and
response units. They compare responses, rankings and selected sets. They
never compare coefficients across families.

Use
[`comparison_objective()`](https://FAkohoue.github.io/DesiredGainR/reference/comparison_objective.md)
for a cross-family analysis. It fixes one desired response target, one
aggregate-weight vector, one optional quadratic utility, and one genetic
covariance matrix. Every fitted method is then evaluated against that
same objective.

[`compare_dg_and_qgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/compare_dg_and_qgsi.md)
accepts one
[`run_dgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_dgsi.md)
result and one
[`run_qgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_qgsi.md)
result. It compares candidate scores, ranks, and selected sets. Its
`decision_summary` reports selected counts, overlap and Jaccard
similarity. The comparison is descriptive. A desired response objective
and a quadratic economic objective are different biological statements.
Rank agreement cannot decide which statement is appropriate.

A QGSI with non-zero curvature retains its documented linear-regression
approximation. The special case with \\\mathbf{W}=0\\ is a linear
genomic index and has the exact linear response. A quadratic common
utility is evaluated exactly from common validation values. This
distinction appears in `Expected_response_basis` and must remain visible
when methods are ranked.

### 6.2 Conditions for a fair comparison

A comprehensive comparison requires common conditions. Otherwise, a
numerical difference can reflect the setup rather than the method.

Use the same:

1.  candidates and trait records;
2.  trait directions and units;
3.  one common genetic covariance for the objective traits;
4.  selected proportion or selection intensity;
5.  aggregate-merit definition when comparing \\R\_{HI}\\ or \\\Delta
    H\\;
6.  training, validation, and test partitions;
7.  simulation scenarios, replicate seeds, and computational budget.

The comparison also needs several outcomes. No single statistic is
sufficient.

| Outcome | Question answered |
|----|----|
| Expected response for every trait | What genetic change is predicted? |
| Target attainment | How much of each requested gain is reached? |
| Feasibility | Can the complete target be reached at the planned intensity? |
| \\R\_{HI}\\ | How closely does the index rank one fixed aggregate merit? |
| \\\Delta H\\ | How much change in that fixed merit is expected? |
| Relative efficiency (RE) | How much response in the declared main trait remains relative to direct selection? |
| Index heritability | How repeatable is the composite score under the fitted covariance model? |
| Mahalanobis alignment | Does the response follow the requested multivariate direction? |
| Selected-set overlap | Do two methods choose the same candidates? |
| Rank correlation | Do two methods order all candidates similarly? |
| Selection probability | Does candidate choice remain stable after prediction uncertainty is propagated? |
| Information-deletion efficiency | Which records contribute useful selection information? |
| Multi-cycle response and diversity | Does the apparent advantage persist across breeding cycles? |

The coefficient of variation of the index, \\CV_I\\, is retained for
reproduction of published analyses. It depends on the score origin and
becomes undefined after centring. Therefore, it should never decide the
preferred method.

### 6.3 Evaluate each family with the criteria it can support

| Family or strategy | Primary criteria | Additional evidence |
|----|----|----|
| Smith-Hazel and base | Expected response by trait, \\R\_{HI}\\, \\\Delta H\\, relative efficiency, and index accuracy | Weight sensitivity, covariance uncertainty, and validation |
| Closed-form desired-gain index | Expected response by trait, feasibility, target attainment, Satoh’s \\\beta\\, and Mahalanobis alignment | Covariance uncertainty and selected-set stability |
| Desired-gain index with iterative search | The same desired-gain criteria, plus selected differential, holdout objective, replicate stability, and comparison with the closed-form starting index | Independent-cohort validation |
| Restricted index | Expected response, restriction residual, loss in common aggregate merit, and Satoh projection where relevant | Sensitivity to the restriction and penalty |
| Mulamba-Mock rank sum | Rank stability, selected-set overlap, and observed differential | Later-cohort validation or simulation |
| Independent culling | Passing rate, failed limits, selected count, and observed differential | Sensitivity to each threshold |
| Elston multiplicative index | Passing rate, margins among survivors, selected set, and observed differential | Threshold sensitivity and validation |
| Within-cohort sequential screen | Attrition at each screen, selected set, and observed differential | Screen-order sensitivity and operational cost |
| Generalised index | Expected response in the objective traits, accuracy for the economic form, and information-deletion efficiency | Prediction-error propagation and validation |
| Quadratic genomic selection index | Linear and quadratic score contributions, observed genomic estimated breeding value differential, approximate expected gains, and selection agreement with a linear benchmark | Out-of-sample utility, true-value accuracy in simulation, and multi-cycle response |

The absence of one criterion is not evidence of poor quality. For
example, independent culling has no linear-index accuracy. Its purpose
is to enforce limits. Conversely, a high \\R\_{HI}\\ cannot rescue a
method that violates a mandatory disease threshold.

## 7. A comprehensive comparison analysis

### 7.1 State the objective

The example seeks gains in grain yield (GY), ears per plant (EPP), and
lower values for plant height (PHT), anthesis date (AD),
anthesis-silking interval (ASI), and grey leaf spot score (GLS). The
target is expressed in genetic standard deviations.

``` r
desired_gains_sd <- c(
  GY = 1.0,
  PHT = 0.4,
  AD = 0.6,
  ASI = 0.5,
  EPP = 0.4,
  GLS = 0.6
)
```

First test feasibility. This separates an unrealistic target from a poor
method.

``` r
feasibility <- gain_feasibility(
  desired_gains_sd,
  dgr_G,
  dgr_P,
  n_candidates = nrow(dgr_candidates),
  n_select = 20L,
  lower_is_better = lower_is_better,
  gain_units = "genetic_sd"
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

### 7.2 Put every method on one explicit scale

The following transformation expresses candidates and covariance
matrices in genetic-standard-deviation units.
[`selection_index()`](https://FAkohoue.github.io/DesiredGainR/reference/selection_index.md)
then receives `scale_traits = FALSE`. Hence, the desired-gain vector and
every reported response use the same units.

``` r
genetic_sd <- sqrt(diag(dgr_G))[traits]
inverse_sd <- diag(1 / genetic_sd, nrow = length(traits))
dimnames(inverse_sd) <- list(traits, traits)

G_sd <- inverse_sd %*% dgr_G %*% inverse_sd
P_sd <- inverse_sd %*% dgr_P %*% inverse_sd

candidates_sd <- dgr_candidates
centred <- sweep(
  as.matrix(dgr_candidates[, traits]),
  2L,
  colMeans(dgr_candidates[, traits]),
  "-"
)
candidates_sd[, traits] <- sweep(centred, 2L, genetic_sd, "/")
```

The economic weights below are implied by the desired-gain direction.
This creates one common merit definition for the mathematical
comparison.

``` r
common_weights <- implied_economic_weights(
  desired_gains_sd,
  G_sd,
  P_sd,
  lower_is_better = lower_is_better,
  gain_units = "trait"
)
round(common_weights, 3)
#>     GY    PHT     AD    ASI    EPP    GLS 
#> 10.887  0.362  5.808 -7.928 -1.515  0.576 
#> attr(,"provenance")
#> [1] "Implied by the supplied desired gains through w = G^-1 P G^-1 d; not an independently estimated economic value. Expressed in the favourable-direction space, so a positive weight favours movement in the breeder-defined direction."
```

Negative implied weights can occur. They are mathematically valid under
correlated response. The desired response remains the biological
statement to interpret.

### 7.3 Fit objective-matched methods and benchmarks

``` r
common_arguments <- list(
  values = candidates_sd,
  trait_cols = traits,
  id_col = "GenoID",
  G = G_sd,
  P = P_sd,
  lower_is_better = lower_is_better,
  center_traits = FALSE,
  scale_traits = FALSE,
  n_select = 20L,
  main_trait = "GY"
)

smith_hazel <- do.call(selection_index, c(
  common_arguments,
  list(method = "smith_hazel", economic_weights = common_weights)
))

pesek_baker <- do.call(selection_index, c(
  common_arguments,
  list(
    method = "pesek_baker",
    desired_gains = desired_gains_sd,
    aggregate_weights = common_weights
  )
))

base_index <- do.call(selection_index, c(
  common_arguments,
  list(method = "base", economic_weights = common_weights)
))

rank_sum <- do.call(selection_index, c(
  common_arguments,
  list(method = "mulamba_mock")
))
```

The first two methods should agree because the economic weights were
derived from the desired gains. This is a mathematical verification. The
base and rank-sum methods provide simpler benchmarks.

Threshold methods require limits rather than a response vector. These
limits use the units of `candidates_sd`. A maximum is supplied for each
trait listed in `lower_is_better`.

``` r
limits_sd <- c(
  GY = -0.50,
  PHT = 0.75,
  AD = 0.75,
  ASI = 0.75,
  EPP = -0.75,
  GLS = 0.75
)

culling <- do.call(selection_index, c(
  common_arguments,
  list(
    method = "independent_culling",
    culling_thresholds = limits_sd
  )
))

elston <- do.call(selection_index, c(
  common_arguments,
  list(method = "elston", culling_thresholds = limits_sd)
))

sequential_screen <- do.call(selection_index, c(
  common_arguments,
  list(method = "tandem", tandem_order = c("GY", "GLS", "AD"))
))
```

### 7.4 Assemble the comparison

``` r
comparison <- compare_selection_methods(
  list(
    Smith_Hazel = smith_hazel,
    Pesek_Baker = pesek_baker,
    Base = base_index,
    Rank_sum = rank_sum,
    Culling = culling,
    Elston = elston,
    Sequential_screen = sequential_screen
  ),
  target_gains = desired_gains_sd
)
comparison
#> <desiredgainr_method_comparison>
#>   Methods: Smith_Hazel, Pesek_Baker, Base, Rank_sum, Culling, Elston, Sequential_screen 
#>   Response units: trait 
#>   All recorded comparison conditions are satisfied.
#>             Method              Family N_selected Common_merit_response
#>        Smith_Hazel         smith_hazel         20                    NA
#>        Pesek_Baker         pesek_baker         20                    NA
#>               Base                base         20                    NA
#>           Rank_sum        mulamba_mock         20                    NA
#>            Culling independent_culling         20                    NA
#>             Elston              elston         20                    NA
#>  Sequential_screen              tandem         20                    NA
#>  Worst_expected_attainment Mahalanobis_alignment Validation_utility_response
#>                  0.5470366               1.00000                          NA
#>                  0.5470366               1.00000                          NA
#>                 -0.1915515               0.83361                          NA
#>                         NA                    NA                          NA
#>                         NA                    NA                          NA
#>                         NA                    NA                          NA
#>                         NA                    NA                          NA
```

Read the fairness table first.

``` r
comparison$fairness
#>                                             Condition Satisfied
#> 1                Same candidates and objective traits      TRUE
#> 2                    Same favourable trait directions      TRUE
#> 3                               Common response units      TRUE
#> 4                                Same number selected      TRUE
#> 5                            Same selection intensity      TRUE
#> 6 One common genetic covariance for response geometry      TRUE
#> 7             Common genetic covariance is invertible      TRUE
#> 8        Validation identifiers aligned when supplied      TRUE
#>                                                                        Interpretation
#> 1                                             Required and enforced by this function.
#> 2                                             Required and enforced by this function.
#> 3                                              Responses are reported in trait units.
#> 4                            Hard culling can retain fewer candidates than requested.
#> 5                        Model-based responses require one common selection pressure.
#> 6                Alignment uses fitted genetic covariance ; disagreement is reported.
#> 7 Mahalanobis alignment and Satoh projection require an invertible covariance matrix.
#> 8 No validation data were supplied. External validation evidence remains unavailable.
```

A culling method may select fewer than 20 candidates. That result is
useful. It shows that the limits and the requested selection quota
conflict.

### 7.5 Compare biological response

``` r
comparison$responses[
  , c(
    "Method", "Trait", "Expected_response",
    "Observed_differential", "Expected_attainment"
  )
]
#>                Method  Trait Expected_response Observed_differential
#>                <char> <char>             <num>                 <num>
#>  1:       Smith_Hazel     GY        0.54703663            2.71839267
#>  2:       Smith_Hazel    PHT        0.21881465            0.51969887
#>  3:       Smith_Hazel     AD        0.32822198            0.64456380
#>  4:       Smith_Hazel    ASI        0.27351832           -0.30315333
#>  5:       Smith_Hazel    EPP        0.21881465            1.05538000
#>  6:       Smith_Hazel    GLS        0.32822198            0.69842471
#>  7:       Pesek_Baker     GY        0.54703663            2.71839267
#>  8:       Pesek_Baker    PHT        0.21881465            0.51969887
#>  9:       Pesek_Baker     AD        0.32822198            0.64456380
#> 10:       Pesek_Baker    ASI        0.27351832           -0.30315333
#> 11:       Pesek_Baker    EPP        0.21881465            1.05538000
#> 12:       Pesek_Baker    GLS        0.32822198            0.69842471
#> 13:              Base     GY        0.35864060            2.68933267
#> 14:              Base    PHT        0.04922368            0.30407054
#> 15:              Base     AD        0.09127695            0.62602380
#> 16:              Base    ASI       -0.09577576           -0.58177000
#> 17:              Base    EPP        0.06607095            0.77743000
#> 18:              Base    GLS        0.09854367            0.13686588
#> 19:          Rank_sum     GY                NA            1.91381933
#> 20:          Rank_sum    PHT                NA            0.66151429
#> 21:          Rank_sum     AD                NA            0.80376780
#> 22:          Rank_sum    ASI                NA            1.14847167
#> 23:          Rank_sum    EPP                NA            1.24483000
#> 24:          Rank_sum    GLS                NA            1.10650706
#> 25:           Culling     GY                NA            1.48661267
#> 26:           Culling    PHT                NA            0.59683721
#> 27:           Culling     AD                NA            0.60975580
#> 28:           Culling    ASI                NA            1.07222167
#> 29:           Culling    EPP                NA            0.60343000
#> 30:           Culling    GLS                NA            0.62050118
#> 31:            Elston     GY                NA            1.37939933
#> 32:            Elston    PHT                NA            0.69501637
#> 33:            Elston     AD                NA            0.50001980
#> 34:            Elston    ASI                NA            1.07734667
#> 35:            Elston    EPP                NA            0.67248000
#> 36:            Elston    GLS                NA            0.84341294
#> 37: Sequential_screen     GY                NA            1.66673933
#> 38: Sequential_screen    PHT                NA           -0.03648821
#> 39: Sequential_screen     AD                NA            0.91060980
#> 40: Sequential_screen    ASI                NA            0.13530500
#> 41: Sequential_screen    EPP                NA            0.55378000
#> 42: Sequential_screen    GLS                NA            1.55563059
#>                Method  Trait Expected_response Observed_differential
#>                <char> <char>             <num>                 <num>
#>     Expected_attainment
#>                   <num>
#>  1:           0.5470366
#>  2:           0.5470366
#>  3:           0.5470366
#>  4:           0.5470366
#>  5:           0.5470366
#>  6:           0.5470366
#>  7:           0.5470366
#>  8:           0.5470366
#>  9:           0.5470366
#> 10:           0.5470366
#> 11:           0.5470366
#> 12:           0.5470366
#> 13:           0.3586406
#> 14:           0.1230592
#> 15:           0.1521282
#> 16:          -0.1915515
#> 17:           0.1651774
#> 18:           0.1642395
#> 19:                  NA
#> 20:                  NA
#> 21:                  NA
#> 22:                  NA
#> 23:                  NA
#> 24:                  NA
#> 25:                  NA
#> 26:                  NA
#> 27:                  NA
#> 28:                  NA
#> 29:                  NA
#> 30:                  NA
#> 31:                  NA
#> 32:                  NA
#> 33:                  NA
#> 34:                  NA
#> 35:                  NA
#> 36:                  NA
#> 37:                  NA
#> 38:                  NA
#> 39:                  NA
#> 40:                  NA
#> 41:                  NA
#> 42:                  NA
#>     Expected_attainment
#>                   <num>
```

Expected response and observed differential answer different questions.
Expected response predicts transmitted change under the covariance
model. Observed differential describes the candidates selected from this
data set. Rank, threshold, and sequential-screen methods retain the
second quantity because a closed-form expected response is unavailable
here.

Target attainment equals expected response divided by the requested
gain. A value of 1 reaches the target for that trait. A value of 0.70
reaches 70%. Negative values move in an unfavourable direction. Read the
worst trait before the mean because one severe failure can be hidden by
strong gains elsewhere.

### 7.6 Compare summary criteria

``` r
comparison$summary[
  , c(
    "Method", "N_selected", "R_HI", "Delta_H", "RE",
    "Index_heritability", "Worst_expected_attainment",
    "Mahalanobis_alignment", "Mahalanobis_residual"
  )
]
#>               Method N_selected      R_HI  Delta_H        RE Index_heritability
#>               <char>      <int>     <num>    <num>     <num>              <num>
#> 1:       Smith_Hazel         20 0.4117023 5.630275 0.5268772          0.2439162
#> 2:       Pesek_Baker         20 0.4117023 5.630275 0.5268772          0.2439162
#> 3:              Base         20 0.3779319 5.168444 0.3454239          0.1428325
#> 4:          Rank_sum         20        NA       NA        NA                 NA
#> 5:           Culling         20        NA       NA        NA                 NA
#> 6:            Elston         20        NA       NA        NA                 NA
#> 7: Sequential_screen         20        NA       NA        NA                 NA
#>    Worst_expected_attainment Mahalanobis_alignment Mahalanobis_residual
#>                        <num>                 <num>                <num>
#> 1:                 0.5470366               1.00000         2.087841e-15
#> 2:                 0.5470366               1.00000         1.064591e-15
#> 3:                -0.1915515               0.83361         3.663562e-01
#> 4:                        NA                    NA                   NA
#> 5:                        NA                    NA                   NA
#> 6:                        NA                    NA                   NA
#> 7:                        NA                    NA                   NA
```

Use each quantity for its defined purpose.

1.  Read `R_HI` and `Delta_H` only where the same aggregate merit was
    supplied.
2.  Read `RE` as efficiency for GY alone.
3.  Read index heritability as repeatability of the composite score.
4.  Read worst attainment as protection against hidden trait failure.
5.  Read Mahalanobis alignment as agreement with the complete desired
    direction.
6.  Read the residual as departure from the requested proportions.

### 7.7 Compare rankings and selected candidates

``` r
round(comparison$rank_correlation, 2)
#>                   Smith_Hazel Pesek_Baker Base Rank_sum Culling Elston
#> Smith_Hazel              1.00        1.00 0.92     0.63    0.31   0.31
#> Pesek_Baker              1.00        1.00 0.92     0.63    0.31   0.31
#> Base                     0.92        0.92 1.00     0.33    0.17   0.17
#> Rank_sum                 0.63        0.63 0.33     1.00    0.43   0.44
#> Culling                  0.31        0.31 0.17     0.43    1.00   1.00
#> Elston                   0.31        0.31 0.17     0.44    1.00   1.00
#> Sequential_screen        0.41        0.41 0.35     0.40    0.24   0.24
#>                   Sequential_screen
#> Smith_Hazel                    0.41
#> Pesek_Baker                    0.41
#> Base                           0.35
#> Rank_sum                       0.40
#> Culling                        0.24
#> Elston                         0.24
#> Sequential_screen              1.00
round(comparison$selected_jaccard, 2)
#>                   Smith_Hazel Pesek_Baker Base Rank_sum Culling Elston
#> Smith_Hazel              1.00        1.00 0.74     0.33    0.18   0.18
#> Pesek_Baker              1.00        1.00 0.74     0.33    0.18   0.18
#> Base                     0.74        0.74 1.00     0.21    0.11   0.11
#> Rank_sum                 0.33        0.33 0.21     1.00    0.21   0.25
#> Culling                  0.18        0.18 0.11     0.21    1.00   0.74
#> Elston                   0.18        0.18 0.11     0.25    0.74   1.00
#> Sequential_screen        0.25        0.25 0.21     0.33    0.21   0.21
#>                   Sequential_screen
#> Smith_Hazel                    0.25
#> Pesek_Baker                    0.25
#> Base                           0.21
#> Rank_sum                       0.33
#> Culling                        0.21
#> Elston                         0.21
#> Sequential_screen              1.00
comparison$selected_overlap
#>                   Smith_Hazel Pesek_Baker Base Rank_sum Culling Elston
#> Smith_Hazel                20          20   17       10       6      6
#> Pesek_Baker                20          20   17       10       6      6
#> Base                       17          17   20        7       4      4
#> Rank_sum                   10          10    7       20       7      8
#> Culling                     6           6    4        7      20     17
#> Elston                      6           6    4        8      17     20
#> Sequential_screen           8           8    7       10       7      7
#>                   Sequential_screen
#> Smith_Hazel                       8
#> Pesek_Baker                       8
#> Base                              7
#> Rank_sum                         10
#> Culling                           7
#> Elston                            7
#> Sequential_screen                20
```

Rank correlation uses every candidate. Jaccard similarity uses the
selected sets. A high rank correlation can coexist with a different
parent list near the selection boundary. Therefore, report both.

### 7.8 Compare the closed-form and iteratively optimised desired-gain index

The iterative procedure remains part of the desired-gain family.
Therefore, compare it first with its closed-form starting index. Use the
same candidates, covariance matrices, trait directions, and selected
count.

[`run_dgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_dgsi.md)
expresses its target in candidate standard-deviation units. The
closed-form example above uses genetic standard deviations. The
conversion must be explicit.

``` r
candidate_sd <- vapply(
  dgr_candidates[, traits], stats::sd, numeric(1L)
)
target_original_units <- desired_gains_sd * genetic_sd
target_candidate_sd <- target_original_units / candidate_sd

iterative_dgsi <- run_dgsi(
  init_data = dgr_candidates["GenoID"],
  cand_data = dgr_candidates,
  trait_cols = traits,
  dg = target_candidate_sd,
  G = dgr_G,
  P = dgr_P,
  lower_is_better = lower_is_better,
  n_select = 20L,
  n_iter = 500L,
  n_rep = 10L,
  seed = 2026L
)
```

Construct one objective and include the iterative result in the same
call.

``` r
comparison_target <- comparison_objective(
  desired_gains = desired_gains_sd,
  aggregate_weights = economic_weights,
  G = dgr_G,
  gain_units = "genetic_sd"
)

cross_family <- compare_selection_methods(
  list(
    Pesek_Baker = pesek_baker,
    Iterative_DGSI = iterative_dgsi
  ),
  objective = comparison_target,
  validation_data = dgr_candidates
)
```

The direct calculations below show the components returned by the
unified comparison.

``` r
closed_rank <- setNames(
  pesek_baker$ranking$rank,
  pesek_baker$ranking$id
)
iterative_rank <- setNames(
  iterative_dgsi$ranked_geno$Rank,
  iterative_dgsi$ranked_geno$GenoID
)

rank_correlation <- cor(
  closed_rank[names(closed_rank)],
  iterative_rank[names(closed_rank)],
  method = "spearman"
)
selected_overlap <- length(intersect(
  pesek_baker$selected$id,
  iterative_dgsi$selected_geno$GenoID
))

response_genetic_sd <- rbind(
  closed_form = pesek_baker$evaluation$expected_response,
  iterative =
    iterative_dgsi$theoretical_response$original_units / genetic_sd
)
attainment <- sweep(
  response_genetic_sd,
  2L,
  desired_gains_sd,
  "/"
)
```

Read four results: (i) transmitted response, (ii) target attainment,
(iii) rank correlation, and (iv) selected overlap. Then inspect the
iterative holdout objective and replicate stability. Improvement on the
fitted candidate set alone is insufficient.

### 7.9 Compare a linear and quadratic genomic index

A quadratic genomic selection index (QGSI) should have a linear genomic
selection index (LGSI) benchmark. Use the same genomic estimated
breeding values (GEBVs), linear economic weights, genomic covariance,
transformation, and selected count. Set the quadratic matrix to zero for
the linear benchmark.

``` r
q_traits <- c("GY", "PHT", "GLS")
linear_values <- c(GY = 1.0, PHT = 0.2, GLS = 0.4)
Gamma <- stats::cov(dgr_gebv[, q_traits])

W_zero <- matrix(
  0, length(q_traits), length(q_traits),
  dimnames = list(q_traits, q_traits)
)
W_curved <- W_zero
diag(W_curved) <- c(-0.05, -0.02, -0.04)

lgsi <- run_qgsi(
  init_data = dgr_gebv["GenoID"],
  gebv_data = dgr_gebv,
  trait_cols = q_traits,
  linear_weights = linear_values,
  W = W_zero,
  Gamma = Gamma,
  lower_is_better = c("PHT", "GLS"),
  n_select = 20L
)
qgsi <- run_qgsi(
  init_data = dgr_gebv["GenoID"],
  gebv_data = dgr_gebv,
  trait_cols = q_traits,
  linear_weights = linear_values,
  W = W_curved,
  Gamma = Gamma,
  lower_is_better = c("PHT", "GLS"),
  n_select = 20L
)

linear_rank <- setNames(lgsi$ranked_geno$Rank, lgsi$ranked_geno$GenoID)
quadratic_rank <- setNames(qgsi$ranked_geno$Rank, qgsi$ranked_geno$GenoID)
cor(linear_rank, quadratic_rank, method = "spearman")
length(intersect(
  lgsi$selected_geno$GenoID,
  qgsi$selected_geno$GenoID
))
lgsi$expected_gain_per_trait
qgsi$expected_gain_per_trait
qgsi$component_summary
```

The numerical weights above illustrate the comparison mechanics. They
are not a breeding recommendation. A quadratic analysis requires
biological or economic justification for every squared and cross-product
term.

Evaluate the two fits using: (i) per-trait response, (ii) rank and
selected-set agreement, (iii) candidate-specific linear and quadratic
contributions, and (iv) out-of-sample nonlinear utility. In simulation,
also compare accuracy against true genetic merit and multi-cycle
response. The linear model is preferred when curvature does not improve
validated utility or change the decision reliably.

### 7.10 Compare every fitted family under one objective

The final comparison uses one declared yardstick. `common_G` is the
genetic covariance in original trait coordinates. `validation_values`
contains one row per fitted candidate and the same objective traits.
These values can be independent predictions, simulated true genetic
values, or later measurements. Their provenance determines the strength
of the evidence.

``` r
common_objective <- comparison_objective(
  desired_gains = desired_gains_sd,
  aggregate_weights = common_aggregate_weights,
  W = common_quadratic_utility,
  G = common_G,
  gain_units = "genetic_sd"
)

all_methods <- compare_selection_methods(
  list(
    Smith_Hazel = smith_hazel,
    Restricted = restricted_fit,
    General = general_fit,
    Desired_gain = pesek_baker,
    Iterative_DGSI = iterative_dgsi,
    Linear_genomic = lgsi,
    Quadratic_genomic = qgsi
  ),
  objective = common_objective,
  validation_data = validation_values
)

all_methods$fairness
all_methods$responses
all_methods$validation_responses
all_methods$validation_utility
all_methods$summary
all_methods$rank_correlation
all_methods$selected_jaccard
```

The `R_HI`, `Delta_H`, and `Index_accuracy` columns retain native family
criteria. They can describe different merits. `Common_merit_response` is
the response in the linear component of the declared common utility. Use
it for a cross-family conclusion about that component.
`Validation_utility_response` evaluates the complete utility on the
supplied validation values. A non-zero quadratic term requires these
values because mean trait responses do not determine the mean of a
quadratic function.

## 8. Add uncertainty and validation

The point comparison above is the beginning of the analysis. A breeding
recommendation also needs stability evidence.

### 8.1 Covariance and weight uncertainty

Use
[`index_uncertainty()`](https://FAkohoue.github.io/DesiredGainR/reference/index_uncertainty.md)
to resample the genetic and phenotypic covariance matrices. Use
[`weight_sensitivity()`](https://FAkohoue.github.io/DesiredGainR/reference/weight_sensitivity.md)
when economic weights are uncertain. Report coefficient intervals,
response intervals, rank stability, and selected-set stability.

### 8.2 Candidate prediction uncertainty

Genomic estimated breeding values (GEBVs) have prediction error. When
candidate-specific prediction error covariance (PEC) is available, use
[`candidate_score_uncertainty()`](https://FAkohoue.github.io/DesiredGainR/reference/candidate_score_uncertainty.md).
Report each candidate’s score interval and selection probability. Give
particular attention to candidates near the selection boundary.

### 8.3 Value of each information source

Use
[`index_information_efficiency()`](https://FAkohoue.github.io/DesiredGainR/reference/index_information_efficiency.md)
for a Smith-Hazel or general economic index. It applies Cunningham’s
deletion-efficiency principle. A trait with little information value can
then be removed from future phenotyping plans if the saving is
operationally important.

### 8.4 Out-of-sample validation

Fit the index with training data. Fix its coefficients and
transformations. Apply
[`predict()`](https://rdrr.io/r/stats/predict.html) to a later cohort or
a held-out set. Compare predicted responses with realised breeding
values, progeny means, or later-stage trial performance. This separates
genuine prediction from reuse of the fitting data.

### 8.5 Multi-cycle validation

Use
[`simulate_selection_cycles()`](https://FAkohoue.github.io/DesiredGainR/reference/simulate_selection_cycles.md)
when the decision concerns several cycles. Compare cumulative response,
joint target attainment, genetic variance, and diversity. Use matched
scenarios and replicate seeds. Read [Multi-cycle
simulation](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-simulation.md)
for calibration and stress testing.

## 9. Choosing the preferred method

A method is preferred when five conditions are satisfied.

1.  It represents the breeder’s actual objective.
2.  Its required inputs have adequate quality.
3.  Its expected response respects every critical trait.
4.  Its ranking remains stable under credible uncertainty.
5.  Its advantage persists in validation or relevant multi-cycle
    scenarios.

The method with the largest single statistic can fail this standard. For
example, a high aggregate-merit response can hide failure of a disease
limit. A stable rank can also follow the wrong biological direction.
Therefore, the final recommendation should report the complete decision
panel.

Close the analysis with a practical statement. Name the selected method.
State why it matches the objective. Report the expected response for
every trait. State the main uncertainty. Then describe how the decision
will be validated in the breeding programme.

## 10. Package coverage and regression tests

DesiredGainR maintains an explicit contract for its public interface.
The contract currently covers 43 exported functions. Every export must
have a manual alias and at least one behavioural test file. Adding an
exported function without both entries causes the test suite to fail.

The family-level tests cover the complete method set described in this
vignette. They fit six classical index families and two operational
comparators on one labelled data set. They also fit all five restricted
families and both general-index objectives. The tests check the method
returned, the selected candidates, trait responses, restriction
residuals, target attainment, rank agreement, and selected-set
agreement. Separate tests cover DGSI, QGSI, uncertainty, feasibility,
simulation, covariance handling, and data validation.

Contributors should run the following commands after changing code or
documentation.

``` r
devtools::document()
devtools::test()
devtools::check()
```

[`devtools::document()`](https://devtools.r-lib.org/reference/document.html)
rebuilds the manual files from the documentation beside each R function.
[`devtools::test()`](https://devtools.r-lib.org/reference/test.html)
checks scientific and programming contracts.
[`devtools::check()`](https://devtools.r-lib.org/reference/check.html)
checks the complete source package in a clean build environment.

## 11. References

- Beavis W, Lamkey K, Mahama AA, Suza W (2023). Multiple Trait
  Selection. In Suza WP and Lamkey KR, editors, *Quantitative Genetics
  for Plant Breeding*. Iowa State University Digital Press.
- Cunningham EP (1969). The relative efficiencies of selection indexes.
  *Acta Agriculturae Scandinavica* **19**:45-48.
- Elston RC (1963). A weight-free index for the purpose of ranking or
  selection with respect to several traits at a time. *Biometrics*
  **19**:85-97.
- Kempthorne O, Nordskog AW (1959). Restricted selection indices.
  *Biometrics* **15**:10-19.
- Mulamba NN, Mock JJ (1978). Improvement of yield potential of the Eto
  Blanco maize population by breeding for plant traits. *Egyptian
  Journal of Genetics and Cytology* **7**:40-51.
- Pesek J, Baker RJ (1969). Desired improvement in relation to selection
  indices. *Canadian Journal of Plant Science* **49**:803-804.
- Satoh M (2024). Characteristics of restricted selection indices and
  geometrical interpretation of restricted breeding values. *Journal of
  Animal Breeding and Genetics* **141**:353-363.
- Smith HF (1936). A discriminant function for plant selection. *Annals
  of Eugenics* **7**:240-250.
- Yamada Y, Yokouchi K, Nishida A (1975). Selection index when genetic
  gains of individual traits are of primary concern. *Japanese Journal
  of Genetics* **50**:33-41.
