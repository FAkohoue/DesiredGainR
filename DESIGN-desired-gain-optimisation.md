# DesiredGainR — design for economic-weight guidance and simulation-optimised desired gains

Status: design proposal, no code written. Date 2026-07-30.
Target: 0.4.0 (exact tools), 0.5.0 (simulation engine), 0.6.0 (optimiser).

Decisions taken: AlphaSimR is the ground-truth simulation engine; all four
ranking modes are in scope; all three mating systems (self-pollinated,
cross-pollinated, clonal) are in scope.

---

## 0. The two results everything else hangs off

Both are exact, cost nothing, and need no simulation. They should ship first.

**Duality between economic weights and desired gains.**
Smith-Hazel is `b = P⁻¹Gw`; Pesek-Baker is `b = P⁻¹G(GP⁻¹G)⁻¹d`. The two
coincide when

```
w = G⁻¹ P G⁻¹ d          (weights implied by desired gains)
d = G P⁻¹ G w            (desired gains implied by weights)
```

**Achievable-gain ellipsoid.** For any index, `R = i·Gb/√(b'Pb)`, so every
attainable response vector satisfies `R'G⁻¹PG⁻¹R = i²`. Pesek-Baker always
lands on it, at

```
R = i · d / √(d' G⁻¹ P G⁻¹ d)
i_required = √(d' G⁻¹ P G⁻¹ d) = √(d' w_implied)
```

Two consequences to state prominently in the documentation:

1. Classical Pesek-Baker honours only the **direction** of `d`. Scaling every
   desired gain by a constant changes nothing. The current `run_dgsi()` API
   invites the opposite belief.
2. Absolute targets have an exact feasibility test. `i_required` converts to a
   selected proportion; if that proportion is below `1/N`, the target is
   unreachable this cycle regardless of index choice.

---

## 1. Objective 1 — guiding economic weights

Four mechanisms, none of which invent economics. Each returns provenance in the
style the package already uses.

### 1.1 `implied_economic_weights()` / `implied_desired_gains()`

The duality above. The reverse direction is the more valuable one in practice:
a breeder proposes weights, and the package reports the gain vector those
weights actually request. Most weight-specification errors become visible
immediately.

### 1.2 `economic_weights_from_profit()`

Accept a user-supplied profit or gross-margin function of trait means and
differentiate numerically at the current population means:
`w_j = ∂π/∂μ_j`. This is the only principled definition of an economic weight.
Return the gradient, the evaluation point, a finite-difference convergence
check, and a warning when the function is materially non-linear over the
plausible response range (in which case a single weight vector is a local
approximation and should be re-derived each cycle).

### 1.3 `elicit_economic_weights()`

Marginal rates of substitution are weight ratios directly: the answer to "how
much yield would you give up for one unit less disease?" is `w_dis/w_yld`.
`p−1` trade-offs pin the vector up to scale. Over-determine with additional
comparisons, check transitivity, and report an inconsistency index. Also
support swing weighting (rank traits by the value of moving each from its worst
to best plausible level, then ratio-scale).

Output must record which questions were asked and what was answered, so the
weight vector is auditable rather than asserted.

### 1.4 `weight_sensitivity()`

Weights are estimates. Report:

- Correlation between the index under `w` and under perturbed `w`.
- The region of weight space over which the selected set is unchanged —
  usually large, which is reassuring, and when it is small that is the finding.
- Which pairwise weight ratios the decision is actually sensitive to. Most are
  not, and knowing which few matter focuses the breeder's effort.

### 1.5 Scale discipline

Always report `w_j` alongside `w_j · σ_gj` — importance per unit of *available
genetic variation*. Raw-scale weights are uninterpretable across traits and are
where most specification errors originate.

---

## 2. Objective 2 — simulation-optimised desired gains

### 2.1 Why AlphaSimR must be the objective, not a validator

The infinitesimal deterministic recursion (Bulmer variance reduction plus
drift) is fast but assumes many loci of small effect, no epistasis, and LD
dynamics that are only approximately right. Where those assumptions fail —
major QTL, strong dominance in clonal crops, small `N_e`, tight LD — the
optimal desired-gain direction can differ from what the recursion predicts.

Therefore the deterministic model must **never gate, filter, or pre-screen**
candidate directions. If it did, the answer would be the deterministic model's
answer with an AlphaSimR rubber stamp, and a true optimum outside the
recursion's blind spot would be discarded before AlphaSimR ever saw it.

### 2.2 Architecture: surrogate-assisted search with AlphaSimR as ground truth

```
                 ┌──────────────────────────────┐
   d-direction   │  AlphaSimR cycle simulation  │   noisy multi-trait
   (unit sphere) │  (ground truth objective)    │   cumulative gain
       │         └──────────────┬───────────────┘        │
       │                        │                        │
       │                 evaluations (d, R)               │
       │                        ▼                        │
       │         ┌──────────────────────────────┐        │
       └─────────┤  Gaussian-process surrogate  │◄───────┘
       proposals │  prior mean = deterministic  │
                 │  recursion (optional)        │
                 └──────────────┬───────────────┘
                                │ acquisition
                                ▼
                        next d to evaluate
```

The deterministic recursion enters **only as the GP's prior mean function**. It
accelerates early exploration when there is little data, and is progressively
overridden as AlphaSimR evaluations accumulate. Specifically:

- It cannot permanently exclude any region: the acquisition function retains an
  exploration term over the whole sphere.
- Its influence is a fitted scalar, so a systematically wrong prior is
  down-weighted automatically.
- `surrogate_prior = "none"` disables it entirely, giving a constant-mean GP,
  for users who do not want it in the loop at all.
- The package reports the discrepancy between prior prediction and AlphaSimR
  outcome at every evaluated point. A systematic gap is itself a finding — it
  says the infinitesimal model does not describe this population, which is
  worth telling the breeder.

This is the only structure I can see that gets a tractable search without
subordinating AlphaSimR to the cheap model.

### 2.3 Search space

Only the direction of `d` matters, so the domain is the unit sphere in `p`
dimensions — `p−1` free parameters. For the 3–8 traits breeders use in
practice this is small enough for Sobol initialisation plus GP refinement.

Normalising `d` also removes the magnitude confusion at its root and lets the
existing DGSI optimiser be rebuilt on a sound footing (this supersedes audit
§3.2 and the `dg`-seeding defect).

### 2.4 Noise handling

AlphaSimR is stochastic; the objective is noisy. Three requirements:

- **Common random numbers.** Use identical founder genomes and identical
  seed streams across candidate `d` values. This removes most of the variance
  from *comparisons* between directions and is worth more than doubling the
  replicate count.
- **Nugget in the GP.** Model replicate noise explicitly.
- **Pareto frontier from posterior means, not raw draws.** Computing a frontier
  from raw replicate outcomes produces a frontier made of lucky runs. This is a
  common and serious error.

### 2.5 Four ranking modes, one engine

All four requested modes are the same GP with a different acquisition function:

| Mode | Objective | Acquisition |
|---|---|---|
| `"pareto"` | non-dominated cumulative gain vectors | expected hypervolume improvement |
| `"economic"` | max `w'R_cum` | expected improvement |
| `"target"` | min distance to stated absolute targets | expected improvement |
| `"constrained"` | max focal trait s.t. `R_j ≥ floor_j` | constrained / feasibility-weighted EI |

`"pareto"` should be the default. It is the only mode that does not require the
breeder to supply the thing objective 1 says they cannot supply, and choosing a
point on a frontier is a far easier judgement than stating weights — it is
weight elicitation by revealed preference.

Never return a single "best `d`". Return the **stability region**: all
directions whose posterior outcome is within a stated tolerance of the
optimum. The optimum is conditional on an estimated `G` and an assumed genetic
architecture, and false precision here would be worse than no answer.

### 2.6 Mating-system cycle templates

Three templates, each a distinct AlphaSimR recurrent-selection loop.

**Self-pollinated / inbred lines** (wheat, rice, bean, cowpea). Cross selected
parents → F1 → selfing generations or `makeDH()` → line evaluation → index →
recycle. Additive variance among fixed lines is roughly `2σ²_A`; recombination
restores selection-induced disequilibrium slowly.

**Cross-pollinated / recurrent selection** (maize, pearl millet, sorghum).
`randCross()` among selected parents each cycle. Classical Bulmer applies: half
the disequilibrium decays per generation of random mating.

**Clonal** (cassava, sweetpotato, banana, potato). Cross → seedling nursery →
multi-stage clonal evaluation → select clones as parents. The decision variable
is **total genetic value** `gv()` (additive + dominance), not breeding value.
This has a direct consequence for the index itself: `G` supplied to the index
must be the **genotypic** covariance, not the additive covariance, and the
package should validate and warn on this rather than accept either silently.

Each template exposes cycle length, stage sizes, per-stage heritability,
`N_e`, and the number of parents recycled.

### 2.7 The largest scientific risk: founder calibration

AlphaSimR at full power requires a specified genetic architecture — QTL number
and effect distribution, dominance, LD structure, genetic correlations, founder
allele frequencies. If those are assumed rather than calibrated, the "optimal
desired gains" are optimal for a fictional population, and the whole feature
becomes an elaborate way of restating the assumptions.

Proposed mitigation: a `calibrate_founders()` step that tunes the simulated
founder population to match the breeder's observed data on

- the estimated genetic covariance matrix `G`,
- observed heritabilities per trait,
- LD decay from real marker data,
- allele frequency spectrum,

and reports match quality per criterion, refusing to proceed silently when the
match is poor. **Open question for you: is real marker and phenotype data
available for calibration in the target programmes, or must the first release
run on assumed architectures?** This changes how the feature should be framed
and how strongly its output should be caveated.

### 2.8 Index re-estimation within the simulation

Each simulated cycle can either (a) rebuild the index from `G` and `P`
re-estimated from that cycle's simulated data, or (b) reuse fixed coefficients.
(a) is what breeders actually do and captures the compounding of estimation
error over cycles; (b) isolates the effect of `d` alone. Recommend (a) as
default with (b) available, and report both when they diverge — a large
divergence means the recommendation is sensitive to estimation error, which the
breeder needs to know.

---

## 3. Proposed API surface

All additive. No existing argument renamed or removed, so `hapblockr` is
unaffected.

```r
# Exact tools — 0.4.0, no new dependencies
implied_economic_weights(d, G, P)
implied_desired_gains(w, G, P)
gain_feasibility(d, G, P, n_candidates, n_select)   # i_required, closest feasible
economic_weights_from_profit(profit_fn, means, ...)
elicit_economic_weights(traits, method = c("tradeoff", "swing"), ...)
weight_sensitivity(w, G, P, values, ...)

# Simulation — 0.5.0, AlphaSimR in Suggests
founder_population(mating_system = c("self", "outcross", "clonal"), ...)
calibrate_founders(founders, observed_G, observed_h2, marker_data = NULL, ...)
simulate_selection_cycles(founders, d, n_cycles, ...)

# Optimiser — 0.6.0
optimize_desired_gains(
  founders, n_cycles,
  mode = c("pareto", "economic", "target", "constrained"),
  budget, surrogate_prior = c("deterministic", "none"),
  parallel = TRUE, checkpoint = NULL, ...
)

# Existing entry point gains a variant switch, default unchanged
run_dgsi(..., variant = c("classical", "optimized"))
```

## 4. Dependencies and CRAN constraints

- `AlphaSimR` in **Suggests**, never Imports. Every example, test and vignette
  chunk that touches it wrapped in `skip_if_not_installed()`.
- Gaussian process: hand-rolled Matérn GP (~150 lines) rather than a new
  dependency, keeping the current lean profile.
- Parallelism over replicates and candidates via `future`/`future.apply`
  (Suggests), with a sequential fallback.
- The optimiser must never run during `R CMD check`. Vignettes covering it must
  be precomputed.
- Checkpoint to disk and resume: a realistic budget is hours, and an
  interrupted run must not lose everything.

## 5. Sequencing

Ship §0 and §1 first, as 0.4.0. They deliver most of objective 1 and the
feasibility half of objective 2 exactly, with no dependencies, no simulation,
and no runtime concerns. Do not gate that value behind the simulation build.

0.5.0: founder construction, calibration, three cycle templates, and
`simulate_selection_cycles()` for a *user-supplied* `d`. Useful on its own —
a breeder can compare a handful of candidate strategies by hand.

0.6.0: the surrogate-assisted optimiser and the four ranking modes.

## 6. Interaction with the open audit items

Compounding an uncertain `G` over `T` cycles amplifies its error, so audit §3.6
(propagating covariance uncertainty) stops being a nice-to-have and becomes a
precondition for §2 being honest. It should land in 0.4.0 alongside the exact
tools, not after the simulator.

Audit §3.2 (no closed-form expected response reported by `run_dgsi()`) is
subsumed by §0 here and should be implemented as part of it.
