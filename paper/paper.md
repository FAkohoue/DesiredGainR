---
title: "DesiredGainR: auditable multi-trait selection indices from breeder intent"
tags:
  - R
  - plant breeding
  - quantitative genetics
  - selection index
  - multi-trait selection
authors:
  - name: Félicien Akohoué
    orcid: 0000-0002-2160-0182
    affiliation: 1
affiliations:
  - name: CGIAR
    index: 1
date: 31 July 2026
bibliography: paper.bib
---

# Summary

A selection index combines several traits into one score so that candidates can
be ranked. The classical economic construction requires a weight for every
trait: how much a unit of grain yield is worth relative to a unit of plant
height. Most breeding programmes cannot supply these. What they can supply is
an intent — lift yield without lengthening the cycle, hold disease resistance
where it is — and an intent is not a weight vector.

`DesiredGainR` accepts desired gains as a direct statement of breeding intent.
It also maps desired gains and aggregate weights algebraically under a stated
covariance model. An implied aggregate-weight vector reproduces an index
direction. It is not an independently estimated economic value. The package
reports whether a stated objective is attainable before an index is built. It
implements the supported selection methods on common data and reports how far
a recommendation depends on estimated covariance matrices.

# Statement of need

Existing R packages compute selection indices from supplied weights. They do
not help with the step that most often goes wrong, which is arriving at the
weights, and they do not report the diagnostics needed to know whether the
answer is stable.

Three problems recur in practice.

**Desired gains are usually unattainable as stated.** Response is bounded by
selection intensity and by available genetic variation. A breeder asking for a
half-day reduction in anthesis-silking interval alongside a 10 per cent yield
increase may be asking for a selection intensity no programme can run.
`gain_feasibility()` reports the required intensity and the selected fraction
it implies, exactly and without simulation, so an infeasible target is caught
before it is pursued.

**The genetic covariance matrix is an estimate, and every reported gain treats
it as known.** Multi-trait restricted maximum likelihood typically estimates
more variance components than the trial has independent families to determine.
`index_uncertainty()` resamples the covariance and refits, and for the
desired-gain families the resulting measure has an exact reference point.
Because $\mathbf{G}\mathbf{b} = \mathbf{d}$ for the Pesek–Baker index, the
model-based response follows the requested proportions when the same covariance
model defines the fit and the evaluation. Departure after covariance
resampling measures sensitivity to covariance estimation. The observed
differential among a finite candidate set remains a separate quantity.

**Many constrained linear objectives have closed-form solutions.** The
Kempthorne–Nordskog, Tallis, Mallard and Harville indices constrain theoretical
response under an aggregate-merit objective. `restricted_index()` implements
these methods. `run_dgsi()` solves a different problem. It searches
desired-gain directions to improve the selected differential in a finite
candidate set.

# Design

The package is organised as layers, and their maturity is stated explicitly in
`MATURITY.md` rather than left to be inferred.

The objective layer (`implied_economic_weights()`, `implied_desired_gains()`,
`gain_feasibility()`, `weight_sensitivity()`, `retrospective_weights()`) is the
distinctive contribution and has no direct equivalent elsewhere.
`retrospective_weights()` in particular recovers the weight vector implied by a
selection decision that has already been made, which is often the only way to
discover what a programme's objective has actually been.

The index layer implements Smith–Hazel, base, Pesek–Baker, Yamada,
Mulamba–Mock and Elston through `selection_index()`. The same function includes
independent culling and a within-cohort sequential screen as operational
comparators. The restricted families are available through
`restricted_index()`. `run_dgsi()` applies iterative search to the established
desired-gain index. `run_qgsi()` fits the quadratic genomic index of
Cerón-Rojas et al.

The simulation layer runs a programme forward through `AlphaSimR` to compare
desired-gain directions over several cycles, because a direction maximising
first-cycle response may exhaust the variance that response depends on.

Multi-environment problems are handled by expansion:
`expand_environments()` builds the separable covariance
$\mathrm{Cov}(g_{je}, g_{kf}) = \mathbf{G}_{jk}\mathbf{C}_{ef}$ over
trait-environment combinations, so the whole index layer applies unchanged.

# Auditability

The package is written on the premise that a selection decision must be
defensible after the fact, which imposes requirements beyond correctness.

Covariance provenance is recorded rather than assumed: whether `P` was supplied
or estimated, its numerical rank, and whether `P - G` is a valid residual
covariance. Inadmissible inputs are refused with an explanation rather than
producing a heritability above one. Where a covariance matrix has been repaired
to make it invertible, the repair and its magnitude are reported alongside any
result derived from it.

Simulation results carry a fingerprint of every input that determined them,
including package versions, and a stored checkpoint that does not match the
current problem stops rather than silently resuming.

`vignette("DesiredGainR-reproduction")` checks the package against tables
published by Rahimi and Debnath [@rahimi2023], computed independently in SAS.

# Empirical evidence

`vignette("DesiredGainR-empirical-validation")` adds evidence from six real
breeding-programme sources rather than treating simulation as its own external
validation. CNA6 rice and a tropical maize haploid-inducer population provide
observed recurrent-cycle trends [@bartholome2023; @fritscheneto2023]. The INIA
Uruguay rice database provides 23 years of stage-forward validation
[@rebollo2023]. A desired-gain index fitted through 2008 and frozen before
scoring 2009--2016 first-stage (E1) cohorts had an area under the receiver
operating characteristic curve (AUC) of 0.773, compared with
0.756 for yield alone; a hierarchical year-and-candidate bootstrap placed the
difference above zero. Genomes-to-Fields maize provides a non-overlapping
year-window transport check [@lopezcruz2023], and the CIMMYT rapid-cycle wheat
experiment provides a published predicted-versus-realised response anchor
[@dreisigacker2023].
The CIMMYT tropical maize rapid-cycle analysis uses the phenotype, marker and
sample-linkage deposit associated with Zhang, Pérez-Rodríguez, Burgueño, Olsen,
Jannink, Buckler, Atlin, Boddupalli, Vargas, San Vicente and Crossa
[@zhang2017data]. The associated article was authored by Zhang,
Pérez-Rodríguez, Burgueño, Olsen, Buckler, Atlin, Boddupalli, Vargas, San
Vicente and Crossa [@zhang2017].

The evidence boundary is part of the result, not only the discussion. These
datasets test observed cycle response, temporal transport and agreement with
historical advancement. They do not identify the causal optimum of a
desired-gain vector, because comparable populations were not prospectively
assigned to alternative vectors. Population-driven recommendations therefore
remain conditional decision support pending such a trial.

# Acknowledgements

Development benefited from detailed technical review by colleagues at CGIAR.

# References
