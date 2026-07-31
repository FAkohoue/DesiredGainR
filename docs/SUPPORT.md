# Support

Use the [GitHub issue
tracker](https://github.com/FAkohoue/DesiredGainR/issues) for
reproducible bugs and documentation problems. Use a discussion forum or
your institutional quantitative-genetics support channel for study
design, breeding-policy choices, and interpretation specific to
confidential programme data.

## Before requesting support

Record the package version,
[`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html), input
dimensions, random seed, the method selected, and the complete warning
or error message. Replace real identifiers and values with a synthetic
example wherever possible.

## Questions the documentation already answers

Several recurring questions are covered in depth, and the relevant
vignette will usually be faster than an issue.

| Question | Where it is answered |
|----|----|
| How do I state a breeding objective? | [`vignette("DesiredGainR-objective")`](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-objective.md) |
| What shape must my inputs take? | [`vignette("DesiredGainR-introduction")`](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-introduction.md) |
| How do I run this from start to finish? | [`vignette("DesiredGainR-pipeline")`](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-pipeline.md) |
| Why are some implied weights negative? | [`vignette("DesiredGainR-objective")`](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-objective.md) |
| Why is relative efficiency below one? | [`?evaluate_index`](https://FAkohoue.github.io/DesiredGainR/reference/evaluate_index.md) |
| Why was my genotype dosage refused? | [`vignette("DesiredGainR-simulation")`](https://FAkohoue.github.io/DesiredGainR/articles/DesiredGainR-simulation.md) |
| Where do I get a genetic covariance matrix? | [`?estimate_genetic_covariance`](https://FAkohoue.github.io/DesiredGainR/reference/estimate_genetic_covariance.md) |

A plain-language guide for readers who approve a selection decision
without running R is available through
[`open_desiredgain_guide()`](https://FAkohoue.github.io/DesiredGainR/reference/open_desiredgain_guide.md).

## What is outside the scope of support

DesiredGainR begins after the genetic evaluation. Questions about
fitting the upstream multi-trait model, analysing field trials, or
estimating breeding values belong to the software that performs those
steps. Questions about parent selection, mate allocation and crossing
plans belong to a mating-design tool.
