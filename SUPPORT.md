# Support

Use the [GitHub issue tracker](https://github.com/FAkohoue/DesiredGainR/issues)
for reproducible bugs and documentation problems. Use a discussion forum or
your institutional quantitative-genetics support channel for study design,
breeding-policy choices, and interpretation specific to confidential programme
data.

## Before requesting support

Record the package version, `sessionInfo()`, input dimensions, random seed,
the method selected, and the complete warning or error message. Replace real
identifiers and values with a synthetic example wherever possible.

## Questions the documentation already answers

Several recurring questions are covered in depth, and the relevant vignette
will usually be faster than an issue.

| Question | Where it is answered |
|---|---|
| How do I state a breeding objective? | `vignette("DesiredGainR-objective")` |
| What shape must my inputs take? | `vignette("DesiredGainR-introduction")` |
| How do I run this from start to finish? | `vignette("DesiredGainR-pipeline")` |
| Why are some implied weights negative? | `vignette("DesiredGainR-objective")` |
| Why is relative efficiency below one? | `?evaluate_index` |
| Why was my genotype dosage refused? | `vignette("DesiredGainR-simulation")` |
| Where do I get a genetic covariance matrix? | `?estimate_genetic_covariance` |

A plain-language guide for readers who approve a selection decision without
running R is available through `open_desiredgain_guide()`.

## What is outside the scope of support

DesiredGainR begins after the genetic evaluation. Questions about fitting the
upstream multi-trait model, analysing field trials, or estimating breeding
values belong to the software that performs those steps. Questions about parent
selection, mate allocation and crossing plans belong to a mating-design tool.
