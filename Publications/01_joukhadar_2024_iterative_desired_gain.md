# Joukhadar et al. (2024): Optimising desired gain indices to maximise selection response

## Bibliographic information

- **Authors:** Reem Joukhadar, Yongjun Li, Rebecca Thistlethwaite, Kerrie L. Forrest, Josquin F. Tibbits, Richard Trethowan, and Matthew J. Hayden
- **Journal:** *Frontiers in Plant Science* 15:1337388
- **Year:** 2024
- **DOI:** 10.3389/fpls.2024.1337388
- **Document type:** Original research
- **Source file:** `fpls-15-1337388 (4).pdf`

## Central problem

Desired-gain selection indices are attractive when reliable economic values are unavailable. The breeder specifies a vector of desired changes, and the index coefficients are calculated to move all traits in the preferred direction. However, the classical desired-gain formulation does not necessarily maximise total selection response, and an arbitrarily chosen desired-gain vector may produce a realised response that differs substantially from the breeder's target.

The paper addresses this limitation by treating the desired-gain vector as an optimisable input. Instead of accepting a single breeder-defined vector, the proposed algorithm iteratively samples candidate desired-gain values and retains those that make the realised response of the selected population approach a user-defined response target.

## Genetic material and phenotypic information

The study used two related but independent wheat populations:

1. A **reference population of 3,331 bread wheat lines**, including Australian checks, breeding lines, and material derived from diverse germplasm sources.
2. A **validation population of 3,005 doubled-haploid lines** derived from crosses involving parents selected from the reference population.

The reference population was evaluated in 30 irrigated field trials between 2014 and 2020 at several Australian locations. Trials included optimal and late sowing treatments. Seven traits were included:

- grain yield (YLD),
- thousand-kernel weight (TKW),
- protein content (Prot),
- screening percentage (Screen),
- stem rust resistance (SR),
- leaf rust resistance (LR), and
- yellow/stripe rust resistance (YR).

The number of trials and observations differed among traits. Grain yield had 30 trials and 11,053 records, whereas each rust trait was evaluated in five trials. Genetic correlations between optimal and late sowing were high for the non-rust traits, ranging from 0.89 for yield to 0.98 for protein, suggesting limited reranking between sowing treatments.

## Genomic prediction

The reference population was genotyped with the 90K wheat SNP array. After quality control and imputation, genomic estimated breeding values were obtained with BayesR. The model assigned marker effects to a mixture of four normal distributions representing null, small, medium, and large effects. Prediction accuracy was evaluated with repeated random cross-validation.

Mean prediction accuracies were:

| Trait | Prediction accuracy |
|---|---:|
| YLD | 0.29 |
| Screen | 0.30 |
| Prot | 0.32 |
| TKW | 0.39 |
| SR | 0.40 |
| LR | 0.42 |
| YR | 0.47 |

The rust traits had the highest average prediction accuracy, while grain yield had the lowest. These values indicate that the index was developed from moderately noisy, but practically relevant, GEBVs rather than from perfectly predicted genetic values.

## Genetic parameters

Narrow-sense heritability ranged from 0.21 for grain yield to 0.59 for yellow rust. Most pairwise genetic correlations were small, ranging approximately from -0.21 to 0.28. The major exception was the strong negative correlation between TKW and screening percentage, estimated at -0.61. This relationship is biologically consistent because small kernels contribute more strongly to the screened fraction.

The low correlations among most traits simplified simultaneous improvement, whereas the TKW-Screen relationship created a strong correlated response that the index needed to accommodate.

## Classical desired-gain index

The paper used the Yamada et al. desired-gain formulation:

\[
\mathbf{b}=\mathbf{P}^{-1}\mathbf{G}
\left(\mathbf{G}\mathbf{P}^{-1}\mathbf{G}\right)^{-1}\mathbf{d},
\]

where:

- \(\mathbf{b}\) is the vector of index coefficients,
- \(\mathbf{P}\) is the variance-covariance matrix of the selection criteria,
- \(\mathbf{G}\) is the genetic variance-covariance matrix, and
- \(\mathbf{d}\) is the desired-gain vector.

The difficulty is that a breeder-supplied \(\mathbf{d}\) does not directly equal the realised response among selected candidates. The relationship depends on the covariance structure, selection intensity, and available distribution of candidates.

## Iterative optimisation method

The proposed method defines a user-specified target response \(\mathbf{d_g}\), expressed in population standard-deviation units. For each iteration, the algorithm:

1. samples a candidate desired-gain vector \(\mathbf{d}\),
2. calculates the corresponding index coefficients \(\mathbf{b}\),
3. ranks the candidate population,
4. selects the highest-ranking individuals,
5. calculates the realised response for each trait, and
6. evaluates the discrepancy between the realised and targeted responses.

A goodness-of-fit function is maximised when the realised response equals the target. The overall penalty is averaged across traits, and a newly sampled vector is accepted when it reduces the discrepancy. The analysis was run for 1,000 iterations.

The covariance matrices were estimated from the 3,331-line reference population and applied to the 3,005 doubled-haploid candidates. The top 100 candidates, approximately 3.3% of the candidate population, were selected. Each scenario was replicated 20 times. The mean correlation among replicate index solutions was approximately 0.95, indicating high repeatability of the optimisation.

## Selection scenarios

### INDEX1: equal desired response

The target was +0.5 standard deviations for YLD, TKW, and Prot and -0.5 standard deviations for Screen and the three rust scores.

Realised responses in the independent candidate population were:

| Trait | Target | Iterative response |
|---|---:|---:|
| YLD | 0.50 | 0.83 |
| TKW | 0.50 | 1.02 |
| Prot | 0.50 | 0.66 |
| Screen | -0.50 | -0.46 |
| SR | -0.50 | -0.56 |
| LR | -0.50 | -0.61 |
| YR | -0.50 | -0.62 |

The realised responses were all in the required direction. The absolute response averaged 0.68 standard deviations. The conventional non-iterative method produced a much more uneven response, ranging from approximately 0.12 to 1.69 standard deviations in absolute value.

### INDEX2: yield-dominant response

The target was +2 standard deviations for YLD, +0.5 for TKW and Prot, and -0.5 for Screen and rust scores.

| Trait | Target | Iterative response in candidates | Non-iterative response |
|---|---:|---:|---:|
| YLD | 2.00 | 1.71 | 2.01 |
| TKW | 0.50 | 0.46 | 0.74 |
| Prot | 0.50 | 0.51 | -0.35 |
| Screen | -0.50 | -0.50 | -0.08 |
| SR | -0.50 | -0.54 | -0.21 |
| LR | -0.50 | -0.54 | 0.13 |
| YR | -0.50 | -0.54 | -0.45 |

The iterative method maintained the desired direction for all traits. In contrast, the classical non-iterative solution produced an unfavourable decrease in protein and an unfavourable increase in leaf-rust score. This scenario provides the clearest evidence that matching a target response is more informative than simply assigning a large nominal desired gain to yield.

### INDEX3: unconstrained maximum-response scenario

The target was set at an intentionally unattainable magnitude of +4 standard deviations for favourable traits and -4 for unfavourable traits. This forced the algorithm to search for the strongest feasible response in the desired directions.

The independent candidate population showed responses of 1.05 for YLD, 0.64 for TKW, 1.10 for Prot, -1.24 for Screen, -1.09 for SR, -1.06 for LR, and -1.11 for YR. The mean absolute response was approximately 1.04 standard deviations.

Because the requested response was impossible at the applied selection intensity, the algorithm effectively behaved as a directional response maximiser. This demonstrates that the method can also be used when the breeder wants the greatest feasible improvement rather than a tightly constrained response.

## Main interpretation

The key contribution is the separation of two concepts that are often conflated:

- the **desired-gain input vector** used mathematically to derive index coefficients, and
- the **realised selection response** observed after candidates are ranked and selected.

The breeder's real decision concerns the latter. The iterative method searches for an input vector that produces the requested realised response. Consequently, the derived index coefficients may appear counterintuitive and should not be interpreted as direct biological or economic importance weights.

## Strengths

- Large, multi-environment reference population.
- Independent validation in 3,005 doubled-haploid candidates.
- Explicit use of realised response rather than nominal input weights.
- Ability to handle equal, trait-dominant, and maximum-response objectives.
- Replicated optimisation demonstrating stable solutions.
- Application to GEBVs, making the method compatible with genomic selection pipelines.

## Limitations and cautions

- The study did not directly prove that the iterative index maximises correlation with true aggregate genetic merit because true economic weights and net merit were undefined.
- Validation was across related generations, but not across several recurrent breeding cycles.
- The covariance matrices were estimated in the reference population and assumed to remain relevant in the candidate population.
- Prediction accuracy was moderate, particularly for grain yield, and uncertainty in GEBVs was not explicitly propagated into the optimisation.
- An infeasible target is converted into a best-feasible response, but the resulting trade-off depends on the implemented penalty function.

## Practical relevance

This method is particularly useful when:

- economic weights are unavailable or politically difficult to define;
- the breeder can express objectives as desired changes in standard-deviation units;
- several quantitative traits must be improved simultaneously;
- some traits require controlled response while others can respond freely; and
- GEBVs or adjusted means are already available for all candidates.

The paper provides a strong methodological extension of desired-gain selection. Its practical message is that breeders should optimise and validate the realised response of an index rather than relying on a single arbitrary desired-gain vector.
