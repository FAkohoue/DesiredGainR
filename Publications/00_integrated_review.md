# Detailed analytical review of seven publications and technical documents on selection indices

## Scope

This document provides a detailed analytical rewrite of seven uploaded PDFs. It is not a verbatim transcription. The publications and presentations were read as complementary sources addressing the theory, computation, implementation, and empirical evaluation of multi-trait selection in plant breeding.

The collection contains:

1. an iterative desired-gain method validated with genomic breeding values in wheat;
2. a broad review of breeding simulation;
3. two practical presentations on index implementation in maize and rice programs;
4. a concise implementation guideline;
5. an R and SAS coding paper for classical indices; and
6. an empirical rice recurrent-selection study.

## Documents reviewed

| No. | Document | Type | Crop or scope | Primary contribution |
|---:|---|---|---|---|
| 1 | Joukhadar et al. (2024) | Original research | Bread wheat | Iterative optimisation of desired-gain indices |
| 2 | Vieira et al. (2025) | Review | General plant breeding | Simulation-based optimisation of selection strategies |
| 3 | Covarrubias-Pazaran, CIMMYT maize | Presentation | Maize | Retrospective reconstruction of a selection procedure |
| 4 | Covarrubias-Pazaran, IRRI programs | Presentation | Rice | Regional index implementation for yield, zinc, genes, and family balance |
| 5 | Practical implementation guideline | Guideline | General plant breeding | Stage-by-stage operational protocol |
| 6 | Rahimi and Debnath (2023) | Methods/software article | Maize example; general application | R and SAS code for optimum, base, and desired-gain indices |
| 7 | Guimarães et al. (2021) | Empirical research | Lowland rice | Comparison of direct selection and three indices |

## Common conceptual framework

All documents address the same fundamental problem: a breeding candidate is valuable because of several traits, but the traits differ in scale, heritability, reliability, economic importance, and genetic correlation. Selecting the best candidate for each trait separately does not necessarily identify the best candidate for total merit.

A general linear selection index is

\[
I=\mathbf{b}'\mathbf{x},
\]

where \(\mathbf{x}\) is a vector of observed, adjusted, or genomically predicted selection criteria and \(\mathbf{b}\) is a vector of coefficients. The main index families in the uploaded documents are:

### Optimum or Smith-Hazel index

\[
\mathbf{b}=\mathbf{P}^{-1}\mathbf{G}\mathbf{a}.
\]

The economic-value vector \(\mathbf{a}\) defines aggregate merit. The covariance matrices translate economic importance into coefficients that account for scale, reliability, and correlation.

### Base index

\[
\mathbf{b}=\mathbf{a}.
\]

This is simple but does not adjust weights for variance or covariance.

### Desired-gain index

The breeder specifies a desired response vector rather than economic values. The Joukhadar paper uses the Yamada formulation and then iteratively changes its desired-gain input until realised response approaches the breeder's target.

### Retrospective index

\[
\mathbf{b}=\mathbf{P}^{-1}\mathbf{s},
\]

where \(\mathbf{s}\) contains historical selection differentials. This reconstructs the linear preference implicit in past selections.

### Rank-sum index

\[
I_i=\sum_j r_{ij}.
\]

Candidates are ranked trait by trait in the preferred direction and ranks are summed. This avoids economic weights and covariance matrices but discards information about the magnitude and reliability of differences.

## Major synthesis

### Economic values and index coefficients are not the same

A recurring message is that the breeder should define the objective, not manually interpret raw coefficients. Economic values, desired gains, and realised response describe the objective. Index coefficients are derived from the covariance structure. Strongly correlated traits can generate small, large, or counterintuitive coefficients even when both traits are important.

### Hard constraints and compensatory traits should be separated

The maize and rice presentations show that operational selection combines culling and ranking. Mandatory resistance genes, unacceptable maturity, excessive plant height, or minimum quality thresholds are often non-compensatory. They should remain hard constraints. Yield, zinc, adaptation, and other quantitative traits can then be combined in an index if compensation is acceptable.

### Realised response is the final criterion

Joukhadar et al. show that a nominal desired-gain vector does not guarantee a matching realised response. Guimarães et al. demonstrate that an index with theoretically meaningful weights can still perform poorly if the weight vector is badly chosen. Rahimi and Debnath show that indices with similar aggregate diagnostics may rank candidates similarly, while another desired-gain formulation can produce a completely different ranking. Therefore, each index should be evaluated by per-trait selection differentials, aggregate response, candidate overlap, and prospective validation.

### Standardisation improves communication

The practical guideline recommends expressing desired response in standard-deviation units. This enables the breeder to say, for example, "increase yield by one standard deviation while decreasing disease score by half a standard deviation." Standardisation does not remove covariance effects, but it makes the requested objective interpretable across traits measured in different units.

### Retrospective weights are a starting point, not an optimum

The CIMMYT and IRRI presentations derive initial coefficients from historical selections. This is useful for translating breeder judgement into an auditable rule. However, the resulting index reproduces past decisions and may reproduce past biases. The index must then be adjusted against the target product profile and validated using future candidates.

### Genomic selection changes the unit of evidence, not the objective

Joukhadar et al. apply the index to GEBVs. The same logic can be applied to BLUEs, BLUPs, GEBVs, or other predictions, provided their covariance structure and reliability are appropriate. The simulation review emphasises that genomic selection is especially valuable when it reduces cycle time and improves selection of low-heritability traits. However, training sets must be updated and remain related to future candidates.

### Diversity must be included in the decision system

The Vieira review and the family-balanced rice presentation show that unrestricted index selection can concentrate contributions in a small number of families. Family quotas, optimal contribution selection, coancestry penalties, or within-family selection can preserve diversity. These constraints generally reduce immediate index response but can improve long-term gain.

## Recommended operational workflow

A robust implementation based on the combined evidence is:

1. **Define the breeding objective.** Specify the target product profile, target population of environments, selection stage, and planning horizon.
2. **Classify traits.** Separate mandatory constraints from compensatory quantitative objectives.
3. **Assemble selection criteria.** Obtain adjusted means, BLUPs, or GEBVs for every candidate and trait.
4. **Orient and standardise traits.** Ensure that positive values consistently represent improvement and express desired response in standard-deviation units.
5. **Estimate covariance matrices.** Prefer additive genetic and prediction-error-aware estimates when possible. Verify matrix ordering and numerical conditioning.
6. **Choose the index family.** Use an optimum index when meaningful economic values exist; desired-gain optimisation when response targets are clearer; retrospective weights to initialise implementation; or rank-sum when covariance and weights are unreliable.
7. **Apply diversity constraints.** Use family balance, coancestry restrictions, or optimal contribution when parent selection affects long-term diversity.
8. **Evaluate response.** Report per-trait selection differentials, aggregate merit, overlap with current selections, and sensitivity to weights.
9. **Validate prospectively.** Test selected parents or lines in progeny and later-stage trials.
10. **Recalibrate each cycle.** Update means, variances, correlations, product-profile priorities, and training data.

## Overall conclusions

The uploaded documents collectively show that selection indices are not simply weighted sums. Their usefulness depends on the definition of the breeding objective, the quality of the underlying predictions, the covariance structure, the direction of selection, and the constraints imposed on diversity and mandatory traits.

The strongest methodological advance in the collection is the iterative desired-gain method of Joukhadar et al., because it optimises the relationship between a breeder-defined target and realised response. The strongest operational contribution is the Covarrubias-Pazaran workflow, which begins by translating the current selection process into an explicit algorithm. The Guimarães rice study provides the clearest empirical warning about economic weights: a simple rank-sum index can outperform more formal indices when their weights are poorly specified. The Vieira review places all these decisions in a long-term framework in which gain, time, cost, GxE, and diversity must be optimised jointly.

---


## Joukhadar et al. (2024): Optimising desired gain indices to maximise selection response

### Bibliographic information

- **Authors:** Reem Joukhadar, Yongjun Li, Rebecca Thistlethwaite, Kerrie L. Forrest, Josquin F. Tibbits, Richard Trethowan, and Matthew J. Hayden
- **Journal:** *Frontiers in Plant Science* 15:1337388
- **Year:** 2024
- **DOI:** 10.3389/fpls.2024.1337388
- **Document type:** Original research
- **Source file:** `fpls-15-1337388 (4).pdf`

### Central problem

Desired-gain selection indices are attractive when reliable economic values are unavailable. The breeder specifies a vector of desired changes, and the index coefficients are calculated to move all traits in the preferred direction. However, the classical desired-gain formulation does not necessarily maximise total selection response, and an arbitrarily chosen desired-gain vector may produce a realised response that differs substantially from the breeder's target.

The paper addresses this limitation by treating the desired-gain vector as an optimisable input. Instead of accepting a single breeder-defined vector, the proposed algorithm iteratively samples candidate desired-gain values and retains those that make the realised response of the selected population approach a user-defined response target.

### Genetic material and phenotypic information

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

### Genomic prediction

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

### Genetic parameters

Narrow-sense heritability ranged from 0.21 for grain yield to 0.59 for yellow rust. Most pairwise genetic correlations were small, ranging approximately from -0.21 to 0.28. The major exception was the strong negative correlation between TKW and screening percentage, estimated at -0.61. This relationship is biologically consistent because small kernels contribute more strongly to the screened fraction.

The low correlations among most traits simplified simultaneous improvement, whereas the TKW-Screen relationship created a strong correlated response that the index needed to accommodate.

### Classical desired-gain index

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

### Iterative optimisation method

The proposed method defines a user-specified target response \(\mathbf{d_g}\), expressed in population standard-deviation units. For each iteration, the algorithm:

1. samples a candidate desired-gain vector \(\mathbf{d}\),
2. calculates the corresponding index coefficients \(\mathbf{b}\),
3. ranks the candidate population,
4. selects the highest-ranking individuals,
5. calculates the realised response for each trait, and
6. evaluates the discrepancy between the realised and targeted responses.

A goodness-of-fit function is maximised when the realised response equals the target. The overall penalty is averaged across traits, and a newly sampled vector is accepted when it reduces the discrepancy. The analysis was run for 1,000 iterations.

The covariance matrices were estimated from the 3,331-line reference population and applied to the 3,005 doubled-haploid candidates. The top 100 candidates, approximately 3.3% of the candidate population, were selected. Each scenario was replicated 20 times. The mean correlation among replicate index solutions was approximately 0.95, indicating high repeatability of the optimisation.

### Selection scenarios

#### INDEX1: equal desired response

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

#### INDEX2: yield-dominant response

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

#### INDEX3: unconstrained maximum-response scenario

The target was set at an intentionally unattainable magnitude of +4 standard deviations for favourable traits and -4 for unfavourable traits. This forced the algorithm to search for the strongest feasible response in the desired directions.

The independent candidate population showed responses of 1.05 for YLD, 0.64 for TKW, 1.10 for Prot, -1.24 for Screen, -1.09 for SR, -1.06 for LR, and -1.11 for YR. The mean absolute response was approximately 1.04 standard deviations.

Because the requested response was impossible at the applied selection intensity, the algorithm effectively behaved as a directional response maximiser. This demonstrates that the method can also be used when the breeder wants the greatest feasible improvement rather than a tightly constrained response.

### Main interpretation

The key contribution is the separation of two concepts that are often conflated:

- the **desired-gain input vector** used mathematically to derive index coefficients, and
- the **realised selection response** observed after candidates are ranked and selected.

The breeder's real decision concerns the latter. The iterative method searches for an input vector that produces the requested realised response. Consequently, the derived index coefficients may appear counterintuitive and should not be interpreted as direct biological or economic importance weights.

### Strengths

- Large, multi-environment reference population.
- Independent validation in 3,005 doubled-haploid candidates.
- Explicit use of realised response rather than nominal input weights.
- Ability to handle equal, trait-dominant, and maximum-response objectives.
- Replicated optimisation demonstrating stable solutions.
- Application to GEBVs, making the method compatible with genomic selection pipelines.

### Limitations and cautions

- The study did not directly prove that the iterative index maximises correlation with true aggregate genetic merit because true economic weights and net merit were undefined.
- Validation was across related generations, but not across several recurrent breeding cycles.
- The covariance matrices were estimated in the reference population and assumed to remain relevant in the candidate population.
- Prediction accuracy was moderate, particularly for grain yield, and uncertainty in GEBVs was not explicitly propagated into the optimisation.
- An infeasible target is converted into a best-feasible response, but the resulting trade-off depends on the implemented penalty function.

### Practical relevance

This method is particularly useful when:

- economic weights are unavailable or politically difficult to define;
- the breeder can express objectives as desired changes in standard-deviation units;
- several quantitative traits must be improved simultaneously;
- some traits require controlled response while others can respond freely; and
- GEBVs or adjusted means are already available for all candidates.

The paper provides a strong methodological extension of desired-gain selection. Its practical message is that breeders should optimise and validate the realised response of an index rather than relying on a single arbitrary desired-gain vector.

---

## Vieira, Nogueira, and Fritsche-Neto (2025): Optimizing quantitative-trait selection using simulation

### Bibliographic information

- **Authors:** Rafael Augusto Vieira, Ana Paula Oliveira Nogueira, and Roberto Fritsche-Neto
- **Journal:** *Frontiers in Plant Science* 16:1495662
- **Year:** 2025
- **DOI:** 10.3389/fpls.2025.1495662
- **Document type:** Narrative review
- **Source file:** `fpls-16-1495662.pdf`

### Purpose and scope

This review synthesises simulation studies that evaluated the optimisation of quantitative-trait selection in plant breeding. The central premise is that simulation provides an intermediate layer between quantitative-genetic theory and expensive field implementation. A simulation can represent genetic architecture, crossing, recombination, population size, selection intensity, heritability, genotype-by-environment interaction, genotyping cost, and breeding-cycle duration before a breeding program commits resources.

The review is broad. It covers classical phenotypic selection, pre-breeding, marker-assisted selection, QTL introgression, genomic selection, multi-trait selection, diversity management, inbreeding and heterosis, cost allocation, genotype-by-environment interaction, and software for breeding simulation.

### Role of simulation in plant breeding

The authors distinguish deterministic and stochastic approaches. Deterministic models use equations based on heritability, selection intensity, and accuracy to predict expected response. These models are efficient for broad comparisons but cannot fully represent crossing, recombination, segregation, linkage, and recurrent selection. Stochastic simulations generate individual genotypes and phenotypes and are therefore better suited to reproducing complete breeding pipelines.

Simulation is presented as a decision-support tool rather than a replacement for empirical testing. Its main functions are to:

- compare breeding strategies under controlled assumptions;
- quantify trade-offs among gain, diversity, cost, and time;
- identify bottlenecks before changing a production pipeline;
- test sensitivity to genetic architecture and heritability;
- optimise population size, number of parents, and selection timing; and
- evaluate long-term consequences that cannot be observed rapidly in field experiments.

The main limitation is model dependence. A simulation can be precise but wrong if its assumptions about QTL number, effect distribution, recombination, GxE, cost, or selection behaviour do not represent the actual program.

### Diversity, pre-breeding, and parental management

A recurring conclusion across simulation studies is that maximising short-term mean performance can accelerate the loss of genetic diversity. Fewer parents and intense selection can increase immediate gain, whereas larger parental sets and controlled coancestry preserve long-term potential.

The review recommends combining strategies rather than relying on a single selection rule. Overlapping cohorts, strategic recycling, optimal contribution methods, and the introduction of diverse germplasm can improve long-term gain. The optimal number of parents depends on the planning horizon:

- fewer parents can increase short-term gain;
- more parents and more crosses are preferable when long-term improvement and risk control are important.

Pre-breeding simulations indicate that exotic or landrace material should be introduced with deliberate background recovery. Multiple conversion versions and residual donor segments must be managed because aggressive recovery can discard useful diversity, whereas insufficient recovery can reduce agronomic adaptation.

### Phenotypic selection and classical breeding

Simulation studies of self-pollinated crops show that the value of early selection depends on trait heritability, genetic architecture, and the reliability of early-generation phenotypes. Early selection can be efficient for highly heritable traits but can create bias for low-heritability traits and traits strongly affected by GxE.

The review discusses bulk, pedigree, single-seed descent, doubled-haploid, and backcross schemes. No scheme is universally optimal. Population size, number of families, progeny per family, and generation at which selection is initiated must be aligned with the breeding objective. Larger populations generally increase the probability of recovering favourable recombinants, but the benefit is conditional on adequate germplasm, a clear target product profile, and sufficient resources to measure the population accurately.

For introgression, simulations support combining foreground selection for the target locus with background selection for recipient-genome recovery. The number and spacing of markers should be sufficient to track linkage around the target, but merely increasing marker density does not guarantee efficiency if population size and recombination are limiting.

### Marker-assisted selection and QTL management

Marker-assisted selection provides the largest advantage when:

- one or a few loci explain a meaningful proportion of variation;
- marker-QTL linkage is strong and stable;
- the selected alleles have relatively large effects;
- the population is sufficiently large to generate recombination; and
- foreground and background selection are balanced.

Simulations commonly show rapid early gains from MAS, followed by diminishing advantage after several cycles. This occurs because large-effect alleles are fixed rapidly, while residual polygenic variation is not fully captured. MAS is therefore more suitable for disease resistance, quality loci, and major adaptation genes than for highly polygenic traits such as yield.

The review emphasises that low heritability does not directly invalidate MAS, because the selection is marker based. However, low heritability reduces the reliability of QTL discovery and effect estimation, thereby weakening the marker-trait relationship on which MAS depends.

### Genomic selection

Genomic selection is presented as the major simulation focus of modern breeding. Its principal advantage is not only prediction accuracy but the possibility of shortening the breeding cycle. Genetic gain per year depends on selection intensity, accuracy, additive genetic variation, and cycle length. A moderate-accuracy model can therefore outperform a high-accuracy phenotypic scheme if it enables much faster recycling.

Major conclusions include:

- genomic selection is especially advantageous for low-heritability and highly polygenic traits;
- rapid-cycle recurrent genomic selection can produce high gain per unit time;
- prediction accuracy declines as candidates become genetically distant from the training population;
- training sets require regular updating, regardless of the assumed genetic architecture;
- including recent and high-performing candidates in the training set can improve relevance, but diversity must be retained;
- model choice depends on architecture: Bayesian models can perform well when fewer loci have moderate or large effects, whereas BLUP-based models are generally robust for many-QTL traits;
- RR-BLUP is flexible across many conditions and remains a practical baseline; and
- relatedness between training and target populations is often more influential than marker density beyond a sufficient threshold.

The review stresses the distinction between short-term and long-term evaluation. A strategy that gives the highest gain in the first few cycles can exhaust useful variance and perform poorly over a longer horizon.

### Multi-trait selection

Multi-trait genomic prediction is particularly useful when a low-heritability target trait is genetically correlated with a higher-heritability secondary trait. Under favourable correlation structures, the secondary trait contributes information and increases prediction accuracy for the target.

The value of multi-trait analysis depends on:

- the magnitude and direction of genetic correlation;
- the heritability of each trait;
- overlap or imbalance in trait records;
- the stage at which each trait is measured; and
- whether the correlation remains stable across environments and cycles.

Selection indices are discussed as a means to translate multi-trait predictions into a single decision criterion. Simulations are useful because they can compare index weights, desired gains, constraints, and culling thresholds while measuring both total response and undesirable correlated changes.

### Inbreeding, heterosis, and epistasis

The review shows that breeding optimisation must account for non-additive processes when they influence the target product. In hybrid breeding, heterotic-group structure and specific combining ability can alter the optimal choice of parents. In recurrent selection, aggressive genomic selection can increase inbreeding and reduce long-term gain. Optimal contribution and mate-allocation methods can control coancestry while retaining high-value alleles.

Epistasis may improve prediction in some simulated architectures, but its practical benefit is inconsistent and often depends on sample size, relatedness, and the stability of interactions across generations. Additive models remain robust for many operational decisions because additive effects determine transmissible response.

### Cost analysis

Simulation enables the breeder to optimise resource allocation rather than accuracy alone. The review describes examples in which marker-assisted selection reduced cost, and genomic selection increased gain when genotyping was substituted for expensive or slow phenotyping.

Low-density genotyping followed by imputation is highlighted as an important cost-saving strategy. One cited study reported cost reductions of up to 87% with minimal loss in prediction accuracy, and as few as 50 segregating markers per genome were sufficient under the evaluated conditions. These numbers are scenario dependent, but the general conclusion is robust: once marker density is sufficient to recover relationships and linkage information, population size and training relevance can be more valuable than further density increases.

### Genotype-by-environment interaction

GxE can bias estimated breeding values and reduce transferability of selection decisions. Simulation studies indicate that environment-specific models, multi-environment genomic prediction, reaction-norm approaches, and environment-specific selection indices can improve response when target environments differ substantially.

The review notes that genotype-by-year interaction is still represented inadequately in many simulations. This is important because year-to-year climatic variation may be less predictable than location differences. Future simulations should represent changing frequencies of drought, heat, excess rainfall, and disease pressure rather than assuming stationary environments.

### Software tools

The review lists several simulation platforms:

- **AlphaSimR:** flexible simulation of genetic architectures, populations, crossing, and selection;
- **ADAM:** structured populations and inheritance tracking;
- **synbreed:** genomic prediction and SNP-based breeding analysis;
- **MoBPS:** modular simulation of livestock and plant breeding programs, including genomic schemes;
- **QU-LINE:** quantitative-trait and breeding-line simulation, including additive and epistatic effects; and
- **PedigreeSim:** meiosis, recombination, and pedigree-based inheritance in structured populations.

Software selection should follow the biological and operational question. A highly detailed simulator is not automatically preferable when the necessary parameters cannot be estimated reliably.

### Integrated recommendations

The review's practical recommendations can be condensed as follows:

1. Define the target product profile, target population of environments, planning horizon, and resource constraint before optimising the selection method.
2. Match selection timing to trait heritability and measurement quality.
3. Use markers for major loci and genomic prediction for diffuse polygenic variation.
4. Update genomic training sets routinely and maintain relatedness to future candidates.
5. Use multi-trait models when a low-heritability target is correlated with an informative secondary trait.
6. Increase population size only when evaluation quality and germplasm diversity are adequate.
7. Balance short-term gain against inbreeding and loss of rare favourable alleles.
8. Model GxE explicitly when breeding values are expected to differ across target environments.
9. Optimise gain per year and gain per unit cost, not prediction accuracy alone.
10. Validate simulation recommendations with retrospective and prospective data from the actual breeding program.

### Strengths

- Wide coverage of classical and genomic breeding strategies.
- Translation of simulation results into operational recommendations.
- Explicit attention to genetic diversity and long-term response.
- Integration of genetic, temporal, and economic dimensions.
- Recognition that training-set management and cycle length are central to genomic selection.

### Limitations

- The paper is a broad narrative review rather than a formal meta-analysis.
- Results originate from simulations with heterogeneous assumptions, crops, population structures, and cost models.
- Some cited conclusions are highly parameter dependent and should not be transferred directly to another breeding program.
- The review summarises many studies but does not provide a unified quantitative benchmark for comparing strategies.
- Software descriptions are concise and do not assess reproducibility, computational scale, or ease of integration with operational databases.

### Main contribution

The review demonstrates that breeding optimisation is a multi-dimensional problem. The best strategy cannot be identified from heritability, accuracy, or gain alone. Simulation is most valuable when it represents the actual pipeline and jointly measures response, cycle time, diversity, GxE, and cost.

---

## Covarrubias-Pazaran: Bringing a selection index into the CIMMYT maize programs

### Document information

- **Author:** Giovanny Eduardo Covarrubias-Pazaran
- **Affiliation shown:** Excellence in Breeding Platform / CGIAR
- **Document type:** Technical presentation, 21 pages
- **Source file:** `Giovanny Selection index in Maize.pdf`

### Purpose

The presentation documents a practical process for replacing parts of a conventional multi-step maize selection procedure with a formal selection index. It does not begin with a theoretical index and search for a use case. Instead, it first reconstructs the breeder's existing decisions, expresses them as an algorithm, derives retrospective weights from observed selection differentials, and compares the resulting index with the current procedure.

The proposed workflow is:

1. review what is known about selection methods and genetic gain;
2. translate the current selection procedure into explicit algorithmic steps;
3. identify the steps that can be replaced by an index;
4. construct an index;
5. compare selections, selection differentials, and total merit; and
6. refine the index.

### Selection index versus culling and tandem selection

Pages 3 and 4 compare three strategies for two traits:

- **Selection index:** candidates are ranked by total merit, \(I=w_1t_1+w_2t_2\).
- **Independent culling:** candidates must pass a threshold for each trait.
- **Tandem-type selection:** candidates are selected sequentially for different traits.

The simulated examples illustrate that selecting extreme individuals for one trait can produce parents with poor performance for another trait. The selection index produced a more balanced increase in total merit. The presentation reports illustrative gains of approximately 70% versus 50% for Trait 1, 66% versus 38% for Trait 2, and 57% versus 53% for total merit, depending on the compared method.

The biological message is that the best parent is not necessarily the most extreme genotype for any single trait. A useful parent should combine acceptable performance across the complete target profile.

### Reconstructing the S3-to-S4 procedure

The selection step from Stage 3 to Stage 4 is represented as a sequence of reduction and selection operations. Each operation includes:

- selection unit,
- operation type,
- trait,
- comparison value,
- direction of selection, and
- whether the value is used as a culling rule.

The example procedure uses performance under several environments or management conditions, abbreviated as OPTY, MDY, LNY, and MLNY, together with anthesis date and grain moisture. Some traits, such as Striga and MLN, were noted as ignored in the illustrated index replacement.

This algorithmic representation is important because it distinguishes two kinds of decisions:

1. **non-compensatory constraints**, for which failure cannot be offset by superiority in another trait; and
2. **compensatory quantitative objectives**, for which an index can rank total merit.

### Retrospective index weights

The presentation uses retrospective weights based on the relationship

\[
\mathbf{b}=\mathbf{P}^{-1}\mathbf{s},
\]

where \(\mathbf{s}\) is the vector of realised selection differentials from the breeder's historical selections and \(\mathbf{P}\) is the phenotypic variance-covariance matrix of the selection criteria.

For the S3-to-S4 example, the reported coefficients were approximately:

| Trait | Retrospective coefficient |
|---|---:|
| Optimum-environment performance | 1.4254 |
| Managed-drought performance | 0.0718 |
| Low-nitrogen performance | 0.1703 |
| Anthesis date | 0.1206 |

These coefficients reconstruct the implicit linear preference expressed by the current selections. They should not be interpreted as universal biological or economic values. They depend on the observed covariance matrix and the breeder's past decisions.

### Comparison with the current procedure

The presentation compares selections from the existing "Yoseph" approach with the index. Many selected materials were common, but each method also selected unique families or lines. This is expected even when total merit is similar because small changes in coefficients can alter rankings near the selection threshold.

The supplementary comparison on page 21 reports the following selection differentials:

| Method | Opt. | MD | LN | AD | MLN | Striga | Grain moisture | Total merit |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Current approach | 1.09 | 0.35 | 0.48 | -0.39 | 0.07 | 0.30 | -0.31 | 1.61 |
| Index | 1.07 | 0.41 | 0.49 | 0.53 | 0.04 | 0.23 | 0.29 | 1.70 |

The index slightly increased total merit in this comparison. However, the sign and desirability of anthesis date and grain moisture must be aligned explicitly before deployment because the table alone does not indicate whether positive or negative values represent improvement.

### Recycling-stage index

A second application focuses on recycling parents using general combining ability under optimum and drought conditions. The index is calculated among families and then within families. Traits include grain yield, anthesis date, anthesis-silking interval, ear-related traits, lodging or breakage traits, ears per plant, and moisture under different conditions.

The derived coefficients show strong positive emphasis on GCA for grain yield and substantial negative emphasis on some undesirable traits. The presentation reports nearly identical total merit for the current approach and the index, approximately 3.1327 and 3.1345, respectively. Nevertheless, the candidate lists differ, demonstrating that an index can reproduce overall selection intensity while changing which lines are retained.

### Conclusions from the presentation

The author concludes that:

- a traditional selection procedure can be translated into an index;
- an index provides consistency across cycles and users;
- initial weights can be obtained retrospectively from historical decisions;
- the index should be adopted and subsequently refined; and
- future indices should be defined by **program, stage, and product profile**.

The proposed next step is to obtain data from Stage 1 and Stage 2, particularly for recycling, across programs such as Latin America, East Africa, and South Asia.

### Strengths

- Starts from the real breeding process rather than an abstract index.
- Makes hidden selection rules explicit and auditable.
- Separates hard culling constraints from compensatory traits.
- Uses historical selections to derive an initial index.
- Compares both total merit and overlap in selected candidates.
- Recognises the need for program-, stage-, and product-profile-specific indices.

### Limitations and implementation cautions

- The presentation is not a peer-reviewed study and provides limited information on sample size, uncertainty, and statistical validation.
- Trait abbreviations and the direction of improvement are not defined consistently on all pages.
- Some selection differentials are presented without standard errors or confidence intervals.
- Retrospective weights reproduce past behaviour; they do not demonstrate that the past behaviour was optimal.
- A linear index cannot safely replace mandatory disease, quality, or maturity thresholds unless those constraints are retained separately.
- The covariance matrix and trait scaling must be updated when the population, target environment, or stage changes.

### Operational interpretation

The main value of the presentation is methodological. It shows how a breeder's informal sequence of filters can be converted into a reproducible decision system. The recommended practice is a hybrid system: use explicit culling for non-negotiable requirements and use a selection index for traits where compensation is biologically and commercially acceptable.

---

## Covarrubias-Pazaran: Bringing a selection index into the IRRI programs

### Document information

- **Author:** Giovanny E. Covarrubias-Pazaran
- **Document type:** Technical presentation, 13 pages
- **Source file:** `Giovanny Selection index in Rice.pdf`

### Objective

The presentation describes the conversion of a rice parent-selection procedure into a formal index. The example combines yield, zinc concentration, major disease-resistance genes, plant height, and flowering time. The approach is retrospective: first recreate the breeder's existing selections, then calculate an index that reflects those decisions, compare outcomes, and refine the weights.

The stated steps are:

1. identify the target breeding pipeline;
2. identify the stage;
3. specify whether the purpose is parent selection or product advancement;
4. translate the current procedure into an algorithm;
5. calculate retrospective weights; and
6. fine-tune the weights.

The approach slide contains "Hard White" market-segment labels that appear inconsistent with the rice-specific traits and IRRI context. These labels likely originated from a reused template and should be corrected before the presentation is used as formal documentation.

### Existing selection algorithm

The current procedure is represented as 13 steps. The selection unit is usually the designation, with a final parentage-level operation. The algorithm includes:

- preliminary reduction based on relative yield and zinc information;
- selection for resistance alleles at **Xa21, Xa5, Pi54, and Pita**;
- minimum or top-ranking requirements for yield breeding value and zinc breeding value;
- upper limits for plant height and flowering time; and
- random retention of three entries per parentage to maintain family representation.

The initial reduction retained 891 candidates. Subsequent steps selected smaller sets associated with specific resistance genes and high yield or zinc. A final parentage step retained 18 selections.

This representation clarifies that the process combines:

- hard constraints or categorical requirements for major genes,
- quantitative ranking for yield and zinc,
- upper limits for plant height and flowering, and
- family-balance rules.

Only the compensatory quantitative parts are natural candidates for replacement by a linear index.

### Critique of extreme-value selection

Pages 4-6 argue that the current procedure gives excessive emphasis to single-trait transgressive individuals. Separate selection for the highest yield and highest zinc can retain candidates that are exceptional for one trait but mediocre for the remaining target profile.

The presentation contrasts:

- index selection for total merit,
- independent culling, and
- tandem-type selection.

The intended message is that parent selection should favour balanced total merit rather than independent extremes, except when a trait is a strict requirement.

### Retrospective weights

Historical selections were used to derive regional coefficients for yield, zinc, and resistance-gene indicators. Reported coefficients were:

| Trait | Bangladesh | ESA | India | Philippines |
|---|---:|---:|---:|---:|
| YLD_BV | 2.905 | 4.764 | 3.371 | 2.015 |
| ZNC_BV | 0.457 | 0.629 | 0.460 | 1.199 |
| Xa21 indicator | -3.148 | 1.003 | -3.533 | -0.844 |
| Xa5 indicator | 0.751 | -0.436 | 0.915 | -0.492 |
| Pita indicator | -0.496 | -1.075 | -1.153 | -1.358 |
| Pi54 indicator | -1.051 | -5.222 | -2.889 | -1.348 |

The different signs and magnitudes show that historical selection priorities varied strongly among regions. They also illustrate why coefficients should not be interpreted directly as biological importance. Coding of resistance alleles, trait orientation, covariance, and the previous sequence of culling steps all affect the retrospective coefficients.

### Comparison of selection methods

For the Philippines example, the presentation compares the breeder's selections, an algorithmic recreation, and index selections.

| Method | Final gain |
|---|---:|
| Josh selections | 1.919 |
| Algorithmic recreation | 1.722 |
| Index selections | 3.361 |

The index produced a markedly larger reported total gain, mainly through increased zinc response while maintaining favourable changes in the resistance indicators.

Across all regions, the reported final gains were:

| Region | Historical selection | Index selection |
|---|---:|---:|
| Bangladesh | 2.477 | 5.457 |
| ESA | 3.481 | 7.018 |
| India | 3.178 | 5.617 |
| Philippines | 1.919 | 3.361 |

These values are internal merit scores rather than realised genetic gains in a subsequent generation. They demonstrate consistency with the index objective but require prospective validation before being interpreted as achieved breeding progress.

### Fine-tuning the objective

The presentation emphasises that retrospective weights are only a starting point. For example, the Philippines yield coefficient was doubled. Under the revised objective, reported final gain increased from 2.230 for the historical selection to 3.772 for the revised index.

This step is essential because a retrospective index answers the question, "What linear rule best approximates the previous decisions?" It does not answer, "What decisions should the program make in the future?" Fine-tuning converts the reconstructed historical preference into an explicit prospective breeding objective.

### Family-balanced index selection

Page 12 compares crosses selected by the breeder and the index and then applies selection both across and within families. This protects family representation while still ranking individuals by total merit.

With family balancing, the index retained a smaller but consistent advantage:

| Region | Historical final gain | Family-balanced index final gain |
|---|---:|---:|
| Bangladesh | 2.477 | 2.833 |
| ESA | 3.481 | 4.038 |
| India | 3.178 | 3.468 |
| Philippines | 1.997 | 2.069 |

The reduced difference is expected because family constraints lower the freedom of the index to select only the highest-scoring individuals. However, family balance can preserve diversity and reduce the risk of overusing a small number of crosses.

### Strengths

- Documents the existing selection process in a reproducible form.
- Integrates major genes, quantitative breeding values, maturity, plant stature, and family representation.
- Provides regional rather than universal coefficients.
- Distinguishes reconstruction of past decisions from fine-tuning future objectives.
- Demonstrates both unrestricted and family-balanced index selection.

### Limitations and cautions

- The presentation does not provide the full covariance matrices, sample sizes, or uncertainty of the reported gains.
- The "final gain" is an index-derived score and not observed response in progeny or later-stage trials.
- Major-gene indicators are more safely treated as constraints when resistance is mandatory.
- Negative or positive coefficients cannot be interpreted without confirming allele coding and trait direction.
- The approach slide contains pipeline labels that appear unrelated to rice.
- The historical selections may include non-linear judgement, pedigree information, and operational knowledge that a linear index cannot recover completely.

### Practical interpretation

The presentation provides a useful template for implementing an index in a real breeding program:

1. formalise the current rules;
2. retain hard requirements as culling constraints;
3. calculate retrospective coefficients for quantitative criteria;
4. compare candidate overlap and response;
5. modify the desired emphasis explicitly; and
6. apply family constraints when diversity and cross representation are important.

Its main contribution is not a universal rice index but a procedure for converting breeder judgement into a transparent and refinable decision rule.

---

## Covarrubias-Pazaran (2021): Practical implementation of selection indices

### Document information

- **Author:** Giovanny Covarrubias-Pazaran
- **Date:** 18 June 2021
- **Document type:** Four-page practical guideline
- **Source file:** `Guideline_Practical Implementation of Selection Indices.pdf`

### Purpose

The guideline provides a concise sequence for introducing selection indices into an operational plant-breeding pipeline. It assumes that the program has:

- a clearly defined target market or product profile;
- a breeding scheme organised into stages;
- a manageable set of selection traits at each stage; and
- access to biometric support.

The document focuses on implementation rather than on a full derivation of selection-index theory.

### Step 1: Map traits across the breeding pipeline

The program should create a stage-by-trait table recording:

- the stage at which each trait is used;
- the surrogate used to represent merit, such as visual score or BLUP; and
- the number of environments supporting the estimate.

The example spans crossing block, nursery, Stage 1, Stage 2, Stage 3, and on-farm testing. This map identifies where information is available and where the evidence supporting selection becomes more reliable.

### Step 2: Select the index stage

A specific stage should be selected for index implementation. The chosen stage should have enough candidates to benefit from ranking and enough information to estimate merit for all important traits. Introducing an index too early can be ineffective if key traits have not been measured; introducing it too late can reduce the opportunity to save resources.

### Step 3: Close trait-information gaps

All traits included in the index must have a quantitative selection criterion at the chosen stage. The guideline recommends modifying the breeding scheme when a trait is represented only by visual judgement or is not measured at all.

This is an important operational requirement. A multi-trait index cannot correct for missing trait information unless a formal prediction model is used to supply it.

### Step 4: Assemble adjusted means

For each candidate and trait, the program should obtain one adjusted value, such as a BLUE, BLUP, or other mixed-model estimate. The data matrix should have candidates in rows and traits in columns.

The use of adjusted means reduces confounding with field design, environment, and incomplete testing. However, the guideline does not discuss prediction-error variance. In a more advanced implementation, differences in reliability among candidates and traits should be considered.

### Step 5: Standardise traits

Adjusted values should be standardised so that traits are expressed in standard-deviation units. Standardisation has two main benefits:

- traits measured in different units become directly comparable; and
- desired responses can be expressed as the number of standard deviations from the population mean.

For a standardised trait, a desired differential of 1 represents selection toward a mean one standard deviation above the population mean, subject to selection intensity and correlations.

### Step 6: Specify desired differentials and derive weights

The guideline recommends beginning with a desired differential of one standard deviation for each trait. The desired differential is the breeder's decision. The resulting index coefficients are mathematical quantities derived from the covariance structure and should not be interpreted directly.

When traits are correlated, a coefficient can differ substantially from the desired differential and may even have a counterintuitive sign. This does not necessarily mean that the index opposes the trait. It can reflect redundancy or compensation among correlated criteria.

The guideline also proposes a retrospective initialisation:

\[
\mathbf{b}=\mathbf{P}^{-1}\mathbf{s},
\]

where \(\mathbf{s}\) contains the historical selection differentials and \(\mathbf{P}\) is the covariance matrix of the selection criteria. Historical differentials can be converted to standard-deviation units and used as the starting objective.

### Feasibility assessment

The guideline recommends using a tool such as DESIRE to determine whether the requested response is feasible given:

- genetic correlations,
- selection intensity, and
- the number of traits.

A breeder may request one standard deviation of gain for every trait, but antagonistic correlations and limited selection intensity can make this impossible. The desired response should therefore be treated as a scenario to be tested rather than as a guaranteed result.

### Covariance matrix used in practice

The document notes that the true genetic covariance matrix is often unavailable. It proposes using the covariance matrix of across-environment adjusted means as a practical surrogate. This is operationally convenient, but it should be recognised as an approximation. The covariance of adjusted means includes estimation error and may differ from the additive genetic covariance relevant to response in the next generation.

### Step 7: Calculate the index and select candidates

The derived coefficients are applied to the standardised trait values:

\[
I_i=\sum_j b_j z_{ij},
\]

where \(I_i\) is the index score for candidate \(i\), \(b_j\) is the coefficient for trait \(j\), and \(z_{ij}\) is the candidate's standardised adjusted value.

Candidates are ranked by the index. The relationship between the mixed model used to obtain adjusted means and the index should be checked to avoid duplicating adjustments or weighting the same information twice.

### Step 8: Recalibrate every selection season

The index should be recalculated and fine-tuned each season. Population means, variances, correlations, product profiles, and operational constraints change over time. A static index can therefore become misaligned with the current breeding objective.

The original document contains a minor numbering inconsistency in the final instruction, referring to repeating Steps 7 and 8. The intended meaning is clear: apply the index, inspect the outcome, and update it in subsequent cycles.

### Main principles

The guideline establishes several important principles:

- Map the complete breeding process before calculating weights.
- Use quantitative, adjusted trait values at the selected stage.
- Standardise traits before expressing desired responses.
- Treat desired differentials as breeder decisions and coefficients as derived quantities.
- Evaluate feasibility under the observed covariance structure.
- Re-estimate the index when the population or objective changes.

### Strengths

- Short and operationally focused.
- Compatible with stage-gate breeding pipelines.
- Encourages standardisation and mixed-model adjusted means.
- Correctly warns against interpreting coefficients as direct importance weights.
- Includes both prospective desired response and retrospective reconstruction.
- Recognises that an index must be recalibrated.

### Limitations and recommended extensions

A complete implementation should add:

- explicit trait direction and units;
- hard culling thresholds for non-compensatory requirements;
- treatment of missing values and unequal reliabilities;
- validation of realised selection differentials;
- sensitivity analysis for covariance-matrix uncertainty;
- constraints on inbreeding or family representation; and
- prospective validation in later-stage or next-cycle material.

### Practical value

The guideline is a useful minimum protocol for moving from informal multi-trait selection to a transparent index. Its most important message is that the breeder should define the desired response, not manually interpret or tune raw index coefficients without considering covariance.

---

## Rahimi and Debnath (2023): SAS and R codes for optimum, base, and desired-gain indices

### Bibliographic information

- **Authors:** Mehdi Rahimi and Sandip Debnath
- **Journal:** *Scientific Reports* 13:18977
- **Year:** 2023
- **DOI:** 10.1038/s41598-023-46368-6
- **Document type:** Methods and software application article
- **Source file:** `Rahimi et al 2023.pdf`

### Objective

The paper provides relatively simple SAS and R implementations for three selection-index methods:

1. optimum or Smith-Hazel-type index,
2. base index, and
3. Pešek-Baker desired-gain index.

The codes calculate index coefficients, candidate scores, expected gain, and several criteria for comparing alternative indices and economic-weight definitions. The intended audience is plant and animal breeders who have trait measurements and variance-covariance matrices but do not use specialised selection-index software.

### Index formulations

#### Optimum index

The linear index is

\[
I=\sum_i b_iX_i,
\]

with coefficients

\[
\mathbf{b}=\mathbf{P}^{-1}\mathbf{G}\mathbf{a},
\]

where \(\mathbf{P}\) is the phenotypic variance-covariance matrix, \(\mathbf{G}\) is the genetic variance-covariance matrix, and \(\mathbf{a}\) is the vector of economic values.

This formulation seeks the linear combination of phenotypic criteria that is maximally correlated with the aggregate breeding objective \(H=\mathbf{a}'\mathbf{g}\).

#### Base index

The base index uses the economic values directly as coefficients:

\[
I=\sum_i a_iX_i, \qquad \mathbf{b}=\mathbf{a}.
\]

It does not require estimates of \(\mathbf{P}\) and \(\mathbf{G}\). The method is simpler but does not adjust coefficients for differences in variance, reliability, or correlation among traits.

#### Pešek-Baker index

The paper implements the desired-gain form as

\[
\mathbf{b}=\mathbf{G}^{-1}\mathbf{d},
\]

where \(\mathbf{d}\) is a vector of desired genetic gains. In the example, \(\mathbf{d}\) was set equal to the genetic standard deviation of each trait, although the code allows the user to enter another desired-gain vector.

### Evaluation criteria

The codes calculate several diagnostics:

- **\(R_{HI}\):** correlation between index and aggregate breeding value;
- **\(\Delta H\):** expected gain in aggregate merit;
- **\(\Delta_j\):** expected response for each trait;
- **RE:** relative efficiency compared with direct selection for the main trait, grain yield in the example; and
- **\(CV_I\):** phenotypic coefficient of variation of the index.

The use of multiple criteria is valuable because an index may rank candidates consistently but still produce poor response for the primary trait or an undesirable correlated change.

### Input and implementation

The software uses four main inputs:

- **X:** phenotype matrix, candidates by traits;
- **P:** phenotypic variance-covariance matrix;
- **G:** genotypic variance-covariance matrix; and
- **a:** economic-value vectors.

The SAS implementation uses PROC IML, while the R version reads the inputs from files and writes results to multiple Excel sheets. Parameters such as number of genotypes, number of traits, selection intensity, main-trait variance, heritability, and main-trait position must be specified.

The paper assumes a selection intensity of 10%, corresponding to a standardised selection differential of approximately 1.76 in the worked example.

### Worked maize example

The example used 28 maize inbred lines evaluated for seven traits in a randomised complete block design with three replications. The traits were:

- plant height,
- number of grains,
- number of rows,
- row length,
- leaf length,
- 100-grain weight, and
- yield.

Three economic-weight schemes were compared for the optimum and base indices:

1. equal weights of 1;
2. correlations of each trait with yield; and
3. standardised coefficients from stepwise regression on yield.

The Pešek-Baker desired gains were the square roots of the diagonal elements of the genetic covariance matrix.

### Reported index equations

For equal economic weights, the paper reports:

\[
I_{opt}=0.963x_1+2.150x_2+1.119x_3-3.851x_4+1.610x_5+0.976x_6-1.411x_7,
\]

\[
I_{base}=x_1+x_2+x_3+x_4+x_5+x_6+x_7,
\]

and

\[
I_{PB}=0.011x_1+5.249x_2+0.552x_3-3.951x_4+3.547x_5+0.063x_6-9.757x_7.
\]

These equations are dataset specific and should not be transferred to another population.

### Main numerical results

Under equal weights:

| Index | \(R_{HI}\) | \(\Delta H\) | RE |
|---|---:|---:|---:|
| Optimum | 0.9887 | 125.7762 | 0.5504 |
| Base | 0.9885 | 128.6944 | 0.5555 |
| Pešek-Baker | 0.0018 | 9.2231 | 0.1986 |

The optimum and base indices were almost identical in candidate ranking. Their score correlation was 0.99979. Correlations of the Pešek-Baker scores with the base and optimum indices were only approximately 0.315 and 0.325.

The top five candidates differed accordingly. The optimum index selected genotypes 2, 4, 19, 22, and 6, whereas the base index selected 1, 4, 19, 22, and 6. The Pešek-Baker index selected a substantially different group: 3, 17, 8, 16, and 13.

The base index with regression coefficients as economic weights had the strongest values for several comparison criteria and was identified by the authors as the preferred index for this dataset.

### Interpretation of relative efficiency

All RE values were below 1. Therefore, none of the multi-trait indices produced a larger expected response for grain yield alone than direct selection for yield. This is not a failure of an index if the purpose is simultaneous improvement of several traits, but it is an important result. A multi-trait index trades some direct yield response for changes in other components of aggregate merit.

The article's conclusion that the base index was preferable applies to the selected objective and data. It should not be interpreted as a general superiority of the base index over the optimum index.

### Strengths

- Provides transparent R and SAS code rather than a black-box application.
- Accepts any number of candidates and traits.
- Compares several economic-weight strategies.
- Produces per-trait response and aggregate diagnostics.
- Includes a complete worked example.
- Stores R outputs in structured Excel sheets.

### Limitations and cautions

- The empirical demonstration is limited to 28 maize lines and seven traits.
- Broad-sense genotypic variance is used in parts of the analysis, whereas response to selection in a recurrent breeding program is primarily determined by transmissible additive effects.
- Correlations and stepwise-regression coefficients are statistical weights, not necessarily economic values or product-profile utilities.
- Stepwise regression can be unstable under collinearity and small sample size.
- The near equality of base and optimum indices is specific to this covariance structure and weight vector.
- Expected gains were not validated in progeny or a subsequent cycle.
- The very low \(R_{HI}\) for the implemented Pešek-Baker example suggests that the chosen desired-gain vector was poorly aligned with the defined aggregate objective.
- Input matrices must be positive definite and consistently ordered; the paper gives limited discussion of numerical diagnostics when inversion is unstable.

### Practical use

The code is useful for teaching, prototyping, and comparing alternative objectives. A robust operational application should:

1. estimate additive genetic and residual covariance matrices with an appropriate mixed model;
2. verify matrix ordering and conditioning;
3. define economic values or desired gains from a target product profile;
4. standardise traits when units differ substantially;
5. evaluate per-trait response and aggregate response;
6. test sensitivity to alternative weights; and
7. validate the selected index prospectively.

The main contribution is accessibility. The paper enables breeders to calculate and compare selection indices in R or SAS with familiar tabular inputs.

---

## Guimarães et al. (2021): Index selection in a rice recurrent-selection population

### Bibliographic information

- **Authors:** Paulo Henrique Ramos Guimarães, Patrícia Guimarães Santos Melo, Antônio Carlos Centeno Cordeiro, Paula Pereira Torga, Paulo Hideo Nakano Rangel, and Adriano Pereira de Castro
- **Journal:** *Euphytica* 217:95
- **Year:** 2021
- **DOI:** 10.1007/s10681-021-02819-7
- **Document type:** Empirical multi-environment selection study
- **Source file:** `s10681-021-02819-7.pdf`

### Objectives

The study evaluated a recurrent-selection population of lowland rice to:

1. estimate genetic parameters for yield, maturity, plant stature, and disease traits;
2. compare direct and indirect selection with three multi-trait indices; and
3. assess the sensitivity of Smith-Hazel and Tai indices to alternative economic weights.

The practical question was whether simultaneous selection could improve grain yield while reducing plant height, days to flowering, and disease incidence.

### Genetic material and environments

The material comprised:

- 198 \(S_{0:2}\) progenies from the second cycle of the CNA 12 lowland-rice recurrent-selection population; and
- four commercial cultivars used as controls.

CNA 12 was created as a broad-based population with resistance to *Magnaporthe oryzae*. Trials were conducted in three Brazilian locations across two seasons: Goianira, Formoso do Araguaia, and Cantá. The environments represented contrasting tropical and equatorial conditions.

An augmented block design was used. Most progenies were unreplicated within an environment, while controls were repeated; one Goianira trial included two progeny replicates.

### Traits

Six traits were analysed:

- grain yield (GY, kg ha\(^{-1}\)),
- plant height (PH, cm),
- days to flowering (DF),
- panicle blast (PB),
- leaf scald (Ls), and
- grain discoloration (Gd).

Disease traits were scored visually from 0, indicating no incidence, to 9, indicating severe incidence. For the breeding objective, positive response was desired for GY and negative response for PH, DF, PB, Ls, and Gd.

### Statistical analysis

A joint analysis of variance included environment, block within environment, genotype type, genotype within type, type-by-environment interaction, and genotype-by-environment interaction. Progenies were random and controls fixed.

The study estimated:

- genetic, environmental, and phenotypic variances;
- broad-sense heritability;
- genotypic coefficient of variation;
- b-variation index; and
- contribution of each source of variation to total sum of squares.

A MANOVA was used to obtain genetic covariance information for selection indices.

### Selection methods

The top 25% of progenies were selected with:

1. **Direct selection (DS)** for each trait separately.
2. **Indirect selection (IS)** based on grain yield.
3. **Smith-Hazel index (SH)** with coefficients \(\mathbf{b}=\mathbf{P}^{-1}\mathbf{G}\mathbf{a}\).
4. **Tai desired-gain index (TA)** with coefficients derived from the inverse genetic covariance matrix and a desired-gain vector.
5. **Mulamba-Mock rank-sum index (MM)**, calculated by ranking candidates for each trait in the desired direction and summing ranks.

For SH and TA, three economic-weight definitions were tested:

- genetic standard deviation;
- b-variation index; and
- a random vector. The reported random weights were -1, 0, 1, 1, 1, and 1 for GY, PH, DF, PB, Ls, and Gd, respectively, according to the coding used in the paper.

The MM index was applied without economic weights.

### Genetic variability

The progenies showed significant variation for all traits and significant genotype-by-environment interaction. Environmental effects explained a large proportion of variation for GY, DF, Ls, and Gd.

Broad-sense heritability estimates were:

| Trait | Heritability (%) |
|---|---:|
| DF | 93.95 |
| Gd | 78.24 |
| PH | 72.81 |
| Ls | 53.87 |
| GY | 48.42 |
| PB | 40.65 |

The high values for flowering, grain discoloration, and plant height suggested favourable conditions for phenotypic selection. Grain yield and panicle blast had moderate heritability and stronger environmental influence.

Genotypic coefficients of variation ranged from 3.79% to 22.22%. Grain discoloration had the highest CVg, followed by panicle blast and leaf scald. Plant height and flowering had lower CVg but high heritability.

### Direct and indirect response

Direct selection produced the largest response for each individual trait:

| Trait | Direct-selection gain (%) |
|---|---:|
| GY | +7.12 |
| PH | -4.30 |
| DF | -7.50 |
| PB | -12.83 |
| Ls | -15.66 |
| Gd | -24.60 |

However, indirect selection for grain yield produced small undesirable increases in plant height and flowering time. This illustrates the limitation of single-trait selection when correlated traits are part of the target product profile.

### Mulamba-Mock index

The MM index produced:

| Trait | Predicted gain (%) |
|---|---:|
| GY | +4.65 |
| PH | -0.50 |
| DF | -0.73 |
| PB | -9.12 |
| Ls | -6.58 |
| Gd | -18.81 |

Although its yield gain was smaller than direct selection, MM moved all six traits in the desired direction. The selected progenies had higher average yield and lower averages for plant height, flowering, and all disease scores than the base population.

The authors considered this the best overall balance among the evaluated methods.

### Smith-Hazel index

With genetic-standard-deviation weights, SH produced approximately:

- GY +6.98%,
- DF -0.34%,
- PB -6.20%,
- Ls -4.53%, and
- Gd -7.60%.

With b-variation weights, SH produced:

- GY +6.94%,
- DF -0.19%,
- PB -6.58%,
- Ls -5.07%, and
- Gd -8.17%.

Neither configuration produced meaningful improvement in plant height. Random weights performed poorly and generated unfavourable changes for most traits.

### Tai index

Tai with genetic-standard-deviation weights produced responses similar to SH for several traits, including approximately +6.93% for yield, -8.48% for grain discoloration, -3.98% for leaf scald, and a small decrease in flowering. The direction reported for panicle blast in the article requires careful reading because the text lists a positive value while discussing desirable disease reduction.

Tai with b-variation or random weights produced weak or unfavourable response for several traits. The result reinforces the sensitivity of parametric selection indices to the chosen economic or desired-gain vector.

### Main conclusion

The CNA 12 population retained sufficient diversity for continued improvement. Direct selection maximised response for an individual trait but could alter correlated traits unfavourably. The MM rank-sum index gave the most balanced multi-trait response and was easier to implement because it did not require economic weights or inversion of covariance matrices.

The study therefore recommended MM for simultaneous selection in this population.

### Strengths

- Empirical evaluation in multiple environments.
- Inclusion of grain yield, adaptation, plant architecture, and three diseases.
- Direct comparison of several index families.
- Quantification of sensitivity to alternative economic weights.
- Clear presentation of predicted per-trait responses.
- Evaluation in a recurrent-selection population with practical breeding relevance.

### Limitations and cautions

- The response values are predicted gains, not realised response in the next cycle.
- Broad-sense heritability includes non-additive variation and may overstate transmissible response.
- The augmented design and unreplicated progenies increase dependence on model assumptions.
- Strong GxE indicates that a single across-environment index may not represent all target environments equally.
- The MM index ignores covariance, reliability, and the magnitude of differences after ranks are assigned.
- Economic-weight comparisons are sensitive to coding and trait direction; random weights were not biologically calibrated.
- The study selected 25%, so relative performance may change under stronger or weaker selection intensity.

### Practical relevance

The paper supports a pragmatic hierarchy:

- use direct selection when one trait is overwhelmingly important and correlated effects are acceptable;
- use Smith-Hazel or desired-gain indices when reliable covariance matrices and meaningful objectives are available; and
- use a rank-sum index when weights are unavailable and a balanced directional response is more important than theoretical optimality.

For operational use, the MM index should still be validated for realised gain, family diversity, and stability across environments.

---
