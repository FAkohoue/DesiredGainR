# Vieira, Nogueira, and Fritsche-Neto (2025): Optimizing quantitative-trait selection using simulation

## Bibliographic information

- **Authors:** Rafael Augusto Vieira, Ana Paula Oliveira Nogueira, and Roberto Fritsche-Neto
- **Journal:** *Frontiers in Plant Science* 16:1495662
- **Year:** 2025
- **DOI:** 10.3389/fpls.2025.1495662
- **Document type:** Narrative review
- **Source file:** `fpls-16-1495662.pdf`

## Purpose and scope

This review synthesises simulation studies that evaluated the optimisation of quantitative-trait selection in plant breeding. The central premise is that simulation provides an intermediate layer between quantitative-genetic theory and expensive field implementation. A simulation can represent genetic architecture, crossing, recombination, population size, selection intensity, heritability, genotype-by-environment interaction, genotyping cost, and breeding-cycle duration before a breeding program commits resources.

The review is broad. It covers classical phenotypic selection, pre-breeding, marker-assisted selection, QTL introgression, genomic selection, multi-trait selection, diversity management, inbreeding and heterosis, cost allocation, genotype-by-environment interaction, and software for breeding simulation.

## Role of simulation in plant breeding

The authors distinguish deterministic and stochastic approaches. Deterministic models use equations based on heritability, selection intensity, and accuracy to predict expected response. These models are efficient for broad comparisons but cannot fully represent crossing, recombination, segregation, linkage, and recurrent selection. Stochastic simulations generate individual genotypes and phenotypes and are therefore better suited to reproducing complete breeding pipelines.

Simulation is presented as a decision-support tool rather than a replacement for empirical testing. Its main functions are to:

- compare breeding strategies under controlled assumptions;
- quantify trade-offs among gain, diversity, cost, and time;
- identify bottlenecks before changing a production pipeline;
- test sensitivity to genetic architecture and heritability;
- optimise population size, number of parents, and selection timing; and
- evaluate long-term consequences that cannot be observed rapidly in field experiments.

The main limitation is model dependence. A simulation can be precise but wrong if its assumptions about QTL number, effect distribution, recombination, GxE, cost, or selection behaviour do not represent the actual program.

## Diversity, pre-breeding, and parental management

A recurring conclusion across simulation studies is that maximising short-term mean performance can accelerate the loss of genetic diversity. Fewer parents and intense selection can increase immediate gain, whereas larger parental sets and controlled coancestry preserve long-term potential.

The review recommends combining strategies rather than relying on a single selection rule. Overlapping cohorts, strategic recycling, optimal contribution methods, and the introduction of diverse germplasm can improve long-term gain. The optimal number of parents depends on the planning horizon:

- fewer parents can increase short-term gain;
- more parents and more crosses are preferable when long-term improvement and risk control are important.

Pre-breeding simulations indicate that exotic or landrace material should be introduced with deliberate background recovery. Multiple conversion versions and residual donor segments must be managed because aggressive recovery can discard useful diversity, whereas insufficient recovery can reduce agronomic adaptation.

## Phenotypic selection and classical breeding

Simulation studies of self-pollinated crops show that the value of early selection depends on trait heritability, genetic architecture, and the reliability of early-generation phenotypes. Early selection can be efficient for highly heritable traits but can create bias for low-heritability traits and traits strongly affected by GxE.

The review discusses bulk, pedigree, single-seed descent, doubled-haploid, and backcross schemes. No scheme is universally optimal. Population size, number of families, progeny per family, and generation at which selection is initiated must be aligned with the breeding objective. Larger populations generally increase the probability of recovering favourable recombinants, but the benefit is conditional on adequate germplasm, a clear target product profile, and sufficient resources to measure the population accurately.

For introgression, simulations support combining foreground selection for the target locus with background selection for recipient-genome recovery. The number and spacing of markers should be sufficient to track linkage around the target, but merely increasing marker density does not guarantee efficiency if population size and recombination are limiting.

## Marker-assisted selection and QTL management

Marker-assisted selection provides the largest advantage when:

- one or a few loci explain a meaningful proportion of variation;
- marker-QTL linkage is strong and stable;
- the selected alleles have relatively large effects;
- the population is sufficiently large to generate recombination; and
- foreground and background selection are balanced.

Simulations commonly show rapid early gains from MAS, followed by diminishing advantage after several cycles. This occurs because large-effect alleles are fixed rapidly, while residual polygenic variation is not fully captured. MAS is therefore more suitable for disease resistance, quality loci, and major adaptation genes than for highly polygenic traits such as yield.

The review emphasises that low heritability does not directly invalidate MAS, because the selection is marker based. However, low heritability reduces the reliability of QTL discovery and effect estimation, thereby weakening the marker-trait relationship on which MAS depends.

## Genomic selection

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

## Multi-trait selection

Multi-trait genomic prediction is particularly useful when a low-heritability target trait is genetically correlated with a higher-heritability secondary trait. Under favourable correlation structures, the secondary trait contributes information and increases prediction accuracy for the target.

The value of multi-trait analysis depends on:

- the magnitude and direction of genetic correlation;
- the heritability of each trait;
- overlap or imbalance in trait records;
- the stage at which each trait is measured; and
- whether the correlation remains stable across environments and cycles.

Selection indices are discussed as a means to translate multi-trait predictions into a single decision criterion. Simulations are useful because they can compare index weights, desired gains, constraints, and culling thresholds while measuring both total response and undesirable correlated changes.

## Inbreeding, heterosis, and epistasis

The review shows that breeding optimisation must account for non-additive processes when they influence the target product. In hybrid breeding, heterotic-group structure and specific combining ability can alter the optimal choice of parents. In recurrent selection, aggressive genomic selection can increase inbreeding and reduce long-term gain. Optimal contribution and mate-allocation methods can control coancestry while retaining high-value alleles.

Epistasis may improve prediction in some simulated architectures, but its practical benefit is inconsistent and often depends on sample size, relatedness, and the stability of interactions across generations. Additive models remain robust for many operational decisions because additive effects determine transmissible response.

## Cost analysis

Simulation enables the breeder to optimise resource allocation rather than accuracy alone. The review describes examples in which marker-assisted selection reduced cost, and genomic selection increased gain when genotyping was substituted for expensive or slow phenotyping.

Low-density genotyping followed by imputation is highlighted as an important cost-saving strategy. One cited study reported cost reductions of up to 87% with minimal loss in prediction accuracy, and as few as 50 segregating markers per genome were sufficient under the evaluated conditions. These numbers are scenario dependent, but the general conclusion is robust: once marker density is sufficient to recover relationships and linkage information, population size and training relevance can be more valuable than further density increases.

## Genotype-by-environment interaction

GxE can bias estimated breeding values and reduce transferability of selection decisions. Simulation studies indicate that environment-specific models, multi-environment genomic prediction, reaction-norm approaches, and environment-specific selection indices can improve response when target environments differ substantially.

The review notes that genotype-by-year interaction is still represented inadequately in many simulations. This is important because year-to-year climatic variation may be less predictable than location differences. Future simulations should represent changing frequencies of drought, heat, excess rainfall, and disease pressure rather than assuming stationary environments.

## Software tools

The review lists several simulation platforms:

- **AlphaSimR:** flexible simulation of genetic architectures, populations, crossing, and selection;
- **ADAM:** structured populations and inheritance tracking;
- **synbreed:** genomic prediction and SNP-based breeding analysis;
- **MoBPS:** modular simulation of livestock and plant breeding programs, including genomic schemes;
- **QU-LINE:** quantitative-trait and breeding-line simulation, including additive and epistatic effects; and
- **PedigreeSim:** meiosis, recombination, and pedigree-based inheritance in structured populations.

Software selection should follow the biological and operational question. A highly detailed simulator is not automatically preferable when the necessary parameters cannot be estimated reliably.

## Integrated recommendations

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

## Strengths

- Wide coverage of classical and genomic breeding strategies.
- Translation of simulation results into operational recommendations.
- Explicit attention to genetic diversity and long-term response.
- Integration of genetic, temporal, and economic dimensions.
- Recognition that training-set management and cycle length are central to genomic selection.

## Limitations

- The paper is a broad narrative review rather than a formal meta-analysis.
- Results originate from simulations with heterogeneous assumptions, crops, population structures, and cost models.
- Some cited conclusions are highly parameter dependent and should not be transferred directly to another breeding program.
- The review summarises many studies but does not provide a unified quantitative benchmark for comparing strategies.
- Software descriptions are concise and do not assess reproducibility, computational scale, or ease of integration with operational databases.

## Main contribution

The review demonstrates that breeding optimisation is a multi-dimensional problem. The best strategy cannot be identified from heritability, accuracy, or gain alone. Simulation is most valuable when it represents the actual pipeline and jointly measures response, cycle time, diversity, GxE, and cost.
