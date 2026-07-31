# Covarrubias-Pazaran: Bringing a selection index into the CIMMYT maize programs

## Document information

- **Author:** Giovanny Eduardo Covarrubias-Pazaran
- **Affiliation shown:** Excellence in Breeding Platform / CGIAR
- **Document type:** Technical presentation, 21 pages
- **Source file:** `Giovanny Selection index in Maize.pdf`

## Purpose

The presentation documents a practical process for replacing parts of a conventional multi-step maize selection procedure with a formal selection index. It does not begin with a theoretical index and search for a use case. Instead, it first reconstructs the breeder's existing decisions, expresses them as an algorithm, derives retrospective weights from observed selection differentials, and compares the resulting index with the current procedure.

The proposed workflow is:

1. review what is known about selection methods and genetic gain;
2. translate the current selection procedure into explicit algorithmic steps;
3. identify the steps that can be replaced by an index;
4. construct an index;
5. compare selections, selection differentials, and total merit; and
6. refine the index.

## Selection index versus culling and tandem selection

Pages 3 and 4 compare three strategies for two traits:

- **Selection index:** candidates are ranked by total merit, \(I=w_1t_1+w_2t_2\).
- **Independent culling:** candidates must pass a threshold for each trait.
- **Tandem-type selection:** candidates are selected sequentially for different traits.

The simulated examples illustrate that selecting extreme individuals for one trait can produce parents with poor performance for another trait. The selection index produced a more balanced increase in total merit. The presentation reports illustrative gains of approximately 70% versus 50% for Trait 1, 66% versus 38% for Trait 2, and 57% versus 53% for total merit, depending on the compared method.

The biological message is that the best parent is not necessarily the most extreme genotype for any single trait. A useful parent should combine acceptable performance across the complete target profile.

## Reconstructing the S3-to-S4 procedure

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

## Retrospective index weights

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

## Comparison with the current procedure

The presentation compares selections from the existing "Yoseph" approach with the index. Many selected materials were common, but each method also selected unique families or lines. This is expected even when total merit is similar because small changes in coefficients can alter rankings near the selection threshold.

The supplementary comparison on page 21 reports the following selection differentials:

| Method | Opt. | MD | LN | AD | MLN | Striga | Grain moisture | Total merit |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Current approach | 1.09 | 0.35 | 0.48 | -0.39 | 0.07 | 0.30 | -0.31 | 1.61 |
| Index | 1.07 | 0.41 | 0.49 | 0.53 | 0.04 | 0.23 | 0.29 | 1.70 |

The index slightly increased total merit in this comparison. However, the sign and desirability of anthesis date and grain moisture must be aligned explicitly before deployment because the table alone does not indicate whether positive or negative values represent improvement.

## Recycling-stage index

A second application focuses on recycling parents using general combining ability under optimum and drought conditions. The index is calculated among families and then within families. Traits include grain yield, anthesis date, anthesis-silking interval, ear-related traits, lodging or breakage traits, ears per plant, and moisture under different conditions.

The derived coefficients show strong positive emphasis on GCA for grain yield and substantial negative emphasis on some undesirable traits. The presentation reports nearly identical total merit for the current approach and the index, approximately 3.1327 and 3.1345, respectively. Nevertheless, the candidate lists differ, demonstrating that an index can reproduce overall selection intensity while changing which lines are retained.

## Conclusions from the presentation

The author concludes that:

- a traditional selection procedure can be translated into an index;
- an index provides consistency across cycles and users;
- initial weights can be obtained retrospectively from historical decisions;
- the index should be adopted and subsequently refined; and
- future indices should be defined by **program, stage, and product profile**.

The proposed next step is to obtain data from Stage 1 and Stage 2, particularly for recycling, across programs such as Latin America, East Africa, and South Asia.

## Strengths

- Starts from the real breeding process rather than an abstract index.
- Makes hidden selection rules explicit and auditable.
- Separates hard culling constraints from compensatory traits.
- Uses historical selections to derive an initial index.
- Compares both total merit and overlap in selected candidates.
- Recognises the need for program-, stage-, and product-profile-specific indices.

## Limitations and implementation cautions

- The presentation is not a peer-reviewed study and provides limited information on sample size, uncertainty, and statistical validation.
- Trait abbreviations and the direction of improvement are not defined consistently on all pages.
- Some selection differentials are presented without standard errors or confidence intervals.
- Retrospective weights reproduce past behaviour; they do not demonstrate that the past behaviour was optimal.
- A linear index cannot safely replace mandatory disease, quality, or maturity thresholds unless those constraints are retained separately.
- The covariance matrix and trait scaling must be updated when the population, target environment, or stage changes.

## Operational interpretation

The main value of the presentation is methodological. It shows how a breeder's informal sequence of filters can be converted into a reproducible decision system. The recommended practice is a hybrid system: use explicit culling for non-negotiable requirements and use a selection index for traits where compensation is biologically and commercially acceptable.
