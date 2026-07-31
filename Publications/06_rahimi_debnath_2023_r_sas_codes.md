# Rahimi and Debnath (2023): SAS and R codes for optimum, base, and desired-gain indices

## Bibliographic information

- **Authors:** Mehdi Rahimi and Sandip Debnath
- **Journal:** *Scientific Reports* 13:18977
- **Year:** 2023
- **DOI:** 10.1038/s41598-023-46368-6
- **Document type:** Methods and software application article
- **Source file:** `Rahimi et al 2023.pdf`

## Objective

The paper provides relatively simple SAS and R implementations for three selection-index methods:

1. optimum or Smith-Hazel-type index,
2. base index, and
3. Pešek-Baker desired-gain index.

The codes calculate index coefficients, candidate scores, expected gain, and several criteria for comparing alternative indices and economic-weight definitions. The intended audience is plant and animal breeders who have trait measurements and variance-covariance matrices but do not use specialised selection-index software.

## Index formulations

### Optimum index

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

### Base index

The base index uses the economic values directly as coefficients:

\[
I=\sum_i a_iX_i, \qquad \mathbf{b}=\mathbf{a}.
\]

It does not require estimates of \(\mathbf{P}\) and \(\mathbf{G}\). The method is simpler but does not adjust coefficients for differences in variance, reliability, or correlation among traits.

### Pešek-Baker index

The paper implements the desired-gain form as

\[
\mathbf{b}=\mathbf{G}^{-1}\mathbf{d},
\]

where \(\mathbf{d}\) is a vector of desired genetic gains. In the example, \(\mathbf{d}\) was set equal to the genetic standard deviation of each trait, although the code allows the user to enter another desired-gain vector.

## Evaluation criteria

The codes calculate several diagnostics:

- **\(R_{HI}\):** correlation between index and aggregate breeding value;
- **\(\Delta H\):** expected gain in aggregate merit;
- **\(\Delta_j\):** expected response for each trait;
- **RE:** relative efficiency compared with direct selection for the main trait, grain yield in the example; and
- **\(CV_I\):** phenotypic coefficient of variation of the index.

The use of multiple criteria is valuable because an index may rank candidates consistently but still produce poor response for the primary trait or an undesirable correlated change.

## Input and implementation

The software uses four main inputs:

- **X:** phenotype matrix, candidates by traits;
- **P:** phenotypic variance-covariance matrix;
- **G:** genotypic variance-covariance matrix; and
- **a:** economic-value vectors.

The SAS implementation uses PROC IML, while the R version reads the inputs from files and writes results to multiple Excel sheets. Parameters such as number of genotypes, number of traits, selection intensity, main-trait variance, heritability, and main-trait position must be specified.

The paper assumes a selection intensity of 10%, corresponding to a standardised selection differential of approximately 1.76 in the worked example.

## Worked maize example

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

## Reported index equations

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

## Main numerical results

Under equal weights:

| Index | \(R_{HI}\) | \(\Delta H\) | RE |
|---|---:|---:|---:|
| Optimum | 0.9887 | 125.7762 | 0.5504 |
| Base | 0.9885 | 128.6944 | 0.5555 |
| Pešek-Baker | 0.0018 | 9.2231 | 0.1986 |

The optimum and base indices were almost identical in candidate ranking. Their score correlation was 0.99979. Correlations of the Pešek-Baker scores with the base and optimum indices were only approximately 0.315 and 0.325.

The top five candidates differed accordingly. The optimum index selected genotypes 2, 4, 19, 22, and 6, whereas the base index selected 1, 4, 19, 22, and 6. The Pešek-Baker index selected a substantially different group: 3, 17, 8, 16, and 13.

The base index with regression coefficients as economic weights had the strongest values for several comparison criteria and was identified by the authors as the preferred index for this dataset.

## Interpretation of relative efficiency

All RE values were below 1. Therefore, none of the multi-trait indices produced a larger expected response for grain yield alone than direct selection for yield. This is not a failure of an index if the purpose is simultaneous improvement of several traits, but it is an important result. A multi-trait index trades some direct yield response for changes in other components of aggregate merit.

The article's conclusion that the base index was preferable applies to the selected objective and data. It should not be interpreted as a general superiority of the base index over the optimum index.

## Strengths

- Provides transparent R and SAS code rather than a black-box application.
- Accepts any number of candidates and traits.
- Compares several economic-weight strategies.
- Produces per-trait response and aggregate diagnostics.
- Includes a complete worked example.
- Stores R outputs in structured Excel sheets.

## Limitations and cautions

- The empirical demonstration is limited to 28 maize lines and seven traits.
- Broad-sense genotypic variance is used in parts of the analysis, whereas response to selection in a recurrent breeding program is primarily determined by transmissible additive effects.
- Correlations and stepwise-regression coefficients are statistical weights, not necessarily economic values or product-profile utilities.
- Stepwise regression can be unstable under collinearity and small sample size.
- The near equality of base and optimum indices is specific to this covariance structure and weight vector.
- Expected gains were not validated in progeny or a subsequent cycle.
- The very low \(R_{HI}\) for the implemented Pešek-Baker example suggests that the chosen desired-gain vector was poorly aligned with the defined aggregate objective.
- Input matrices must be positive definite and consistently ordered; the paper gives limited discussion of numerical diagnostics when inversion is unstable.

## Practical use

The code is useful for teaching, prototyping, and comparing alternative objectives. A robust operational application should:

1. estimate additive genetic and residual covariance matrices with an appropriate mixed model;
2. verify matrix ordering and conditioning;
3. define economic values or desired gains from a target product profile;
4. standardise traits when units differ substantially;
5. evaluate per-trait response and aggregate response;
6. test sensitivity to alternative weights; and
7. validate the selected index prospectively.

The main contribution is accessibility. The paper enables breeders to calculate and compare selection indices in R or SAS with familiar tabular inputs.
