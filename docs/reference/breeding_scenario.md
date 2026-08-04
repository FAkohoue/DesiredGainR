# Describe a breeding simulation scenario

A scenario records the biological assumptions that define a simulation.
Keeping these assumptions in one object reduces long function calls. It
also makes comparisons auditable.

## Usage

``` r
breeding_scenario(
  label,
  programme = c("self", "outcross", "clonal", "testcross"),
  architecture = list(),
  evaluation = list(),
  environments = list()
)
```

## Arguments

- label:

  Short scenario name.

- programme:

  One of self, outcross, clonal, or testcross.

- architecture:

  Named list. Supported fields include qtl_per_chromosome,
  markers_per_chromosome, effect_distribution, qtl_shape,
  dominance_degree, and dominance_variance.

- evaluation:

  Named list describing phenotype, GEBV, or GCA information.

- environments:

  Named list describing the target environments and their genetic
  correlation structure.

## Value

An object of class desiredgainr_scenario.

## Details

The programme can represent self-pollinated, outcrossing, clonal, or
testcross evaluation. Testcross describes the source of general
combining ability information. Parent allocation and cross design remain
with HapBlockR.

## References

Beavis WD, Mahama AA, Suza W (2023). Simulation Modeling. In
*Quantitative Genetics for Plant Breeding*. Iowa State University
Digital Press.
