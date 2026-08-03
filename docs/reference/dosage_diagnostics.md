# Summarise heterozygosity and missingness in a dosage matrix

No universally appropriate residual-heterozygosity threshold exists for
self-pollinated material, because the level depends on the generation,
mating history, crop, population type, genotyping error rate and
quality-control procedure. It must therefore be measured from each
dataset rather than assumed. Run this function before
[`haplotypes_from_inbred_dosage()`](https://FAkohoue.github.io/DesiredGainR/reference/haplotypes_from_inbred_dosage.md)
and inspect the per-individual and per-marker distributions, since a low
overall rate can still conceal a small number of badly affected
individuals or markers.

## Usage

``` r
dosage_diagnostics(dosage)
```

## Arguments

- dosage:

  Variant-by-individual matrix coded 0, 1, 2 with `NA` for missing
  calls.

## Value

An object of class `desiredgainr_dosage_diagnostics` giving the overall
heterozygous and missing rates together with per-individual and
per-marker tables.

## See also

[`haplotypes_from_inbred_dosage()`](https://FAkohoue.github.io/DesiredGainR/reference/haplotypes_from_inbred_dosage.md)

## Examples

``` r
dosage <- matrix(
  c(0, 2, 1, 0, 2, 0, 2, NA, 2),
  nrow = 3,
  dimnames = list(c("v1", "v2", "v3"), c("l1", "l2", "l3"))
)
dosage_diagnostics(dosage)
#> <desiredgainr_dosage_diagnostics>
#>   3 variants x 3 individuals
#>   Heterozygous calls: 12.500%   Missing calls: 11.111%
#>   Per-individual heterozygosity, maximum 33.333% (l1)
#>   Per-marker heterozygosity, maximum 33.333%
```
