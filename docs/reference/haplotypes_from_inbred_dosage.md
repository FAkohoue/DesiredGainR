# Derive phased haplotypes from inbred-line dosage

Phase is required for simulation because it determines which favourable
alleles segregate together. However, in a diploid inbred line a dosage
of 0 or 2 is already unambiguous: both homologues carry the same allele.
Therefore dosage from doubled haploids, recombinant inbred lines or
advanced selfed generations can be converted without external phasing,
and only the heterozygous and missing calls require a decision.

## Usage

``` r
haplotypes_from_inbred_dosage(
  dosage,
  heterozygous_policy = c("error", "drop_variant", "drop_individual", "mask"),
  missing_policy = c("error", "drop_variant", "drop_individual", "mask")
)
```

## Arguments

- dosage:

  Variant-by-individual matrix coded 0, 1, 2 with `NA` for missing
  calls.

- heterozygous_policy:

  Treatment of dosage 1. `"error"` refuses to proceed, `"drop_variant"`
  removes affected variants, `"drop_individual"` removes affected
  individuals, and `"mask"` sets the calls to `NA` for resolution by
  [`founder_haplotypes()`](https://FAkohoue.github.io/DesiredGainR/reference/founder_haplotypes.md).

- missing_policy:

  Treatment of missing calls, with the same options.

## Value

A list with `hap1` and `hap2` suitable for
[`founder_haplotypes()`](https://FAkohoue.github.io/DesiredGainR/reference/founder_haplotypes.md),
carrying a `conversion` attribute that records the observed rates and
everything removed or masked.

## Details

Hence this helper serves self-pollinated diploid programmes, such as
wheat, rice and common bean, that hold dosage rather than phased data.
It must not be used for outcrossing or clonal material, where
heterozygosity is pervasive. Where heterozygosity is appreciable, or
where cis and trans relationships matter, phase the data externally
instead.

No threshold is imposed on residual heterozygosity, because no
universally appropriate level exists. The rates are measured, reported
and returned instead, and heterozygous calls are never resolved
silently: an explicit policy is always required. Inspect
[`dosage_diagnostics()`](https://FAkohoue.github.io/DesiredGainR/reference/dosage_diagnostics.md)
first, since a low overall rate can conceal a few badly affected
individuals or markers, and the cheaper of `"drop_variant"` and
`"drop_individual"` depends on which.

## See also

[`dosage_diagnostics()`](https://FAkohoue.github.io/DesiredGainR/reference/dosage_diagnostics.md),
[`founder_haplotypes()`](https://FAkohoue.github.io/DesiredGainR/reference/founder_haplotypes.md)

## Examples

``` r
dosage <- matrix(
  c(0, 2, 2, 0, 2, 0),
  nrow = 3,
  dimnames = list(c("v1", "v2", "v3"), c("line1", "line2"))
)
haplotypes_from_inbred_dosage(dosage)
#> $hap1
#>    line1 line2
#> v1     0     0
#> v2     1     1
#> v3     1     0
#> 
#> $hap2
#>    line1 line2
#> v1     0     0
#> v2     1     1
#> v3     1     0
#> 
#> attr(,"conversion")
#> attr(,"conversion")$diagnostics
#> <desiredgainr_dosage_diagnostics>
#>   3 variants x 2 individuals
#>   Heterozygous calls: 0.000%   Missing calls: 0.000%
#>   Per-individual heterozygosity, maximum 0.000% (line1)
#>   Per-marker heterozygosity, maximum 0.000%
#> 
#> attr(,"conversion")$heterozygous_policy
#> [1] "error"
#> 
#> attr(,"conversion")$missing_policy
#> [1] "error"
#> 
#> attr(,"conversion")$n_heterozygous_calls
#> [1] 0
#> 
#> attr(,"conversion")$n_missing_calls
#> [1] 0
#> 
#> attr(,"conversion")$removed_variants
#> character(0)
#> 
#> attr(,"conversion")$removed_individuals
#> character(0)
#> 
#> attr(,"conversion")$basis
#> [1] "Phase was determined from homozygous dosage only. No heterozygous call was assigned to a homologue, because dosage does not identify which alleles are in cis across loci."
#> 
```
