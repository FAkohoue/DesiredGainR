# Validate and package phased founder haplotypes

Multi-cycle simulation in DesiredGainR is calibrated to the breeder's
own germplasm rather than to an assumed genetic architecture. Hence the
founder population is always built from real phased marker data.

## Usage

``` r
founder_haplotypes(
  hap1 = NULL,
  hap2 = NULL,
  map,
  missing_policy = c("error", "drop_variant", "drop_individual"),
  variant_col = "variant_id",
  chromosome_col = "chromosome",
  position_morgan_col = "position_morgan",
  position_bp_col = "position_bp",
  recombination_rate = 1e-08,
  haplotypes = NULL
)
```

## Arguments

- hap1, hap2:

  Variant-by-individual homologue matrices coded 0, 1 or `NA`.

- map:

  Data frame describing the variants, requiring the variant identifier,
  the chromosome, and either a genetic or a physical position.

- missing_policy:

  Treatment of missing calls, which AlphaSimR cannot accept: `"error"`,
  `"drop_variant"` or `"drop_individual"`.

- variant_col, chromosome_col:

  Column names in `map`.

- position_morgan_col:

  Column giving genetic position in Morgans. When absent,
  `position_bp_col` and `recombination_rate` are used instead.

- position_bp_col:

  Column giving physical position in base pairs.

- recombination_rate:

  Morgans per base pair, used only when genetic positions are absent.
  The default corresponds to one centimorgan per megabase.

- haplotypes:

  Optional list of homologue matrices. Only a list of length two is
  accepted; see Details.

## Value

An object of class `desiredgainr_founders` holding the per-chromosome
genetic map and interleaved haplotype matrices required by AlphaSimR,
together with the individual and variant identifiers and a record of any
calls removed.

## Input shape

`hap1` and `hap2` are variant-by-individual matrices coded 0, 1 or `NA`,
giving the allele carried on the first and second homologue. Row names
must be variant identifiers and column names individual identifiers, and
both matrices must carry identical dimnames in identical order. This is
the shape returned by
[`HapBlockR::read_phased_vcf()`](https://FAkohoue.github.io/HapBlockR/reference/read_phased_vcf.html),
whose `hap1` and `hap2` can be passed directly, with its `snp_info`
supplying `map`.

## Why a haplotype is coded 0 or 1, and dosage is not accepted

A haplotype is one chromosome copy, so at a biallelic locus it carries a
single allele rather than an allele count. Genotype dosage coded 0, 1, 2
is the sum of the homologues; it records how many alternative alleles an
individual carries but not which copy carries them, and therefore does
not identify which alleles lie in cis across loci. Since that
co-occurrence is the linkage disequilibrium the simulation exists to
represent, dosage is rejected. Where the material is highly inbred
diploid, use
[`haplotypes_from_inbred_dosage()`](https://FAkohoue.github.io/DesiredGainR/reference/haplotypes_from_inbred_dosage.md),
which derives phase from homozygous calls alone.

## Diploid only

The contract is deliberately diploid and biallelic. Extending `hap1` and
`hap2` to further homologues would not make the simulation polyploid,
because autopolyploid and allopolyploid meiosis require an explicit
pairing and recombination model, including the treatment of multivalents
and double reduction, which differs by crop. Supplying a polyploid
dataset therefore raises an error rather than silently applying a
diploid meiotic model.

## Per-marker haplotypes, not block haplotype states

Two haplotype representations are common in haplotype-aware pipelines,
and only the first belongs here.

- Per-marker phased haplotypes:

  One row per marker, coded 0 or 1, stating the allele carried on a
  given chromosome copy. A genome simulator models recombination between
  individual loci and therefore needs their separate positions.

- Block haplotype states:

  One entry per linkage-disequilibrium block, recording which distinct
  haplotype an individual carries in that block. These are the correct
  representation for a haplotype relationship matrix, block association
  testing and local breeding values, but they have aggregated across the
  constituent markers and no longer provide a chromosome-wide homologue
  representation or a marker-level recombination map.

Read a phased variant call format file with a reader that preserves
phase. A general genotype reader may return dosage even for a phased
file, in which case the phase is discarded on import.

## See also

[`haplotypes_from_inbred_dosage()`](https://FAkohoue.github.io/DesiredGainR/reference/haplotypes_from_inbred_dosage.md),
[`founder_population()`](https://FAkohoue.github.io/DesiredGainR/reference/founder_population.md)

## Examples

``` r
variants <- c("v1", "v2", "v3", "v4")
individuals <- c("a", "b")
hap1 <- matrix(c(0, 1, 0, 1, 1, 0, 1, 0), nrow = 4,
               dimnames = list(variants, individuals))
hap2 <- matrix(c(0, 0, 1, 1, 1, 1, 0, 0), nrow = 4,
               dimnames = list(variants, individuals))
map <- data.frame(
  variant_id = variants,
  chromosome = c(1, 1, 2, 2),
  position_bp = c(1e6, 5e6, 1e6, 4e6)
)
founder_haplotypes(hap1, hap2, map)
#> <desiredgainr_founders>
#>   2 individuals, 4 variants, 2 chromosomes, diploid
#>   Map: converted from physical position at 1e-08 Morgan per base pair 
#>   Missing calls before resolution: 0.000% (policy: error)
```
