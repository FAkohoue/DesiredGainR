# Draft message to HapBlockR users

Subject: Haplotype format required for the DesiredGainR simulation layer

------------------------------------------------------------------------

Dear colleagues,

Specifying economic weights remains the main obstacle to routine use of
selection indices, and I am therefore extending DesiredGainR with a
multi-cycle simulation layer that identifies which desired-gain weights
deliver the greatest genetic gain over several breeding cycles. The
founder population for that simulation is built from real marker data
rather than from a simulated genome, so that its linkage disequilibrium
(LD) and allele-frequency structure are those of the actual breeding
programme. Since several of you already run HapBlockR, I would like to
confirm which haplotype output you hold before I fix the input contract.

HapBlockR produces two distinct haplotype representations. The first is
per-marker phased haplotypes from `read_phased_vcf()`, returned as
`hap1` and `hap2` matrices coded 0 or 1, with variants in rows and
individuals in columns. The second is block-level haplotype states from
`extract_haplotypes()` or `infer_block_haplotypes()`, which record which
distinct haplotype an individual carries within each LD block.

Only the first can drive a genome simulator. Recombination is modelled
between individual loci, whereas the block states have already
aggregated across the markers they contain, so their within-block
positions are no longer available. Hence the simulation requires
`read_phased_vcf()` output specifically.

My questions are therefore:

1.  Which route do you currently use to obtain haplotypes:
    `phase_with_beagle()` followed by `read_phased_vcf()`,
    `extract_haplotypes()` or `infer_block_haplotypes()` for block
    states, or `read_geno()` returning dosage only?

2.  Do you retain the per-marker `hap1` and `hap2` matrices, or only the
    block-level products derived from them?

3.  If your genotypes are unphased, is external phasing with Beagle
    practical in your pipeline, or would you prefer to supply dosage
    coded 0, 1, 2?

4.  Which crop and ploidy are you working with? I currently support any
    ploidy, but I would like to know whether the tetraploid and
    hexaploid cases are actually needed.

5.  For self-pollinated material, what residual heterozygosity do you
    typically observe? DesiredGainR can derive phase from dosage for
    inbred lines, because a dosage of 0 or 2 is unambiguous, but only
    when heterozygous calls are rare.

If it is convenient, the following four lines answer questions 1, 2 and
4 in one step:

``` r
phased <- HapBlockR::read_phased_vcf("your_phased_file.vcf")
dim(phased$hap1)
sort(unique(as.vector(phased$hap1)))
phased$hap1[1:5, 1:3]
```

Your answers will determine two decisions: (i) whether DesiredGainR
requires phased input in every case or also accepts dosage for inbred
material, and (ii) whether the diploid interface is sufficient or
polyploid input should be exposed as the default. A short or partial
reply is entirely welcome.

With best regards,

Félicien Akohoue
