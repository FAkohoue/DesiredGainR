# Finding C1 of AUDIT-2026-07-31.md.
#
# The previous diversity metric recomputed allele frequencies from the
# population being measured. Centring on a sample's own frequencies forces
# every marker column to sum to zero, so the whole relationship matrix sums to
# zero and the mean off-diagonal element equals -1/(n-1) whatever has happened
# to the germplasm. The metric was therefore a function of the number of
# parents and almost nothing else, and it was used as the diversity axis of the
# Pareto frontier in optimize_desired_gains(include_diversity = TRUE).
#
# These tests fix the base frequencies at the founders. The first is a pure
# algebra check that needs no simulation; the remainder need AlphaSimR.

test_that("self-centring destroys the signal it was meant to carry", {
  # This reproduces the defect directly, so the reason for the fix stays
  # visible in the test suite rather than only in the audit document.
  self_centred_mean <- function(geno) {
    frequency <- colMeans(geno) / 2
    keep <- frequency > 0 & frequency < 1
    geno <- geno[, keep, drop = FALSE]
    frequency <- frequency[keep]
    centred <- sweep(geno, 2L, 2 * frequency, "-")
    relationship <- tcrossprod(centred) / (2 * sum(frequency * (1 - frequency)))
    mean(relationship[upper.tri(relationship)])
  }

  # The exact identity. Centring on the sample's own frequencies forces every
  # column to sum to zero, so the whole matrix sums to zero and
  #   sum(off-diagonal) = -sum(diagonal),
  # giving mean(off-diagonal) = -mean(diagonal)/(n - 1).
  # It is NOT -1/(n-1): that would require the mean diagonal to be exactly 1,
  # which holds only when the measured population's heterozygosity happens to
  # equal the scaling constant's. Asserting the stronger form was wrong.
  self_centred_diagonal <- function(geno) {
    frequency <- colMeans(geno) / 2
    keep <- frequency > 0 & frequency < 1
    geno <- geno[, keep, drop = FALSE]
    frequency <- frequency[keep]
    centred <- sweep(geno, 2L, 2 * frequency, "-")
    relationship <- tcrossprod(centred) / (2 * sum(frequency * (1 - frequency)))
    mean(diag(relationship))
  }

  set.seed(4L)
  n <- 30L
  diverse <- matrix(stats::rbinom(n * 200L, 2L, 0.5), nrow = n)
  # A population reduced to a handful of distinct genotypes, each replicated:
  # unambiguously far less diverse than the first.
  founders <- matrix(stats::rbinom(3L * 200L, 2L, 0.5), nrow = 3L)
  narrow <- founders[rep(seq_len(3L), length.out = n), , drop = FALSE]

  for (geno in list(diverse, narrow)) {
    expect_equal(
      self_centred_mean(geno),
      -self_centred_diagonal(geno) / (n - 1),
      tolerance = 1e-8
    )
  }
  # Both are pinned by n alone, to within the difference in mean diagonal, so
  # the metric cannot separate two populations of obviously different
  # diversity.
  expect_lt(abs(self_centred_mean(diverse) - self_centred_mean(narrow)), 0.02)
})

test_that("a fixed base separates populations that self-centring cannot", {
  set.seed(5L)
  n <- 30L
  base_geno <- matrix(stats::rbinom(60L * 200L, 2L, 0.5), nrow = 60L)
  base_frequency <- colMeans(base_geno) / 2
  keep <- base_frequency > 0 & base_frequency < 1
  base_frequency <- base_frequency[keep]
  scaling <- 2 * sum(base_frequency * (1 - base_frequency))

  fixed_base_mean <- function(geno) {
    geno <- geno[, keep, drop = FALSE]
    centred <- sweep(geno, 2L, 2 * base_frequency, "-")
    relationship <- tcrossprod(centred) / scaling
    mean(relationship[upper.tri(relationship)])
  }

  diverse <- base_geno[sample.int(60L, n), , drop = FALSE]
  founders <- base_geno[sample.int(60L, 3L), , drop = FALSE]
  narrow <- founders[rep(seq_len(3L), length.out = n), , drop = FALSE]

  expect_gt(fixed_base_mean(narrow), fixed_base_mean(diverse))
})

phased_founders <- function(n_ind = 60L, n_var_per_chr = 50L, n_chr = 2L,
                            seed = 11L) {
  set.seed(seed)
  n_var <- n_var_per_chr * n_chr
  variant_id <- sprintf("v%04d", seq_len(n_var))
  individual_id <- sprintf("f%03d", seq_len(n_ind))
  draw <- function() {
    matrix(
      stats::rbinom(n_var * n_ind, 1L, 0.5),
      nrow = n_var,
      dimnames = list(variant_id, individual_id)
    )
  }
  map <- data.frame(
    variant_id = variant_id,
    chromosome = rep(seq_len(n_chr), each = n_var_per_chr),
    position_bp = rep(seq_len(n_var_per_chr) * 1e6, times = n_chr),
    stringsAsFactors = FALSE
  )
  founder_haplotypes(draw(), draw(), map)
}

test_that("relatedness against a fixed base accumulates under selection", {
  skip_on_cran()
  skip_if_not_installed("AlphaSimR")

  traits <- c("YLD", "DIS")
  G <- matrix(
    c(1.00, -0.25, -0.25, 0.64), 2L,
    dimnames = list(traits, traits)
  )
  setup <- founder_population(
    phased_founders(),
    G = G, h2 = c(YLD = 0.4, DIS = 0.4),
    n_qtl_per_chromosome = 20L, seed = 3L
  )
  SP <- setup$SP
  base <- DesiredGainR:::.dgr_base_frequency(setup$founder_pop, SP)

  relationship <- numeric(6L)
  inbreeding <- numeric(6L)
  current <- setup$founder_pop
  for (cycle in seq_len(6L)) {
    measure <- DesiredGainR:::.dgr_mean_relationship(current, SP, base)
    relationship[cycle] <- measure$relationship
    inbreeding[cycle] <- measure$inbreeding
    # Deliberately severe: six parents recycled each cycle.
    parents <- AlphaSimR::selectInd(
      current,
      nInd = min(6L, AlphaSimR::nInd(current)),
      use = "pheno", trait = 1L, simParam = SP
    )
    current <- AlphaSimR::randCross(parents, nCrosses = 60L, simParam = SP)
  }

  # The defective metric returned -1/(n-1) at every cycle, so it could not
  # produce an increase at any tolerance. A fixed base must.
  expect_gt(relationship[6L], relationship[1L])
  expect_gt(stats::cor(seq_len(6L), relationship), 0.8)
  expect_gt(inbreeding[6L], inbreeding[1L])
  # It must also not be pinned to a function of the sample size alone.
  expect_gt(abs(relationship[6L] + 1 / (AlphaSimR::nInd(current) - 1)), 1e-3)
})

test_that("simulate_selection_cycles() reports accumulating inbreeding", {
  skip_on_cran()
  skip_if_not_installed("AlphaSimR")

  traits <- c("YLD", "DIS")
  G <- matrix(
    c(1.00, -0.25, -0.25, 0.64), 2L,
    dimnames = list(traits, traits)
  )
  setup <- founder_population(
    phased_founders(),
    G = G, h2 = c(YLD = 0.4, DIS = 0.4),
    n_qtl_per_chromosome = 20L, seed = 3L
  )
  simulation <- simulate_selection_cycles(
    setup,
    desired_gains = c(YLD = 1, DIS = 0.5),
    n_cycles = 5L,
    mating_system = "outcross",
    n_parents = 8L,
    n_crosses = 40L,
    n_progeny_per_cross = 4L,
    seed = 7L
  )

  table <- simulation$cycle_table
  expect_true("mean_inbreeding" %in% names(table))
  by_cycle <- unique(table[, c("cycle", "mean_inbreeding")])
  data.table::setorderv(by_cycle, "cycle")
  inbreeding <- by_cycle$mean_inbreeding
  expect_gt(inbreeding[length(inbreeding)], inbreeding[1L])
  expect_gt(stats::cor(by_cycle$cycle, inbreeding), 0.8)
  expect_true(is.finite(simulation$final_inbreeding))
  expect_false(is.null(simulation$diversity_basis))
  # Effective size derived from the change in inbreeding must be positive and
  # of a plausible order for eight recycled parents.
  effective <- table$effective_size[is.finite(table$effective_size)]
  expect_true(length(effective) > 0L)
  expect_true(all(effective > 0))
})
