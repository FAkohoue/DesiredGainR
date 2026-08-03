# Tests for the founder-haplotype contract and the multi-cycle simulation.
#
# The AlphaSimR-dependent tests are skipped when the package is absent, so the
# suite remains runnable in a minimal environment.

make_phased <- function(n_ind = 24L, n_var_per_chr = 30L, n_chr = 2L,
                        seed = 7L) {
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
  list(hap1 = draw(), hap2 = draw(), map = map)
}

test_that("founder_haplotypes interleaves homologues correctly", {
  phased <- make_phased(n_ind = 4L, n_var_per_chr = 5L, n_chr = 1L)
  founders <- founder_haplotypes(
    phased$hap1, phased$hap2, phased$map
  )
  expect_s3_class(founders, "desiredgainr_founders")
  expect_equal(founders$n_individuals, 4L)
  expect_equal(founders$n_chromosomes, 1L)

  block <- founders$haplotypes[[1L]]
  expect_equal(nrow(block), 8L)
  expect_equal(ncol(block), 5L)
  # Odd rows are the first homologue, even rows the second.
  expect_equal(as.integer(block[1L, ]), as.integer(phased$hap1[, 1L]))
  expect_equal(as.integer(block[2L, ]), as.integer(phased$hap2[, 1L]))
  expect_equal(as.integer(block[7L, ]), as.integer(phased$hap1[, 4L]))
  expect_equal(as.integer(block[8L, ]), as.integer(phased$hap2[, 4L]))
})

test_that("the genetic map starts at zero and increases strictly", {
  phased <- make_phased(n_ind = 4L, n_var_per_chr = 6L, n_chr = 2L)
  founders <- founder_haplotypes(phased$hap1, phased$hap2, phased$map)
  for (chromosome_map in founders$gen_map) {
    expect_equal(min(chromosome_map), 0)
    expect_true(all(diff(chromosome_map) > 0))
  }
})

test_that("dosage input and mismatched dimnames are rejected", {
  phased <- make_phased(n_ind = 4L, n_var_per_chr = 5L, n_chr = 1L)
  dosage <- phased$hap1 + phased$hap2
  expect_error(
    founder_haplotypes(dosage, phased$hap2, phased$map),
    "only 0, 1 and NA"
  )
  # The message must point the user at the inbred conversion route.
  expect_error(
    founder_haplotypes(dosage, phased$hap2, phased$map),
    "haplotypes_from_inbred_dosage"
  )

  shuffled <- phased$hap2
  colnames(shuffled) <- rev(colnames(shuffled))
  expect_error(
    founder_haplotypes(phased$hap1, shuffled, phased$map),
    "identical dimnames"
  )
})

test_that("polyploid input is refused with an explanatory error", {
  phased <- make_phased(n_ind = 3L, n_var_per_chr = 4L, n_chr = 1L)
  extra <- make_phased(n_ind = 3L, n_var_per_chr = 4L, n_chr = 1L, seed = 8L)
  tetraploid <- list(
    h1 = phased$hap1, h2 = phased$hap2,
    h3 = extra$hap1, h4 = extra$hap2
  )
  expect_error(
    founder_haplotypes(map = phased$map, haplotypes = tetraploid),
    "Only diploid material is supported"
  )
  # The message must explain why, rather than merely stating the restriction.
  expect_error(
    founder_haplotypes(map = phased$map, haplotypes = tetraploid),
    "double\\s+reduction"
  )
})

test_that("missing haplotype calls are resolved under an explicit policy", {
  phased <- make_phased(n_ind = 5L, n_var_per_chr = 6L, n_chr = 1L)
  phased$hap1[2L, 3L] <- NA_integer_

  expect_error(
    founder_haplotypes(phased$hap1, phased$hap2, phased$map),
    "missing haplotype call"
  )

  dropped <- founder_haplotypes(
    phased$hap1, phased$hap2, phased$map,
    missing_policy = "drop_variant"
  )
  expect_equal(dropped$n_variants, 5L)
  expect_equal(dropped$missing_data$removed_variants, "v0002")
  expect_equal(dropped$missing_data$n_calls_affected, 1L)

  by_individual <- founder_haplotypes(
    phased$hap1, phased$hap2, phased$map,
    missing_policy = "drop_individual"
  )
  expect_equal(by_individual$n_individuals, 4L)
  expect_equal(by_individual$missing_data$removed_individuals, "f003")
})

test_that("dosage diagnostics summarise heterozygosity and missingness", {
  dosage <- matrix(
    c(
      0L, 2L, 1L,
      0L, 2L, 0L,
      2L, NA_integer_, 2L
    ),
    nrow = 3L,
    dimnames = list(c("v1", "v2", "v3"), c("l1", "l2", "l3"))
  )
  diagnostics <- dosage_diagnostics(dosage)
  expect_s3_class(diagnostics, "desiredgainr_dosage_diagnostics")
  expect_equal(diagnostics$n_calls, 8L)
  expect_equal(diagnostics$overall_heterozygosity, 1 / 8)
  expect_equal(diagnostics$overall_missing, 1 / 9)
  expect_equal(nrow(diagnostics$by_individual), 3L)
  expect_equal(nrow(diagnostics$by_marker), 3L)
  # The heterozygous call sits in v3 of l1.
  expect_equal(
    diagnostics$by_marker[variant == "v3", heterozygosity], 1 / 3
  )
})

test_that("inbred dosage converts losslessly at homozygous calls", {
  dosage <- matrix(
    c(0L, 2L, 2L, 0L, 2L, 0L),
    nrow = 3L,
    dimnames = list(c("v1", "v2", "v3"), c("line1", "line2"))
  )
  converted <- haplotypes_from_inbred_dosage(dosage)
  expect_equal(converted$hap1, converted$hap2)
  expect_equal(converted$hap1 + converted$hap2, dosage)
  conversion <- attr(converted, "conversion")
  expect_equal(conversion$n_heterozygous_calls, 0L)
  expect_equal(conversion$diagnostics$overall_heterozygosity, 0)
})

test_that("heterozygous calls are never resolved silently", {
  dosage <- matrix(
    c(0L, 2L, 1L, 0L, 2L, 0L),
    nrow = 3L,
    dimnames = list(c("v1", "v2", "v3"), c("line1", "line2"))
  )
  expect_error(
    haplotypes_from_inbred_dosage(dosage),
    "heterozygous call"
  )
  # No policy assigns a heterozygote to a homologue; the available options
  # remove or mask it.
  expect_false("random" %in% eval(
    formals(haplotypes_from_inbred_dosage)$heterozygous_policy
  ))

  dropped <- haplotypes_from_inbred_dosage(
    dosage,
    heterozygous_policy = "drop_variant"
  )
  expect_equal(nrow(dropped$hap1), 2L)
  expect_equal(attr(dropped, "conversion")$removed_variants, "v3")

  masked <- haplotypes_from_inbred_dosage(
    dosage,
    heterozygous_policy = "mask"
  )
  expect_true(is.na(masked$hap1["v3", "line1"]))
  expect_true(is.na(masked$hap2["v3", "line1"]))
})

test_that("no arbitrary heterozygosity threshold is imposed", {
  set.seed(4)
  dosage <- matrix(
    sample(c(0L, 1L, 2L), 200L, replace = TRUE),
    nrow = 20L,
    dimnames = list(paste0("v", 1:20), paste0("i", 1:10))
  )
  # Heavily heterozygous material is reported, not silently refused on a
  # threshold that has no universal value.
  diagnostics <- dosage_diagnostics(dosage)
  expect_gt(diagnostics$overall_heterozygosity, 0.2)
  converted <- haplotypes_from_inbred_dosage(
    dosage,
    heterozygous_policy = "mask"
  )
  expect_equal(
    attr(converted, "conversion")$diagnostics$overall_heterozygosity,
    diagnostics$overall_heterozygosity
  )
})

test_that("masked dosage flows into founder_haplotypes under a policy", {
  set.seed(6)
  n_var <- 12L
  dosage <- matrix(
    sample(c(0L, 2L), n_var * 4L, replace = TRUE),
    nrow = n_var,
    dimnames = list(sprintf("v%02d", seq_len(n_var)), paste0("line", 1:4))
  )
  dosage[3L, 2L] <- 1L
  converted <- haplotypes_from_inbred_dosage(
    dosage,
    heterozygous_policy = "mask"
  )
  map <- data.frame(
    variant_id = rownames(dosage),
    chromosome = rep(1:2, each = n_var / 2L),
    position_bp = rep(seq_len(n_var / 2L) * 1e6, times = 2L),
    stringsAsFactors = FALSE
  )
  founders <- founder_haplotypes(
    converted$hap1, converted$hap2, map,
    missing_policy = "drop_variant"
  )
  expect_equal(founders$n_variants, n_var - 1L)
  expect_equal(founders$missing_data$removed_variants, "v03")
})

test_that("converted inbred dosage feeds founder_haplotypes directly", {
  set.seed(5)
  n_var <- 12L
  dosage <- matrix(
    sample(c(0L, 2L), n_var * 4L, replace = TRUE),
    nrow = n_var,
    dimnames = list(sprintf("v%02d", seq_len(n_var)), paste0("line", 1:4))
  )
  converted <- haplotypes_from_inbred_dosage(dosage)
  map <- data.frame(
    variant_id = rownames(dosage),
    chromosome = rep(1:2, each = n_var / 2L),
    position_bp = rep(seq_len(n_var / 2L) * 1e6, times = 2L),
    stringsAsFactors = FALSE
  )
  founders <- founder_haplotypes(
    converted$hap1, converted$hap2, map
  )
  expect_equal(founders$ploidy, 2L)
  expect_equal(founders$n_individuals, 4L)
})

test_that("an uninformative map column name gives an actionable error", {
  phased <- make_phased(n_ind = 4L, n_var_per_chr = 5L, n_chr = 1L)
  renamed <- phased$map
  names(renamed)[names(renamed) == "variant_id"] <- "SNP"
  expect_error(
    founder_haplotypes(phased$hap1, phased$hap2, renamed),
    "read_phased_vcf"
  )
  # Naming the columns explicitly resolves it.
  founders <- founder_haplotypes(
    phased$hap1, phased$hap2, renamed,
    variant_col = "SNP"
  )
  expect_equal(founders$n_variants, 5L)
})

test_that("a missing variant in the map is reported", {
  phased <- make_phased(n_ind = 4L, n_var_per_chr = 5L, n_chr = 1L)
  short_map <- phased$map[-1L, , drop = FALSE]
  expect_error(
    founder_haplotypes(phased$hap1, phased$hap2, short_map),
    "missing 1 variant"
  )
})

test_that("a single-variant chromosome is rejected", {
  phased <- make_phased(n_ind = 4L, n_var_per_chr = 5L, n_chr = 1L)
  phased$map$chromosome[1L] <- 99L
  expect_error(
    founder_haplotypes(phased$hap1, phased$hap2, phased$map),
    "fewer than two variants"
  )
})

test_that("founder_population calibrates trait variance to G", {
  skip_if_not_installed("AlphaSimR")
  skip_on_cran()
  phased <- make_phased(n_ind = 60L, n_var_per_chr = 60L, n_chr = 2L)
  founders <- founder_haplotypes(phased$hap1, phased$hap2, phased$map)
  traits <- c("YLD", "DIS")
  G <- matrix(
    c(1.00, -0.30, -0.30, 0.64), 2,
    dimnames = list(traits, traits)
  )
  setup <- founder_population(
    founders,
    G = G, h2 = c(YLD = 0.3, DIS = 0.5),
    n_qtl_per_chromosome = 25L, seed = 3L
  )
  expect_s3_class(setup, "desiredgainr_sim_setup")
  expect_equal(setup$trait_cols, traits)

  observed <- stats::var(AlphaSimR::bv(setup$founder_pop, simParam = setup$SP))
  # Calibration is stochastic at this founder size, so the tolerance is loose;
  # the check is that the requested scale was honoured, not reproduced exactly.
  expect_equal(observed[1L, 1L], G[1L, 1L], tolerance = 0.5)
  expect_lt(observed[1L, 2L], 0)
})

test_that("founder_population restores the caller's RNG state", {
  skip_if_not_installed("AlphaSimR")
  skip_on_cran()
  phased <- make_phased(n_ind = 20L, n_var_per_chr = 20L, n_chr = 2L)
  founders <- founder_haplotypes(phased$hap1, phased$hap2, phased$map)
  traits <- c("A", "B")
  G <- diag(c(1, 1))
  dimnames(G) <- list(traits, traits)

  set.seed(99)
  before <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  invisible(founder_population(
    founders,
    G = G, h2 = 0.4, n_qtl_per_chromosome = 10L, seed = 5L
  ))
  expect_identical(
    get(".Random.seed", envir = globalenv(), inherits = FALSE), before
  )
})

test_that("simulation records gain and diversity across cycles", {
  skip_if_not_installed("AlphaSimR")
  skip_on_cran()
  phased <- make_phased(n_ind = 60L, n_var_per_chr = 50L, n_chr = 2L)
  founders <- founder_haplotypes(phased$hap1, phased$hap2, phased$map)
  traits <- c("YLD", "DIS")
  G <- matrix(
    c(1.00, -0.25, -0.25, 0.64), 2,
    dimnames = list(traits, traits)
  )
  setup <- founder_population(
    founders,
    G = G, h2 = c(YLD = 0.4, DIS = 0.4),
    n_qtl_per_chromosome = 20L, seed = 3L
  )

  simulation <- simulate_selection_cycles(
    setup,
    desired_gains = c(YLD = 1, DIS = 1),
    n_cycles = 3L, mating_system = "outcross",
    n_parents = 10L, n_crosses = 20L, n_progeny_per_cross = 5L,
    lower_is_better = "DIS", seed = 11L
  )
  expect_s3_class(simulation, "desiredgainr_simulation")
  # n_cycles + 1 rows per trait: cycle 0 is the founder population before any
  # selection, and cycles 1 to n_cycles are transmitted selection responses.
  # Earlier versions recorded n_cycles rows because they measured the
  # candidates of each cycle rather than their selected descendants.
  expect_equal(nrow(simulation$cycle_table), 4L * 2L)
  expect_equal(sort(unique(simulation$cycle_table$cycle)), 0:3)

  founder_rows <- simulation$cycle_table[cycle == 0L]
  expect_true(all(founder_rows$cumulative_gain == 0))
  expect_true(all(is.na(founder_rows$selection_intensity)))

  # Selecting for higher YLD and lower DIS must move both in the intended
  # direction of the favourable space.
  expect_gt(simulation$cumulative_gain[["YLD"]], 0)
  expect_lt(simulation$cumulative_gain[["DIS"]], 0)

  # Gain must accumulate: the last cycle must be further from the founders than
  # the first. This is the robust directional claim at this population size.
  yield_by_cycle <- simulation$cycle_table[
    trait == "YLD"
  ][order(cycle), cumulative_gain]
  expect_length(yield_by_cycle, 4L)
  expect_equal(yield_by_cycle[1L], 0)
  # Index position 2 is cycle 1, the first transmitted response.
  expect_gt(yield_by_cycle[4L], yield_by_cycle[2L])

  # Genetic variance is recorded but its trajectory is deliberately not
  # asserted. The Bulmer reduction at ten parents from one hundred progeny is
  # partly restored by recombination each cycle, and the residual change is
  # smaller than the sampling error of a variance estimated from one hundred
  # individuals over forty quantitative trait loci. Detecting it reliably needs
  # a far larger simulation than a test suite should run.
  variance_by_cycle <- simulation$cycle_table[
    trait == "YLD", genetic_variance
  ]
  expect_length(variance_by_cycle, 4L)
  expect_true(all(is.finite(variance_by_cycle) & variance_by_cycle > 0))
})

test_that("the clonal system requires a dominance setup", {
  skip_if_not_installed("AlphaSimR")
  skip_on_cran()
  phased <- make_phased(n_ind = 30L, n_var_per_chr = 20L, n_chr = 2L)
  founders <- founder_haplotypes(phased$hap1, phased$hap2, phased$map)
  traits <- c("A", "B")
  G <- diag(c(1, 1))
  dimnames(G) <- list(traits, traits)
  additive_only <- founder_population(
    founders,
    G = G, h2 = 0.4, n_qtl_per_chromosome = 10L, seed = 2L
  )
  expect_error(
    simulate_selection_cycles(
      additive_only, c(A = 1, B = 1),
      n_cycles = 1L,
      mating_system = "clonal"
    ),
    "requires a setup built with"
  )
})

test_that("desired-gain direction changes the realised gain profile", {
  skip_if_not_installed("AlphaSimR")
  skip_on_cran()
  phased <- make_phased(n_ind = 60L, n_var_per_chr = 50L, n_chr = 2L)
  founders <- founder_haplotypes(phased$hap1, phased$hap2, phased$map)
  traits <- c("YLD", "DIS")
  G <- matrix(
    c(1.00, -0.40, -0.40, 0.64), 2,
    dimnames = list(traits, traits)
  )
  setup <- founder_population(
    founders,
    G = G, h2 = c(YLD = 0.4, DIS = 0.4),
    n_qtl_per_chromosome = 20L, seed = 3L
  )
  # A single stochastic run cannot support a comparison between directions, so
  # each is averaged over several independent seeds.
  run <- function(d, seeds = c(11L, 12L, 13L, 14L)) {
    gains <- vapply(seeds, function(s) {
      simulate_selection_cycles(
        setup, d,
        n_cycles = 3L, mating_system = "outcross",
        n_parents = 10L, n_crosses = 20L, n_progeny_per_cross = 5L,
        lower_is_better = "DIS", seed = s
      )$cumulative_gain[c("YLD", "DIS")]
    }, numeric(2L))
    stats::setNames(rowMeans(gains), c("YLD", "DIS"))
  }
  balanced <- run(c(YLD = 1, DIS = 1))
  yield_led <- run(c(YLD = 4, DIS = 1))

  # The desired-gain index makes the expected response proportional to the
  # desired-gain vector, so the property it controls is the ratio between
  # traits rather than the absolute gain in either. Comparing ratios also
  # cancels the shared variation in overall response magnitude, which is what
  # makes this comparison detectable at a testable population size.
  ratio <- function(gain) abs(gain[["YLD"]]) / abs(gain[["DIS"]])
  expect_gt(ratio(yield_led), ratio(balanced))

  # Both directions must still improve both traits.
  expect_gt(balanced[["YLD"]], 0)
  expect_lt(balanced[["DIS"]], 0)
  expect_gt(yield_led[["YLD"]], 0)
  expect_lt(yield_led[["DIS"]], 0)
})

test_that("RR-BLUP selection is cross-fitted and reports held-out accuracy", {
  skip_if_not_installed("AlphaSimR")
  skip_on_cran()

  phased <- make_phased(
    n_ind = 40L, n_var_per_chr = 30L, n_chr = 2L, seed = 801L
  )
  founders <- founder_haplotypes(phased$hap1, phased$hap2, phased$map)
  traits <- c("YLD", "DIS")
  G <- matrix(c(1, -0.15, -0.15, 0.7), 2,
    dimnames = list(traits, traits)
  )
  setup <- founder_population(
    founders,
    G = G, h2 = 0.5,
    n_qtl_per_chromosome = 8L,
    n_markers_per_chromosome = 10L,
    seed = 802L
  )
  run <- function() {
    simulate_selection_cycles(
      setup,
      desired_gains = c(YLD = 1, DIS = 0.5),
      lower_is_better = "DIS", mating_system = "outcross",
      n_cycles = 1L, n_parents = 8L, n_crosses = 10L,
      n_progeny_per_cross = 4L,
      prediction = list(method = "rrblup", folds = 4L),
      seed = 803L
    )
  }
  first <- run()
  second <- run()
  expect_equal(first$cumulative_gain, second$cumulative_gain)
  expect_identical(first$prediction$method, "rrblup")
  accuracy <- first$cycle_table[cycle == 1L, prediction_accuracy]
  expect_length(accuracy, length(traits))
  expect_true(all(is.finite(accuracy)))
  expect_true(all(accuracy >= -1 & accuracy <= 1))
  calibration <- first$cycle_table[
    cycle == 1L, prediction_calibration_slope
  ]
  expect_true(all(is.finite(calibration)))
  expect_true(first$provenance$simparam_isolated)
})
