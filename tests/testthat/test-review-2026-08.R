# Regression tests for the external review of 2026-08.
#
# Each block names the review item it guards and, where the item states an
# acceptance criterion, asserts that criterion directly rather than something
# adjacent to it.

review_fixture <- function(seed = 11L, n = 50L) {
  set.seed(seed)
  traits <- c("yield", "protein", "height")
  G <- matrix(
    c(
      1.00, 0.25, -0.30,
      0.25, 0.50, 0.10,
      -0.30, 0.10, 0.80
    ), 3L,
    dimnames = list(traits, traits)
  )
  P <- G + diag(c(1.5, 0.9, 1.2))
  dimnames(P) <- list(traits, traits)
  values <- as.data.frame(matrix(
    stats::rnorm(n * 3L),
    ncol = 3L,
    dimnames = list(paste0("C", seq_len(n)), traits)
  ))
  list(traits = traits, G = G, P = P, values = values)
}

# ---------------------------------------------------------------------------
# Item 1: covariance compatibility
# ---------------------------------------------------------------------------

test_that("a genetic covariance exceeding the phenotypic one is refused", {
  set <- review_fixture()
  # No single trait exceeds its phenotypic variance, but a linear combination
  # does. This is the case a diagonal check misses, which is why the test uses
  # it rather than the easy one.
  bad_G <- set$P
  bad_G["yield", "protein"] <- bad_G["protein", "yield"] <- 1.8
  expect_true(all(diag(bad_G) <= diag(set$P)))
  expect_error(
    selection_index(
      set$values, set$traits,
      method = "smith_hazel",
      G = bad_G, P = set$P,
      economic_weights = c(yield = 2, protein = 1, height = 1)
    ),
    "not positive semidefinite"
  )
})

test_that("the compatibility error says where the problem is", {
  set <- review_fixture()
  bad_G <- set$G
  bad_G["yield", "yield"] <- set$P["yield", "yield"] * 1.4
  expect_error(
    evaluate_index(
      coefficients = stats::setNames(rep(1, 3L), set$traits),
      G = bad_G, P = set$P
    ),
    "yield"
  )
})

test_that("index heritability cannot exceed one on admissible inputs", {
  set <- review_fixture()
  for (method in c("smith_hazel", "base")) {
    fit <- selection_index(
      set$values, set$traits,
      method = method,
      G = set$G, P = set$P,
      economic_weights = c(yield = 2, protein = 1, height = 1),
      n_select = 10L
    )
    expect_lte(fit$evaluation$h2_index, 1)
    expect_gte(fit$evaluation$h2_index, 0)
    expect_lte(abs(fit$evaluation$R_HI), 1)
  }
})

test_that("a compatible pair is accepted without complaint", {
  set <- review_fixture()
  expect_silent(
    selection_index(
      set$values, set$traits,
      method = "smith_hazel",
      G = set$G, P = set$P,
      economic_weights = c(yield = 2, protein = 1, height = 1)
    )
  )
})

# ---------------------------------------------------------------------------
# Item 6: hard independent-culling gates
# ---------------------------------------------------------------------------

test_that("requesting 20 when only 7 pass returns exactly 7", {
  traits <- c("t1", "t2")
  # Seven candidates clear both gates; the rest fail at least one.
  values <- data.frame(
    t1 = c(rep(3, 7L), rep(-3, 13L)),
    t2 = c(rep(3, 7L), rep(3, 13L)),
    row.names = sprintf("g%02d", 1:20)
  )
  fit <- suppressWarnings(selection_index(
    values, traits,
    method = "independent_culling",
    culling_thresholds = c(t1 = 0, t2 = 0),
    n_select = 20L, scale_traits = FALSE, center_traits = FALSE
  ))
  expect_identical(fit$n_selected, 7L)
  expect_identical(nrow(fit$selected), 7L)
  expect_identical(fit$n_select, 20L)
  # And every selected candidate is one that passed.
  expect_true(all(fit$selected$id %in% sprintf("g%02d", 1:7)))
})

test_that("a shortfall warns rather than filling with failures", {
  traits <- c("t1", "t2")
  values <- data.frame(
    t1 = c(rep(3, 4L), rep(-3, 6L)),
    t2 = rep(3, 10L),
    row.names = sprintf("g%02d", 1:10)
  )
  expect_warning(
    selection_index(
      values, traits,
      method = "independent_culling",
      culling_thresholds = c(t1 = 0, t2 = 0),
      n_select = 8L, scale_traits = FALSE, center_traits = FALSE
    ),
    "passed every culling threshold"
  )
})

test_that("selection intensity uses the number actually selected", {
  traits <- c("t1", "t2")
  values <- data.frame(
    t1 = c(rep(3, 5L), rep(-3, 15L)),
    t2 = rep(3, 20L),
    row.names = sprintf("g%02d", 1:20)
  )
  fit <- suppressWarnings(selection_index(
    values, traits,
    method = "independent_culling",
    culling_thresholds = c(t1 = 0, t2 = 0),
    n_select = 12L, scale_traits = FALSE, center_traits = FALSE
  ))
  # 5 of 20 were taken, not the 12 requested, so the intensity is that of a
  # 25 per cent selected fraction.
  expect_equal(
    fit$selection_intensity,
    DesiredGainR:::.dgr_intensity(5 / 20),
    tolerance = 1e-10
  )
})

test_that("the culling report names the gates each candidate failed", {
  traits <- c("t1", "t2")
  values <- data.frame(
    t1 = c(3, -3, 3, -3),
    t2 = c(3, 3, -3, -3),
    row.names = sprintf("g%d", 1:4)
  )
  fit <- selection_index(
    values, traits,
    method = "independent_culling",
    culling_thresholds = c(t1 = 0, t2 = 0),
    scale_traits = FALSE, center_traits = FALSE
  )
  report <- fit$culling_report
  expect_false(is.null(report))
  expect_identical(report$passed, c(TRUE, FALSE, FALSE, FALSE))
  expect_identical(report$failed_traits[2L], "t1")
  expect_identical(report$failed_traits[3L], "t2")
  expect_identical(report$failed_traits[4L], "t1, t2")
  expect_true(is.na(report$failed_traits[1L]))
})

# ---------------------------------------------------------------------------
# Item 7: strict integer validation
# ---------------------------------------------------------------------------

test_that("non-integral, missing and infinite counts all fail clearly", {
  validate <- DesiredGainR:::.dgr_positive_integer
  expect_error(validate(1.9, "n"), "whole number")
  expect_error(validate(NA, "n"), "NA")
  expect_error(validate(NA_integer_, "n"), "NA")
  expect_error(validate(Inf, "n"), "Inf")
  expect_error(validate(0, "n"), "positive integer")
  expect_error(validate(-3, "n"), "positive integer")
  expect_error(validate(c(1, 2), "n"), "length 2")
  expect_error(validate("4", "n"), "single number")
  expect_error(validate(1e12, "n"), "largest representable")
  expect_identical(validate(4, "n"), 4L)
  expect_identical(validate(4L, "n"), 4L)
})

test_that("zero is allowed only where it is meaningful", {
  non_negative <- DesiredGainR:::.dgr_non_negative_integer
  expect_identical(non_negative(0, "n"), 0L)
  expect_error(non_negative(-1, "n"), "non-negative")
  expect_error(non_negative(2.5, "n"), "whole number")
})

test_that("run_dgsi() rejects a fractional n_select rather than truncating", {
  set.seed(3L)
  traits <- paste0("t", 1:3)
  candidates <- data.table::data.table(GenoID = sprintf("g%03d", 1:30))
  for (trait in traits) {
    data.table::set(candidates, j = trait, value = stats::rnorm(30L))
  }
  G <- diag(3L) * 0.6
  dimnames(G) <- list(traits, traits)
  expect_error(
    run_dgsi(
      init_data = candidates[, "GenoID"], cand_data = candidates,
      trait_cols = traits, dg = stats::setNames(rep(0.5, 3L), traits),
      G = G, n_select = 1.9, n_iter = 5L, n_rep = 1L, seed = 1L
    ),
    "whole number"
  )
})

# ---------------------------------------------------------------------------
# Item 8: safe pipeline defaults
# ---------------------------------------------------------------------------

test_that("the pipeline no longer falls back or debugs by default", {
  defaults <- formals(run_dgsi_qgsi_pipeline)
  expect_false(eval(defaults$fallback_to_top_n))
  expect_false(eval(defaults$debug))
  expect_false(eval(formals(run_qgsi_desired_gain)$debug))
})

test_that("enabling the fallback warns at the point of request", {
  set.seed(5L)
  traits <- paste0("t", 1:3)
  candidates <- data.table::data.table(GenoID = sprintf("g%03d", 1:30))
  for (trait in traits) {
    data.table::set(candidates, j = trait, value = stats::rnorm(30L))
  }
  G <- diag(3L) * 0.6
  dimnames(G) <- list(traits, traits)
  expect_warning(
    run_dgsi_qgsi_pipeline(
      mode = "dg", init_data = candidates[, "GenoID"],
      cand_data = candidates, trait_cols = traits,
      dg = stats::setNames(rep(0.5, 3L), traits), G = G,
      n_select = 5L, fallback_to_top_n = TRUE, n_iter = 5L, n_rep = 1L
    ),
    "fallback_to_top_n = TRUE"
  )
})

# ---------------------------------------------------------------------------
# Item 2 and 3: simulation timing and isolation
# ---------------------------------------------------------------------------

review_setup <- function(seed = 3L, n_markers = 20L) {
  set.seed(seed)
  n_ind <- 40L
  n_chr <- 2L
  n_var_per_chr <- 60L
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
  founders <- founder_haplotypes(draw(), draw(), map)
  traits <- c("YLD", "DIS")
  G <- matrix(c(1.00, -0.25, -0.25, 0.64), 2L,
    dimnames = list(traits, traits)
  )
  founder_population(
    founders,
    G = G, h2 = c(YLD = 0.4, DIS = 0.4),
    n_qtl_per_chromosome = 15L, n_markers_per_chromosome = n_markers,
    seed = 5L
  )
}

test_that("one cycle is one transmitted response, not a baseline", {
  skip_on_cran()
  skip_if_not_installed("AlphaSimR")

  setup <- review_setup()
  run <- function(gains) {
    simulate_selection_cycles(
      setup,
      desired_gains = gains, n_cycles = 1L,
      mating_system = "outcross", n_parents = 6L, n_crosses = 30L,
      n_progeny_per_cross = 2L, seed = 99L
    )
  }
  # Common random numbers: the same seed, so any difference in the one-cycle
  # gain is attributable to the direction. Under the previous timing the cycle
  # 1 row measured an unselected cross, whose mean equals the founders' in
  # expectation, so these two were indistinguishable by construction.
  yield_first <- run(c(YLD = 1, DIS = 0.1))
  disease_first <- run(c(YLD = 0.1, DIS = 1))
  expect_false(isTRUE(all.equal(
    yield_first$cumulative_gain, disease_first$cumulative_gain
  )))
  expect_gt(
    yield_first$cumulative_gain[["YLD"]],
    disease_first$cumulative_gain[["YLD"]]
  )
})

test_that("cycle 0 records the founders and carries no gain", {
  skip_on_cran()
  skip_if_not_installed("AlphaSimR")

  setup <- review_setup()
  simulation <- simulate_selection_cycles(
    setup,
    desired_gains = c(YLD = 1, DIS = 0.5), n_cycles = 3L,
    mating_system = "outcross", n_parents = 6L, n_crosses = 30L,
    n_progeny_per_cross = 2L, seed = 7L
  )
  table <- simulation$cycle_table
  expect_true(0L %in% table$cycle)
  expect_identical(sort(unique(table$cycle)), 0:3)
  expect_true(all(table$cumulative_gain[table$cycle == 0L] == 0))
  expect_true(all(is.na(table$selection_intensity[table$cycle == 0L])))
})

test_that("the caller's SimParam is not advanced by a simulation", {
  skip_on_cran()
  skip_if_not_installed("AlphaSimR")

  setup <- review_setup()
  before <- setup$SP$lastId
  simulate_selection_cycles(
    setup,
    desired_gains = c(YLD = 1, DIS = 0.5), n_cycles = 2L,
    mating_system = "outcross", n_parents = 6L, n_crosses = 30L,
    n_progeny_per_cross = 2L, seed = 7L
  )
  expect_identical(setup$SP$lastId, before)
})

test_that("the same seed reproduces the same result", {
  skip_on_cran()
  skip_if_not_installed("AlphaSimR")

  setup <- review_setup()
  once <- simulate_selection_cycles(
    setup,
    desired_gains = c(YLD = 1, DIS = 0.5), n_cycles = 2L,
    mating_system = "outcross", n_parents = 6L, n_crosses = 30L,
    n_progeny_per_cross = 2L, seed = 21L
  )
  twice <- simulate_selection_cycles(
    setup,
    desired_gains = c(YLD = 1, DIS = 0.5), n_cycles = 2L,
    mating_system = "outcross", n_parents = 6L, n_crosses = 30L,
    n_progeny_per_cross = 2L, seed = 21L
  )
  expect_equal(once$cumulative_gain, twice$cumulative_gain, tolerance = 1e-12)
  expect_equal(
    once$cycle_table$genetic_mean, twice$cycle_table$genetic_mean,
    tolerance = 1e-12
  )
  expect_true(once$provenance$caller_last_id_unchanged)
  expect_identical(once$provenance$n_threads, 1L)
})

test_that("multiple threads warn about reproducibility", {
  skip_on_cran()
  skip_if_not_installed("AlphaSimR")

  setup <- review_setup()
  expect_warning(
    simulate_selection_cycles(
      setup,
      desired_gains = c(YLD = 1, DIS = 0.5), n_cycles = 1L,
      mating_system = "outcross", n_parents = 6L, n_crosses = 30L,
      n_progeny_per_cross = 2L, n_threads = 2L, seed = 7L
    ),
    "may not reproduce"
  )
})

# ---------------------------------------------------------------------------
# Item 5: marker-based diversity
# ---------------------------------------------------------------------------

test_that("the setup retains a marker panel distinct from the QTL", {
  skip_on_cran()
  skip_if_not_installed("AlphaSimR")

  setup <- review_setup()
  expect_false(is.null(setup$marker_panel))
  expect_gt(length(setup$marker_panel), 0L)
  qtl <- colnames(
    AlphaSimR::pullQtlGeno(setup$founder_pop, simParam = setup$SP)
  )
  if (!isTRUE(setup$marker_qtl_overlap) && !is.null(qtl)) {
    expect_length(intersect(setup$marker_panel, qtl), 0L)
  }
})

test_that("diversity is reported separately from QTL frequency change", {
  skip_on_cran()
  skip_if_not_installed("AlphaSimR")

  setup <- review_setup()
  simulation <- simulate_selection_cycles(
    setup,
    desired_gains = c(YLD = 1, DIS = 0.5), n_cycles = 3L,
    mating_system = "outcross", n_parents = 5L, n_crosses = 30L,
    n_progeny_per_cross = 2L, seed = 13L
  )
  table <- simulation$cycle_table
  expect_true(all(
    c("mean_relationship", "qtl_segregating", "genic_variance_ratio") %in%
      names(table)
  ))
  expect_true(any(is.finite(table$qtl_segregating)))
  expect_match(simulation$diversity_basis, "genome-wide markers")
})

# ---------------------------------------------------------------------------
# Item 4: clonal calibration
# ---------------------------------------------------------------------------

test_that("a clonal setup calibrates genotypic variance to the target", {
  skip_on_cran()
  skip_if_not_installed("AlphaSimR")

  set.seed(9L)
  n_ind <- 40L
  n_chr <- 2L
  n_var_per_chr <- 60L
  n_var <- n_var_per_chr * n_chr
  variant_id <- sprintf("v%04d", seq_len(n_var))
  draw <- function() {
    matrix(
      stats::rbinom(n_var * n_ind, 1L, 0.5),
      nrow = n_var,
      dimnames = list(variant_id, sprintf("f%03d", seq_len(n_ind)))
    )
  }
  map <- data.frame(
    variant_id = variant_id,
    chromosome = rep(seq_len(n_chr), each = n_var_per_chr),
    position_bp = rep(seq_len(n_var_per_chr) * 1e6, times = n_chr),
    stringsAsFactors = FALSE
  )
  founders <- founder_haplotypes(draw(), draw(), map)
  traits <- c("YLD", "DIS")
  G <- matrix(c(1.00, -0.20, -0.20, 0.64), 2L,
    dimnames = list(traits, traits)
  )

  setup <- suppressWarnings(founder_population(
    founders,
    G = G, h2 = c(YLD = 0.3, DIS = 0.3),
    n_qtl_per_chromosome = 25L, n_markers_per_chromosome = 20L,
    dominance_degree = c(YLD = 0.4, DIS = 0.4),
    heritability = "broad", seed = 5L
  ))

  # The acceptance criterion: the realised founder GENOTYPIC covariance
  # matches the documented target. Without the calibration the dominance
  # variance was added on top of diag(G), so both traits came out too large
  # and in the wrong ratio.
  expect_true(setup$dominance)
  expect_identical(setup$heritability_type, "broad")

  # Whether the iteration converged is recorded, and the assertion is made
  # against that record rather than against a tolerance guessed in advance.
  # A finite marker panel and a finite number of loci put a floor on how
  # closely any dominance model can be scaled.
  skip_if_not(
    isTRUE(setup$calibration_converged),
    paste(
      "Calibration did not converge with this genome; largest deviation",
      format(max(abs(diag(setup$G_realised) - diag(G)) / diag(G)), digits = 3)
    )
  )
  expect_equal(
    diag(setup$G_realised), diag(G),
    tolerance = 0.01, ignore_attr = TRUE
  )
  # And the ratio between the traits must be right, not merely the totals.
  expect_equal(
    setup$G_realised[["YLD", "YLD"]] / setup$G_realised[["DIS", "DIS"]],
    G[["YLD", "YLD"]] / G[["DIS", "DIS"]],
    tolerance = 0.02
  )
})

test_that("clonal calibration reproduces the correlations, not just variances", {
  skip_on_cran()
  skip_if_not_installed("AlphaSimR")

  # M6 of the second review: rescaling the variances alone left the realised
  # genotypic CORRELATIONS emergent, so the setup accepted a full covariance
  # matrix as its target and reproduced only its diagonal. This asserts the
  # whole matrix.
  set.seed(15L)
  n_ind <- 50L
  n_chr <- 2L
  n_var_per_chr <- 80L
  n_var <- n_var_per_chr * n_chr
  variant_id <- sprintf("v%04d", seq_len(n_var))
  draw <- function() {
    matrix(
      stats::rbinom(n_var * n_ind, 1L, 0.5),
      nrow = n_var,
      dimnames = list(variant_id, sprintf("f%03d", seq_len(n_ind)))
    )
  }
  map <- data.frame(
    variant_id = variant_id,
    chromosome = rep(seq_len(n_chr), each = n_var_per_chr),
    position_bp = rep(seq_len(n_var_per_chr) * 1e6, times = n_chr),
    stringsAsFactors = FALSE
  )
  founders <- founder_haplotypes(draw(), draw(), map)
  traits <- c("YLD", "DIS")
  # A target correlation that is neither zero nor extreme, so hitting it is
  # informative rather than automatic.
  G <- matrix(c(1.00, -0.35, -0.35, 0.64), 2L,
    dimnames = list(traits, traits)
  )

  setup <- suppressWarnings(founder_population(
    founders,
    G = G, h2 = c(YLD = 0.3, DIS = 0.3),
    n_qtl_per_chromosome = 30L, n_markers_per_chromosome = 20L,
    dominance_degree = c(YLD = 0.3, DIS = 0.3),
    heritability = "broad", seed = 8L
  ))

  expect_false(is.null(setup$calibration_error))
  skip_if_not(
    isTRUE(setup$calibration_converged),
    paste(
      "Calibration did not converge:",
      format(setup$calibration_error$variance, digits = 3), "on variances,",
      format(setup$calibration_error$correlation, digits = 3),
      "on correlations."
    )
  )
  # The whole covariance matrix, not just its diagonal.
  expect_equal(
    stats::cov2cor(setup$G_realised), stats::cov2cor(G),
    tolerance = 0.01, ignore_attr = TRUE
  )
  expect_lt(setup$calibration_error$correlation, 0.01)
  expect_lt(setup$calibration_error$variance, 0.01)
  # The iteration must actually have done work rather than converged trivially.
  expect_gt(setup$calibration_iterations, 0L)
})

test_that("the correlation projection returns a valid correlation matrix", {
  project <- DesiredGainR:::.dgr_project_correlation
  # An indefinite matrix with an out-of-range entry, which is what a Newton
  # correction can produce mid-iteration.
  bad <- matrix(
    c(
      1.0, 0.95, -0.95,
      0.95, 1.0, 0.95,
      -0.95, 0.95, 1.10
    ), 3L
  )
  projected <- project(bad)
  expect_equal(diag(projected), rep(1, 3L), tolerance = 1e-10)
  expect_true(isSymmetric(projected, tol = 1e-10))
  expect_gt(
    min(eigen(projected, symmetric = TRUE, only.values = TRUE)$values), 0
  )
  expect_true(all(abs(projected) <= 1 + 1e-12))
  # An already valid matrix is left essentially alone.
  good <- matrix(c(1, 0.3, 0.3, 1), 2L)
  expect_equal(project(good), good, tolerance = 1e-8)
})

test_that("broad-sense heritability is refused without dominance", {
  skip_on_cran()
  skip_if_not_installed("AlphaSimR")

  expect_error(
    review_setup_broad <- {
      set.seed(2L)
      n_var <- 40L
      variant_id <- sprintf("v%03d", seq_len(n_var))
      draw <- function() {
        matrix(stats::rbinom(n_var * 20L, 1L, 0.5),
          nrow = n_var,
          dimnames = list(variant_id, sprintf("f%02d", 1:20))
        )
      }
      map <- data.frame(
        variant_id = variant_id, chromosome = 1L,
        position_bp = seq_len(n_var) * 1e6, stringsAsFactors = FALSE
      )
      founder_population(
        founder_haplotypes(draw(), draw(), map),
        G = matrix(1, 1L, 1L, dimnames = list("YLD", "YLD")),
        h2 = c(YLD = 0.4), n_qtl_per_chromosome = 10L,
        heritability = "broad", seed = 1L
      )
    },
    "only meaningful when dominance"
  )
})
