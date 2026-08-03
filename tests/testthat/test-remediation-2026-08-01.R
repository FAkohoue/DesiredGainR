test_that("DGSI validation observations never enter replicate fitting", {
  set.seed(101)
  traits <- c("yield", "quality")
  candidates <- data.frame(
    GenoID = sprintf("g%03d", 1:80),
    yield = rnorm(80), quality = rnorm(80)
  )
  validation_a <- data.frame(yield = rnorm(30), quality = rnorm(30))
  validation_b <- transform(validation_a,
    yield = -20 * yield,
    quality = 10 + quality
  )
  G <- matrix(c(0.6, 0.1, 0.1, 0.5), 2,
    dimnames = list(traits, traits)
  )
  P <- G + diag(c(0.8, 0.7))
  dimnames(P) <- list(traits, traits)

  fit <- function(validation) {
    run_dgsi(
      candidates["GenoID"], candidates, traits,
      dg = c(yield = 0.5, quality = 0.3), G = G, P = P,
      validation_data = validation, n_select = 12L, n_iter = 30L,
      n_rep = 4L, seed = 11L
    )
  }
  a <- fit(validation_a)
  b <- fit(validation_b)

  coefficients_a <- lapply(a$all_reps, `[[`, "b")
  coefficients_b <- lapply(b$all_reps, `[[`, "b")
  expect_equal(coefficients_a, coefficients_b, tolerance = 0)
  expect_identical(a$optimism$selection_rule, "external validation data")
  expect_identical(a$optimism$holdout$training_rows, seq_len(nrow(candidates)))
})

test_that("estimated-P compatibility is invariant to trait units", {
  n <- 60L
  x1 <- as.numeric(scale(seq_len(n), center = TRUE, scale = TRUE)) * 100
  x2 <- as.numeric(scale(rep(c(-1, 1), length.out = n),
    center = TRUE, scale = TRUE
  ))
  candidates <- data.frame(
    GenoID = sprintf("g%03d", seq_len(n)),
    large = x1, small = x2
  )
  traits <- c("large", "small")
  G <- diag(c(10000, 2))
  dimnames(G) <- list(traits, traits)

  fit <- function(frame, covariance) {
    suppressWarnings(run_dgsi(
      frame["GenoID"], frame, traits,
      dg = c(large = 0.3, small = 0.3),
      G = covariance, n_select = 10L, n_iter = 2L, n_rep = 1L,
      allow_incompatible_estimated_P = TRUE, seed = 3L
    ))
  }
  raw <- fit(candidates, G)
  rescaled_candidates <- transform(candidates, large = large / 100)
  transform_matrix <- diag(c(1 / 100, 1))
  rescaled_G <- transform_matrix %*% G %*% transform_matrix
  dimnames(rescaled_G) <- list(traits, traits)
  rescaled <- fit(rescaled_candidates, rescaled_G)

  expect_identical(
    raw$covariance_provenance$compatibility$status,
    rescaled$covariance_provenance$compatibility$status
  )
  expect_equal(
    raw$covariance_provenance$compatibility$
      smallest_standardised_residual_eigenvalue,
    rescaled$covariance_provenance$compatibility$
      smallest_standardised_residual_eigenvalue,
    tolerance = 1e-10
  )
})

test_that("optimizer fingerprints cover complete founders and objectives", {
  skip_if_not_installed("AlphaSimR")
  setup <- list(
    G_target = diag(2), G_realised = diag(2), h2 = c(0.4, 0.4),
    heritability_type = "narrow",
    founders = list(
      haplotypes = matrix(0L, 100, 100),
      gen_map = data.frame(chr = 1L, pos = seq_len(100))
    ),
    n_qtl_per_chromosome = 10L, marker_panel = NULL,
    dominance = FALSE, seed = 1L
  )
  fingerprint <- function(objective, local_setup = setup) {
    DesiredGainR:::.dgr_optimiser_fingerprint(
      local_setup, c("A", "B"), c("A", "B"), list(n_cycles = 2L),
      n_replicates = 2L, seed = 4L, include_diversity = FALSE,
      non_negative = TRUE, mode = "target",
      objective_parameters = list(target_gains = objective),
      search_parameters = list(n_candidates = 20L)
    )
  }
  first <- fingerprint(c(A = 1, B = 1))
  changed_setup <- setup
  changed_setup$founders$haplotypes[100, 100] <- 1L
  expect_false(identical(
    first$founder_haplotypes,
    fingerprint(c(A = 1, B = 1), changed_setup)$founder_haplotypes
  ))
  expect_false(identical(
    first$objective_parameters,
    fingerprint(c(A = 2, B = 1))$objective_parameters
  ))
})

test_that("optimizer pool exhaustion is rejected before simulation", {
  skip_if_not_installed("AlphaSimR")
  fake_setup <- structure(
    list(trait_cols = c("A", "B")),
    class = c("desiredgainr_sim_setup", "list")
  )
  expect_error(
    optimize_desired_gains(
      fake_setup,
      budget = 4L, n_initial = 2L, n_candidates = 1L,
      include_diversity = FALSE, verbose = FALSE
    ),
    "n_initial \\+ n_candidates"
  )
})

test_that("multi-environment stability contrasts feed restricted_index", {
  traits <- c("yield", "protein")
  environments <- c("wet", "dry")
  G <- matrix(c(1, 0.2, 0.2, 0.5), 2,
    dimnames = list(traits, traits)
  )
  P <- G + diag(c(1.2, 0.8))
  dimnames(P) <- list(traits, traits)
  expanded <- expand_environments(
    traits, environments, G, P,
    environment_correlation = 0.5,
    economic_weights = c(yield = 2, protein = 1),
    include_stability = TRUE
  )
  set.seed(8)
  values <- as.data.frame(matrix(
    rnorm(60 * length(expanded$labels)),
    ncol = length(expanded$labels),
    dimnames = list(sprintf("g%03d", 1:60), expanded$labels)
  ))
  fit <- restricted_index(
    values, expanded$labels,
    method = "kempthorne_nordskog",
    G = expanded$G, P = expanded$P,
    economic_weights = expanded$economic_weights,
    constraint_matrix = expanded$stability$constraint_matrix,
    n_select = 12L
  )
  expect_lt(fit$constraint$largest_violation, 1e-8)
})

test_that("widen_environments rejects duplicate genotype-environment keys", {
  duplicated <- data.frame(
    genotype = c("g1", "g1"), environment = c("wet", "wet"),
    yield = c(1, 2)
  )
  expect_error(
    widen_environments(duplicated, "genotype", "environment", "yield"),
    "duplicate records"
  )
})

test_that("covariance propagation uses residual draws end to end", {
  skip_if_not_installed("AlphaSimR")
  set.seed(41)
  n_variant <- 30L
  n_founder <- 16L
  variant_id <- sprintf("v%03d", seq_len(n_variant))
  founder_id <- sprintf("f%02d", seq_len(n_founder))
  draw_haplotype <- function() {
    matrix(
      rbinom(n_variant * n_founder, 1L, 0.5),
      nrow = n_variant,
      dimnames = list(variant_id, founder_id)
    )
  }
  map <- data.frame(
    variant_id = variant_id,
    chromosome = rep(1:2, each = n_variant / 2L),
    position_bp = rep(seq_len(n_variant / 2L) * 1e6, 2L)
  )
  founders <- founder_haplotypes(draw_haplotype(), draw_haplotype(), map)
  traits <- c("yield", "quality")
  G <- matrix(c(1, 0.2, 0.2, 0.6), 2,
    dimnames = list(traits, traits)
  )
  P <- G + matrix(c(1.0, 0.25, 0.25, 0.8), 2,
    dimnames = list(traits, traits)
  )
  setup <- founder_population(
    founders, G,
    h2 = 0.4, n_qtl_per_chromosome = 5L, seed = 7L
  )
  directions <- rbind(c(1, 0), c(0, 1))
  colnames(directions) <- traits

  propagated <- propagate_covariance_uncertainty(
    setup, directions,
    genetic_df = 20L, P = P, residual_df = 20L,
    n_covariance_draws = 2L, n_replicates = 2L, n_cycles = 1L,
    rank_weights = c(yield = 1, quality = 1), seed = 13L,
    mating_system = "outcross", n_parents = 6L, n_crosses = 6L,
    n_progeny_per_cross = 3L
  )
  expect_s3_class(propagated, "desiredgainr_covariance_uncertainty")
  expect_equal(dim(propagated$point_estimate), c(2L, 2L))
  expect_equal(dim(propagated$covariance_draw_mean), c(2L, 2L))
  expect_true(all(is.finite(propagated$variance_components$monte_carlo)))
  expect_identical(
    propagated$rank_churn$note,
    "Computed from the user-declared scalarisation weights."
  )
})
