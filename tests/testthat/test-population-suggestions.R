test_that("the exact active-set solver satisfies known convex optima", {
  identity_solution <- DesiredGainR:::.dgr_quadratic_lower_bound(
    diag(2), c(1, 0.5)
  )
  expect_equal(identity_solution$solution, c(1, 0.5), tolerance = 1e-12)
  expect_equal(identity_solution$objective, 1.25, tolerance = 1e-12)
  expect_true(identity_solution$globally_optimal)

  # With negative cross-curvature, raising the second coordinate above its
  # zero lower bound reduces the quadratic cost. The exact KKT solution is
  # (1, 1/2), not the superficial corner (1, 0).
  B <- matrix(c(2, -1, -1, 2), 2)
  correlated <- DesiredGainR:::.dgr_quadratic_lower_bound(B, c(1, 0))
  expect_equal(correlated$solution, c(1, 0.5), tolerance = 1e-10)
  expect_equal(correlated$objective, 1.5, tolerance = 1e-10)
  expect_lt(correlated$kkt_residual, 1e-10)
})

test_that("population scaling uses the breeder's target covariance", {
  traits <- c("A", "B")
  target <- diag(c(4, 9))
  realised <- diag(c(1, 16))
  dimnames(target) <- dimnames(realised) <- list(traits, traits)
  setup <- list(
    trait_cols = traits, G_target = target, G_realised = realised,
    residual_covariance = diag(c(2, 3)), dominance = FALSE,
    heritability_type = "narrow", h2 = c(A = 0.5, B = 0.5)
  )
  expect_equal(
    DesiredGainR:::.dgr_population_genetic_covariance(setup), target
  )
  expect_equal(
    DesiredGainR:::.dgr_setup_phenotypic_covariance(setup, target),
    target + setup$residual_covariance
  )

  setup$residual_covariance <- NULL
  setup$dominance <- TRUE
  expect_error(
    DesiredGainR:::.dgr_setup_phenotypic_covariance(setup, target),
    "cannot be recovered exactly"
  )
})

test_that("realised covariance is a calibration diagnostic, not a unit switch", {
  traits <- c("A", "B")
  target <- diag(c(1, 4))
  realised <- matrix(c(1.02, 0.02, 0.02, 3.9), 2)
  dimnames(target) <- dimnames(realised) <- list(traits, traits)
  diagnostic <- DesiredGainR:::.dgr_population_model_diagnostic(list(
    trait_cols = traits, G_target = target, G_realised = realised
  ))
  expect_identical(diagnostic$status, "adequate")
  expect_true(diagnostic$calibrated)

  bad <- realised
  bad[2, 2] <- 8
  diagnostic <- DesiredGainR:::.dgr_population_model_diagnostic(list(
    trait_cols = traits, G_target = target, G_realised = bad
  ))
  expect_identical(diagnostic$status, "inadequate")
})

test_that("support decisions use simultaneous lower and upper bounds", {
  expect_identical(
    DesiredGainR:::.dgr_support_status(
      list(
        joint_exact_simultaneous_lower = 0.82,
        joint_exact_simultaneous_upper = 0.95
      ),
      0.8
    ),
    "supported"
  )
  expect_identical(
    DesiredGainR:::.dgr_support_status(
      list(
        joint_exact_simultaneous_lower = 0.3,
        joint_exact_simultaneous_upper = 0.7
      ),
      0.8
    ),
    "not_supported"
  )
  expect_identical(
    DesiredGainR:::.dgr_support_status(
      list(
        joint_exact_simultaneous_lower = 0.7,
        joint_exact_simultaneous_upper = 0.9
      ),
      0.8
    ),
    "uncertain"
  )
})

test_that("exact probability and median bounds are conservative", {
  expect_equal(
    DesiredGainR:::.dgr_simultaneous_binomial_lower(0, 50, 0.05, 10),
    0
  )
  expect_equal(
    DesiredGainR:::.dgr_simultaneous_binomial_lower(50, 50, 0.05, 10),
    stats::qbeta(0.005, 50, 1),
    tolerance = 1e-14
  )
  expect_equal(
    DesiredGainR:::.dgr_simultaneous_binomial_upper(0, 50, 0.05, 10),
    stats::qbeta(0.995, 1, 50),
    tolerance = 1e-14
  )

  x <- seq_len(100)
  bound <- DesiredGainR:::.dgr_simultaneous_median_lower(x, 0.05, 5)
  k <- as.integer(bound)
  expect_lte(stats::pbinom(k - 1L, 100, 0.5), 0.01)
  if (k < 100L) {
    expect_gt(stats::pbinom(k, 100, 0.5), 0.01)
  }
})

test_that("population summaries honour trait-specific minima and direction", {
  traits <- c("Yield", "Disease")
  directions <- rbind(c(1, 0), c(0, 1))
  colnames(directions) <- traits
  replicates <- list(
    rbind(Yield = c(1.2, 1.1, 0.9, 1.3), Disease = c(-0.6, -0.7, -0.8, -0.4)),
    rbind(Yield = c(0.3, 0.4, 0.2, 0.1), Disease = c(-1.2, -1.1, -1.3, -1.0))
  )
  result <- DesiredGainR:::.dgr_population_candidate_summary(
    directions, replicates, traits,
    minimum_gains = c(Yield = 1, Disease = 0.5),
    scale = c(Yield = 1, Disease = 1),
    direction = c(Yield = 1, Disease = -1),
    level = 0.95, seed = 10
  )$table

  expect_equal(result$joint_probability, c(2.5 / 5, 0.5 / 5))
  expect_equal(result$minimum_probability_Yield, c(3.5 / 5, 0.5 / 5))
  expect_equal(result$minimum_probability_Disease, c(3.5 / 5, 4.5 / 5))
  expect_equal(result$mean_gain_Yield, c(1.125, 0.25))
  expect_gt(
    result$median_worst_trait_gain[1],
    result$median_worst_trait_gain[2]
  )
})

test_that("suggest_desired_gains validates the compact public interface", {
  skip_if_not_installed("AlphaSimR")
  fake_setup <- structure(
    list(trait_cols = c("A", "B")),
    class = c("desiredgainr_sim_setup", "list")
  )
  expect_error(
    suggest_desired_gains(fake_setup, c(A = 1, B = -0.2)),
    "non-negative"
  )
  expect_error(
    suggest_desired_gains(
      fake_setup, c(A = 1, B = 0.2),
      control = list(replictes = 20)
    ),
    "Unknown control field"
  )
  expect_error(
    suggest_desired_gains(
      fake_setup, c(A = 1, B = 0.2),
      programme = list(parent_count = 20)
    ),
    "Unknown programme field"
  )
})

test_that("population-driven recommendation runs end to end", {
  skip_if_not_installed("AlphaSimR")
  skip_on_cran()

  set.seed(71)
  n_var <- 40L
  n_founders <- 24L
  variants <- sprintf("v%03d", seq_len(n_var))
  individuals <- sprintf("f%02d", seq_len(n_founders))
  draw <- function() {
    matrix(
      stats::rbinom(n_var * n_founders, 1L, 0.5),
      nrow = n_var, dimnames = list(variants, individuals)
    )
  }
  map <- data.frame(
    variant_id = variants,
    chromosome = rep(1:2, each = n_var / 2L),
    position_bp = rep(seq_len(n_var / 2L) * 1e6, 2L)
  )
  founders <- founder_haplotypes(draw(), draw(), map)
  traits <- c("YLD", "DIS")
  G <- matrix(c(1, -0.2, -0.2, 0.64), 2,
    dimnames = list(traits, traits)
  )
  setup <- founder_population(
    founders,
    G = G, h2 = 0.5,
    n_qtl_per_chromosome = 8L, seed = 72L
  )

  suppressWarnings(
    result <- suggest_desired_gains(
      setup,
      minimum_gains = c(YLD = 0.2, DIS = 0.1),
      lower_is_better = "DIS",
      programme = list(
        mating_system = "outcross", n_parents = 6L,
        n_crosses = 8L, n_progeny_per_cross = 3L
      ),
      control = list(
        n_cycles = 1L, budget = 3L, n_initial = 3L,
        n_replicates = 2L, n_candidates = 30L,
        screening_replicates = 3L, confirmation_replicates = 3L,
        confirmation_finalists = 1L, search_starts = 2L,
        uncertainty = list(
          architecture_draws = 2L, covariance_draws = 2L,
          genetic_df = 20L, residual_df = 40L
        ),
        verbose = FALSE, seed = 73L
      )
    )
  )
  expect_s3_class(result, "desiredgainr_gain_suggestion")
  expect_equal(result$minimum_gains, c(YLD = 0.2, DIS = 0.1))
  expect_true(result$analytical_feasibility$maximum_common_gain > 0)
  expect_lt(result$analytical_feasibility$minimum_kkt_residual, 1e-7)
  expect_equal(
    sqrt(sum(result$minimum_recommendation$desired_gain_direction^2)),
    1,
    tolerance = 1e-10
  )
  expect_true("screening" %in% names(result$minimum_recommendation))
  expect_false("discovery" %in% names(result$minimum_recommendation))
  expect_true(all(c(
    "joint_probability", "joint_exact_simultaneous_lower",
    "median_worst_exact_simultaneous_lower"
  ) %in% names(result$candidate_results)))
  expect_true(nrow(result$confirmation) %in% 1:2)
  expect_true(result$minimum_recommendation$decision_status %in%
    c("supported", "uncertain", "not_supported"))
  expect_true(result$search_stability$status %in% c("resolved", "unresolved"))
  expect_identical(result$outer_validation$architecture$status, "assessed")
  expect_identical(result$outer_validation$covariance$status, "assessed")
})
