# Tests for the Gaussian-process surrogate and the desired-gain search.
#
# The surrogate is tested on deterministic analytic functions, which is where
# its correctness can actually be established. The search itself is tested
# against AlphaSimR only in a small smoke test, because a full run is far too
# slow for a test suite.

test_that("the Halton sequence spreads points through the unit cube", {
  points <- DesiredGainR:::.dgr_halton(64L, 3L)
  expect_equal(dim(points), c(64L, 3L))
  expect_true(all(points > 0 & points < 1))
  # A quasi-random sequence should cover each margin far more evenly than the
  # extremes of an independent uniform sample would.
  for (j in seq_len(3L)) {
    counts <- table(cut(points[, j], breaks = seq(0, 1, by = 0.25)))
    expect_true(all(counts >= 10L))
  }
})

test_that("directions lie on the unit sphere and honour the orthant", {
  cube <- DesiredGainR:::.dgr_halton(50L, 4L)
  non_negative <- DesiredGainR:::.dgr_directions_from_cube(cube, TRUE)
  expect_equal(sqrt(rowSums(non_negative^2)), rep(1, 50L), tolerance = 1e-10)
  expect_true(all(non_negative >= 0))

  unconstrained <- DesiredGainR:::.dgr_directions_from_cube(cube, FALSE)
  expect_equal(sqrt(rowSums(unconstrained^2)), rep(1, 50L), tolerance = 1e-10)
  expect_true(any(unconstrained < 0))
})

test_that("interval-constrained directions stay inside the admissible cone", {
  intervals <- rbind(
    A = c(lower = 0.5, upper = 1.5),
    B = c(lower = 0.25, upper = 0.75),
    C = c(lower = 0.1, upper = 0.8)
  )
  cube <- DesiredGainR:::.dgr_halton(80L, 3L)
  directions <- DesiredGainR:::.dgr_directions_from_intervals(cube, intervals)

  expect_equal(sqrt(rowSums(directions^2)), rep(1, 80L), tolerance = 1e-10)
  # Each normalised ray must intersect the original interval box. The lower
  # and upper bounds imply an admissible interval for its common multiplier.
  feasible <- apply(directions, 1L, function(direction) {
    multiplier_lower <- max(intervals[, "lower"] / direction)
    multiplier_upper <- min(intervals[, "upper"] / direction)
    multiplier_lower <= multiplier_upper + 1e-12
  })
  expect_true(all(feasible))
})

test_that("interval validation catches empty and fixed search spaces", {
  traits <- c("A", "B")
  expect_error(
    DesiredGainR:::.dgr_validate_gain_intervals(
      rbind(A = c(lower = 0, upper = 0), B = c(lower = 0, upper = 0)),
      traits, TRUE
    ),
    "every trait to zero"
  )
  expect_error(
    DesiredGainR:::.dgr_validate_gain_intervals(
      rbind(A = c(lower = 1, upper = 1), B = c(lower = 0.5, upper = 0.5)),
      traits, TRUE
    ),
    "single direction"
  )
  expect_error(
    DesiredGainR:::.dgr_validate_gain_intervals(
      rbind(A = c(lower = -1, upper = 1), B = c(lower = 0, upper = 1)),
      traits, TRUE
    ),
    "non_negative = FALSE"
  )
})

test_that("interval probabilities are exact Jeffreys-binomial summaries", {
  traits <- c("Yield", "Disease")
  intervals <- rbind(
    Yield = c(lower = 0.5, upper = 1.5),
    Disease = c(lower = 0.25, upper = 1.0)
  )
  # Rows are raw cumulative responses; disease is lower-is-better.
  always <- rbind(
    Yield = c(0.7, 0.8, 0.9, 1.0),
    Disease = c(-0.4, -0.5, -0.6, -0.7)
  )
  half <- rbind(
    Yield = c(0.7, 0.8, 0.2, 0.3),
    Disease = c(-0.4, -0.5, -0.6, -0.7)
  )
  set.seed(100)
  result <- DesiredGainR:::.dgr_interval_attainment(
    list(always, half), traits, intervals,
    outcome_scale = c(Yield = 1, Disease = 1),
    favourable_direction = c(Yield = 1, Disease = -1),
    level = 0.95
  )

  expect_equal(result$joint_probability, c(4.5 / 5, 2.5 / 5))
  expect_equal(
    result$joint_lower[1], stats::qbeta(0.025, 4.5, 0.5),
    tolerance = 1e-12
  )
  expect_equal(
    result$per_trait_probability[1, ],
    c(Yield = 0.9, Disease = 0.9)
  )
  expect_equal(
    result$per_trait_probability[2, ],
    c(Yield = 0.5, Disease = 0.9)
  )
  expect_equal(
    sum(result$bootstrap_selection_frequency), 1,
    tolerance = 1e-12
  )
})

test_that("non-dominated filtering is correct on a known frontier", {
  objectives <- rbind(
    c(1, 3),
    c(2, 2),
    c(3, 1),
    c(1, 1), # dominated by every other row
    c(2, 1) # dominated by c(3, 1) and c(2, 2)
  )
  keep <- DesiredGainR:::.dgr_non_dominated(objectives)
  expect_equal(keep, c(TRUE, TRUE, TRUE, FALSE, FALSE))
})

test_that("the Gaussian process interpolates a smooth function", {
  set.seed(11)
  X <- DesiredGainR:::.dgr_halton(60L, 2L)
  truth <- function(x) sin(3 * x[, 1L]) + 0.5 * x[, 2L]^2
  y <- truth(X)

  model <- DesiredGainR:::.dgr_gp_fit(X, y)
  fitted <- DesiredGainR:::.dgr_gp_predict(model, X)
  # With a negligible nugget the posterior mean must reproduce the training
  # data closely, and the posterior standard deviation must be small there.
  expect_lt(max(abs(fitted$mean - y)), 0.05)
  expect_lt(max(fitted$sd), 0.2)

  X_new <- DesiredGainR:::.dgr_halton(25L, 2L, skip = 300L)
  predicted <- DesiredGainR:::.dgr_gp_predict(model, X_new)
  expect_lt(
    sqrt(mean((predicted$mean - truth(X_new))^2)),
    0.1
  )
  expect_true(all(predicted$sd >= 0))
})

test_that("the Gaussian process reverts to the prior away from the data", {
  set.seed(12)
  X <- matrix(stats::runif(40L, 0, 0.2), ncol = 2L)
  y <- X[, 1L] + X[, 2L]
  model <- DesiredGainR:::.dgr_gp_fit(X, y)

  near <- DesiredGainR:::.dgr_gp_predict(model, matrix(c(0.1, 0.1), nrow = 1L))
  far <- DesiredGainR:::.dgr_gp_predict(model, matrix(c(5, 5), nrow = 1L))
  # Uncertainty must grow with distance from the observations, which is what
  # keeps the acquisition function exploring.
  expect_gt(far$sd, near$sd)
})

test_that("expected improvement behaves at the boundaries", {
  # No uncertainty and no improvement gives no expected improvement.
  expect_equal(
    DesiredGainR:::.dgr_expected_improvement(1, 0, best = 2, xi = 0), 0
  )
  # Uncertainty alone creates expected improvement even at the incumbent.
  expect_gt(
    DesiredGainR:::.dgr_expected_improvement(2, 1, best = 2, xi = 0), 0
  )
  # Improvement increases monotonically with the posterior mean.
  values <- DesiredGainR:::.dgr_expected_improvement(
    c(1, 2, 3), rep(0.5, 3L),
    best = 1, xi = 0
  )
  expect_true(all(diff(values) > 0))
})

test_that("feasibility probability tracks the floor", {
  expect_equal(
    DesiredGainR:::.dgr_probability_feasible(1, 1, floor = 1), 0.5,
    tolerance = 1e-8
  )
  expect_gt(
    DesiredGainR:::.dgr_probability_feasible(3, 1, floor = 1),
    DesiredGainR:::.dgr_probability_feasible(0, 1, floor = 1)
  )
  # With no uncertainty the probability collapses to an indicator.
  expect_equal(
    DesiredGainR:::.dgr_probability_feasible(c(0, 2), c(0, 0), floor = 1),
    c(0, 1)
  )
})

test_that("expected improvement finds the optimum of an analytic function", {
  skip_on_cran()
  set.seed(13)
  # A deterministic two-dimensional problem with its maximum at (0.3, 0.7).
  objective <- function(x) {
    -((x[, 1L] - 0.3)^2 + (x[, 2L] - 0.7)^2)
  }
  X <- DesiredGainR:::.dgr_halton(12L, 2L)
  y <- objective(X)
  pool <- DesiredGainR:::.dgr_halton(1500L, 2L, skip = 700L)

  for (iteration in seq_len(20L)) {
    model <- DesiredGainR:::.dgr_gp_fit(X, y)
    prediction <- DesiredGainR:::.dgr_gp_predict(model, pool)
    acquisition <- DesiredGainR:::.dgr_expected_improvement(
      prediction$mean, prediction$sd, max(y)
    )
    chosen <- which.max(acquisition)
    X <- rbind(X, pool[chosen, , drop = FALSE])
    y <- c(y, objective(pool[chosen, , drop = FALSE]))
  }
  best <- X[which.max(y), ]
  expect_lt(sqrt(sum((best - c(0.3, 0.7))^2)), 0.1)
})

test_that("the search runs end to end and returns a frontier", {
  skip_if_not_installed("AlphaSimR")
  skip_on_cran()

  set.seed(21)
  n_var <- 60L
  variant_id <- sprintf("v%03d", seq_len(n_var))
  individual_id <- sprintf("f%02d", seq_len(40L))
  draw <- function() {
    matrix(
      stats::rbinom(n_var * 40L, 1L, 0.5),
      nrow = n_var,
      dimnames = list(variant_id, individual_id)
    )
  }
  map <- data.frame(
    variant_id = variant_id,
    chromosome = rep(1:2, each = n_var / 2L),
    position_bp = rep(seq_len(n_var / 2L) * 1e6, times = 2L),
    stringsAsFactors = FALSE
  )
  founders <- founder_haplotypes(draw(), draw(), map)
  traits <- c("YLD", "DIS")
  G <- matrix(
    c(1.00, -0.30, -0.30, 0.64), 2,
    dimnames = list(traits, traits)
  )
  setup <- founder_population(
    founders,
    G = G, h2 = 0.4, n_qtl_per_chromosome = 15L, seed = 3L
  )
  intervals <- define_desired_gain_intervals(
    lower = c(YLD = 0.4, DIS = -1.0),
    upper = c(YLD = 1.5, DIS = -0.2),
    lower_is_better = "DIS",
    gain_units = "genetic_sd"
  )

  optimisation <- optimize_desired_gains(
    setup,
    n_cycles = 2L, mode = "pareto",
    budget = 8L, n_initial = 6L, n_replicates = 1L,
    n_candidates = 200L, include_diversity = FALSE,
    desired_gain_intervals = intervals,
    mating_system = "outcross", n_parents = 8L, n_crosses = 12L,
    n_progeny_per_cross = 4L, lower_is_better = "DIS",
    verbose = FALSE, seed = 5L
  )

  expect_s3_class(optimisation, "desiredgainr_optimisation")
  expect_identical(optimisation$desired_gain_units, "genetic_sd")
  expect_equal(
    optimisation$simulation_directions,
    sweep(optimisation$directions, 2L, sqrt(diag(G)), "*")
  )
  expect_equal(nrow(optimisation$directions), 8L)
  expect_equal(sqrt(rowSums(optimisation$directions^2)), rep(1, 8L),
    tolerance = 1e-8
  )
  expect_true(nrow(optimisation$pareto_set) >= 1L)
  expect_true(all(
    c("observed_YLD", "posterior_YLD", "pareto_optimal") %in%
      names(optimisation$results)
  ))
  # The Pareto flag must be derived from the posterior columns.
  favourable_posterior <- optimisation$posterior_objectives
  favourable_posterior[, "DIS"] <- -favourable_posterior[, "DIS"]
  expect_equal(
    optimisation$results$pareto_optimal,
    DesiredGainR:::.dgr_non_dominated(favourable_posterior)
  )

  expect_warning(
    interval_result <- optimize_desired_gains(
      setup,
      n_cycles = 1L, mode = "interval",
      budget = 4L, n_initial = 3L, n_replicates = 2L,
      n_candidates = 80L, desired_gain_intervals = intervals,
      mating_system = "outcross", n_parents = 8L, n_crosses = 12L,
      n_progeny_per_cross = 4L, lower_is_better = "DIS",
      verbose = FALSE, seed = 6L
    ),
    "fewer than 20 replicates"
  )
  expect_s3_class(interval_result, "desiredgainr_optimisation")
  expect_true(all(interval_result$recommendations$
    joint_high_gain_probability > 0))
  expect_true(all(interval_result$recommendations$
    joint_high_gain_probability < 1))
  expect_equal(length(interval_result$recommended_probability), 3L)
  expect_identical(interval_result$mode, "interval")
})

test_that("scalar modes validate their required arguments", {
  skip_if_not_installed("AlphaSimR")
  fake_setup <- structure(
    list(trait_cols = c("A", "B")),
    class = c("desiredgainr_sim_setup", "list")
  )
  expect_error(
    optimize_desired_gains(fake_setup, mode = "economic", budget = 2L),
    "economic_weights"
  )
  expect_error(
    optimize_desired_gains(fake_setup, mode = "target", budget = 2L),
    "target_gains"
  )
  expect_error(
    optimize_desired_gains(fake_setup, mode = "interval", budget = 2L),
    "requires desired_gain_intervals"
  )
  horizon_intervals <- define_desired_gain_intervals(
    lower = c(A = 0.1, B = 0.1),
    upper = c(A = 1, B = 1),
    horizon_cycles = 3
  )
  expect_error(
    optimize_desired_gains(
      fake_setup,
      mode = "interval", budget = 2L, n_cycles = 2L,
      desired_gain_intervals = horizon_intervals
    ),
    "defined for 3 cycles"
  )
  expect_error(
    optimize_desired_gains(
      fake_setup,
      mode = "constrained", budget = 2L, focal_trait = "Z"
    ),
    "focal_trait must name"
  )
})

test_that("covariance propagation rejects legacy optimisation coordinates", {
  skip_if_not_installed("AlphaSimR")
  fake_setup <- structure(
    list(trait_cols = c("A", "B")),
    class = c("desiredgainr_sim_setup", "list")
  )
  legacy <- structure(
    list(directions = matrix(c(1, 0), nrow = 1L)),
    class = c("desiredgainr_optimisation", "list")
  )
  expect_error(
    propagate_covariance_uncertainty(
      fake_setup, legacy,
      genetic_df = 20,
      n_covariance_draws = 2, n_replicates = 2
    ),
    "does not contain trait-unit simulation directions"
  )
})
