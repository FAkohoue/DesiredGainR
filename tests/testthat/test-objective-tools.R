# Property tests for the objective-setting layer. Each formula is checked
# against an independent construction rather than restated.

test_traits <- c("t1", "t2", "t3")

test_G <- matrix(
  c(
    0.60, 0.12, -0.08,
    0.12, 0.45, 0.05,
    -0.08, 0.05, 0.70
  ),
  3, dimnames = list(test_traits, test_traits)
)

test_P <- matrix(
  c(
    1.20, 0.25, -0.10,
    0.25, 0.90, 0.08,
    -0.10, 0.08, 1.35
  ),
  3, dimnames = list(test_traits, test_traits)
)

test_that("orientation is declared rather than signed by hand", {
  d <- c(t1 = 0.5, t2 = 0.3, t3 = 0.2)

  # Declaring t2 as lower-is-better must give the same answer as orienting the
  # covariance matrices by hand and signing the desired gain.
  D <- diag(c(1, -1, 1))
  dimnames(D) <- list(test_traits, test_traits)
  by_hand <- implied_economic_weights(
    c(t1 = 0.5, t2 = -0.3, t3 = 0.2),
    G = test_G, P = test_P
  )
  declared <- implied_economic_weights(
    d, G = test_G, P = test_P, lower_is_better = "t2"
  )
  # The hand-oriented result is in the raw direction; flipping t2 back aligns
  # the two.
  expect_equal(
    as.numeric(declared), as.numeric(by_hand * c(1, -1, 1)),
    tolerance = 1e-8
  )
})

test_that("the translation round trip preserves the units it was given", {
  d <- c(t1 = 0.5, t2 = 0.3, t3 = 0.2)
  for (units in c("trait", "genetic_sd", "phenotypic_sd")) {
    w <- implied_economic_weights(
      d, test_G, test_P, lower_is_better = "t2", gain_units = units
    )
    back <- implied_desired_gains(
      w, test_G, test_P, lower_is_better = "t2", gain_units = units
    )
    expect_equal(as.numeric(back), as.numeric(d), tolerance = 1e-8,
                 info = units)
  }
})

test_that("feasibility respects the declared trait directions", {
  # Asking to improve a trait, having declared that improvement means a
  # decrease, must not be the same request as asking it to increase.
  improving <- gain_feasibility(
    c(t1 = 0.4, t2 = 0.4, t3 = 0.3), test_G, test_P,
    n_candidates = 500, n_select = 50, lower_is_better = "t2"
  )
  raw <- gain_feasibility(
    c(t1 = 0.4, t2 = 0.4, t3 = 0.3), test_G, test_P,
    n_candidates = 500, n_select = 50
  )
  expect_false(isTRUE(all.equal(
    improving$required_intensity, raw$required_intensity
  )))
})

test_that("a required proportion below one candidate is not rounded up", {
  # A target can require a proportion that is attainable in principle but
  # smaller than a single candidate in the population to hand. Reporting a
  # required count of one would overstate what that population can deliver.
  d <- c(t1 = 1, t2 = 1, t3 = 1)
  probe <- gain_feasibility(
    d, test_G, test_P, n_candidates = 1000, n_select = 100
  )
  skip_if(
    !is.finite(probe$required_proportion),
    "The required intensity exceeds what normal truncation can deliver."
  )
  # Choose a population small enough that the required proportion falls below
  # one candidate.
  n_small <- max(2L, floor(0.5 / probe$required_proportion))
  tight <- gain_feasibility(
    d, test_G, test_P,
    n_candidates = n_small, n_select = max(1L, floor(0.1 * n_small))
  )
  expect_lt(tight$required_proportion * n_small, 1)
  expect_true(is.na(tight$required_n_select))
  expect_output(print(tight), "fewer than one")
})

test_that("an unreachable intensity is reported as such, not as a proportion", {
  # Beyond the reach of normal truncation there is no proportion at all, and
  # the report must say so rather than invent one.
  extreme <- gain_feasibility(
    c(t1 = 30, t2 = 30, t3 = 30), test_G, test_P,
    n_candidates = 100, n_select = 10
  )
  expect_true(is.na(extreme$required_proportion))
  expect_true(is.na(extreme$required_n_select))
  expect_false(extreme$feasible_in_this_population)
  expect_output(print(extreme), "unattainable under normal truncation")
})

test_that("the weight and desired-gain translations are mutual inverses", {
  d <- c(t1 = 0.5, t2 = 0.3, t3 = -0.2)
  w <- implied_economic_weights(d, test_G, test_P)
  back <- implied_desired_gains(w, test_G, test_P)
  expect_equal(as.numeric(back), as.numeric(d), tolerance = 1e-8)

  w0 <- c(t1 = 1.0, t2 = 0.4, t3 = -0.6)
  d0 <- implied_desired_gains(w0, test_G, test_P)
  expect_equal(
    as.numeric(implied_economic_weights(d0, test_G, test_P)),
    as.numeric(w0),
    tolerance = 1e-8
  )
})

test_that("implied weights reproduce the Smith-Hazel coefficients", {
  d <- c(t1 = 0.5, t2 = 0.3, t3 = 0.2)
  w <- implied_economic_weights(d, test_G, test_P)

  P_inv <- solve(test_P)
  smith_hazel <- as.numeric(P_inv %*% test_G %*% w)
  middle <- test_G %*% P_inv %*% test_G
  yamada <- as.numeric(P_inv %*% test_G %*% solve(middle, d))

  expect_equal(smith_hazel, yamada, tolerance = 1e-8)
})

test_that("required intensity matches the achievable-response ellipsoid", {
  d <- c(t1 = 0.4, t2 = 0.25, t3 = 0.30)
  feasibility <- gain_feasibility(
    d, test_G, test_P, n_candidates = 1000, n_select = 100
  )

  G_inv <- solve(test_G)
  expected <- sqrt(as.numeric(crossprod(d, G_inv %*% test_P %*% G_inv %*% d)))
  expect_equal(feasibility$required_intensity, expected, tolerance = 1e-8)

  # Any Yamada index lands exactly on the ellipsoid R' G^-1 P G^-1 R = i^2.
  P_inv <- solve(test_P)
  middle <- test_G %*% P_inv %*% test_G
  b <- as.numeric(P_inv %*% test_G %*% solve(middle, d))
  intensity <- feasibility$planned_intensity
  response <- intensity * as.numeric(test_G %*% b) /
    sqrt(as.numeric(crossprod(b, test_P %*% b)))
  expect_equal(
    as.numeric(crossprod(response, G_inv %*% test_P %*% G_inv %*% response)),
    intensity^2,
    tolerance = 1e-6
  )
  expect_equal(response, as.numeric(feasibility$attainable_response),
               tolerance = 1e-6)
})

test_that("feasibility is invariant to rescaling the desired-gain direction", {
  d <- c(t1 = 0.4, t2 = 0.25, t3 = 0.30)
  one <- gain_feasibility(d, test_G, test_P, n_candidates = 500, n_select = 50)
  ten <- gain_feasibility(
    10 * d, test_G, test_P, n_candidates = 500, n_select = 50
  )
  # The requested magnitude changes the required intensity tenfold, but the
  # attainable direction is unchanged.
  expect_equal(ten$required_intensity, 10 * one$required_intensity,
               tolerance = 1e-8)
  expect_equal(
    as.numeric(ten$attainable_response),
    as.numeric(one$attainable_response),
    tolerance = 1e-8
  )
})

test_that("an unreachable target is reported as infeasible", {
  extreme <- c(t1 = 40, t2 = 40, t3 = 40)
  feasibility <- gain_feasibility(
    extreme, test_G, test_P, n_candidates = 300, n_select = 30
  )
  expect_false(feasibility$feasible_at_planned_intensity)
  expect_false(feasibility$feasible_in_this_population)
  expect_lt(feasibility$attainable_fraction, 1)
})

test_that("required proportion inverts the selection-intensity relation", {
  d <- c(t1 = 0.30, t2 = 0.20, t3 = 0.25)
  feasibility <- gain_feasibility(
    d, test_G, test_P, n_candidates = 2000, n_select = 200
  )
  p <- feasibility$required_proportion
  skip_if(!is.finite(p))
  recovered <- stats::dnorm(stats::qnorm(1 - p)) / p
  expect_equal(recovered, feasibility$required_intensity, tolerance = 1e-6)
})

test_that("retrospective weights recover a known linear preference", {
  set.seed(404)
  n <- 400L
  population <- matrix(
    stats::rnorm(n * 3), ncol = 3, dimnames = list(NULL, test_traits)
  )
  truth <- c(t1 = 1.0, t2 = 0.5, t3 = -0.3)
  score <- as.numeric(population %*% truth)
  selected <- population[order(-score)[seq_len(40L)], , drop = FALSE]

  recovered <- retrospective_weights(
    selected, population, test_traits
  )
  # b = P^-1 s recovers the generating weights up to a positive scalar,
  # because the differential is proportional to P times the true weights.
  scaled_truth <- truth / sqrt(sum(truth^2))
  scaled_recovered <- recovered$coefficients /
    sqrt(sum(recovered$coefficients^2))
  expect_gt(stats::cor(scaled_truth, scaled_recovered), 0.98)
  expect_equal(recovered$n_selected, 40L)
  expect_equal(recovered$selected_proportion, 0.1)
})

test_that("effective weights flag a dominating trait when asked", {
  b <- c(t1 = 1, t2 = 1, t3 = 1)
  lopsided <- diag(c(0.001, 0.001, 100))
  dimnames(lopsided) <- list(test_traits, test_traits)
  expect_warning(
    table <- effective_weights(b, lopsided, warn = TRUE),
    "effective weight"
  )
  expect_equal(nrow(table), 3L)
  expect_equal(sum(table$Genetic_share), 1, tolerance = 1e-8)
  expect_equal(which.max(table$Genetic_share), 3L)
})

test_that("weight sensitivity is high under small perturbation", {
  set.seed(77)
  values <- matrix(
    stats::rnorm(120 * 3), ncol = 3,
    dimnames = list(paste0("c", 1:120), test_traits)
  )
  w <- c(t1 = 1, t2 = 0.5, t3 = 0.25)

  tight <- weight_sensitivity(
    w, values, test_G, test_P, n_select = 12,
    relative_sd = 0.02, n_draws = 60L, seed = 1
  )
  loose <- weight_sensitivity(
    w, values, test_G, test_P, n_select = 12,
    relative_sd = 1.20, n_draws = 60L, seed = 1
  )
  expect_gt(tight$stability_proportion, loose$stability_proportion)
  expect_gt(stats::median(tight$agreement), stats::median(loose$agreement))
  expect_true(all(tight$agreement >= 0 & tight$agreement <= 1))
})

test_that("weight sensitivity restores the caller's RNG state", {
  set.seed(11)
  before <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  values <- matrix(
    stats::rnorm(60 * 3), ncol = 3,
    dimnames = list(paste0("c", 1:60), test_traits)
  )
  set.seed(11)
  invisible(weight_sensitivity(
    c(t1 = 1, t2 = 1, t3 = 1), values, test_G, test_P,
    n_select = 6, n_draws = 20L
  ))
  expect_identical(
    get(".Random.seed", envir = globalenv(), inherits = FALSE), before
  )
})

test_that("ill-conditioned matrices are reported and warned about", {
  well_conditioned <- diag(2)
  dimnames(well_conditioned) <- list(c("a", "b"), c("a", "b"))

  # Eigenvalues 2 - 1e-10 and 1e-10, so the condition number is about 2e10 and
  # the inversion has lost most of the available precision.
  nearly_singular <- matrix(
    c(1, 1 - 1e-10, 1 - 1e-10, 1), 2,
    dimnames = list(c("a", "b"), c("a", "b"))
  )
  diagnostics <- matrix_diagnostics(nearly_singular, "test")
  expect_lt(diagnostics$reciprocal_condition, 1e-8)
  expect_gt(diagnostics$condition_number, 1e8)

  expect_warning(
    implied_economic_weights(
      c(a = 1, b = 1), nearly_singular, well_conditioned
    ),
    "ill-conditioned"
  )

  # A merely correlated matrix must stay quiet, or the warning becomes noise.
  correlated <- matrix(
    c(1, 0.9, 0.9, 1), 2,
    dimnames = list(c("a", "b"), c("a", "b"))
  )
  expect_gt(matrix_diagnostics(correlated)$reciprocal_condition, 1e-8)
  expect_silent(
    implied_economic_weights(c(a = 1, b = 1), correlated, well_conditioned)
  )
})

test_that("effective weights survive a degenerate covariance matrix", {
  # A rank-deficient G can drive the effective weights to zero or to a
  # non-finite value. A diagnostic must degrade gracefully rather than break
  # the fit that called it.
  singular <- matrix(0, 2L, 2L, dimnames = list(c("a", "b"), c("a", "b")))
  table <- effective_weights(c(a = 1, b = 1), singular)
  expect_equal(nrow(table), 2L)
  expect_true(all(is.na(table$Genetic_share)))
  expect_silent(effective_weights(c(a = 0, b = 0), diag(2) * 1))
})

test_that("the dominance threshold adapts to the number of traits", {
  # With two traits a 60% share is an unremarkable emphasis and must stay
  # quiet; the default threshold is halfway between an equal share and total
  # dominance.
  G2 <- diag(2)
  dimnames(G2) <- list(c("a", "b"), c("a", "b"))
  expect_silent(
    effective_weights(c(a = 1.5, b = 1), G2, warn = TRUE)
  )
  expect_warning(
    effective_weights(c(a = 20, b = 1), G2, warn = TRUE),
    "effective weight"
  )
  # Warnings are off by default, because concentration can reflect a
  # deliberately asymmetric objective rather than mismatched scales.
  expect_silent(effective_weights(c(a = 20, b = 1), G2))
})
