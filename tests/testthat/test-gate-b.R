# Regression tests for the Gate B findings (M2, M3, M7, M8, M12) and the
# expert-reference additions (Q2, Q5, Q7).

gate_b_fixture <- function(seed = 23L, n = 80L) {
  set.seed(seed)
  traits <- c("yield", "height", "protein")
  G <- matrix(
    c(
      1.0, 0.4, 0.2,
      0.4, 0.8, 0.1,
      0.2, 0.1, 0.5
    ), 3L,
    dimnames = list(traits, traits)
  )
  P <- G + diag(c(1.2, 1.0, 0.9))
  dimnames(P) <- list(traits, traits)
  values <- as.data.frame(matrix(
    stats::rnorm(n * 3L),
    ncol = 3L,
    dimnames = list(sprintf("g%03d", seq_len(n)), traits)
  ))
  list(traits = traits, G = G, P = P, values = values)
}

# ---------------------------------------------------------------------------
# Q2: restricted and proportional-gain indices
# ---------------------------------------------------------------------------

test_that("the restricted index holds the restricted trait at zero", {
  fixture <- gate_b_fixture()
  fit <- restricted_index(
    fixture$values, fixture$traits,
    method = "kempthorne_nordskog",
    G = fixture$G, P = fixture$P,
    economic_weights = c(yield = 2, height = 0, protein = 1),
    restricted_traits = "height", n_select = 15L
  )
  # The defining property: expected response in the restricted trait is zero.
  expect_equal(
    unname(fit$constraint$achieved_response[["height"]]), 0,
    tolerance = 1e-8
  )
  expect_lt(fit$constraint$largest_violation, 1e-8)
  expect_equal(
    unname(fit$evaluation$expected_response[["height"]]), 0,
    tolerance = 1e-6
  )
})

test_that("restricting a trait costs merit relative to Smith-Hazel", {
  # A constraint cannot increase the correlation with net merit; if it appears
  # to, the restriction has been applied to the wrong quantity.
  fixture <- gate_b_fixture()
  weights <- c(yield = 2, height = 0, protein = 1)
  unrestricted <- selection_index(
    fixture$values, fixture$traits,
    method = "smith_hazel",
    G = fixture$G, P = fixture$P, economic_weights = weights, n_select = 15L
  )
  restricted <- restricted_index(
    fixture$values, fixture$traits,
    method = "kempthorne_nordskog",
    G = fixture$G, P = fixture$P, economic_weights = weights,
    restricted_traits = "height", n_select = 15L
  )
  expect_lte(
    restricted$evaluation$R_HI, unrestricted$evaluation$R_HI + 1e-8
  )
})

test_that("the projection form agrees with Kempthorne-Nordskog", {
  fixture <- gate_b_fixture()
  weights <- c(yield = 2, height = 0, protein = 1)
  arguments <- list(
    values = fixture$values, trait_cols = fixture$traits,
    G = fixture$G, P = fixture$P, economic_weights = weights,
    restricted_traits = "height", n_select = 15L
  )
  a <- do.call(
    restricted_index,
    c(arguments, list(method = "kempthorne_nordskog"))
  )
  b <- do.call(
    restricted_index,
    c(arguments, list(method = "restricted_smith_hazel"))
  )
  expect_equal(coef(a), coef(b), tolerance = 1e-10)
})

test_that("Tallis holds the responses in the requested proportions", {
  fixture <- gate_b_fixture()
  fit <- restricted_index(
    fixture$values, fixture$traits,
    method = "tallis",
    G = fixture$G, P = fixture$P,
    economic_weights = c(yield = 2, height = 1, protein = 1),
    target_gains = c(yield = 2, protein = 1), n_select = 15L
  )
  response <- fit$evaluation$expected_response
  # Only the ratio is constrained, not the magnitude.
  expect_equal(
    response[["yield"]] / response[["protein"]], 2,
    tolerance = 1e-6
  )
})

test_that("Mallard enforces published response proportions, not magnitudes", {
  fixture <- gate_b_fixture()
  fit <- restricted_index(
    fixture$values, fixture$traits,
    method = "mallard",
    G = fixture$G, P = fixture$P,
    economic_weights = c(yield = 2, height = 1, protein = 1),
    target_gains = c(yield = 0.4, height = 0.1), n_select = 15L
  )
  achieved <- fit$constraint$achieved_response
  expect_equal(
    unname(achieved[["yield"]] / achieved[["height"]]), 4,
    tolerance = 1e-6
  )
  expect_lt(fit$constraint$largest_violation, 1e-8)
})

test_that("Harville is the published proportionality formulation", {
  fixture <- gate_b_fixture()
  arguments <- list(
    values = fixture$values, trait_cols = fixture$traits,
    G = fixture$G, P = fixture$P,
    economic_weights = c(yield = 2, height = 1, protein = 1),
    target_gains = c(yield = 2, protein = 1), n_select = 15L
  )
  tallis <- do.call(
    restricted_index,
    c(arguments, list(method = "tallis"))
  )
  harville <- do.call(
    restricted_index, c(arguments, list(method = "harville", penalty = Inf))
  )
  expect_equal(coef(tallis), coef(harville), tolerance = 1e-10)
  expect_equal(
    harville$evaluation$expected_response[["yield"]] /
      harville$evaluation$expected_response[["protein"]],
    2,
    tolerance = 1e-6
  )

  expect_error(
    do.call(
      restricted_index, c(arguments, list(method = "harville", penalty = 1))
    ),
    "not Harville's published"
  )
})

test_that("restricted_index() validates its constraints", {
  fixture <- gate_b_fixture()
  base <- list(
    values = fixture$values, trait_cols = fixture$traits,
    G = fixture$G, P = fixture$P,
    economic_weights = c(yield = 2, height = 1, protein = 1)
  )
  expect_error(
    do.call(
      restricted_index,
      c(base, list(method = "kempthorne_nordskog"))
    ),
    "restricted_traits or constraint_matrix is required"
  )
  expect_error(
    do.call(restricted_index, c(base, list(
      method = "kempthorne_nordskog", restricted_traits = fixture$traits
    ))),
    "At least one trait must be left free"
  )
  expect_error(
    do.call(restricted_index, c(base, list(
      method = "kempthorne_nordskog", restricted_traits = "absent"
    ))),
    "Unknown constrained traits"
  )
  expect_error(
    do.call(restricted_index, c(base, list(method = "tallis"))),
    "target_gains is required"
  )
  expect_error(
    do.call(restricted_index, c(base, list(
      method = "tallis", target_gains = c(yield = 1, height = 0)
    ))),
    "zero target"
  )
})

# ---------------------------------------------------------------------------
# Q7: fitted-object ergonomics
# ---------------------------------------------------------------------------

test_that("coef() and summary() work on a fitted index", {
  fixture <- gate_b_fixture()
  fit <- selection_index(
    fixture$values, fixture$traits,
    method = "smith_hazel",
    G = fixture$G, P = fixture$P,
    economic_weights = c(yield = 2, height = 1, protein = 1), n_select = 10L
  )
  expect_identical(coef(fit), fit$coefficients)
  summarised <- summary(fit)
  expect_s3_class(summarised, "desiredgainr_index_summary")
  expect_identical(nrow(summarised$coefficients), 3L)
  expect_true("h2_index" %in% summarised$criteria$Criterion)
  expect_output(print(summarised), "desiredgainr_index_summary")
})

# ---------------------------------------------------------------------------
# Q5: upstream adapters
# ---------------------------------------------------------------------------

test_that("import_covariance() reorders traits rather than assuming order", {
  traits <- c("yield", "protein")
  supplied <- matrix(
    c(0.5, 0.2, 0.2, 1.0), 2L,
    dimnames = list(c("protein", "yield"), c("protein", "yield"))
  )
  imported <- import_covariance(supplied, traits, source = "sommer")
  # The yield variance is 1.0 in the supplied matrix and must stay with yield.
  expect_equal(imported$covariance[["yield", "yield"]], 1.0)
  expect_equal(imported$covariance[["protein", "protein"]], 0.5)
  expect_identical(rownames(imported$covariance), traits)
})

test_that("a correlation matrix is rescaled rather than used as covariance", {
  traits <- c("yield", "protein")
  correlation <- matrix(c(1, 0.3, 0.3, 1), 2L,
    dimnames = list(traits, traits)
  )
  imported <- import_covariance(
    correlation, traits,
    source = "matrix",
    is_correlation = TRUE, variances = c(yield = 4, protein = 1)
  )
  expect_equal(diag(imported$covariance), c(yield = 4, protein = 1))
  expect_equal(imported$covariance[["yield", "protein"]], 0.3 * 2)
  # Passing it unflagged must at least warn.
  expect_warning(
    import_covariance(correlation, traits, source = "matrix"),
    "correlation matrix"
  )
})

test_that("named variance components are unpacked into a matrix", {
  traits <- c("yield", "protein")
  components <- c(
    "yield:yield" = 1.0, "yield:protein" = 0.2, "protein:protein" = 0.5
  )
  imported <- import_covariance(components, traits, source = "asreml")
  expect_equal(imported$covariance[["yield", "protein"]], 0.2)
  expect_equal(imported$covariance[["protein", "yield"]], 0.2)
})

test_that("an incomplete set of components is refused", {
  traits <- c("yield", "protein")
  expect_error(
    import_covariance(c("yield:yield" = 1), traits, source = "asreml"),
    "do not cover every trait pair"
  )
  expect_error(
    import_covariance(matrix(1:4, 2L), c("a", "b", "c"), source = "matrix"),
    "does not match"
  )
})

test_that("an imported G incompatible with P is refused", {
  traits <- c("yield", "protein")
  G <- matrix(c(2, 0, 0, 0.5), 2L, dimnames = list(traits, traits))
  P <- matrix(c(1, 0, 0, 1), 2L, dimnames = list(traits, traits))
  expect_error(
    import_covariance(G, traits, source = "matrix", P = P),
    "not positive semidefinite"
  )
})

# ---------------------------------------------------------------------------
# M7: the estimated-P exception is explicit
# ---------------------------------------------------------------------------

dgsi_case <- function(n = 40L, p = 3L, seed = 31L) {
  set.seed(seed)
  traits <- paste0("t", seq_len(p))
  candidates <- data.table::data.table(GenoID = sprintf("g%03d", seq_len(n)))
  for (trait in traits) {
    data.table::set(candidates, j = trait, value = stats::rnorm(n))
  }
  G <- diag(p) * 0.5
  dimnames(G) <- list(traits, traits)
  list(
    traits = traits, G = G, candidates = candidates,
    dg = stats::setNames(rep(0.5, p), traits)
  )
}

test_that("compatibility status is recorded whether or not P was supplied", {
  case <- dgsi_case()
  fit <- run_dgsi(
    init_data = case$candidates[, "GenoID"], cand_data = case$candidates,
    trait_cols = case$traits, dg = case$dg, G = case$G,
    n_select = 8L, n_iter = 5L, n_rep = 1L, seed = 1L
  )
  compatibility <- fit$covariance_provenance$compatibility
  expect_true(is.finite(compatibility$smallest_residual_eigenvalue))
  expect_true(compatibility$status %in%
    c("admissible", "within sampling allowance"))
  expect_false(compatibility$override)
  expect_true(fit$covariance_provenance$P_was_estimated)
})

test_that("a badly inadmissible estimated P stops unless acknowledged", {
  case <- dgsi_case()
  # A genetic variance far above the phenotypic one cannot be sampling error.
  case$G <- diag(3L) * 20
  dimnames(case$G) <- list(case$traits, case$traits)
  expect_error(
    run_dgsi(
      init_data = case$candidates[, "GenoID"], cand_data = case$candidates,
      trait_cols = case$traits, dg = case$dg, G = case$G,
      n_select = 8L, n_iter = 5L, n_rep = 1L, seed = 1L
    ),
    "allow_incompatible_estimated_P"
  )
  fit <- suppressWarnings(run_dgsi(
    init_data = case$candidates[, "GenoID"], cand_data = case$candidates,
    trait_cols = case$traits, dg = case$dg, G = case$G,
    n_select = 8L, n_iter = 5L, n_rep = 1L, seed = 1L,
    allow_incompatible_estimated_P = TRUE
  ))
  expect_true(fit$covariance_provenance$compatibility$override)
  expect_identical(
    fit$covariance_provenance$compatibility$status, "inadmissible"
  )
})

# ---------------------------------------------------------------------------
# M8: replicate selection no longer uses the training objective
# ---------------------------------------------------------------------------

test_that("the winning replicate is chosen on held-out candidates", {
  case <- dgsi_case(n = 120L)
  fit <- run_dgsi(
    init_data = case$candidates[, "GenoID"], cand_data = case$candidates,
    trait_cols = case$traits, dg = case$dg, G = case$G,
    n_select = 10L, n_iter = 20L, n_rep = 4L, seed = 2L
  )
  expect_identical(fit$optimism$selection_rule, "internal pre-fit holdout")
  expect_false(is.null(fit$optimism$holdout))
  expect_gt(fit$optimism$holdout$n_holdout, 0L)
  # The chosen objective need no longer be the minimum, which is the point:
  # the minimum was a biased estimate.
  expect_true(is.finite(fit$optimism$chosen_objective))
  expect_true(is.finite(fit$optimism$chosen_training_objective))
})

test_that("the training rule remains available and is labelled", {
  case <- dgsi_case(n = 120L)
  fit <- run_dgsi(
    init_data = case$candidates[, "GenoID"], cand_data = case$candidates,
    trait_cols = case$traits, dg = case$dg, G = case$G,
    n_select = 10L, n_iter = 20L, n_rep = 4L, seed = 2L,
    replicate_selection = "training"
  )
  expect_identical(fit$optimism$selection_rule, "training")
  expect_equal(
    fit$optimism$chosen_objective, fit$optimism$minimum_objective
  )
  expect_match(fit$optimism$note, "biased downward")
})

# ---------------------------------------------------------------------------
# M3: the surrogate accepts known per-point noise
# ---------------------------------------------------------------------------

test_that("a heteroskedastic nugget is used when noise is supplied", {
  set.seed(41L)
  X <- matrix(stats::runif(40L * 3L), ncol = 3L)
  X <- X / sqrt(rowSums(X^2))
  y <- as.numeric(X %*% c(1, -1, 0.5))
  noise <- rep(c(1e-6, 1e-1), length.out = nrow(X))
  model <- DesiredGainR:::.dgr_gp_fit(X, y, noise_variance = noise)
  expect_true(model$heteroskedastic)
  quiet <- DesiredGainR:::.dgr_gp_fit(X, y)
  expect_false(quiet$heteroskedastic)
  # Both must still interpolate the training points reasonably.
  prediction <- DesiredGainR:::.dgr_gp_predict(model, X)
  expect_equal(as.numeric(prediction$mean), y, tolerance = 0.2)
})
