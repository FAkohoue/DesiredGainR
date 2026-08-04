test_that("general information index reduces to classical equations", {
  traits <- c("yield", "quality")
  G <- matrix(c(1, 0.2, 0.2, 0.6), 2,
    dimnames = list(traits, traits)
  )
  P <- G + diag(c(0.8, 0.7))
  dimnames(P) <- list(traits, traits)
  set.seed(4)
  values <- as.data.frame(matrix(
    stats::rnorm(80), ncol = 2,
    dimnames = list(paste0("g", 1:40), traits)
  ))
  model <- selection_information(values, P = P, C = G, G = G)
  economic <- c(yield = 2, quality = 1)
  desired <- c(yield = 1, quality = 0.5)

  fit_economic <- generalized_index(model, economic, "economic", n_select = 8)
  expect_equal(
    unname(fit_economic$coefficients),
    as.numeric(solve(P, G %*% economic)),
    tolerance = 1e-10
  )

  fit_desired <- generalized_index(
    model, desired, "desired_gain", n_select = 8
  )
  expected_direction <- crossprod(G, fit_desired$coefficients)
  expect_equal(
    as.numeric(expected_direction / expected_direction[1]),
    unname(desired / desired[1]),
    tolerance = 1e-8
  )

  oriented <- generalized_index(
    model, desired, "desired_gain", n_select = 8,
    lower_is_better = "quality"
  )
  signed <- generalized_index(
    model, c(yield = 1, quality = -0.5), "desired_gain", n_select = 8
  )
  expect_equal(oriented$coefficients, signed$coefficients, tolerance = 1e-10)
  expect_lt(oriented$expected_response[["quality"]], 0)
})

test_that("rectangular information model predicts a smaller objective", {
  information <- c("record", "family", "marker")
  objective <- c("yield", "disease")
  P <- diag(c(2, 1.5, 1.2))
  dimnames(P) <- list(information, information)
  C <- matrix(c(
    0.8, -0.2,
    0.5, -0.1,
    0.6, -0.3
  ), 3, 2, byrow = TRUE, dimnames = list(information, objective))
  G <- diag(c(1, 0.7))
  dimnames(G) <- list(objective, objective)
  joint <- rbind(cbind(P, C), cbind(t(C), G))
  expect_true(min(eigen(joint, symmetric = TRUE)$values) > 0)
  set.seed(5)
  values <- as.data.frame(matrix(
    stats::rnorm(90), ncol = 3,
    dimnames = list(paste0("g", 1:30), information)
  ))
  model <- selection_information(values, P, C, G)
  fit <- generalized_index(
    model, c(yield = 1, disease = -0.4), "desired_gain", n_select = 6
  )
  expect_named(fit$coefficients, information)
  expect_named(fit$expected_response, objective)
  expect_equal(nrow(fit$selected), 6)
})

test_that("invalid cross-covariance is rejected", {
  variables <- c("x1", "x2")
  P <- G <- diag(2)
  dimnames(P) <- dimnames(G) <- list(variables, variables)
  C <- matrix(2, 2, 2, dimnames = list(variables, variables))
  values <- as.data.frame(matrix(1:20, ncol = 2))
  names(values) <- variables
  expect_error(
    selection_information(values, P, C, G),
    "valid joint covariance"
  )
})

test_that("candidate prediction errors propagate to selection probabilities", {
  traits <- c("t1", "t2")
  G <- diag(c(1, 0.5))
  P <- diag(c(2, 1))
  dimnames(G) <- dimnames(P) <- list(traits, traits)
  values <- data.frame(
    t1 = c(2, 1, 0),
    t2 = c(1, 1, 1),
    row.names = c("a", "b", "c")
  )
  fit <- selection_index(
    values, traits, method = "smith_hazel",
    G = G, P = P, economic_weights = c(t1 = 1, t2 = 1),
    scale_traits = FALSE, n_select = 1
  )
  low <- diag(c(0.001, 0.001))
  dimnames(low) <- list(traits, traits)
  uncertainty <- candidate_score_uncertainty(
    fit, low, n_draws = 500, seed = 9
  )
  expect_gt(
    uncertainty$summary[id == "a", selection_probability],
    0.95
  )
  expect_true(all(uncertainty$summary$score_se > 0))
})

test_that("Satoh projection gives proportional restricted breeding values", {
  traits <- c("t1", "t2", "t3")
  G <- matrix(c(
    1, 0.2, 0.1,
    0.2, 0.8, 0.15,
    0.1, 0.15, 0.6
  ), 3, dimnames = list(traits, traits))
  values <- data.frame(
    t1 = c(1, -0.2),
    t2 = c(0.4, 0.7),
    t3 = c(-0.1, 0.5),
    row.names = c("a", "b")
  )
  direction <- c(t1 = 3, t2 = 4, t3 = 5)
  fit <- restricted_breeding_values(values, G, direction = direction)
  expected <- fit$beta * matrix(direction, 2, 3, byrow = TRUE)
  expect_equal(unname(fit$restricted), unname(expected), tolerance = 1e-10)
  expect_equal(fit$largest_violation, 0, tolerance = 1e-10)

  response <- evaluate_restricted_response(
    c(t1 = 0.3, t2 = 0.4, t3 = 0.5), direction, G
  )
  expect_equal(response$beta, 0.1, tolerance = 1e-10)
  expect_equal(response$mahalanobis_residual, 0, tolerance = 1e-10)

  decreasing <- restricted_breeding_values(
    values, G, direction = direction, lower_is_better = "t3"
  )
  expect_true(all(decreasing$restriction$direction[["t3"]] < 0))
})

test_that("restricted index reports Satoh criterion for a full direction", {
  traits <- c("t1", "t2", "t3")
  G <- matrix(c(
    1, 0.2, 0.1,
    0.2, 0.8, 0.15,
    0.1, 0.15, 0.6
  ), 3, dimnames = list(traits, traits))
  P <- G + diag(3)
  dimnames(P) <- list(traits, traits)
  set.seed(2)
  values <- as.data.frame(matrix(
    stats::rnorm(150), ncol = 3,
    dimnames = list(paste0("g", 1:50), traits)
  ))
  fit <- restricted_index(
    values, traits, method = "harville", G = G, P = P,
    target_gains = c(t1 = 3, t2 = 4, t3 = 5),
    scale_traits = FALSE, n_select = 10
  )
  expect_s3_class(fit$restricted_breeding_values, "desiredgainr_restricted_bv")
  expect_s3_class(
    fit$constraint$satoh_response,
    "desiredgainr_restricted_response"
  )
  expect_lt(fit$constraint$largest_violation, 1e-8)
})

test_that("Cunningham efficiency identifies redundant information", {
  traits <- c("strong", "weak")
  G <- diag(c(1, 0.01))
  P <- diag(c(1.2, 10))
  dimnames(G) <- dimnames(P) <- list(traits, traits)
  values <- data.frame(
    strong = seq(-2, 2, length.out = 30),
    weak = stats::rnorm(30)
  )
  fit <- selection_index(
    values, traits, method = "smith_hazel", G = G, P = P,
    economic_weights = c(strong = 1, weak = 1),
    scale_traits = FALSE, n_select = 5
  )
  efficiency <- index_information_efficiency(fit)
  expect_gt(
    efficiency[Information == "weak", Efficiency_after_deletion],
    0.99
  )
  expect_error(
    index_information_efficiency(selection_index(
      values, traits, method = "base",
      economic_weights = c(strong = 1, weak = 1),
      scale_traits = FALSE, n_select = 5
    )),
    "Smith-Hazel or general economic"
  )
  information <- selection_information(values, P = P, C = G, G = G)
  general <- generalized_index(
    information, c(strong = 1, weak = 1), "economic", n_select = 5
  )
  expect_equal(nrow(index_information_efficiency(general)), 2)
})

test_that("Elston index keeps every floor firm", {
  values <- data.frame(
    yield = c(2, 3, 1),
    quality = c(2, 0.5, 3),
    row.names = c("balanced", "low_quality", "low_yield")
  )
  expect_warning(
    fit <- selection_index(
      values, c("yield", "quality"), method = "elston",
      culling_thresholds = c(yield = 1.5, quality = 1.5),
      center_traits = FALSE, scale_traits = FALSE, n_select = 2
    ),
    "keeps the floors firm"
  )
  expect_equal(fit$n_selected, 1)
  expect_equal(fit$selected$id, "balanced")
  expect_true(fit$culling_report[id == "balanced", passed])
})

test_that("culling limits use original units and respect trait direction", {
  values <- data.frame(
    yield = c(8, 9, 10),
    disease = c(5, 3, 2),
    row.names = c("fails_both", "at_limits", "passes")
  )
  fit <- selection_index(
    values, c("yield", "disease"), method = "independent_culling",
    culling_thresholds = c(yield = 9, disease = 3),
    lower_is_better = "disease",
    center_traits = TRUE, scale_traits = TRUE, n_select = 2
  )
  expect_setequal(fit$selected$id, c("at_limits", "passes"))
  expect_false(fit$culling_report[id == "fails_both", passed])
})

test_that("method comparison checks fairness and reports attainment", {
  traits <- c("yield", "disease")
  values <- data.frame(
    yield = c(-1.2, -0.4, 0.1, 0.5, 0.9, 1.4),
    disease = c(1.1, 0.3, -0.1, -0.6, -0.2, -1.2),
    row.names = paste0("L", 1:6)
  )
  G <- matrix(c(1, -0.2, -0.2, 0.8), 2,
    dimnames = list(traits, traits)
  )
  P <- matrix(c(1.8, -0.2, -0.2, 1.4), 2,
    dimnames = list(traits, traits)
  )
  objective <- c(yield = 1, disease = 0.5)
  fits <- list(
    Smith_Hazel = selection_index(
      values, traits, method = "smith_hazel", G = G, P = P,
      economic_weights = objective, lower_is_better = "disease",
      scale_traits = FALSE, n_select = 2
    ),
    Base = selection_index(
      values, traits, method = "base", G = G, P = P,
      economic_weights = objective, lower_is_better = "disease",
      scale_traits = FALSE, n_select = 2
    )
  )
  comparison <- compare_selection_methods(
    fits, target_gains = c(yield = 0.5, disease = 0.2)
  )
  expect_s3_class(comparison, "desiredgainr_method_comparison")
  expect_equal(nrow(comparison$summary), 2)
  expect_equal(nrow(comparison$responses), 4)
  expect_true(all(comparison$fairness$Satisfied))
  expect_equal(unname(diag(comparison$selected_jaccard)), c(1, 1))
  expect_true(all(is.finite(comparison$summary$Mahalanobis_alignment)))
})

test_that("scenario object keeps architecture assumptions together", {
  scenario <- breeding_scenario(
    "oligogenic",
    programme = "outcross",
    architecture = list(
      qtl_per_chromosome = 5,
      effect_distribution = "gamma",
      qtl_shape = 0.4
    ),
    evaluation = list(type = "gebv")
  )
  expect_s3_class(scenario, "desiredgainr_scenario")
  expect_equal(scenario$architecture$qtl_per_chromosome, 5)
  expect_equal(scenario$evaluation$type, "gebv")
  expect_error(
    breeding_scenario("clone", programme = "clonal"),
    "dominance_degree"
  )
})

test_that("scenario calibration and stress testing form one auditable path", {
  skip_if_not_installed("AlphaSimR")
  skip_on_cran()

  set.seed(17)
  n_individuals <- 24L
  n_variants <- 40L
  variant_names <- sprintf("v%03d", seq_len(n_variants))
  individual_names <- sprintf("f%03d", seq_len(n_individuals))
  draw_haplotype <- function() {
    matrix(
      stats::rbinom(n_variants * n_individuals, 1L, 0.5),
      nrow = n_variants,
      dimnames = list(variant_names, individual_names)
    )
  }
  map <- data.frame(
    variant_id = variant_names,
    chromosome = rep(1:2, each = 20L),
    position_bp = rep(seq_len(20L) * 1e6, 2L)
  )
  founders <- founder_haplotypes(draw_haplotype(), draw_haplotype(), map)
  traits <- c("yield", "disease")
  G <- matrix(c(1, -0.2, -0.2, 0.7), 2,
    dimnames = list(traits, traits)
  )
  scenario <- breeding_scenario(
    "gamma effects",
    programme = "outcross",
    architecture = list(
      qtl_per_chromosome = 5L,
      effect_distribution = "gamma",
      qtl_shape = 0.8
    )
  )
  setup <- suppressWarnings(founder_population(
    founders, G = G, h2 = c(yield = 0.4, disease = 0.4),
    seed = 23L, scenario = scenario
  ))
  expect_equal(setup$qtl_effect_distribution, "gamma")
  calibration <- simulation_calibration(setup)
  expect_s3_class(calibration, "desiredgainr_calibration")
  expect_true(all(c("Check", "Status") %in% names(calibration$checks)))

  stress <- stress_test_desired_gains(
    setups = list(first = setup, repeated = setup),
    desired_gains = c(yield = 1, disease = 0.5),
    options = list(
      simulation = list(
        n_cycles = 1L,
        mating_system = "outcross",
        n_parents = 4L,
        n_crosses = 4L,
        n_progeny_per_cross = 3L,
        lower_is_better = "disease"
      ),
      minimum_gains = c(yield = 0.1, disease = 0.1),
      min_replicates = 2L,
      max_replicates = 2L,
      batch_size = 2L,
      seed = 31L
    )
  )
  expect_s3_class(stress, "desiredgainr_stress_test")
  expect_true(stress$matched_replicate_seeds)
  expect_equal(stress$summary$Mean_utility[1],
    stress$summary$Mean_utility[2], tolerance = 1e-12
  )
  expect_true(all(stress$summary$Regret == 0))
})
