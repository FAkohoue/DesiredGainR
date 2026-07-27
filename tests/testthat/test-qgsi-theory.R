test_that("QGSI scores and contributions reproduce the matrix equation", {
  traits <- c("t1", "t2")
  values <- data.table::data.table(
    GenoID = c("A", "B", "C", "D"),
    t1 = c(-1.5, -0.5, 0.5, 1.5),
    t2 = c(1, -1, -1, 1)
  )
  w <- c(t1 = 0.8, t2 = -0.2)
  W <- matrix(
    c(0.10, 0.04, 0.04, -0.06), 2,
    dimnames = list(traits, traits)
  )
  result <- run_qgsi(
    values[, .(GenoID)], values, traits,
    linear_weights = w, W = W,
    center_traits = FALSE
  )
  X <- as.matrix(values[, ..traits])
  expected_linear <- as.numeric(X %*% w)
  expected_quadratic <- rowSums((X %*% W) * X)
  by_id <- result$ranked_geno[
    match(values$GenoID, result$ranked_geno$GenoID)
  ]

  expect_equal(by_id$LinearPart, expected_linear)
  expect_equal(by_id$QuadraticPart, expected_quadratic)
  expect_equal(by_id$QGSI, expected_linear + expected_quadratic)
  expect_equal(
    rowSums(result$linear_contributions[, -"GenoID"]),
    expected_linear
  )
  expect_equal(
    rowSums(result$quadratic_contributions[, -"GenoID"]),
    expected_quadratic
  )
})

test_that("Gamma and theoretical parameters reproduce published equations", {
  traits <- c("t1", "t2")
  values <- data.table::data.table(
    GenoID = paste0("G", 1:5),
    t1 = c(-2, -1, 0, 1, 2),
    t2 = c(-1, 1, 0, -1, 1)
  )
  w <- c(t1 = 1.2, t2 = 0.4)
  W <- matrix(
    c(0.10, -0.03, -0.03, 0.05), 2,
    dimnames = list(traits, traits)
  )
  Gamma <- matrix(
    c(1.5, 0.25, 0.25, 0.8), 2,
    dimnames = list(traits, traits)
  )
  result <- run_qgsi(
    values[, .(GenoID)], values, traits,
    linear_weights = w, W = W, Gamma = Gamma,
    true_G = Gamma, center_traits = FALSE,
    selection_proportion = 0.4
  )
  tr <- function(M) sum(diag(M))
  expected_mean <- tr(W %*% Gamma)
  expected_linear_variance <- as.numeric(crossprod(w, Gamma %*% w))
  expected_quadratic_variance <- 2 * tr(
    W %*% Gamma %*% W %*% Gamma
  )
  parameters <- result$theoretical_parameters

  expect_equal(parameters$model_expected_index, expected_mean)
  expect_equal(
    parameters$linear_index_variance, expected_linear_variance
  )
  expect_equal(
    parameters$quadratic_index_variance, expected_quadratic_variance
  )
  expect_equal(parameters$squared_index_merit_correlation, 1)
  expect_equal(parameters$mean_squared_prediction_error, 0)
  expect_equal(
    result$expected_gain_per_trait$Expected_Genetic_Gain,
    as.numeric(
      parameters$selection_intensity * Gamma %*% w /
        sqrt(expected_linear_variance)
    )
  )
})

test_that("empirical Gamma uses Equation 19.2 with divisor g", {
  traits <- c("t1", "t2")
  values <- data.table::data.table(
    GenoID = paste0("G", 1:4),
    t1 = c(-3, -1, 1, 3),
    t2 = c(2, -2, -2, 2)
  )
  W <- diag(c(t1 = 0.1, t2 = 0.1))
  dimnames(W) <- list(traits, traits)
  result <- run_qgsi(
    values[, .(GenoID)], values, traits,
    linear_weights = c(t1 = 1, t2 = 1),
    W = W,
    center_traits = TRUE
  )
  X <- scale(as.matrix(values[, ..traits]), center = TRUE, scale = FALSE)
  expect_equal(result$Gamma, crossprod(X) / nrow(X))
  expect_equal(result$covariance_provenance$divisor, nrow(X))
})

test_that("relationship-adjusted Gamma honours genotype names and rank", {
  traits <- c("t1", "t2")
  values <- data.table::data.table(
    GenoID = c("G3", "G1", "G4", "G2"),
    t1 = c(-1, 0, 1, 2),
    t2 = c(2, -1, 0, 1)
  )
  ids <- sort(values$GenoID)
  Phi <- diag(4)
  dimnames(Phi) <- list(ids, ids)
  W <- diag(c(t1 = 0.1, t2 = -0.1))
  dimnames(W) <- list(traits, traits)

  result <- run_qgsi(
    values[, .(GenoID)], values, traits,
    linear_weights = c(t1 = 1, t2 = 0.2),
    W = W,
    relationship_matrix = Phi
  )
  expect_match(result$covariance_provenance$equation, "19.1", fixed = TRUE)
  expect_equal(result$covariance_provenance$relationship_rank, 4)
  expect_equal(
    result$covariance_provenance$inverse,
    "spectral inverse"
  )
})

test_that("QGSI selection is exact and reports observed differential honestly", {
  values <- data.table::data.table(
    GenoID = paste0("G", 1:10),
    t1 = seq(-2, 2, length.out = 10),
    t2 = rep(0, 10)
  )
  W <- diag(c(t1 = 0, t2 = 0))
  dimnames(W) <- list(c("t1", "t2"), c("t1", "t2"))
  result <- run_qgsi(
    values[, .(GenoID)], values, c("t1", "t2"),
    linear_weights = c(t1 = 1, t2 = 0),
    W = W,
    n_select = 3
  )
  expect_equal(sum(result$ranked_geno$Selected), 3)
  expect_equal(result$selection$selection_proportion, 0.3)
  expect_true(
    "Observed_GEBV_differential" %in%
      names(result$observed_selection_differential)
  )
  expect_false(any(grepl(
    "realised",
    names(result$observed_selection_differential),
    ignore.case = TRUE
  )))
})
