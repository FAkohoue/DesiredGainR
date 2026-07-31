# Monte Carlo checks of the closed-form quantities against simulation.
#
# These tests validate the published equations rather than re-asserting the
# implementation. They are the tests that catch a wrong denominator or a
# wrong correlation definition, which equation-echoing tests cannot.

# Draw from MVN(0, Sigma) without depending on MASS.
rmvn <- function(n, Sigma) {
  decomposition <- eigen(Sigma, symmetric = TRUE)
  root <- decomposition$vectors %*%
    diag(sqrt(pmax(decomposition$values, 0)), nrow(Sigma))
  matrix(stats::rnorm(n * nrow(Sigma)), n) %*% t(root)
}

test_that("QGSI index mean and variance match simulation", {
  skip_on_cran()
  set.seed(2026)
  traits <- c("t1", "t2", "t3")
  Gamma <- matrix(
    c(
      1.00, 0.30, -0.20,
      0.30, 0.80, 0.10,
      -0.20, 0.10, 1.40
    ),
    3, dimnames = list(traits, traits)
  )
  w <- c(t1 = 1.0, t2 = 0.5, t3 = -0.4)
  W <- matrix(
    c(
      0.12, 0.03, -0.02,
      0.03, -0.07, 0.01,
      -0.02, 0.01, 0.09
    ),
    3, dimnames = list(traits, traits)
  )

  n <- 400000L
  gamma_draws <- rmvn(n, Gamma)
  colnames(gamma_draws) <- traits
  index <- as.numeric(gamma_draws %*% w) +
    rowSums((gamma_draws %*% W) * gamma_draws)

  tr <- function(M) sum(diag(M))
  theoretical_mean <- tr(W %*% Gamma)
  theoretical_variance <- as.numeric(crossprod(w, Gamma %*% w)) +
    2 * tr(W %*% Gamma %*% W %*% Gamma)

  expect_equal(mean(index), theoretical_mean, tolerance = 0.02)
  expect_equal(stats::var(index), theoretical_variance, tolerance = 0.02)

  # And the package reproduces those closed forms.
  candidates <- data.table::data.table(
    GenoID = paste0("G", seq_len(200)),
    t1 = gamma_draws[seq_len(200), 1],
    t2 = gamma_draws[seq_len(200), 2],
    t3 = gamma_draws[seq_len(200), 3]
  )
  fit <- run_qgsi(
    candidates[, .(GenoID)], candidates, traits,
    linear_weights = w, W = W, Gamma = Gamma,
    center_traits = FALSE, selection_proportion = 0.1
  )
  expect_equal(
    fit$theoretical_parameters$model_expected_index,
    theoretical_mean
  )
  expect_equal(
    fit$theoretical_parameters$total_index_variance,
    theoretical_variance
  )
})

test_that("Cov(gamma, I) equals Gamma w, which the gain formula rests on", {
  skip_on_cran()
  set.seed(11)
  traits <- c("t1", "t2")
  Gamma <- matrix(
    c(1.0, 0.35, 0.35, 0.7), 2,
    dimnames = list(traits, traits)
  )
  w <- c(t1 = 1.0, t2 = 0.6)
  W <- matrix(
    c(0.15, 0.05, 0.05, -0.10), 2,
    dimnames = list(traits, traits)
  )

  n <- 400000L
  gamma_draws <- rmvn(n, Gamma)
  index <- as.numeric(gamma_draws %*% w) +
    rowSums((gamma_draws %*% W) * gamma_draws)

  # The quadratic term contributes nothing to the covariance with gamma,
  # because the third central moments of a zero-mean multivariate normal
  # vanish. This identity is what makes Gamma %*% w the numerator of the
  # per-trait gain, and it is exact rather than approximate.
  observed_covariance <- as.numeric(stats::cov(gamma_draws, index))
  expect_equal(
    observed_covariance, as.numeric(Gamma %*% w), tolerance = 0.03
  )
})

test_that("per-trait gain is a linear-regression approximation, not exact", {
  skip_on_cran()
  set.seed(11)
  traits <- c("t1", "t2")
  Gamma <- matrix(
    c(1.0, 0.35, 0.35, 0.7), 2,
    dimnames = list(traits, traits)
  )
  w <- c(t1 = 1.0, t2 = 0.6)
  W <- matrix(
    c(0.15, 0.05, 0.05, -0.10), 2,
    dimnames = list(traits, traits)
  )

  n <- 400000L
  proportion <- 0.10
  gamma_draws <- rmvn(n, Gamma)
  index <- as.numeric(gamma_draws %*% w) +
    rowSums((gamma_draws %*% W) * gamma_draws)
  keep <- index >= stats::quantile(index, 1 - proportion)
  observed_gain <- colMeans(gamma_draws[keep, , drop = FALSE])

  tr <- function(M) sum(diag(M))
  intensity <- stats::dnorm(stats::qnorm(1 - proportion)) / proportion
  linear_variance <- as.numeric(crossprod(w, Gamma %*% w))
  total_variance <- linear_variance +
    2 * tr(W %*% Gamma %*% W %*% Gamma)

  prediction_total_sd <- as.numeric(
    intensity * Gamma %*% w / sqrt(total_variance)
  )
  prediction_linear_sd <- as.numeric(
    intensity * Gamma %*% w / sqrt(linear_variance)
  )

  # Two separate approximations enter the per-trait gain, and they can be
  # examined one at a time.
  #
  # The first is that the selection differential achieved on the index equals
  # the normal-theory value i * sd(I). The index is a quadratic form and is
  # right skewed, so truncating its upper tail captures more than a normal
  # variate would and the achieved differential is larger.
  observed_differential <- mean(index[keep]) - mean(index)
  normal_differential <- intensity * stats::sd(index)
  expect_gt(observed_differential, normal_differential)

  # The second is that E[gamma | I] is linear in I, which also fails for a
  # quadratic index. Substituting the observed differential corrects only the
  # first approximation and therefore does not improve the prediction: it
  # removes an error that was partly cancelling the second one.
  regression_prediction <- as.numeric(
    Gamma %*% w * observed_differential / stats::var(index)
  )
  expect_gt(
    max(abs(observed_gain - regression_prediction)),
    max(abs(observed_gain - prediction_total_sd))
  )
  # Both remain the same order of accuracy, so neither is exact and the closer
  # agreement of the normal-theory form is a partial cancellation rather than
  # evidence that it is the better model.
  expect_lt(
    max(abs(observed_gain - regression_prediction)) /
      max(abs(observed_gain - prediction_total_sd)),
    5
  )

  # The normal-theory prediction is close but not exact, and the residual
  # error is trait-specific because E[gamma | I] is not linear in I.
  expect_equal(observed_gain, prediction_total_sd, tolerance = 0.15)

  # At this curvature the quadratic variance is only about four per cent of
  # the total, so the two candidate denominators differ by roughly two per
  # cent while the approximation error is five to nine per cent. Simulation
  # therefore cannot arbitrate between them here. The total index standard
  # deviation is used because it is the correct linear-regression denominator,
  # Cov(gamma, I) / sd(I), not because it fits better in this configuration.
  denominator_gap <- max(
    abs(prediction_total_sd - prediction_linear_sd) / abs(prediction_total_sd)
  )
  approximation_error <- max(
    abs(observed_gain - prediction_total_sd) / abs(prediction_total_sd)
  )
  expect_lt(denominator_gap, approximation_error)
})

test_that("squared index-merit correlation matches simulation", {
  skip_on_cran()
  set.seed(77)
  traits <- c("t1", "t2")
  Gamma <- matrix(
    c(1.0, 0.20, 0.20, 0.6), 2,
    dimnames = list(traits, traits)
  )
  true_G <- matrix(
    c(1.8, 0.45, 0.45, 1.1), 2,
    dimnames = list(traits, traits)
  )
  w <- c(t1 = 1.0, t2 = 0.5)
  W <- matrix(
    c(0.10, 0.02, 0.02, -0.06), 2,
    dimnames = list(traits, traits)
  )

  tr <- function(M) sum(diag(M))
  var_i <- as.numeric(crossprod(w, Gamma %*% w)) +
    2 * tr(W %*% Gamma %*% W %*% Gamma)
  var_h <- as.numeric(crossprod(w, true_G %*% w)) +
    2 * tr(W %*% true_G %*% W %*% true_G)
  cov_hi <- as.numeric(crossprod(w, Gamma %*% w)) +
    2 * tr(W %*% Gamma %*% W %*% true_G)

  # A squared correlation is bounded; the pre-0.3.1 variance ratio is not
  # constrained to be, and here the two differ materially.
  general <- cov_hi^2 / (var_i * var_h)
  ratio <- var_i / var_h
  expect_gte(general, 0)
  expect_lte(general, 1)
  expect_false(isTRUE(all.equal(general, ratio)))

  candidates <- data.table::data.table(
    GenoID = paste0("G", 1:50),
    t1 = stats::rnorm(50),
    t2 = stats::rnorm(50)
  )
  fit <- run_qgsi(
    candidates[, .(GenoID)], candidates, traits,
    linear_weights = w, W = W, Gamma = Gamma, true_G = true_G,
    center_traits = FALSE, selection_proportion = 0.2
  )
  expect_equal(
    fit$theoretical_parameters$squared_index_merit_correlation, general
  )
  expect_equal(
    fit$theoretical_parameters$variance_ratio_index_to_merit, ratio
  )
})

test_that("DGSI coefficients solve the Pesek-Baker system", {
  set.seed(5)
  traits <- c("t1", "t2", "t3")
  candidates <- data.table::data.table(
    GenoID = paste0("G", seq_len(60)),
    t1 = stats::rnorm(60),
    t2 = stats::rnorm(60),
    t3 = stats::rnorm(60)
  )
  P <- matrix(
    c(
      1.20, 0.30, 0.10,
      0.30, 0.90, -0.20,
      0.10, -0.20, 1.05
    ),
    3, dimnames = list(traits, traits)
  )
  G <- matrix(
    c(
      0.60, 0.15, 0.05,
      0.15, 0.45, -0.10,
      0.05, -0.10, 0.50
    ),
    3, dimnames = list(traits, traits)
  )

  fit <- run_dgsi(
    candidates[, .(GenoID)], candidates, traits,
    dg = c(t1 = 0.5, t2 = 0.3, t3 = 0.2),
    P = P, G = G,
    n_select = 12, n_iter = 30, n_rep = 2, seed = 5
  )

  # b = P^-1 G (G' P^-1 G)^-1 d must reproduce the optimised d exactly, up to
  # the ridge constants, via G' b = d.
  reconstructed <- as.numeric(crossprod(fit$G, fit$coefficients))
  expect_equal(
    reconstructed, as.numeric(fit$optimised_d),
    tolerance = 1e-5
  )
})

test_that("run_dgsi restores the caller's RNG state", {
  set.seed(123)
  before <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  reference <- stats::runif(3)

  set.seed(123)
  candidates <- data.table::data.table(
    GenoID = paste0("G", 1:20),
    t1 = stats::rnorm(20),
    t2 = stats::rnorm(20)
  )
  set.seed(123)
  G <- diag(2)
  dimnames(G) <- list(c("t1", "t2"), c("t1", "t2"))
  invisible(run_dgsi(
    candidates[, .(GenoID)], candidates, c("t1", "t2"),
    dg = c(t1 = 0.3, t2 = 0.2), G = G,
    n_select = 5, n_iter = 5, n_rep = 2, seed = 999
  ))

  expect_identical(
    get(".Random.seed", envir = globalenv(), inherits = FALSE), before
  )
  expect_equal(stats::runif(3), reference)
})

test_that("DGSI selection is invariant to input row order", {
  set.seed(31)
  traits <- c("t1", "t2")
  # Rounding creates genuine ties in the index score without making the
  # covariance matrix singular.
  candidates <- data.table::data.table(
    GenoID = sprintf("G%02d", 1:30),
    t1 = round(stats::rnorm(30), 1),
    t2 = round(stats::rnorm(30), 1)
  )
  G <- stats::cov(as.matrix(candidates[, ..traits]))

  fit_a <- run_dgsi(
    candidates[, .(GenoID)], candidates, traits,
    dg = c(t1 = 0.4, t2 = 0.3), G = G,
    n_select = 9, n_iter = 15, n_rep = 2, seed = 31
  )
  shuffled <- candidates[order(GenoID, decreasing = TRUE)]
  fit_b <- run_dgsi(
    shuffled[, .(GenoID)], shuffled, traits,
    dg = c(t1 = 0.4, t2 = 0.3), G = G,
    n_select = 9, n_iter = 15, n_rep = 2, seed = 31
  )

  expect_equal(
    sort(fit_a$selected_geno$GenoID),
    sort(fit_b$selected_geno$GenoID)
  )
})
