test_that("prediction covariance is labelled as a prediction estimand", {
  values <- data.frame(
    t1 = c(-1, 0, 1, 2),
    t2 = c(2, 1, 0, -1)
  )
  result <- estimate_genetic_covariance(
    values,
    trait_cols = c("t1", "t2"),
    method = "prediction_covariance"
  )

  expect_s3_class(result, "desiredgainr_covariance_estimate")
  expect_equal(result$G, stats::cov(values))
  expect_match(result$estimand, "predictions")
  expect_false(result$diagnostics$full_cross_trait_pev)
})

test_that("full PEV correction adds variances and covariances", {
  values <- data.frame(
    t1 = c(-1, 0, 1, 2),
    t2 = c(2, 1, 0, -1)
  )
  pev <- matrix(c(0.20, 0.04, 0.04, 0.10), 2)
  dimnames(pev) <- list(c("t1", "t2"), c("t1", "t2"))

  result <- estimate_genetic_covariance(
    values,
    trait_cols = c("t1", "t2"),
    method = "pev_corrected",
    prediction_error_covariance = pev
  )

  expect_equal(result$G, stats::cov(values) + pev)
  expect_equal(result$correction, pev)
  expect_match(result$estimand, "total genetic covariance")
  expect_true(result$diagnostics$full_cross_trait_pev)
})

test_that("genotype-specific PEV matrices are averaged", {
  values <- cbind(t1 = -2:2, t2 = c(-1, 0, 2, 0, -1))
  pev <- array(0, dim = c(2, 2, 5))
  for (i in seq_len(5)) {
    pev[, , i] <- diag(c(i / 100, i / 200))
  }

  result <- estimate_genetic_covariance(
    values,
    method = "pev_corrected",
    prediction_error_covariance = pev
  )

  expect_equal(
    unname(diag(result$correction)),
    c(0.03, 0.015)
  )
})

test_that("SE correction is explicit about uncorrected off-diagonals", {
  values <- data.frame(
    t1 = c(-1, 0, 1, 2),
    t2 = c(2, 1, 0, -1)
  )
  se <- data.frame(
    t1 = c(0.1, 0.2, 0.3, 0.4),
    t2 = c(0.2, 0.2, 0.4, 0.4)
  )

  result <- estimate_genetic_covariance(
    values,
    trait_cols = c("t1", "t2"),
    method = "se_diagonal_corrected",
    prediction_se = se
  )

  expect_equal(diag(result$correction), colMeans(se^2))
  expect_equal(result$correction[1, 2], 0)
  expect_match(result$estimand, "off-diagonals")
})

test_that("relationship adjustment aligns IDs and uses effective rank", {
  values <- cbind(
    t1 = c(-1, 0, 1, 2),
    t2 = c(1, -1, 2, 0)
  )
  ids <- paste0("G", 1:4)
  K <- diag(4)
  dimnames(K) <- list(rev(ids), rev(ids))

  result <- estimate_genetic_covariance(
    values,
    method = "relationship_adjusted",
    relationship_matrix = K,
    ids = ids
  )

  centred <- sweep(values, 2, colMeans(values), "-")
  expect_equal(result$G, crossprod(centred) / 4)
  expect_equal(result$diagnostics$relationship_rank, 4)
  expect_match(result$estimand, "predictions")
})

test_that("invalid covariance approximations are rejected", {
  values <- cbind(t1 = 1:4, t2 = 4:1)
  bad_pev <- matrix(c(1, 3, 3, 1), 2)

  expect_error(
    estimate_genetic_covariance(
      values,
      method = "pev_corrected",
      prediction_error_covariance = bad_pev
    ),
    "positive semidefinite"
  )
  expect_error(
    estimate_genetic_covariance(
      values,
      method = "relationship_adjusted",
      relationship_matrix = diag(4)
    ),
    "ids must contain"
  )
})

test_that("DGSI preserves covariance-estimate provenance", {
  set.seed(51)
  values <- data.table(
    GenoID = paste0("G", 1:20),
    t1 = rnorm(20),
    t2 = rnorm(20)
  )
  estimate <- estimate_genetic_covariance(
    values,
    trait_cols = c("t1", "t2"),
    method = "prediction_covariance"
  )

  result <- run_dgsi(
    init_data = values[, .(GenoID)],
    cand_data = values,
    trait_cols = c("t1", "t2"),
    dg = c(t1 = 0.4, t2 = 0.2),
    G = estimate,
    n_select = 4,
    n_iter = 5,
    n_rep = 2,
    seed = 51
  )

  expect_equal(
    result$covariance_provenance$G$method,
    "prediction_covariance"
  )
  expect_match(
    result$covariance_provenance$G$estimand,
    "predictions"
  )
  expect_equal(result$G, estimate$G)
})
