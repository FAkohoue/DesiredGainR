test_that("QGSI rejects invalid weight and covariance specifications", {
  values <- data.table::data.table(
    GenoID = paste0("G", 1:6),
    t1 = -2:3,
    t2 = c(1, -1, 2, -2, 0, 1)
  )
  asymmetric <- matrix(c(1, 0.2, 0.1, 1), 2)
  dimnames(asymmetric) <- list(c("t1", "t2"), c("t1", "t2"))
  expect_error(
    run_qgsi(
      values[, .(GenoID)], values, c("t1", "t2"),
      linear_weights = c(t1 = 1, t2 = 1),
      W = asymmetric
    ),
    "symmetric"
  )

  W <- diag(2)
  dimnames(W) <- list(c("t1", "t2"), c("t1", "t2"))
  expect_error(
    run_qgsi(
      values[, .(GenoID)], values, c("t1", "t2"),
      linear_weights = c(t1 = 1, t2 = 1),
      W = W,
      Gamma = diag(c(1, -1))
    ),
    "positive semidefinite"
  )
  expect_error(
    run_qgsi(
      values[, .(GenoID)], values, c("t1", "t2"),
      linear_weights = c(t1 = 1, t2 = 1),
      W = W,
      n_select = 2,
      selection_proportion = 0.5
    ),
    "at most one"
  )
})

test_that("legacy desired-gain QGSI wrapper is explicit about semantics", {
  values <- data.table::data.table(
    GenoID = paste0("G", 1:6),
    t1 = -2:3,
    t2 = c(1, -1, 2, -2, 0, 1)
  )
  W <- diag(c(t1 = 0.1, t2 = -0.1))
  dimnames(W) <- list(c("t1", "t2"), c("t1", "t2"))
  # The wrapper raises two distinct warnings, and both are intended. The first
  # states that desired gains are not economic weights; the second, from
  # run_qgsi(), states that the wrapper's center_traits = FALSE default leaves
  # the reported theoretical parameters describing a centred index while the
  # scores are uncentred. Asserting both keeps either from being lost.
  expect_warning(
    expect_warning(
      result <- run_qgsi_desired_gain(
        values[, .(GenoID)],
        values,
        c("t1", "t2"),
        dg = c(t1 = 1, t2 = 0.5),
        W_d = W,
        debug = FALSE
      ),
      "desired gains are not"
    ),
    "centred covariance"
  )
  expect_true("QGSI_DG" %in% names(result$ranked_geno))
})

test_that("desired-gain eligibility mode requires thresholds", {
  values <- data.table::data.table(
    GenoID = paste0("G", 1:8),
    t1 = rnorm(8),
    t2 = rnorm(8)
  )
  G <- stats::cov(as.matrix(values[, .(t1, t2)]))
  expect_error(
    run_dgsi(
      init_data = values[, .(GenoID)],
      cand_data = values,
      trait_cols = c("t1", "t2"),
      dg = c(t1 = 0.3, t2 = 0.2),
      G = G,
      select_mode = "eligible_top_n",
      debug = FALSE
    ),
    "trait_min is required"
  )
})
