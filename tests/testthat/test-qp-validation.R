test_that("generic active-set QP solver matches analytical clipping", {
  H <- diag(c(1, 2, 4, 8))
  f <- c(-3, 4, -2, 16)
  lower <- c(-1, -1, 0, -Inf)
  upper <- c(2, 3, 0.4, 1)
  expected <- pmin(upper, pmax(lower, -f / diag(H)))
  result <- DesiredGainR:::.dgr_box_qp(H, f, lower, upper)
  expect_equal(result$solution, expected, tolerance = 1e-12)
  expect_lt(result$kkt$overall, 1e-12)
  expect_identical(result$status, "solved")
})

test_that("generic QP solver handles fixed and infinite bounds", {
  H <- matrix(c(3, -0.5, 0, -0.5, 2, 0.2, 0, 0.2, 1), 3)
  f <- c(-2, 0.5, -1)
  lower <- c(0.25, -Inf, -1)
  upper <- c(0.25, 0.1, Inf)
  result <- DesiredGainR:::.dgr_box_qp(H, f, lower, upper)
  expect_equal(result$solution[1], 0.25, tolerance = 1e-12)
  expect_lte(result$solution[2], 0.1 + 1e-10)
  expect_lt(result$kkt$overall, 1e-8)
})

test_that("population lower-bound solver is permutation invariant", {
  B <- matrix(c(2, -0.4, 0.2, -0.4, 1.5, -0.1, 0.2, -0.1, 1), 3)
  lower <- c(1, 0.3, 0)
  original <- DesiredGainR:::.dgr_quadratic_lower_bound(B, lower)$solution
  permutation <- c(3, 1, 2)
  transformed <- DesiredGainR:::.dgr_quadratic_lower_bound(
    B[permutation, permutation], lower[permutation]
  )$solution
  expect_equal(transformed[order(permutation)], original, tolerance = 1e-10)
})

test_that("committed independent solver report is complete and passing", {
  report_path <- system.file(
    "validation", "qp-validation-results.csv",
    package = "DesiredGainR"
  )
  if (!nzchar(report_path)) {
    report_path <- testthat::test_path(
      "..", "..", "inst", "validation", "qp-validation-results.csv"
    )
  }
  report <- utils::read.csv(report_path, check.names = FALSE)
  expect_gte(nrow(report), 50L)
  expect_true(all(report$pass))
  expect_setequal(unique(report$class), c(
    "analytical", "well_conditioned", "conditioned",
    "constructed_active", "degenerate_active_set", "invariance"
  ))
  expect_true(all(c(
    "clarabel_version", "osqp_version", "clarabel_dual_difference",
    "clarabel_reordered_solution_difference", "internal_kkt_residual",
    "osqp_unpolished_kkt_residual", "osqp_polished_kkt_residual"
  ) %in% names(report)))
})
