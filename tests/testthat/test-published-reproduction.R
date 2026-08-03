# Item Q4: independent reproduction against published results.
#
# Every other test in this suite checks the package against itself or against
# an equation this package also implements. These check it against numbers
# produced by someone else's software on someone else's data: Rahimi and
# Debnath (2023), Scientific Reports 13:18977, computed in SAS PROC IML and
# published in tables.
#
# Tolerances are stated in each test and reflect the four-significant-figure
# rounding of the published tables, not a tuned fit.

test_that("the published desired gains are the genetic standard deviations", {
  reference <- rahimi_debnath_2023()
  # The article states d = sqrt(vecdiag(G)). This is definitional, so it is a
  # transcription check rather than a scientific one: it catches a typo in the
  # frozen values before the substantive tests run on them.
  expect_equal(
    sqrt(reference$genetic_variances), reference$desired_gains,
    tolerance = 1e-12
  )
  expect_length(reference$traits, 7L)
  expect_named(reference$desired_gains, reference$traits)
})

test_that("published Pesek-Baker gains are exactly proportional to d", {
  # THE reproduction. Pesek-Baker sets b = G^-1 d, so the expected correlated
  # response is
  #
  #   Delta = i * G b / sqrt(b' P b) = i * d / sqrt(b' P b),
  #
  # exactly proportional to d with no approximation. If the published table
  # satisfies this, the authors' independent implementation and this package's
  # understanding of the method agree on the defining property of the family.
  #
  # This is a stronger check than comparing coefficients would be, because it
  # cannot be satisfied by accident: seven ratios must coincide.
  reference <- rahimi_debnath_2023()
  ratio <- reference$pesek_baker_gain / reference$desired_gains

  expect_equal(
    unname(stats::sd(ratio) / mean(ratio)), 0,
    tolerance = 2e-3,
    label = "coefficient of variation of the response-to-gain ratio"
  )
  # The common ratio is i / sqrt(b' P b), a single number for the whole index.
  expect_equal(unname(mean(ratio)), 0.33585, tolerance = 1e-3)
  # And every trait individually, to the precision the table is printed at.
  for (trait in reference$traits) {
    expect_equal(
      unname(ratio[[trait]]), 0.33585,
      tolerance = 5e-3,
      label = paste("response-to-gain ratio for", trait)
    )
  }
})

test_that("our Pesek-Baker implementation reproduces the same property", {
  # The published numbers satisfy the identity. This checks that our
  # implementation does too, on a covariance structure of our own, which is
  # what makes the previous test evidence about this package rather than only
  # about the article.
  set.seed(19L)
  traits <- paste0("t", 1:5)
  A <- matrix(stats::rnorm(25L), 5L)
  G <- crossprod(A) + diag(5L) * 2
  dimnames(G) <- list(traits, traits)
  P <- G + diag(runif(5L, 1, 3))
  dimnames(P) <- list(traits, traits)

  # The article's own choice of d: the genetic standard deviations.
  d <- stats::setNames(sqrt(diag(G)), traits)
  b <- solve(G) %*% d
  response <- as.numeric(G %*% b)
  ratio <- response / d
  expect_equal(
    unname(stats::sd(ratio) / mean(ratio)), 0,
    tolerance = 1e-10
  )

  # And through the exported interface, where the transform must not break it.
  values <- as.data.frame(matrix(
    stats::rnorm(200L),
    ncol = 5L,
    dimnames = list(paste0("g", 1:40), traits)
  ))
  fit <- selection_index(
    values, traits,
    method = "pesek_baker", G = G, P = P,
    desired_gains = d, scale_traits = FALSE, n_select = 10L
  )
  fitted_ratio <- fit$evaluation$expected_response / d
  expect_equal(
    unname(stats::sd(fitted_ratio) / mean(fitted_ratio)), 0,
    tolerance = 1e-8
  )
})

test_that("the published R_HI of 0.0018 is explained, not reproduced", {
  # Rahimi and Debnath report R_HI = 0.0018 for Pesek-Baker against 0.9887 for
  # the optimum index on the same data. That is not a defect in their
  # arithmetic: it follows from using the desired-gain vector d as the
  # aggregate weights a, which is what their code does and what this package
  # did until version 0.5.0.
  #
  # The quantity is then
  #
  #   R_HI = d'd / (sqrt(d' G^-1 P G^-1 d) * sqrt(d' G d))
  #
  # which is not the correlation between the index and net merit under any
  # definition of merit, and collapses when trait scales differ by orders of
  # magnitude. In this dataset they span 0.224 to 45.046, a factor of 200.
  reference <- rahimi_debnath_2023()
  expect_lt(reference$pesek_baker_criteria[["R_HI"]], 0.01)
  expect_gt(reference$optimum_method1_criteria[["R_HI"]], 0.98)

  scale_span <- max(reference$desired_gains) / min(reference$desired_gains)
  expect_gt(scale_span, 100)

  # Desired gains do not define net merit. For an explicit comparison we use
  # one common aggregate objective; the package must never replace it with a
  # different implied objective for each desired-gain direction.
  set.seed(23L)
  traits <- reference$traits
  p <- length(traits)
  sd_vector <- reference$desired_gains
  correlation <- diag(p)
  correlation[lower.tri(correlation)] <-
    correlation[upper.tri(correlation)] <- 0.25
  G <- outer(sd_vector, sd_vector) * correlation
  dimnames(G) <- list(traits, traits)
  P <- G + diag(sd_vector^2 * 1.5)
  dimnames(P) <- list(traits, traits)

  values <- as.data.frame(matrix(
    stats::rnorm(p * 60L, sd = rep(sd_vector, each = 60L)),
    ncol = p,
    dimnames = list(paste0("g", 1:60), traits)
  ))
  fit <- selection_index(
    values, traits,
    method = "pesek_baker", G = G, P = P,
    desired_gains = sd_vector,
    aggregate_weights = reference$economic_weights_method1,
    scale_traits = FALSE, n_select = 12L
  )
  expect_true(is.finite(fit$evaluation$R_HI))
  expect_lte(abs(fit$evaluation$R_HI), 1)
  expect_equal(
    unname(fit$aggregate_weights),
    unname(reference$economic_weights_method1)
  )
})

test_that("the optimum index outperforms Pesek-Baker on relative efficiency", {
  # A qualitative ordering the article reports and any correct implementation
  # must agree with: the optimum index, which uses economic weights, achieves
  # more response in the main trait than a desired-gain index constrained to a
  # fixed proportional direction.
  reference <- rahimi_debnath_2023()
  expect_gt(
    reference$optimum_method1_criteria[["RE"]],
    reference$pesek_baker_criteria[["RE"]]
  )
  expect_gt(
    reference$optimum_method1_gain[["yield"]],
    reference$pesek_baker_gain[["yield"]]
  )
})

test_that("the frozen reference values carry their provenance", {
  reference <- rahimi_debnath_2023()
  expect_match(reference$citation, "Scientific Reports")
  expect_match(reference$citation, "10.1038/s41598-023-46368-6")
  expect_match(reference$note, "external check")
})
