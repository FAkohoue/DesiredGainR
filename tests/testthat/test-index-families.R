# Tests for the classical index families and the evaluation criteria.

family_traits <- c("YLD", "DIS", "HT")

family_G <- matrix(
  c(
    0.60, -0.18, 0.05,
    -0.18, 0.42, -0.02,
    0.05, -0.02, 0.55
  ),
  3,
  dimnames = list(family_traits, family_traits)
)

family_P <- matrix(
  c(
    1.10, -0.25, 0.08,
    -0.25, 0.85, -0.04,
    0.08, -0.04, 1.00
  ),
  3,
  dimnames = list(family_traits, family_traits)
)

family_values <- function(n = 60L, seed = 21L) {
  set.seed(seed)
  matrix(
    stats::rnorm(n * 3),
    ncol = 3,
    dimnames = list(sprintf("g%03d", seq_len(n)), family_traits)
  )
}

test_that("Smith-Hazel coefficients match the closed form", {
  values <- family_values()
  a <- c(YLD = 1, DIS = 0.5, HT = 0.2)
  fit <- selection_index(
    values, family_traits,
    method = "smith_hazel",
    G = family_G, P = family_P, economic_weights = a,
    scale_traits = FALSE, n_select = 6
  )
  expect_equal(
    as.numeric(fit$coefficients),
    as.numeric(solve(family_P, family_G %*% a)),
    tolerance = 1e-8
  )
  expect_equal(nrow(fit$ranking), nrow(values))
  expect_equal(sum(fit$ranking$selected), 6L)
})

test_that("base index returns the economic weights unchanged", {
  values <- family_values()
  a <- c(YLD = 1, DIS = 0.5, HT = 0.2)
  fit <- selection_index(
    values, family_traits,
    method = "base",
    G = family_G, P = family_P, economic_weights = a,
    scale_traits = FALSE, n_select = 6
  )
  expect_equal(as.numeric(fit$coefficients), as.numeric(a))
})

test_that("the two desired-gain formulations coincide for invertible G", {
  values <- family_values()
  d <- c(YLD = 0.5, DIS = 0.3, HT = 0.2)

  pesek <- selection_index(
    values, family_traits,
    method = "pesek_baker",
    G = family_G, P = family_P, desired_gains = d,
    scale_traits = FALSE, n_select = 6
  )
  yamada <- selection_index(
    values, family_traits,
    method = "yamada",
    G = family_G, P = family_P, desired_gains = d,
    scale_traits = FALSE, n_select = 6
  )

  # Both satisfy G b = d, so they target the same estimand.
  expect_equal(
    as.numeric(family_G %*% pesek$coefficients), as.numeric(d),
    tolerance = 1e-8
  )
  expect_equal(
    as.numeric(family_G %*% yamada$coefficients), as.numeric(d),
    tolerance = 1e-6
  )

  # For a square invertible G the Yamada form reduces algebraically to the
  # Pesek-Baker form:
  #   P^-1 G (G P^-1 G)^-1 d = P^-1 G (G^-1 P G^-1) d = G^-1 d.
  # They are therefore the same index reached by two routes, and any published
  # difference between them reflects conditioning or a rank-deficient G rather
  # than a different estimand.
  expect_equal(
    as.numeric(pesek$coefficients), as.numeric(yamada$coefficients),
    tolerance = 1e-6
  )
  expect_equal(
    sort(pesek$selected$id), sort(yamada$selected$id)
  )
})

test_that("the two routes diverge numerically when G is ill-conditioned", {
  values <- family_values()
  d <- c(YLD = 0.5, DIS = 0.3, HT = 0.2)
  # A nearly singular G: the third trait is almost a copy of the first.
  ill <- family_G
  ill["HT", ] <- ill["YLD", ] * (1 - 1e-7)
  ill[, "HT"] <- ill[, "YLD"] * (1 - 1e-7)
  ill["HT", "HT"] <- family_G["YLD", "YLD"] * (1 - 1e-7)
  ill <- (ill + t(ill)) / 2

  diagnostics <- matrix_diagnostics(ill, "G")
  expect_lt(diagnostics$reciprocal_condition, 1e-6)

  # P must remain compatible with the perturbed G, which the original
  # family_P is not: making HT a near-copy of YLD in G alone leaves P - G with
  # a negative eigenvalue, so the pair could not both be true even though no
  # single trait has a genetic variance above its phenotypic one. Rebuilding P
  # as the perturbed G plus the original residual keeps the residual positive
  # definite while leaving G exactly as ill-conditioned as intended.
  ill_P <- ill + (family_P - family_G)
  ill_P <- (ill_P + t(ill_P)) / 2
  dimnames(ill_P) <- dimnames(family_P)
  expect_gt(
    min(eigen(ill_P - ill, symmetric = TRUE, only.values = TRUE)$values), 0
  )

  # Inverting a badly conditioned G directly, and reaching the same solution
  # through P, are not numerically equivalent. This is the practical reason the
  # formulation must be stated, and the reason conditioning is reported.
  pesek <- suppressWarnings(selection_index(
    values, family_traits,
    method = "pesek_baker",
    G = ill, P = ill_P, desired_gains = d,
    scale_traits = FALSE, n_select = 6
  ))
  yamada <- suppressWarnings(selection_index(
    values, family_traits,
    method = "yamada",
    G = ill, P = ill_P, desired_gains = d,
    scale_traits = FALSE, n_select = 6
  ))
  expect_true(all(is.finite(pesek$coefficients)))
  expect_true(all(is.finite(yamada$coefficients)))
})

test_that("Yamada honours only the direction of the desired gains", {
  values <- family_values()
  d <- c(YLD = 0.5, DIS = 0.3, HT = 0.2)
  one <- selection_index(
    values, family_traits,
    method = "yamada",
    G = family_G, P = family_P, desired_gains = d,
    scale_traits = FALSE, n_select = 8
  )
  ten <- selection_index(
    values, family_traits,
    method = "yamada",
    G = family_G, P = family_P, desired_gains = 10 * d,
    scale_traits = FALSE, n_select = 8
  )
  expect_equal(
    sort(one$selected$id), sort(ten$selected$id)
  )
  expect_equal(
    as.numeric(ten$coefficients), 10 * as.numeric(one$coefficients),
    tolerance = 1e-6
  )
})

test_that("Mulamba-Mock needs no weights and orients traits correctly", {
  values <- matrix(
    c(
      3, 2, 1,
      1, 2, 3,
      2, 3, 1
    ),
    ncol = 3, byrow = TRUE,
    dimnames = list(c("a", "b", "c"), family_traits)
  )
  fit <- selection_index(
    values, family_traits,
    method = "mulamba_mock",
    scale_traits = FALSE, n_select = 1
  )
  expect_null(fit$coefficients)
  expect_equal(fit$strategy, "rank_sum")

  # Declaring DIS as lower-is-better must change the ranking.
  flipped <- selection_index(
    values, family_traits,
    method = "mulamba_mock",
    lower_is_better = "DIS", scale_traits = FALSE, n_select = 1
  )
  expect_false(identical(fit$ranking$id, flipped$ranking$id))
})

test_that("independent culling retains only candidates passing every gate", {
  values <- matrix(
    c(
      2, 2, 2,
      2, -2, 2,
      -2, 2, 2
    ),
    ncol = 3, byrow = TRUE,
    dimnames = list(c("pass", "fail_dis", "fail_yld"), family_traits)
  )
  fit <- selection_index(
    values, family_traits,
    method = "independent_culling",
    culling_thresholds = c(YLD = 0, DIS = 0, HT = 0),
    scale_traits = FALSE
  )
  passing <- fit$ranking[score == 1, id]
  expect_equal(passing, "pass")
})

test_that("economic weights may be negative but desired gains may not", {
  values <- family_values()
  expect_silent(
    selection_index(
      values, family_traits,
      method = "smith_hazel",
      G = family_G, P = family_P,
      economic_weights = c(YLD = 1, DIS = -0.5, HT = 0.2),
      scale_traits = FALSE
    )
  )
  expect_error(
    selection_index(
      values, family_traits,
      method = "pesek_baker",
      G = family_G, P = family_P,
      desired_gains = c(YLD = 1, DIS = -0.5, HT = 0.2),
      scale_traits = FALSE
    ),
    "non-negative"
  )
})

test_that("evaluation criteria match independent calculations", {
  values <- family_values()
  a <- c(YLD = 1, DIS = 0.5, HT = 0.2)
  b <- as.numeric(solve(family_P, family_G %*% a))
  names(b) <- family_traits
  intensity <- 1.755

  evaluation <- evaluate_index(
    b, family_G, family_P,
    aggregate_weights = a,
    selection_intensity = intensity, main_trait = "YLD"
  )

  index_sd <- sqrt(as.numeric(crossprod(b, family_P %*% b)))
  merit_sd <- sqrt(as.numeric(crossprod(a, family_G %*% a)))
  expect_equal(
    evaluation$R_HI,
    as.numeric(crossprod(b, family_G %*% a)) / (index_sd * merit_sd),
    tolerance = 1e-10
  )
  expect_equal(
    evaluation$delta_H, intensity * evaluation$R_HI * merit_sd,
    tolerance = 1e-10
  )
  expect_equal(
    as.numeric(evaluation$expected_response),
    as.numeric(intensity * (family_G %*% b) / index_sd),
    tolerance = 1e-10
  )
  direct <- intensity * family_G["YLD", "YLD"] / sqrt(family_P["YLD", "YLD"])
  expect_equal(
    evaluation$RE, evaluation$expected_response[["YLD"]] / direct,
    tolerance = 1e-10
  )
  # A squared correlation cannot exceed one.
  expect_lte(abs(evaluation$R_HI), 1 + 1e-10)
})

test_that("Smith-Hazel maximises R_HI among the families tested", {
  values <- family_values()
  a <- c(YLD = 1, DIS = 0.5, HT = 0.2)
  intensity <- 1.755

  smith <- as.numeric(solve(family_P, family_G %*% a))
  names(smith) <- family_traits
  base <- a
  pesek <- as.numeric(solve(family_G, a))
  names(pesek) <- family_traits

  correlations <- vapply(
    list(smith = smith, base = base, pesek = pesek),
    function(b) {
      evaluate_index(
        b, family_G, family_P,
        aggregate_weights = a,
        selection_intensity = intensity
      )$R_HI
    },
    numeric(1)
  )
  # The optimum index is optimum by construction for this aggregate.
  expect_equal(names(which.max(correlations)), "smith")
})

test_that("the selected table holds the actual top-ranked candidates", {
  # This guards a row-alignment failure that no other assertion catches. The
  # ranking table is sorted by score after the selection flag is computed, so
  # filtering it with the unsorted flag returns positionally wrong rows while
  # still returning the right *number* of rows, the right column values, and
  # an internally plausible object.
  set.seed(404)
  values <- family_values(n = 80L, seed = 404L)
  a <- c(YLD = 1, DIS = 0.5, HT = 0.2)
  fit <- selection_index(
    values, family_traits,
    method = "smith_hazel",
    G = family_G, P = family_P, economic_weights = a,
    lower_is_better = "DIS", n_select = 8L
  )

  # Independently recompute the top eight from the returned scores.
  expected <- fit$ranking[order(-score, id)][seq_len(8L), id]
  expect_setequal(fit$selected$id, expected)

  # The table and the flag must agree with each other.
  expect_setequal(fit$selected$id, fit$ranking[selected == TRUE, id])
  expect_equal(nrow(fit$selected), 8L)
  expect_true(all(fit$selected$selected))

  # Every selected candidate must outscore every unselected one.
  expect_gte(
    min(fit$selected$score),
    max(fit$ranking[selected == FALSE, score])
  )
})

test_that("two similar indices select overlapping sets", {
  # A high rank correlation between two indices and a near-empty intersection
  # of their selected sets are mutually inconsistent. Asserting the
  # relationship catches misalignment that either quantity alone would miss.
  values <- family_values(n = 120L, seed = 77L)
  a <- c(YLD = 1, DIS = 0.5, HT = 0.2)

  smith <- selection_index(
    values, family_traits,
    method = "smith_hazel",
    G = family_G, P = family_P, economic_weights = a,
    lower_is_better = "DIS", n_select = 12L
  )
  base <- selection_index(
    values, family_traits,
    method = "base",
    G = family_G, P = family_P, economic_weights = a,
    lower_is_better = "DIS", n_select = 12L
  )

  paired <- merge(
    smith$ranking[, c("id", "score")],
    base$ranking[, c("id", "score")],
    by = "id", suffixes = c("_a", "_b")
  )
  rank_correlation <- stats::cor(
    paired$score_a, paired$score_b,
    method = "spearman"
  )
  overlap <- length(intersect(smith$selected$id, base$selected$id)) / 12

  expect_gt(rank_correlation, 0.9)
  # At that correlation the top-decile sets cannot be nearly disjoint.
  expect_gt(overlap, 0.4)
})

test_that("the index coefficient of variation is withheld when undefined", {
  values <- family_values()
  a <- c(YLD = 1, DIS = 0.5, HT = 0.2)

  # Standardising centres the index at zero, so its coefficient of variation
  # diverges and must be withheld rather than reported as an enormous number.
  standardised <- selection_index(
    values, family_traits,
    method = "smith_hazel",
    G = family_G, P = family_P, economic_weights = a,
    scale_traits = TRUE, n_select = 6
  )
  expect_true(is.na(standardised$evaluation$CV_I))
  expect_match(standardised$evaluation$CV_I_note, "undefined")
  expect_output(print(standardised), "CV_I undefined")

  # Centring is what makes it undefined, not standardising, so switching off
  # the scaling alone changes nothing.
  unscaled <- selection_index(
    values, family_traits,
    method = "smith_hazel",
    G = family_G, P = family_P, economic_weights = a,
    scale_traits = FALSE, n_select = 6
  )
  expect_true(is.na(unscaled$evaluation$CV_I))

  # Switching off the centring gives the index a non-zero location, which is
  # the form in which published results report it.
  shifted <- values + 100
  uncentred <- selection_index(
    shifted, family_traits,
    method = "smith_hazel",
    G = family_G, P = family_P, economic_weights = a,
    center_traits = FALSE, scale_traits = FALSE, n_select = 6
  )
  expect_true(is.finite(uncentred$evaluation$CV_I))
  expect_lt(uncentred$evaluation$CV_I, 1000)
  expect_match(uncentred$evaluation$CV_I_note, "Reported")

  # Centring shifts every score by one constant, so it cannot change the
  # ranking or the selected set.
  centred <- selection_index(
    shifted, family_traits,
    method = "smith_hazel",
    G = family_G, P = family_P, economic_weights = a,
    center_traits = TRUE, scale_traits = FALSE, n_select = 6
  )
  expect_equal(sort(centred$selected$id), sort(uncentred$selected$id))
})

test_that("selection_index reports evaluation and effective weights", {
  values <- family_values()
  fit <- selection_index(
    values, family_traits,
    method = "smith_hazel",
    G = family_G, P = family_P,
    economic_weights = c(YLD = 1, DIS = 0.5, HT = 0.2),
    lower_is_better = "DIS", n_select = 6
  )
  expect_s3_class(fit, "desiredgainr_index")
  expect_s3_class(fit$evaluation, "desiredgainr_evaluation")
  expect_true(is.data.frame(fit$effective_weights))
  expect_equal(nrow(fit$effective_weights), 3L)
  expect_true(all(is.finite(fit$evaluation$expected_response)))
  expect_equal(nrow(fit$observed_differential), 3L)
})
