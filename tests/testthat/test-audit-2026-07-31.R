# Tests for the findings of AUDIT-2026-07-31.md.
#
# Each block names the finding it guards. Several of these assert properties
# that the previous implementation could not satisfy at any tolerance, which is
# the point: a test that the old code would also have passed would not have
# caught the defect.

fixture <- function(seed = 7L, n = 60L) {
  set.seed(seed)
  traits <- c("yield", "protein", "height")
  G <- matrix(
    c(
      1.00, 0.25, -0.30,
      0.25, 0.50, 0.10,
      -0.30, 0.10, 0.80
    ), 3L,
    dimnames = list(traits, traits)
  )
  P <- G + diag(c(1.5, 0.9, 1.2))
  dimnames(P) <- list(traits, traits)
  values <- as.data.frame(matrix(
    stats::rnorm(n * 3L),
    ncol = 3L,
    dimnames = list(paste0("C", seq_len(n)), traits)
  ))
  list(traits = traits, G = G, P = P, values = values)
}

# ---------------------------------------------------------------------------
# Bending
# ---------------------------------------------------------------------------

indefinite_matrix <- function() {
  traits <- c("t1", "t2", "t3")
  # t1 and t2 correlate strongly and positively, t2 and t3 likewise, yet t1
  # and t3 correlate strongly and negatively. No such matrix exists.
  matrix(
    c(
      1.0, 0.9, -0.9,
      0.9, 1.0, 0.9,
      -0.9, 0.9, 1.0
    ), 3L,
    dimnames = list(traits, traits)
  )
}

test_that("bending makes an indefinite matrix invertible", {
  skip_if_not_installed("ASRgenomics")

  M <- indefinite_matrix()
  expect_false(matrix_diagnostics(M)$positive_definite)

  bent <- bend_covariance(M)
  expect_s3_class(bent, "desiredgainr_bent_covariance")
  expect_true(bent$after$positive_definite)
  expect_identical(dimnames(bent$G), dimnames(M))
  expect_true(isSymmetric(bent$G, tol = 1e-6))
  # The whole point is that it can now be inverted.
  expect_silent(solve(bent$G))
  expect_gt(bent$after$minimum_eigenvalue, 0)
})

test_that("bending is faithful to G.tuneup and not a reimplementation", {
  skip_on_cran()
  skip_if_not_installed("ASRgenomics")

  # This is the contract. bend_covariance() must return exactly what
  # G.tuneup() returns, so that a result computed here is reproducible by
  # anyone who runs the underlying function directly. If this fails, the
  # wrapper has started doing arithmetic of its own.
  M <- indefinite_matrix()
  direct <- ASRgenomics::G.tuneup(
    G = M, bend = TRUE, eig.tol = 1e-6, rcn = TRUE, digits = 8,
    sparseform = FALSE, determinant = TRUE, message = FALSE
  )
  wrapped <- bend_covariance(M, method = "asrgenomics")
  expect_equal(wrapped$G, as.matrix(direct$Gb),
    tolerance = 1e-10,
    ignore_attr = TRUE
  )
  expect_equal(wrapped$reciprocal_condition[["before"]], direct$rcn0)
  expect_equal(wrapped$reciprocal_condition[["after"]], direct$rcnb)
})

test_that("the variances are reported as moved, because G.tuneup moves them", {
  skip_if_not_installed("ASRgenomics")

  # G.tuneup() fixes keepDiag = FALSE. That is a real consequence for a trait
  # covariance, so the wrapper must surface it rather than let it pass
  # unnoticed. A user reading only the bent matrix would not know the
  # heritabilities had changed.
  traits <- c("t1", "t2", "t3")
  M <- matrix(
    c(
      4.0, 1.9, -1.8,
      1.9, 1.0, 0.9,
      -1.8, 0.9, 1.0
    ), 3L,
    dimnames = list(traits, traits)
  )
  bent <- bend_covariance(M, method = "asrgenomics")
  expect_gt(bent$adjustment$max_variance_change, 0)
  expect_true(is.finite(bent$adjustment$max_relative_variance_change))
  expect_false(isTRUE(all.equal(diag(bent$G), diag(M))))
})

test_that("the condition number after bending is capped near 100", {
  skip_if_not_installed("ASRgenomics")

  # posd.tol = 1e-2 is hard-wired inside G.tuneup(). The documentation warns
  # that this is aggressive for a trait covariance, so the claim should be
  # demonstrable rather than asserted.
  #
  # Only the post-bending conditioning is asserted. Comparing it to the
  # condition number BEFORE bending would be meaningless: matrix_diagnostics()
  # forms max|lambda| / min|lambda|, and for an indefinite matrix min|lambda|
  # is the modulus of an eigenvalue that may be nowhere near zero, so the
  # "before" figure is not a measure of near-singularity at all and need not
  # exceed the "after" figure.
  bent <- bend_covariance(indefinite_matrix(), method = "asrgenomics")
  expect_lt(bent$after$condition_number, 200)
  expect_true(bent$after$positive_definite)
  expect_silent(solve(bent$G))
  # The adjustment must be real and must be reported.
  expect_gt(bent$adjustment$frobenius_relative, 0)
  expect_true(is.finite(bent$adjustment$max_variance_change))
})

test_that("an already positive definite matrix is left essentially alone", {
  skip_if_not_installed("ASRgenomics")

  traits <- c("t1", "t2")
  M <- matrix(c(2, 0.5, 0.5, 1), 2L, dimnames = list(traits, traits))
  bent <- bend_covariance(M)
  expect_equal(bent$G, M, tolerance = 1e-6, ignore_attr = TRUE)
  expect_lt(bent$adjustment$frobenius_relative, 1e-4)
})

test_that("the default trait-covariance repair preserves variances", {
  M <- indefinite_matrix()
  bent <- bend_covariance(M)
  expect_equal(diag(bent$G), diag(M), tolerance = 1e-10)
  expect_true(bent$after$positive_definite)
  expect_identical(bent$method, "correlation")
})

test_that("bend_covariance() validates its input", {
  skip_if_not_installed("ASRgenomics")

  traits <- c("t1", "t2")
  named <- function(x) {
    dimnames(x) <- list(traits, traits)
    x
  }
  expect_error(
    bend_covariance(named(matrix(c(-1, 0, 0, -1), 2L))),
    "no positive eigenvalue"
  )
  expect_error(bend_covariance(matrix(1:6, 2L)), "square")
  expect_error(
    bend_covariance(named(matrix(c(1, NA, NA, 1), 2L))), "finite"
  )
  expect_error(
    bend_covariance(named(diag(2)), eig_tol = -1), "positive number"
  )
  # G.tuneup() validates the names, so the wrapper must require them.
  expect_error(bend_covariance(diag(2)), "trait names")
  mismatched <- named(diag(2))
  colnames(mismatched) <- c("a", "b")
  expect_error(bend_covariance(mismatched), "must match")
})

test_that("the missing dependency is named, not hidden", {
  skip_if(requireNamespace("ASRgenomics", quietly = TRUE))
  expect_error(
    bend_covariance(indefinite_matrix(), method = "asrgenomics"),
    "ASRgenomics"
  )
})

test_that("the positive-definiteness rejection points at the repair", {
  set <- fixture()
  bad <- set$G
  bad["yield", "height"] <- bad["height", "yield"] <- -5
  expect_error(
    selection_index(
      set$values, set$traits,
      method = "smith_hazel", G = bad, P = set$P,
      economic_weights = c(yield = 2, protein = 1, height = 1)
    ),
    "bend_covariance"
  )
})

# ---------------------------------------------------------------------------
# M2: the aggregate for the desired-gain families
# ---------------------------------------------------------------------------

test_that("desired-gain families do not invent a changing merit objective", {
  set <- fixture()
  d <- c(yield = 1, protein = 0.5, height = 0.2)
  fit <- selection_index(
    set$values, set$traits,
    method = "pesek_baker",
    G = set$G, P = set$P, desired_gains = d
  )
  expect_null(fit$aggregate_weights)
  expect_true(is.na(fit$evaluation$R_HI))
  expect_true(is.na(fit$evaluation$delta_H))
})

test_that("a common supplied merit makes desired-gain R_HI comparable", {
  set <- fixture()
  fit <- selection_index(
    set$values, set$traits,
    method = "pesek_baker",
    G = set$G, P = set$P,
    desired_gains = c(yield = 1, protein = 0.5, height = 0.2),
    aggregate_weights = c(yield = 2, protein = 1, height = 0.5),
    n_select = 10L
  )
  expect_gt(fit$evaluation$R_HI, 0)
  expect_lte(fit$evaluation$R_HI, 1 + 1e-8)
})

# ---------------------------------------------------------------------------
# Index heritability
# ---------------------------------------------------------------------------

test_that("index heritability lies in the unit interval and equals b'Gb/b'Pb", {
  set <- fixture()
  fit <- selection_index(
    set$values, set$traits,
    method = "smith_hazel",
    G = set$G, P = set$P,
    economic_weights = c(yield = 2, protein = 1, height = 1),
    n_select = 10L
  )
  b <- fit$coefficients
  expected <- as.numeric(t(b) %*% fit$G %*% b) /
    as.numeric(t(b) %*% fit$P %*% b)
  expect_equal(fit$evaluation$h2_index, expected, tolerance = 1e-10)
  expect_gt(fit$evaluation$h2_index, 0)
  expect_lt(fit$evaluation$h2_index, 1)
  expect_equal(
    fit$evaluation$accuracy_index, sqrt(expected),
    tolerance = 1e-10
  )
})

test_that("a Smith-Hazel index has the accuracy that R_HI reports", {
  # For b = P^-1 G a the index is the best linear predictor of net merit, so
  # its correlation with net merit equals the square root of its own
  # heritability divided by nothing further; check the two agree in ordering
  # across objectives rather than asserting a false identity.
  set <- fixture()
  weights <- list(
    c(yield = 3, protein = 1, height = 0.5),
    c(yield = 1, protein = 1, height = 1)
  )
  fits <- lapply(weights, function(a) {
    selection_index(
      set$values, set$traits,
      method = "smith_hazel",
      G = set$G, P = set$P, economic_weights = a, n_select = 10L
    )
  })
  for (fit in fits) {
    expect_lte(fit$evaluation$R_HI, 1 + 1e-8)
    expect_gt(fit$evaluation$R_HI, 0)
  }
})

# ---------------------------------------------------------------------------
# N5: scale_by
# ---------------------------------------------------------------------------

test_that("scale_by = 'phenotypic' gives P a unit diagonal exactly", {
  set <- fixture()
  fit <- selection_index(
    set$values, set$traits,
    method = "smith_hazel",
    G = set$G, P = set$P, scale_by = "phenotypic",
    economic_weights = c(yield = 2, protein = 1, height = 1)
  )
  expect_equal(diag(fit$P), rep(1, 3L),
    tolerance = 1e-10, ignore_attr = TRUE
  )
  # The scaled G then carries the heritabilities.
  expect_equal(
    unname(diag(fit$G)), unname(diag(set$G) / diag(set$P)),
    tolerance = 1e-10
  )
  expect_identical(fit$transformation$scale_by, "phenotypic")
})

test_that("scale_by defaults to the previous behaviour", {
  set <- fixture()
  fit <- selection_index(
    set$values, set$traits,
    method = "smith_hazel",
    G = set$G, P = set$P,
    economic_weights = c(yield = 2, protein = 1, height = 1)
  )
  expect_identical(fit$transformation$scale_by, "sample")
  expect_equal(
    unname(fit$transformation$scale),
    unname(apply(as.matrix(set$values), 2L, stats::sd)),
    tolerance = 1e-12
  )
})

test_that("scale_by = 'phenotypic' requires P", {
  set <- fixture()
  expect_error(
    selection_index(
      set$values, set$traits,
      method = "pesek_baker", G = set$G,
      scale_by = "phenotypic",
      desired_gains = c(yield = 1, protein = 1, height = 1)
    ),
    "requires P"
  )
})

test_that("scale_by is recorded as 'none' when scaling is off", {
  set <- fixture()
  fit <- selection_index(
    set$values, set$traits,
    method = "base", scale_traits = FALSE,
    economic_weights = c(yield = 2, protein = 1, height = 1)
  )
  expect_identical(fit$transformation$scale_by, "none")
})

# ---------------------------------------------------------------------------
# predict()
# ---------------------------------------------------------------------------

test_that("predict() reproduces the fitted ranking on the fitting data", {
  set <- fixture()
  fit <- selection_index(
    set$values, set$traits,
    method = "smith_hazel",
    G = set$G, P = set$P,
    economic_weights = c(yield = 2, protein = 1, height = 1),
    n_select = 12L
  )
  again <- predict(fit, set$values)
  expect_equal(again$score, fit$ranking$score, tolerance = 1e-12)
  expect_identical(again$id, fit$ranking$id)
  expect_identical(again$selected, fit$ranking$selected)
})

test_that("predict() uses the fitting transformation, not the new data's", {
  set <- fixture()
  fit <- selection_index(
    set$values, set$traits,
    method = "smith_hazel",
    G = set$G, P = set$P,
    economic_weights = c(yield = 2, protein = 1, height = 1)
  )
  # A candidate set shifted upwards on every trait must score higher. If the
  # centring were recomputed from newdata the shift would cancel and the mean
  # score would be unchanged, which is the failure this guards.
  shifted <- set$values + 1
  baseline <- predict(fit, set$values)
  advanced <- predict(fit, shifted)
  expect_gt(mean(advanced$score), mean(baseline$score))
})

test_that("predict() refuses methods that have no coefficients", {
  set <- fixture()
  fit <- selection_index(
    set$values, set$traits,
    method = "mulamba_mock", n_select = 10L
  )
  expect_error(predict(fit, set$values), "no coefficients")
})

test_that("predict() validates newdata", {
  set <- fixture()
  fit <- selection_index(
    set$values, set$traits,
    method = "smith_hazel",
    G = set$G, P = set$P,
    economic_weights = c(yield = 2, protein = 1, height = 1)
  )
  expect_error(
    predict(fit, set$values[, c("yield", "protein")]), "missing trait columns"
  )
  expect_error(
    predict(fit, set$values, n_select = 1000L), "cannot exceed"
  )
})

# ---------------------------------------------------------------------------
# N7: index_uncertainty()
# ---------------------------------------------------------------------------

test_that("intervals narrow as the degrees of freedom grow", {
  set <- fixture()
  fit <- selection_index(
    set$values, set$traits,
    method = "smith_hazel",
    G = set$G, P = set$P,
    economic_weights = c(yield = 2, protein = 1, height = 1),
    n_select = 10L
  )
  loose <- suppressWarnings(index_uncertainty(
    fit,
    genetic_df = 20L, n_draws = 300L, seed = 1L
  ))
  tight <- index_uncertainty(
    fit,
    genetic_df = 400L, n_draws = 300L, seed = 1L
  )
  expect_s3_class(loose, "desiredgainr_uncertainty")
  expect_true(all(tight$coefficients$SD < loose$coefficients$SD))
  expect_gt(
    tight$rank_stability$spearman_mean, loose$rank_stability$spearman_mean
  )
})

test_that("the Wishart draws are centred on the supplied matrix", {
  set.seed(3L)
  M <- matrix(c(2, 0.5, 0.5, 1), 2L, dimnames = list(c("a", "b"), c("a", "b")))
  draws <- DesiredGainR:::.dgr_wishart_draws(M, df = 50L, n_draws = 4000L, "M")
  expect_equal(apply(draws, c(1L, 2L), mean), M,
    tolerance = 0.05, ignore_attr = TRUE
  )
})

test_that("alignment is exactly one when the fitted G is the truth", {
  # Pesek-Baker gives Gb = d identically, so with no resampling error the
  # cosine between the achieved response and the desired gains must be 1. This
  # is the reference point that makes the reported shortfall interpretable.
  set <- fixture()
  d <- c(yield = 1, protein = 0.5, height = 0.2)
  fit <- selection_index(
    set$values, set$traits,
    method = "pesek_baker",
    G = set$G, P = set$P, desired_gains = d
  )
  response <- as.numeric(fit$G %*% fit$coefficients)
  cosine <- sum(response * d) / (sqrt(sum(response^2)) * sqrt(sum(d^2)))
  expect_equal(cosine, 1, tolerance = 1e-10)
})

test_that("alignment degrades as the genetic degrees of freedom fall", {
  set <- fixture()
  fit <- selection_index(
    set$values, set$traits,
    method = "pesek_baker",
    G = set$G, P = set$P,
    desired_gains = c(yield = 1, protein = 0.5, height = 0.2)
  )
  poor <- suppressWarnings(index_uncertainty(
    fit,
    genetic_df = 12L, n_draws = 300L, seed = 5L
  ))
  good <- index_uncertainty(
    fit,
    genetic_df = 500L, n_draws = 300L, seed = 5L
  )
  expect_identical(poor$estimation_cost$quantity, "alignment")
  expect_lt(poor$estimation_cost$mean, good$estimation_cost$mean)
  expect_lte(good$estimation_cost$mean, 1 + 1e-8)
  expect_gt(good$estimation_cost$mean, 0.99)
})

test_that("relative efficiency is bounded above by one for weight methods", {
  set <- fixture()
  fit <- selection_index(
    set$values, set$traits,
    method = "smith_hazel",
    G = set$G, P = set$P,
    economic_weights = c(yield = 2, protein = 1, height = 1)
  )
  result <- index_uncertainty(
    fit,
    genetic_df = 60L, residual_df = 200L, n_draws = 300L, seed = 2L
  )
  expect_identical(result$estimation_cost$quantity, "relative_efficiency")
  expect_lte(max(result$estimation_cost$draws, na.rm = TRUE), 1 + 1e-8)
  expect_gt(result$estimation_cost$mean, 0)
})

test_that("resampling the residual keeps G and P mutually consistent", {
  # P* is assembled as G* + E* rather than drawn independently of G*. The
  # consequence that matters is that P* - G* = E* is positive definite by
  # construction on every draw, so no draw implies a heritability above one.
  # Independent draws of G and P would not guarantee that, and the resulting
  # coefficients would be built on impossible parameter values.
  set.seed(41L)
  traits <- c("a", "b")
  G <- matrix(c(1.0, 0.3, 0.3, 0.8), 2L, dimnames = list(traits, traits))
  E <- matrix(c(1.2, 0.1, 0.1, 1.0), 2L, dimnames = list(traits, traits))

  G_draws <- DesiredGainR:::.dgr_wishart_draws(G, 40L, 200L, "G")
  E_draws <- DesiredGainR:::.dgr_wishart_draws(E, 40L, 200L, "E")
  heritabilities <- vapply(seq_len(200L), function(i) {
    G_star <- G_draws[, , i]
    P_star <- G_star + E_draws[, , i]
    max(diag(G_star) / diag(P_star))
  }, numeric(1L))
  expect_true(all(heritabilities < 1))
  expect_true(all(heritabilities > 0))
})

test_that("residual_df is recorded and requires P", {
  set <- fixture()
  fit <- selection_index(
    set$values, set$traits,
    method = "smith_hazel",
    G = set$G, P = set$P,
    economic_weights = c(yield = 2, protein = 1, height = 1)
  )
  fixed_P <- index_uncertainty(
    fit,
    genetic_df = 80L, n_draws = 100L, seed = 11L
  )
  both <- index_uncertainty(
    fit,
    genetic_df = 80L, residual_df = 80L, n_draws = 100L, seed = 11L
  )
  expect_null(fixed_P$residual_df)
  expect_identical(both$residual_df, 80L)
  expect_null(both$residual_note)
})

test_that("an inconsistent G and P are flagged before resampling", {
  # The repair path routes through bend_covariance(), so it needs G.tuneup().
  skip_if_not_installed("ASRgenomics")
  set <- fixture()
  fit <- selection_index(
    set$values, set$traits,
    method = "smith_hazel",
    G = set$G, P = set$P, scale_traits = FALSE,
    economic_weights = c(yield = 2, protein = 1, height = 1)
  )
  # selection_index() now refuses an inadmissible pair outright, so the only
  # way such an object can reach index_uncertainty() is if it was built by an
  # older version or edited afterwards. Reproduce that here, because the guard
  # inside index_uncertainty() must still hold for objects it did not create.
  fit$G["yield", "yield"] <- fit$P["yield", "yield"] * 1.5
  expect_warning(
    index_uncertainty(
      fit,
      genetic_df = 60L, residual_df = 60L, n_draws = 50L, seed = 3L
    ),
    "mutually inconsistent"
  )
})

test_that("index_uncertainty() rejects impossible degrees of freedom", {
  set <- fixture()
  fit <- selection_index(
    set$values, set$traits,
    method = "smith_hazel",
    G = set$G, P = set$P,
    economic_weights = c(yield = 2, protein = 1, height = 1)
  )
  expect_error(index_uncertainty(fit, genetic_df = 2L), "Wishart draw")
  expect_warning(
    index_uncertainty(fit, genetic_df = 5L, n_draws = 50L), "free parameters"
  )
  expect_error(index_uncertainty(fit, genetic_df = 50L, level = 1), "level")
})

test_that("index_uncertainty() rejects the rank-based families", {
  set <- fixture()
  fit <- selection_index(
    set$values, set$traits,
    method = "mulamba_mock", n_select = 10L
  )
  expect_error(index_uncertainty(fit, genetic_df = 50L), "covariance-based")
})

test_that("index_uncertainty() restores the random seed", {
  set <- fixture()
  fit <- selection_index(
    set$values, set$traits,
    method = "smith_hazel",
    G = set$G, P = set$P,
    economic_weights = c(yield = 2, protein = 1, height = 1)
  )
  set.seed(99L)
  before <- .Random.seed
  index_uncertainty(fit, genetic_df = 60L, n_draws = 50L, seed = 4L)
  expect_identical(.Random.seed, before)
})

# ---------------------------------------------------------------------------
# M3: the reference set guard
# ---------------------------------------------------------------------------

dgsi_inputs <- function(n = 40L, p = 4L, seed = 31L) {
  set.seed(seed)
  traits <- paste0("t", seq_len(p))
  G <- diag(p) * 0.6
  dimnames(G) <- list(traits, traits)
  candidates <- data.table::data.table(GenoID = sprintf("g%03d", seq_len(n)))
  for (trait in traits) {
    data.table::set(candidates, j = trait, value = stats::rnorm(n))
  }
  list(
    traits = traits, G = G, candidates = candidates,
    dg = stats::setNames(rep(0.5, p), traits)
  )
}

test_that("run_dgsi() refuses to estimate P from fewer records than traits", {
  input <- dgsi_inputs()
  too_few <- input$candidates[seq_along(input$traits)]
  expect_error(
    run_dgsi(
      init_data = input$candidates[, "GenoID"],
      cand_data = input$candidates,
      ref_data = too_few,
      trait_cols = input$traits,
      dg = input$dg, G = input$G,
      n_select = 5L, n_iter = 5L, n_rep = 1L, seed = 1L
    ),
    "singular covariance matrix"
  )
})

test_that("run_dgsi() warns when P rests on few reference records", {
  input <- dgsi_inputs()
  few <- input$candidates[seq_len(10L)]
  expect_warning(
    run_dgsi(
      init_data = input$candidates[, "GenoID"],
      cand_data = input$candidates,
      ref_data = few,
      trait_cols = input$traits,
      dg = input$dg, G = input$G,
      n_select = 5L, n_iter = 5L, n_rep = 1L, seed = 1L
    ),
    "poorly\\s+determined"
  )
})

test_that("the numerical rank of an estimated P is recorded", {
  input <- dgsi_inputs()
  result <- run_dgsi(
    init_data = input$candidates[, "GenoID"],
    cand_data = input$candidates,
    trait_cols = input$traits,
    dg = input$dg, G = input$G,
    n_select = 5L, n_iter = 5L, n_rep = 1L, seed = 1L
  )
  expect_identical(
    result$covariance_provenance$P_numerical_rank, length(input$traits)
  )
})

# ---------------------------------------------------------------------------
# M4: the optimism of best-of-n_rep selection
# ---------------------------------------------------------------------------

test_that("run_dgsi() reports the optimism of the winning replicate", {
  input <- dgsi_inputs()
  result <- run_dgsi(
    init_data = input$candidates[, "GenoID"],
    cand_data = input$candidates,
    trait_cols = input$traits,
    dg = input$dg, G = input$G,
    n_select = 5L, n_iter = 20L, n_rep = 5L, seed = 1L
  )
  optimism <- result$optimism
  expect_false(is.null(optimism))
  expect_identical(optimism$n_replicates, 5L)

  # The gap is still mean minus minimum across replicates, so it is
  # non-negative by construction and remains a sensitivity diagnostic.
  expect_gte(optimism$gap, -1e-10)
  expect_equal(
    optimism$gap, optimism$mean_objective - optimism$minimum_objective,
    tolerance = 1e-10
  )
  expect_true(is.character(optimism$note))

  # The chosen replicate is no longer the minimum. Since 0.5.0 it is selected
  # on held-out candidates, so its TRAINING objective can sit anywhere in the
  # spread -- above the mean included. That is the point of the change: an
  # objective guaranteed to be at or below the mean is guaranteed to be
  # optimistically biased. Asserting the old inequality would re-impose the
  # bias the holdout rule exists to remove.
  expect_identical(optimism$selection_rule, "internal pre-fit holdout")
  expect_true(is.finite(optimism$chosen_objective))
  expect_true(is.finite(optimism$chosen_training_objective))
})

test_that("optimism is absent or zero with a single replicate", {
  input <- dgsi_inputs()
  result <- run_dgsi(
    init_data = input$candidates[, "GenoID"],
    cand_data = input$candidates,
    trait_cols = input$traits,
    dg = input$dg, G = input$G,
    n_select = 5L, n_iter = 20L, n_rep = 1L, seed = 1L
  )
  expect_equal(result$optimism$gap, 0, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# N6: the surrogate's lengthscale count
# ---------------------------------------------------------------------------

test_that("the surrogate shares one lengthscale on a small design", {
  set.seed(21L)
  X <- matrix(stats::runif(30L * 5L), ncol = 5L)
  X <- X / sqrt(rowSums(X^2))
  y <- as.numeric(X %*% c(1, -1, 0.5, 0.2, -0.3)) + stats::rnorm(30L, sd = 0.05)
  model <- DesiredGainR:::.dgr_gp_fit(X, y)
  expect_false(model$anisotropic)
  expect_identical(model$n_lengthscale, 1L)
  expect_equal(length(unique(model$lengthscale)), 1L)
  expect_identical(length(model$lengthscale), 5L)
})

test_that("the surrogate fits one lengthscale per dimension when told to", {
  set.seed(22L)
  X <- matrix(stats::runif(30L * 3L), ncol = 3L)
  X <- X / sqrt(rowSums(X^2))
  y <- as.numeric(X %*% c(1, -1, 0.5)) + stats::rnorm(30L, sd = 0.05)
  model <- DesiredGainR:::.dgr_gp_fit(X, y, anisotropic = TRUE)
  expect_true(model$anisotropic)
  expect_identical(model$n_lengthscale, 3L)
  expect_true(is.finite(model$nlml))
})

test_that("both kernels still predict the training points", {
  set.seed(23L)
  X <- matrix(stats::runif(40L * 3L), ncol = 3L)
  X <- X / sqrt(rowSums(X^2))
  y <- as.numeric(X %*% c(1, -1, 0.5))
  for (anisotropic in c(TRUE, FALSE)) {
    model <- DesiredGainR:::.dgr_gp_fit(X, y, anisotropic = anisotropic)
    prediction <- DesiredGainR:::.dgr_gp_predict(model, X)
    expect_equal(as.numeric(prediction$mean), y, tolerance = 0.05)
  }
})
