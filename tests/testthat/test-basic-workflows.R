test_that("run_qgsi returns canonical rankings and theory", {
  skip_if_not(
    !inherits(
      try(utils::data("dgr_gebv",
        package = "DesiredGainR",
        envir = environment()
      ), silent = TRUE),
      "try-error"
    ),
    "Example data not yet generated; run data-raw/generate_example_data.R."
  )
  utils::data("dgr_traits", package = "DesiredGainR", envir = environment())
  utils::data("dgr_candidates",
    package = "DesiredGainR",
    envir = environment()
  )

  traits <- dgr_traits$trait
  lower <- dgr_traits$trait[dgr_traits$direction == "decrease"]
  w <- stats::setNames(
    c(1.5, 0.4, 0.6, 0.5, 0.4, 0.6), traits
  )
  W <- diag(c(0.05, -0.03, -0.03, -0.04, 0.02, -0.03))
  dimnames(W) <- list(traits, traits)

  result <- run_qgsi(
    init_data = dgr_candidates[, c("GenoID", "Family")],
    gebv_data = dgr_gebv,
    trait_cols = traits,
    linear_weights = w,
    W = W,
    lower_is_better = lower,
    scale_traits = TRUE,
    n_select = 20
  )

  expect_s3_class(result, "quadratic_genomic_index")
  expect_true(all(
    c("LinearPart", "QuadraticPart", "QGSI", "Rank", "Selected") %in%
      names(result$ranked_geno)
  ))
  expect_equal(sum(result$ranked_geno$Selected), 20)
  expect_equal(nrow(result$ranked_geno), nrow(dgr_gebv))
  expect_match(result$covariance_provenance$equation, "19.2", fixed = TRUE)
  expect_true(is.finite(
    result$theoretical_parameters$total_index_variance
  ))
  expect_false(
    result$theoretical_parameters$accuracy_and_mspe_available
  )
})

test_that("the example data exercise the standardisation warning", {
  skip_if_not(
    !inherits(
      try(utils::data("dgr_candidates",
        package = "DesiredGainR",
        envir = environment()
      ), silent = TRUE),
      "try-error"
    ),
    "Example data not yet generated; run data-raw/generate_example_data.R."
  )
  utils::data("dgr_traits", package = "DesiredGainR", envir = environment())
  utils::data("dgr_G", package = "DesiredGainR", envir = environment())
  utils::data("dgr_P", package = "DesiredGainR", envir = environment())

  traits <- dgr_traits$trait
  lower <- dgr_traits$trait[dgr_traits$direction == "decrease"]
  dg <- stats::setNames(rep(0.5, length(traits)), traits)

  # Trait standard deviations span a factor of about one hundred and twenty, so
  # leaving the traits unstandardised must be flagged.
  expect_warning(
    run_dgsi(
      init_data = dgr_candidates[, c("GenoID", "Family")],
      cand_data = dgr_candidates,
      trait_cols = traits,
      dg = dg, G = dgr_G, P = dgr_P,
      lower_is_better = lower,
      scale_traits = FALSE,
      n_select = 20, n_iter = 10, n_rep = 2, seed = 1
    ),
    "scale_traits = TRUE"
  )

  # Standardising, which is what the CGIAR guideline recommends, is quiet.
  expect_silent(
    run_dgsi(
      init_data = dgr_candidates[, c("GenoID", "Family")],
      cand_data = dgr_candidates,
      trait_cols = traits,
      dg = dg, G = dgr_G, P = dgr_P,
      lower_is_better = lower,
      scale_traits = TRUE,
      n_select = 20, n_iter = 10, n_rep = 2, seed = 1
    )
  )
})

test_that("run_qgsi warns when trait scales let one trait dominate", {
  skip_if_not(
    !inherits(
      try(utils::data("dgr_gebv",
        package = "DesiredGainR",
        envir = environment()
      ), silent = TRUE),
      "try-error"
    ),
    "Example data not yet generated; run data-raw/generate_example_data.R."
  )
  utils::data("dgr_traits", package = "DesiredGainR", envir = environment())

  traits <- dgr_traits$trait
  lower <- dgr_traits$trait[dgr_traits$direction == "decrease"]
  w <- stats::setNames(c(1.0, 0.2, 0.5, 0.4, 0.3, 0.5), traits)
  W <- diag(c(0.05, -0.03, -0.03, -0.04, 0.02, -0.03))
  dimnames(W) <- list(traits, traits)

  # Economic weights multiply the genomic estimated breeding values directly,
  # so on the original scale plant height overwhelms the index despite carrying
  # nearly the smallest weight.
  expect_warning(
    unscaled <- run_qgsi(
      init_data = dgr_gebv["GenoID"], gebv_data = dgr_gebv,
      trait_cols = traits, linear_weights = w, W = W,
      lower_is_better = lower, scale_traits = FALSE, n_select = 20
    ),
    "carries"
  )

  scaled <- run_qgsi(
    init_data = dgr_gebv["GenoID"], gebv_data = dgr_gebv,
    trait_cols = traits, linear_weights = w, W = W,
    lower_is_better = lower, scale_traits = TRUE, n_select = 20
  )

  # The consequence is not cosmetic: unstandardised, the expected gain in grain
  # yield is driven negative even though the index is selecting for more of it.
  unscaled_gain <- unscaled$expected_gain_per_trait$Expected_Genetic_Gain[
    unscaled$expected_gain_per_trait$Trait == "GY"
  ]
  scaled_gain <- scaled$expected_gain_per_trait$Expected_Genetic_Gain[
    scaled$expected_gain_per_trait$Trait == "GY"
  ]
  expect_lt(unscaled_gain, 0)
  expect_gt(scaled_gain, 0)
})

test_that("compare_dg_and_qgsi merges canonical outputs", {
  set.seed(8)
  traits <- c("t1", "t2")
  values <- data.table::data.table(
    GenoID = paste0("G", 1:30),
    t1 = rnorm(30),
    t2 = rnorm(30)
  )
  G <- stats::cov(as.matrix(values[, ..traits]))
  W <- diag(c(t1 = 0.05, t2 = -0.03))
  dimnames(W) <- list(traits, traits)
  dg_result <- run_dgsi(
    values[, .(GenoID)], values, traits,
    dg = c(t1 = 0.5, t2 = 0.2), G = G,
    n_select = 5, n_iter = 20, n_rep = 2, seed = 8
  )
  qgsi_result <- run_qgsi(
    values[, .(GenoID)], values, traits,
    linear_weights = c(t1 = 1, t2 = 0.3), W = W,
    n_select = 5
  )

  comparison <- compare_dg_and_qgsi(
    dg_result, qgsi_result,
    debug = FALSE
  )
  expect_true(all(
    c("DG_SelectionIndex", "DG_Rank", "QGSI", "QGSI_Rank") %in%
      names(comparison$comparison_table)
  ))
  expect_equal(nrow(comparison$comparison_table), nrow(values))
  expect_match(comparison$interpretation, "descriptive", fixed = TRUE)
})

test_that("combined pipeline never reuses dg as QGSI weights", {
  set.seed(9)
  traits <- c("t1", "t2")
  values <- data.table::data.table(
    GenoID = paste0("G", 1:24),
    t1 = rnorm(24),
    t2 = rnorm(24)
  )
  G <- stats::cov(as.matrix(values[, ..traits]))
  W <- diag(c(t1 = 0.05, t2 = -0.03))
  dimnames(W) <- list(traits, traits)

  expect_error(
    run_dgsi_qgsi_pipeline(
      mode = "qgsi",
      init_data = values[, .(GenoID)],
      trait_cols = traits,
      gebv_data = values,
      W = W,
      debug = FALSE
    ),
    "qgsi_linear_weights is required"
  )

  result <- run_dgsi_qgsi_pipeline(
    mode = "both",
    init_data = values[, .(GenoID)],
    trait_cols = traits,
    dg = c(t1 = 0.5, t2 = 0.2),
    cand_data = values,
    G = G,
    n_select = 5,
    n_iter = 20,
    n_rep = 2,
    gebv_data = values,
    qgsi_linear_weights = c(t1 = 1, t2 = 0.3),
    W = W,
    qgsi_n_select = 5,
    debug = FALSE
  )
  expect_true(all(
    c("dg_result", "qgsi_result", "comparison_result") %in% names(result)
  ))
})

test_that("installed demonstration script remains syntactically valid", {
  demo <- system.file(
    "examples", "demo_pipeline.R",
    package = "DesiredGainR"
  )
  expect_true(nzchar(demo))
  expect_silent(parse(file = demo))
})
