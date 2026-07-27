test_that("run_qgsi returns canonical rankings and theory", {
  ext <- system.file("extdata", package = "DesiredGainR")
  pheno <- data.table::fread(file.path(ext, "example_pheno.csv"))
  gebv <- data.table::fread(file.path(ext, "example_gebv.csv"))
  traits <- c("YLD", "MY", "MI", "BL", "NBL", "VHB")
  w <- c(YLD = 1.5, MY = 0.5, MI = 0.5, BL = 1, NBL = 1, VHB = 1)
  W <- diag(
    c(YLD = 0.05, MY = 0.02, MI = 0.02,
      BL = -0.03, NBL = -0.03, VHB = -0.03)
  )
  dimnames(W) <- list(traits, traits)

  result <- run_qgsi(
    init_data = pheno[, .(GenoID, Family)],
    gebv_data = gebv,
    trait_cols = traits,
    linear_weights = w,
    W = W,
    lower_is_better = c("BL", "NBL", "VHB"),
    n_select = 5
  )

  expect_s3_class(result, "quadratic_genomic_index")
  expect_true(all(
    c("LinearPart", "QuadraticPart", "QGSI", "Rank", "Selected") %in%
      names(result$ranked_geno)
  ))
  expect_equal(sum(result$ranked_geno$Selected), 5)
  expect_equal(nrow(result$ranked_geno), nrow(gebv))
  expect_match(result$covariance_provenance$equation, "19.2", fixed = TRUE)
  expect_true(is.finite(
    result$theoretical_parameters$total_index_variance
  ))
  expect_false(
    result$theoretical_parameters$accuracy_and_mspe_available
  )
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
    dg_result, qgsi_result, debug = FALSE
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
    "examples", "demo_pipeline.R", package = "DesiredGainR"
  )
  expect_true(nzchar(demo))
  expect_silent(parse(file = demo))
})
