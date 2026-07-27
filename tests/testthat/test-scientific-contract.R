test_that("DGSI automatically chooses the best replicate", {
  set.seed(19)
  traits <- c("yield", "disease")
  values <- data.table::data.table(
    GenoID = paste0("G", seq_len(40)),
    yield = rnorm(40),
    disease = rnorm(40)
  )
  metadata <- values[, .(GenoID)]
  G <- stats::cov(as.matrix(values[, ..traits]))

  result <- run_dgsi(
    init_data = metadata,
    cand_data = values,
    trait_cols = traits,
    dg = c(yield = 0.5, disease = 0.3),
    G = G,
    lower_is_better = "disease",
    n_select = 8,
    n_iter = 25,
    n_rep = 4,
    seed = 99
  )

  expect_equal(
    result$objective,
    min(result$replicate_diagnostics$Objective)
  )
  expect_equal(sum(result$replicate_diagnostics$Chosen), 1)
  expect_equal(result$best_replicate,
               result$replicate_diagnostics[Chosen == TRUE, Replicate])
  expect_equal(sum(result$ranked_geno$Selected), 8)
})

test_that("eligibility thresholds are followed by index ranking", {
  values <- data.table::data.table(
    GenoID = paste0("G", seq_len(12)),
    t1 = seq(-1.1, 1.1, length.out = 12),
    t2 = rev(seq(-1.1, 1.1, length.out = 12))
  )
  G <- stats::cov(as.matrix(values[, .(t1, t2)]))
  result <- run_dgsi(
    init_data = values[, .(GenoID)],
    cand_data = values,
    trait_cols = c("t1", "t2"),
    dg = c(t1 = 0.4, t2 = 0.4),
    G = G,
    select_mode = "eligible_top_n",
    trait_min = c(t1 = -0.5, t2 = -0.5),
    n_select = 3,
    n_iter = 10,
    n_rep = 2,
    seed = 1
  )

  expect_true(all(result$ranked_geno[Selected == TRUE, Eligible]))
  expect_equal(sum(result$ranked_geno$Selected), 3)
  expect_equal(
    result$ranked_geno[Eligible == TRUE][order(-SelectionIndex)][1:3, GenoID],
    result$ranked_geno[Selected == TRUE][order(-SelectionIndex), GenoID]
  )
})

test_that("missing values are rejected unless policy is explicit", {
  values <- data.table::data.table(
    GenoID = paste0("G", 1:6),
    t1 = c(1, 2, NA, 4, 5, 6),
    t2 = 6:1
  )
  G <- diag(2)
  dimnames(G) <- list(c("t1", "t2"), c("t1", "t2"))
  expect_error(
    run_dgsi(
      values[, .(GenoID)], values, c("t1", "t2"),
      dg = c(t1 = 0.2, t2 = 0.2), G = G,
      n_select = 2, n_iter = 2, n_rep = 1
    ),
    "Missing trait values"
  )
})

test_that("QGSI requires W and returns candidate-specific contributions", {
  values <- data.table::data.table(
    GenoID = paste0("G", 1:5),
    t1 = -2:2,
    t2 = c(2, 1, 0, -1, -2)
  )
  expect_error(
    run_qgsi(
      values[, .(GenoID)], values, c("t1", "t2"),
      linear_weights = c(t1 = 1, t2 = 1)
    ),
    "W is required"
  )

  W <- matrix(c(0.2, -0.1, -0.1, 0.3), 2)
  dimnames(W) <- list(c("t1", "t2"), c("t1", "t2"))
  result <- run_qgsi(
    values[, .(GenoID)], values, c("t1", "t2"),
    linear_weights = c(t1 = 1, t2 = 1),
    W = W
  )
  expect_s3_class(result, "quadratic_genomic_index")
  expect_equal(
    rowSums(result$quadratic_contributions[, -"GenoID"]),
    result$ranked_geno[
      match(
        result$quadratic_contributions$GenoID,
        result$ranked_geno$GenoID
      ),
      QuadraticPart
    ]
  )
})
