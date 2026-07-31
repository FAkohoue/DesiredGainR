# The example data exist to demonstrate specific properties. If they stop
# having those properties the demonstrations become misleading, so the
# properties are asserted rather than assumed.
#
# These tests run only once data-raw/generate_example_data.R has been executed.

skip_if_no_data <- function() {
  skip_if_not(
    exists("dgr_G", where = asNamespace("DesiredGainR"), inherits = FALSE) ||
      requireNamespace("DesiredGainR", quietly = TRUE) &&
        !inherits(try(utils::data("dgr_G", package = "DesiredGainR",
                                  envir = environment()), silent = TRUE),
                  "try-error"),
    "Example data not yet generated; run data-raw/generate_example_data.R."
  )
}

test_that("the covariance matrices are consistent with the declared traits", {
  skip_if_no_data()
  traits <- dgr_traits$trait

  expect_equal(colnames(dgr_G), traits)
  expect_equal(colnames(dgr_P), traits)
  expect_true(isSymmetric(unname(dgr_G)))
  expect_true(isSymmetric(unname(dgr_P)))

  # Heritability must be recoverable from the two matrices, or the metadata and
  # the matrices are telling different stories.
  expect_equal(
    as.numeric(diag(dgr_G) / diag(dgr_P)),
    dgr_traits$heritability,
    tolerance = 1e-8
  )
  # The genetic standard deviations must match the declared ones.
  expect_equal(
    as.numeric(sqrt(diag(dgr_G))), dgr_traits$genetic_sd, tolerance = 1e-8
  )

  # Both must be invertible, since every index inverts at least one of them.
  expect_true(matrix_diagnostics(dgr_G, "G")$positive_definite)
  expect_true(matrix_diagnostics(dgr_P, "P")$positive_definite)
})

test_that("the data carry the antagonism the examples depend on", {
  skip_if_no_data()
  correlation <- stats::cov2cor(dgr_G)

  # Grain yield and anthesis date correlate positively while the objective
  # requires yield up and anthesis date down. Without that tension the
  # feasibility and Pareto demonstrations have nothing to show.
  expect_gt(correlation["GY", "AD"], 0.1)
  expect_equal(dgr_traits$direction[dgr_traits$trait == "GY"], "increase")
  expect_equal(dgr_traits$direction[dgr_traits$trait == "AD"], "decrease")

  # Yield and the anthesis-silking interval are strongly antagonistic in the
  # favourable-direction sense.
  expect_lt(correlation["GY", "ASI"], -0.4)

  # Off-diagonal structure must be real, not decorative.
  off_diagonal <- correlation[upper.tri(correlation)]
  expect_gt(mean(abs(off_diagonal)), 0.15)
})

test_that("trait scales differ enough to exercise the diagnostics", {
  skip_if_no_data()
  # A factor of five triggers the run_dgsi() standardisation warning; the real
  # spread here is far larger, which is the point.
  ratio <- max(dgr_traits$genetic_sd) / min(dgr_traits$genetic_sd)
  expect_gt(ratio, 50)
})

test_that("candidates realise the declared genetic covariance", {
  skip_if_no_data()
  traits <- dgr_traits$trait
  values <- as.matrix(dgr_candidates[, traits, drop = FALSE])

  expect_equal(nrow(dgr_candidates), 200L)
  expect_false(anyDuplicated(dgr_candidates$GenoID) > 0L)
  expect_equal(length(unique(dgr_candidates$Family)), 20L)

  # The phenotypic covariance of the adjusted means should approach dgr_P.
  # Sampling error at two hundred candidates is appreciable, so the tolerance
  # is on the correlation structure rather than the raw covariances.
  observed <- stats::cov2cor(stats::cov(values))
  expected <- stats::cov2cor(dgr_P)
  expect_lt(max(abs(observed - expected)), 0.20)

  # Trait means should sit near the declared programme means.
  expect_equal(
    as.numeric(colMeans(values)), dgr_traits$mean,
    tolerance = 0.05
  )
})

test_that("genomic estimates are shrunken relative to genetic values", {
  skip_if_no_data()
  traits <- dgr_traits$trait
  gebv <- as.matrix(dgr_gebv[, traits, drop = FALSE])

  expect_equal(dgr_gebv$GenoID, dgr_candidates$GenoID)
  # A prediction retains only the reliable fraction of the genetic variance,
  # so its variance must fall below the genetic variance for every trait.
  expect_true(all(diag(stats::cov(gebv)) < diag(dgr_G)))
})

test_that("haplotypes satisfy the founder contract", {
  skip_if_no_data()
  expect_true(all(dgr_hap1 %in% c(0L, 1L)))
  expect_true(all(dgr_hap2 %in% c(0L, 1L)))
  expect_identical(dimnames(dgr_hap1), dimnames(dgr_hap2))
  expect_equal(colnames(dgr_hap1), dgr_candidates$GenoID)
  expect_equal(rownames(dgr_hap1), dgr_map$variant_id)

  # Every marker must be polymorphic, or it carries no information.
  frequency <- rowMeans(dgr_hap1 + dgr_hap2) / 2
  expect_true(all(frequency > 0 & frequency < 1))

  # The map column names must match the founder_haplotypes() defaults, so that
  # it can be passed without further argument.
  expect_true(all(
    c("variant_id", "chromosome", "position_bp") %in% names(dgr_map)
  ))
  founders <- founder_haplotypes(dgr_hap1, dgr_hap2, dgr_map)
  expect_equal(founders$n_individuals, 200L)
  expect_equal(founders$n_chromosomes, 3L)
})

test_that("markers within a block are in linkage disequilibrium", {
  skip_if_no_data()
  dosage <- t(dgr_hap1 + dgr_hap2)
  # Adjacent markers on the same chromosome should be far more correlated than
  # markers on different chromosomes, or the simulation has no useful linkage
  # structure for a recombination model to act on.
  same_chromosome <- dgr_map$chromosome[-1L] == dgr_map$chromosome[-nrow(dgr_map)]
  adjacent <- vapply(which(same_chromosome), function(i) {
    abs(stats::cor(dosage[, i], dosage[, i + 1L]))
  }, numeric(1))

  set.seed(1)
  distant <- replicate(200L, {
    pair <- sample.int(ncol(dosage), 2L)
    if (dgr_map$chromosome[pair[1L]] == dgr_map$chromosome[pair[2L]]) {
      return(NA_real_)
    }
    abs(stats::cor(dosage[, pair[1L]], dosage[, pair[2L]]))
  })
  expect_gt(mean(adjacent, na.rm = TRUE), mean(distant, na.rm = TRUE))
})

test_that("the historical decision supports weight recovery", {
  skip_if_no_data()
  traits <- dgr_traits$trait
  expect_equal(sum(dgr_history$selected), 40L)
  expect_equal(dgr_history$GenoID, dgr_candidates$GenoID)

  recovered <- retrospective_weights(
    selected_values = dgr_candidates[dgr_history$selected, traits],
    population_values = dgr_candidates[, traits],
    trait_cols = traits
  )
  expect_equal(recovered$n_selected, 40L)
  expect_equal(recovered$selected_proportion, 0.2)
  # The recovered differentials must point the right way once oriented: an
  # improving decision raises grain yield and lowers anthesis-silking interval.
  expect_gt(recovered$selection_differential[["GY"]], 0)
  expect_lt(recovered$selection_differential[["ASI"]], 0)
})
