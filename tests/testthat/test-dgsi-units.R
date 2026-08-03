# Finding M1 of the review of 2026-07-31.
#
# run_dgsi() documents `dg` and `realised_response` in candidate
# standard-deviation units. The Yamada coefficient map, however, operates in
# the units of G and P. With scale_traits = TRUE those coincide. With
# scale_traits = FALSE the analysis space is raw trait units, and passing an
# SD-unit vector straight into the solve asks for d_j RAW units per trait
# rather than d_j standard deviations.
#
# The consequence is a misdirected index, not a labelling nicety, and it scales
# with how far apart the trait units are.

units_fixture <- function(seed = 17L, n = 120L, scale_factor = 10) {
  set.seed(seed)
  traits <- c("big", "small")
  # "big" is measured on a scale ten times coarser than "small" but the two
  # are otherwise identical, so any correct method must treat them alike.
  G <- diag(c(scale_factor^2, 1)) * 0.5
  dimnames(G) <- list(traits, traits)
  P <- diag(c(scale_factor^2, 1))
  dimnames(P) <- list(traits, traits)
  values <- data.table::data.table(
    GenoID = sprintf("g%03d", seq_len(n)),
    big = stats::rnorm(n, sd = scale_factor),
    small = stats::rnorm(n, sd = 1)
  )
  list(traits = traits, G = G, P = P, values = values, factor = scale_factor)
}

run_unscaled <- function(fixture, dg, seed = 4L) {
  run_dgsi(
    init_data = fixture$values[, "GenoID"],
    cand_data = fixture$values,
    trait_cols = fixture$traits,
    dg = dg, G = fixture$G, P = fixture$P,
    scale_traits = FALSE,
    n_select = 20L, n_iter = 1L, n_rep = 1L, seed = seed
  )
}

test_that("equal desired gains give equal standardised response", {
  fixture <- units_fixture()
  dg <- c(big = 1, small = 1)
  # n_iter = 1 keeps the stochastic search from compensating, so this tests
  # the coefficient map itself rather than the optimiser's ability to recover
  # from a bad starting parameterisation.
  fit <- suppressWarnings(run_unscaled(fixture, dg))

  # The classical Yamada solution at the supplied desired gains is the search
  # starting point, and is the quantity the units bug corrupted.
  baseline <- fit$non_iterated$realised_response
  expect_equal(
    unname(baseline[["big"]]), unname(baseline[["small"]]),
    tolerance = 0.35
  )
})

test_that("the index does not depend on the units traits are measured in", {
  # The same population, with one trait rescaled. Selected candidates and
  # standardised response must be unchanged: a unit change carries no genetic
  # information. Before the fix the coarse-scaled trait received roughly a
  # tenth of its intended share of the response.
  coarse <- units_fixture(scale_factor = 10)
  fine <- units_fixture(scale_factor = 1)
  dg <- c(big = 1, small = 0.5)

  coarse_fit <- suppressWarnings(run_unscaled(coarse, dg))
  fine_fit <- suppressWarnings(run_unscaled(fine, dg))

  expect_equal(
    coarse_fit$non_iterated$realised_response,
    fine_fit$non_iterated$realised_response,
    tolerance = 1e-6
  )
  # The selection flag in ranked_geno is `Selected`; `Chosen` is the column of
  # replicate_diagnostics that marks the winning search.
  expect_identical(
    sort(coarse_fit$ranked_geno[Selected == TRUE, GenoID]),
    sort(fine_fit$ranked_geno[Selected == TRUE, GenoID])
  )
})

test_that("scaling the traits gives the same answer as not scaling them", {
  # scale_traits only changes the space the solve happens in, so with the
  # units correctly converted the two routes must agree on standardised
  # response. This is the invariance the previous code lacked.
  fixture <- units_fixture()
  dg <- c(big = 1, small = 0.5)
  unscaled <- suppressWarnings(run_unscaled(fixture, dg))
  scaled <- suppressWarnings(run_dgsi(
    init_data = fixture$values[, "GenoID"],
    cand_data = fixture$values,
    trait_cols = fixture$traits,
    dg = dg, G = fixture$G, P = fixture$P,
    scale_traits = TRUE,
    n_select = 20L, n_iter = 1L, n_rep = 1L, seed = 4L
  ))
  expect_equal(
    unscaled$non_iterated$realised_response,
    scaled$non_iterated$realised_response,
    tolerance = 0.05
  )
})

test_that("the theoretical transmitted response is reported and coherent", {
  # Finding M9. The empirical differential among selected candidates and the
  # expected inherited response are different quantities and must not be
  # confused; both are now reported.
  fixture <- units_fixture()
  fit <- suppressWarnings(run_unscaled(fixture, c(big = 1, small = 0.5)))
  theoretical <- fit$theoretical_response
  expect_false(is.null(theoretical))
  expect_named(
    theoretical$analysis_units, fixture$traits,
    ignore.order = FALSE
  )
  expect_true(all(is.finite(theoretical$analysis_units)))
  expect_identical(theoretical$equation, "i * G b / sqrt(b' P b)")

  # Selecting 20 of 120 implies a standardised intensity near 1.4.
  expect_equal(
    theoretical$selection_intensity,
    DesiredGainR:::.dgr_intensity(20 / 120),
    tolerance = 1e-8
  )
  # It must point the same way as the requested gains: both traits favourable.
  expect_true(all(theoretical$standardised > 0))
})
