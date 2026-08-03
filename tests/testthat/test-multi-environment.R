# Item Q3: multi-environment and genotype-by-environment structure.

me_fixture <- function() {
  traits <- c("yield", "protein")
  environments <- c("irrigated", "rainfed", "lowland")
  G <- matrix(c(1.0, 0.2, 0.2, 0.5), 2L, dimnames = list(traits, traits))
  P <- G + diag(c(1.5, 0.9))
  dimnames(P) <- list(traits, traits)
  list(traits = traits, environments = environments, G = G, P = P)
}

test_that("the expansion has the separable Kronecker structure", {
  fixture <- me_fixture()
  expanded <- expand_environments(
    fixture$traits, fixture$environments, fixture$G, fixture$P,
    environment_correlation = 0.6
  )
  p <- length(fixture$traits)
  e <- length(fixture$environments)
  expect_identical(dim(expanded$G), c(p * e, p * e))
  expect_identical(expanded$labels[1L], "yield_irrigated")

  # The defining property: Cov(g_je, g_kf) = G_jk * C_ef.
  expect_equal(
    expanded$G[["yield_irrigated", "yield_rainfed"]],
    fixture$G[["yield", "yield"]] * 0.6
  )
  expect_equal(
    expanded$G[["yield_irrigated", "protein_rainfed"]],
    fixture$G[["yield", "protein"]] * 0.6
  )
  # Same environment: the environmental correlation is 1.
  expect_equal(
    expanded$G[["yield_irrigated", "protein_irrigated"]],
    fixture$G[["yield", "protein"]]
  )
})

test_that("the expanded P and G remain a compatible pair", {
  # Expanding P with the genetic environmental correlation would not guarantee
  # this, which is why the residual is expanded separately and P reassembled.
  fixture <- me_fixture()
  for (correlation in c(0, 0.3, 0.9)) {
    expanded <- expand_environments(
      fixture$traits, fixture$environments, fixture$G, fixture$P,
      environment_correlation = correlation
    )
    smallest <- min(eigen(
      expanded$P - expanded$G,
      symmetric = TRUE, only.values = TRUE
    )$values)
    expect_gte(smallest, -1e-8)
  }
})

test_that("a multi-environment problem feeds the ordinary index layer", {
  fixture <- me_fixture()
  expanded <- expand_environments(
    fixture$traits, fixture$environments, fixture$G, fixture$P,
    environment_correlation = 0.6,
    economic_weights = c(yield = 2, protein = 1),
    environment_weights = c(irrigated = 0.5, rainfed = 0.3, lowland = 0.2)
  )
  set.seed(5L)
  values <- as.data.frame(matrix(
    stats::rnorm(40L * length(expanded$labels)),
    ncol = length(expanded$labels),
    dimnames = list(paste0("g", 1:40), expanded$labels)
  ))
  fit <- selection_index(
    values, expanded$labels,
    method = "smith_hazel",
    G = expanded$G, P = expanded$P,
    economic_weights = expanded$economic_weights, n_select = 8L
  )
  expect_length(coef(fit), length(expanded$labels))
  expect_lte(fit$evaluation$R_HI, 1)
  expect_gt(fit$evaluation$R_HI, 0)
})

test_that("environment weights are normalised and applied per trait", {
  fixture <- me_fixture()
  expanded <- expand_environments(
    fixture$traits, fixture$environments, fixture$G, fixture$P,
    environment_correlation = 0.5,
    economic_weights = c(yield = 2, protein = 1),
    environment_weights = c(irrigated = 6, rainfed = 3, lowland = 1)
  )
  expect_equal(sum(expanded$environment_weights), 1)
  expect_equal(
    expanded$economic_weights[["yield_irrigated"]], 2 * 0.6,
    tolerance = 1e-12
  )
  expect_equal(
    expanded$economic_weights[["protein_lowland"]], 1 * 0.1,
    tolerance = 1e-12
  )
  # Total weight is preserved across the expansion.
  expect_equal(sum(expanded$economic_weights), 3, tolerance = 1e-12)
})

test_that("the interaction share reflects the environmental correlation", {
  fixture <- me_fixture()
  high <- expand_environments(
    fixture$traits, fixture$environments, fixture$G, fixture$P,
    environment_correlation = 0.95
  )
  low <- expand_environments(
    fixture$traits, fixture$environments, fixture$G, fixture$P,
    environment_correlation = 0.1
  )
  expect_lt(high$interaction_share, low$interaction_share)
  expect_match(high$structure, "almost identically")
  expect_match(low$structure, "very differently")
})

test_that("stability contrasts measure departure from the trait mean", {
  fixture <- me_fixture()
  expanded <- expand_environments(
    fixture$traits, fixture$environments, fixture$G, fixture$P,
    environment_correlation = 0.4, include_stability = TRUE
  )
  contrast <- expanded$stability$contrast
  # One contrast per trait per environment, less one per trait for the
  # rank deficiency.
  expect_identical(nrow(contrast), 2L * (3L - 1L))
  # Each row sums to zero: it is a deviation, so a genotype uniformly above
  # its own mean everywhere has zero stability score.
  expect_equal(rowSums(contrast), rep(0, nrow(contrast)), tolerance = 1e-12)
  # And the contrast covariance is a valid covariance matrix.
  expect_gte(
    min(eigen(expanded$stability$G,
      symmetric = TRUE,
      only.values = TRUE
    )$values),
    -1e-8
  )
})

test_that("a perfectly correlated environment set has no stability variance", {
  # With an environmental correlation of exactly 1 there is no interaction, so
  # every genotype performs at its own mean everywhere and the stability
  # contrasts must carry no genetic variance at all.
  fixture <- me_fixture()
  expanded <- expand_environments(
    fixture$traits, fixture$environments, fixture$G, fixture$P,
    environment_correlation = 1, include_stability = TRUE
  )
  expect_lt(max(abs(expanded$stability$G)), 1e-8)
})

test_that("expand_environments() validates its inputs", {
  fixture <- me_fixture()
  expect_error(
    expand_environments(fixture$traits, "only_one", fixture$G, fixture$P, 0.5),
    "at least two unique names"
  )
  expect_error(
    expand_environments(
      fixture$traits, fixture$environments, fixture$G, fixture$P,
      environment_correlation = 1.5
    ),
    "correlation matrix or a single value"
  )
  # A correlation matrix that is not positive semidefinite.
  bad <- matrix(c(1, 0.9, -0.9, 0.9, 1, 0.9, -0.9, 0.9, 1), 3L,
    dimnames = list(fixture$environments, fixture$environments)
  )
  expect_error(
    expand_environments(
      fixture$traits, fixture$environments, fixture$G, fixture$P,
      environment_correlation = bad
    ),
    "negative eigenvalue"
  )
  # And an inadmissible trait-level pair is caught before expansion.
  expect_error(
    expand_environments(
      fixture$traits, fixture$environments, fixture$P, fixture$G, 0.5
    ),
    "not positive semidefinite"
  )
})

test_that("widen_environments() reshapes and handles missing genotypes", {
  long <- expand.grid(
    geno = paste0("g", 1:5), env = c("irrigated", "rainfed"),
    stringsAsFactors = FALSE
  )
  set.seed(2L)
  long$yield <- stats::rnorm(nrow(long))
  long$protein <- stats::rnorm(nrow(long))

  wide <- widen_environments(long, "geno", "env", c("yield", "protein"))
  expect_identical(nrow(wide), 5L)
  expect_true("yield_irrigated" %in% names(wide))
  expect_true("protein_rainfed" %in% names(wide))

  # A genotype missing from one environment.
  incomplete <- long[-1L, ]
  expect_error(
    widen_environments(incomplete, "geno", "env", c("yield", "protein")),
    "not observed in every"
  )
  dropped <- widen_environments(
    incomplete, "geno", "env", c("yield", "protein"),
    missing_policy = "drop"
  )
  expect_identical(nrow(dropped), 4L)
  imputed <- widen_environments(
    incomplete, "geno", "env", c("yield", "protein"),
    missing_policy = "mean_impute"
  )
  expect_identical(nrow(imputed), 5L)
  expect_false(anyNA(imputed))
})

# ---------------------------------------------------------------------------
# M4: covariance draws
# ---------------------------------------------------------------------------

test_that("drawn covariance pairs are admissible by construction", {
  traits <- c("yield", "protein")
  G <- matrix(c(1.0, 0.2, 0.2, 0.5), 2L, dimnames = list(traits, traits))
  P <- matrix(c(2.5, 0.4, 0.4, 1.2), 2L, dimnames = list(traits, traits))
  set.seed(7L)
  draws <- draw_covariance_pairs(
    G, P,
    genetic_df = 30L, residual_df = 100L, n_draws = 50L
  )
  expect_length(draws, 50L)
  admissible <- vapply(draws, function(draw) {
    min(eigen(draw$P - draw$G, symmetric = TRUE, only.values = TRUE)$values)
  }, numeric(1L))
  # Every single draw, not most of them: P* = G* + E* makes this structural.
  expect_true(all(admissible > 0))
  # And the draws are centred on the supplied matrices.
  mean_G <- Reduce(`+`, lapply(draws, function(d) d$G)) / length(draws)
  expect_equal(mean_G, G, tolerance = 0.25, ignore_attr = TRUE)
})

test_that("draws without P carry no residual", {
  traits <- c("a", "b")
  G <- matrix(c(1, 0.1, 0.1, 1), 2L, dimnames = list(traits, traits))
  draws <- draw_covariance_pairs(G, genetic_df = 40L, n_draws = 3L)
  expect_null(draws[[1L]]$P)
  expect_null(draws[[1L]]$E)
  expect_identical(dim(draws[[1L]]$G), c(2L, 2L))
})
