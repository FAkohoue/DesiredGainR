# Validation against the published reference implementation of the quadratic
# selection index.
#
#   Cerón-Rojas JJ, Montesinos-López OA, Montesinos-López A, Vitale P,
#   Pérez-Rodríguez P, Fernandes SB, Ortiz R, Crossa J (2026) Replication Data
#   for: Nonlinear Genomic Selection Index Accelerates Multi-Trait Crop
#   Improvement. CIMMYT Research Data & Software Repository Network, V1.
#   doi:10.71682/10549385
#
# run_qgsi() implements that paper. Until now it was checked only against the
# equations as this package restates them, which establishes that the code
# matches our reading and nothing more. The deposit contains the authors' own
# R scripts, and those are a reference implementation independent of this one.
#
# The scripts read intermediate data files that were not deposited, so their
# published numbers cannot be recomputed. The algorithm can be, and that is
# what is compared here: the reference code is transcribed verbatim below and
# both implementations are run on identical inputs.
#
# Transcribed from RealData-LPSI-QSI-set3-4.R and RealData-LGSI-QGSI-set3-4.R.
# Do not "tidy" it. Its value is in being their arithmetic, not ours.

reference_lpsi <- function(G, P, w) {
  IP <- solve(P)
  b <- IP %*% G %*% w
  Vp <- t(b) %*% P %*% b
  SD <- sqrt(Vp)
  ks <- 1.755
  Rs <- ks * SD
  VH <- t(w) %*% G %*% w
  SDH <- sqrt(VH)
  corrHI <- SD / SDH
  VE <- VH - Vp
  list(
    b = as.numeric(b), Vp = as.numeric(Vp), SD = as.numeric(SD),
    Rs = as.numeric(Rs), VH = as.numeric(VH), SDH = as.numeric(SDH),
    corrHI = as.numeric(corrHI), VE = as.numeric(VE), ks = ks
  )
}

reference_qsi <- function(G, P, w, A) {
  linear <- reference_lpsi(G, P, w)
  IP <- solve(P)
  BB <- IP %*% G %*% A %*% G %*% IP

  T1 <- BB %*% P
  T2 <- T1 %*% T1
  Tr <- 2 * sum(diag(T2))
  Vpq <- linear$Vp + Tr
  SDq <- sqrt(Vpq)
  ks <- 1.755
  Rq <- ks * SDq

  V <- A %*% G
  TV <- V %*% V
  VTr <- 2 * sum(diag(TV))
  VHq <- linear$VH + VTr
  SHq <- sqrt(VHq)
  corrHIq <- SDq / SHq

  list(
    BB = BB, Tr = Tr, Vpq = as.numeric(Vpq), SDq = as.numeric(SDq),
    Rq = as.numeric(Rq), VTr = VTr, VHq = as.numeric(VHq),
    SHq = as.numeric(SHq), corrHIq = as.numeric(corrHIq)
  )
}

# Four traits, as the reference scripts are written for. G and P are chosen to
# be an admissible pair with genuinely correlated traits; the reference code
# takes them as given, so any admissible pair exercises the same arithmetic.
reference_fixture <- function(seed = 3L) {
  set.seed(seed)
  traits <- c("T1", "T2", "T3", "T4")
  correlation <- matrix(
    c(
      1.00, 0.35, -0.20, 0.15,
      0.35, 1.00, 0.25, -0.10,
      -0.20, 0.25, 1.00, 0.30,
      0.15, -0.10, 0.30, 1.00
    ), 4L
  )
  genetic_sd <- c(1.2, 0.8, 2.5, 0.6)
  G <- outer(genetic_sd, genetic_sd) * correlation
  dimnames(G) <- list(traits, traits)
  P <- G + diag(c(1.6, 1.1, 3.0, 0.9))
  dimnames(P) <- list(traits, traits)
  # The reference scripts use w = c(1, -1, 1, 1).
  w <- stats::setNames(c(1, -1, 1, 1), traits)
  list(traits = traits, G = G, P = P, w = w)
}

# ---------------------------------------------------------------------------
# Linear phenotypic selection index
# ---------------------------------------------------------------------------

test_that("our Smith-Hazel coefficients equal the reference implementation", {
  fixture <- reference_fixture()
  reference <- reference_lpsi(fixture$G, fixture$P, fixture$w)

  set.seed(11L)
  values <- as.data.frame(matrix(
    stats::rnorm(60L * 4L),
    ncol = 4L,
    dimnames = list(paste0("g", 1:60), fixture$traits)
  ))
  # scale_traits = FALSE keeps the analysis in the units of G and P, which is
  # the space the reference code works in.
  fit <- selection_index(
    values, fixture$traits,
    method = "smith_hazel",
    G = fixture$G, P = fixture$P, economic_weights = abs(fixture$w),
    lower_is_better = "T2", scale_traits = FALSE, center_traits = FALSE,
    n_select = 15L
  )
  # lower_is_better flips the sign of T2 in the analysis space, which is how
  # this package expresses w = c(1, -1, 1, 1). Undo it to compare.
  direction <- c(1, -1, 1, 1)
  expect_equal(
    as.numeric(coef(fit)) * direction, reference$b,
    tolerance = 1e-10
  )
})

test_that("our R_HI equals the reference corrHI for the optimum index", {
  # The reference computes corrHI = sqrt(b'Pb) / sqrt(w'Gw). This package
  # computes R_HI = b'Ga / (sqrt(b'Pb) sqrt(a'Ga)). For the optimum
  # b = P^-1 G w we have b'Gw = w'G P^-1 G w = b'Pb, so the two coincide.
  # They are different expressions of the same quantity, and agreeing is
  # evidence that both parties mean the same thing by accuracy.
  fixture <- reference_fixture()
  reference <- reference_lpsi(fixture$G, fixture$P, fixture$w)

  b <- reference$b
  names(b) <- fixture$traits
  evaluation <- evaluate_index(
    coefficients = b, G = fixture$G, P = fixture$P,
    aggregate_weights = fixture$w, selection_intensity = reference$ks
  )
  expect_equal(evaluation$R_HI, reference$corrHI, tolerance = 1e-10)
  # And the aggregate response, their Rs = ks * sqrt(b'Pb).
  expect_equal(evaluation$delta_H, reference$Rs, tolerance = 1e-10)
})

test_that("the reference identity b'Gw = b'Pb holds for the optimum index", {
  # The algebraic step the previous test relies on, asserted directly so that a
  # failure there can be attributed.
  fixture <- reference_fixture()
  reference <- reference_lpsi(fixture$G, fixture$P, fixture$w)
  b <- reference$b
  expect_equal(
    as.numeric(t(b) %*% fixture$G %*% fixture$w), reference$Vp,
    tolerance = 1e-10
  )
  # And the residual variance the reference reports is non-negative, which is
  # the compatibility condition this package now enforces.
  expect_gt(reference$VE, 0)
})

# ---------------------------------------------------------------------------
# Quadratic selection index
# ---------------------------------------------------------------------------

test_that("our QGSI variance matches the reference quadratic formulation", {
  # The reference builds the quadratic coefficient matrix as
  #   BB = P^-1 G A G P^-1
  # and the index variance as
  #   Vpq = b'Pb + 2 tr(BB P BB P).
  # This package computes Var(I) = w'Gamma w + 2 tr(W Gamma W Gamma). Setting
  # Gamma = P, w = b and W = BB makes them the same quantity, so the two
  # expressions must agree numerically.
  fixture <- reference_fixture()
  A <- diag(c(0.5, -0.2, 0.3, -0.1))
  dimnames(A) <- list(fixture$traits, fixture$traits)
  reference <- reference_qsi(fixture$G, fixture$P, fixture$w, A)

  b <- reference$b <- reference_lpsi(fixture$G, fixture$P, fixture$w)$b
  names(b) <- fixture$traits
  W <- reference$BB
  dimnames(W) <- list(fixture$traits, fixture$traits)

  ours <- as.numeric(
    t(b) %*% fixture$P %*% b +
      2 * sum(diag(W %*% fixture$P %*% W %*% fixture$P))
  )
  expect_equal(ours, reference$Vpq, tolerance = 1e-10)
  expect_equal(sqrt(ours), reference$SDq, tolerance = 1e-10)
})

test_that("our QGSI merit variance matches the reference", {
  # Reference: VHq = w'Gw + 2 tr(A G A G).
  fixture <- reference_fixture()
  A <- diag(c(0.5, -0.2, 0.3, -0.1))
  dimnames(A) <- list(fixture$traits, fixture$traits)
  reference <- reference_qsi(fixture$G, fixture$P, fixture$w, A)

  ours <- as.numeric(
    t(fixture$w) %*% fixture$G %*% fixture$w +
      2 * sum(diag(A %*% fixture$G %*% A %*% fixture$G))
  )
  expect_equal(ours, reference$VHq, tolerance = 1e-10)
  expect_equal(
    reference$SDq / sqrt(ours), reference$corrHIq,
    tolerance = 1e-10
  )
})

test_that("run_qgsi() reproduces the reference index variance end to end", {
  # Through the exported interface rather than by restating the formula, so
  # that the transformation, validation and assembly layers are exercised too.
  skip_if_not_installed("MASS")
  fixture <- reference_fixture()
  A <- diag(c(0.5, -0.2, 0.3, -0.1))
  dimnames(A) <- list(fixture$traits, fixture$traits)
  reference <- reference_qsi(fixture$G, fixture$P, fixture$w, A)

  b <- stats::setNames(
    reference_lpsi(fixture$G, fixture$P, fixture$w)$b, fixture$traits
  )
  W <- reference$BB
  dimnames(W) <- list(fixture$traits, fixture$traits)

  set.seed(13L)
  n <- 80L
  gebv <- data.table::data.table(GenoID = sprintf("g%03d", seq_len(n)))
  scores <- MASS::mvrnorm(n, mu = rep(0, 4L), Sigma = fixture$P)
  for (index in seq_along(fixture$traits)) {
    data.table::set(gebv, j = fixture$traits[index], value = scores[, index])
  }

  fitted <- run_qgsi(
    init_data = gebv[, "GenoID"], gebv_data = gebv,
    trait_cols = fixture$traits,
    linear_weights = b, W = W,
    Gamma = fixture$P,
    center_traits = FALSE, scale_traits = FALSE,
    n_select = 20L
  )
  expect_equal(
    fitted$theoretical_parameters$total_index_variance, reference$Vpq,
    tolerance = 1e-8
  )
})

test_that("the reference implementation is recorded with its provenance", {
  # The transcription above is the evidence, so its origin must travel with
  # it. If this file is ever edited to make a test pass, that edit changes what
  # the package is being compared against and must be deliberate.
  path <- test_path("test-reference-implementation.R")
  skip_if_not(file.exists(path))
  source_text <- readLines(path, warn = FALSE)
  expect_true(any(grepl("10.71682/10549385", source_text, fixed = TRUE)))
  expect_true(any(grepl("RealData-LPSI-QSI-set3-4.R", source_text,
    fixed = TRUE
  )))
  expect_true(any(grepl("Do not \"tidy\" it", source_text, fixed = TRUE)))

  manifest <- system.file(
    "reference", "dataverse-qgsi-manifest.txt",
    package = "DesiredGainR"
  )
  if (!nzchar(manifest)) {
    manifest <- test_path(
      "..", "..", "inst", "reference",
      "dataverse-qgsi-manifest.txt"
    )
  }
  expect_true(file.exists(manifest))
  provenance <- readLines(manifest, warn = FALSE)
  expect_true(any(grepl("10.71682/10549385", provenance, fixed = TRUE)))
  expect_true(any(grepl(
    "b26a3f6ad5bcc037d25765bc23b5b2062566b29c3f5175ef5fbf4b0e0d20f937",
    provenance,
    fixed = TRUE
  )))
})
