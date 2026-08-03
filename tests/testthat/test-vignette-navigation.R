test_that("every breeder reading-table destination exists and is exposed", {
  expected <- c(
    "DesiredGainR-workflow",
    "DesiredGainR-objective",
    "DesiredGainR-desired-gain-intervals",
    "DesiredGainR-covariance",
    "DesiredGainR-index-families",
    "DesiredGainR-dgsi",
    "DesiredGainR-qgsi",
    "DesiredGainR-simulation",
    "DesiredGainR-empirical-validation",
    "DesiredGainR-interoperation"
  )

  root <- testthat::test_path("..", "..")
  repository_files <- file.path(
    root,
    c(
      "vignettes/DesiredGainR-introduction.Rmd",
      "_pkgdown.yml",
      "README.md"
    )
  )
  testthat::skip_if_not(
    all(file.exists(repository_files)),
    "Repository documentation sources are not installed by R CMD build"
  )
  introduction <- readLines(
    file.path(root, "vignettes", "DesiredGainR-introduction.Rmd"),
    warn = FALSE
  )
  pkgdown <- readLines(file.path(root, "_pkgdown.yml"), warn = FALSE)
  readme <- readLines(file.path(root, "README.md"), warn = FALSE)

  for (article in expected) {
    expect_true(
      file.exists(file.path(root, "vignettes", paste0(article, ".Rmd"))),
      info = paste("Missing vignette source:", article)
    )
    expect_true(
      any(grepl(paste0(article, ".html"), introduction, fixed = TRUE)),
      info = paste("Reading table does not link to:", article)
    )
    expect_true(
      any(grepl(article, pkgdown, fixed = TRUE)),
      info = paste("pkgdown does not expose:", article)
    )
    expect_true(
      any(grepl(paste0(article, ".html"), readme, fixed = TRUE)),
      info = paste("README documentation table does not link to:", article)
    )
  }
})

test_that("DGSI retains trait names used by its print contract", {
  set.seed(104)
  traits <- c("yield", "disease")
  candidates <- data.frame(
    GenoID = paste0("G", seq_len(32)),
    yield = stats::rnorm(32),
    disease = stats::rnorm(32)
  )
  G <- matrix(
    c(0.7, 0.1, 0.1, 0.5), 2,
    dimnames = list(traits, traits)
  )
  P <- matrix(
    c(1.2, 0.1, 0.1, 0.9), 2,
    dimnames = list(traits, traits)
  )

  fit <- run_dgsi(
    init_data = candidates["GenoID"],
    cand_data = candidates,
    trait_cols = traits,
    dg = c(yield = 0.8, disease = 0.4),
    G = G,
    P = P,
    lower_is_better = "disease",
    n_select = 8L,
    n_iter = 5L,
    n_rep = 1L,
    replicate_selection = "training",
    seed = 104L
  )

  expect_identical(fit$trait_cols, traits)
  expect_match(paste(capture.output(print(fit)), collapse = "\n"), "Traits: 2")
})

test_that("vignettes use MathJax-safe vector notation", {
  root <- testthat::test_path("..", "..")
  vignette_dir <- file.path(root, "vignettes")
  testthat::skip_if_not(
    dir.exists(vignette_dir),
    "Repository vignette sources are not installed by R CMD build"
  )
  vignette_files <- list.files(
    vignette_dir,
    pattern = "[.]Rmd$",
    full.names = TRUE
  )
  expect_gt(length(vignette_files), 0L)
  source <- paste(vapply(
    vignette_files,
    function(path) paste(readLines(path, warn = FALSE), collapse = "\n"),
    character(1L)
  ), collapse = "\n")

  expect_false(grepl("\\\\boldsymbol", source))
  expect_false(grepl("\\\\bm[{]", source))
})
