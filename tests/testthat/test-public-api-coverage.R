method_coverage_fixture <- function() {
  traits <- c("yield", "disease", "height")
  values <- data.frame(
    id = paste0("L", seq_len(12)),
    yield = seq(-1.1, 1.1, length.out = 12),
    disease = c(0.9, 0.1, -0.6, 0.5, -0.8, 0.2,
                -1.0, 0.7, -0.3, 0.4, -0.5, -0.9),
    height = c(0.8, -0.4, 0.2, -0.8, 0.6, -0.2,
               0.4, -0.6, 0.1, -0.1, 0.5, -0.5)
  )
  G <- matrix(
    c(
      1.00, -0.15, -0.10,
      -0.15, 0.80, 0.05,
      -0.10, 0.05, 0.60
    ),
    nrow = 3,
    dimnames = list(traits, traits)
  )
  P <- G + diag(c(0.70, 0.60, 0.50))
  dimnames(P) <- list(traits, traits)
  list(
    traits = traits,
    values = values,
    G = G,
    P = P,
    economic = c(yield = 1.0, disease = 0.4, height = 0.2),
    desired = c(yield = 1.0, disease = 0.5, height = 0.3),
    limits = c(yield = -0.5, disease = 0.7, height = 0.6),
    lower = c("disease", "height")
  )
}

fit_all_classical_families <- function(set) {
  common <- list(
    values = set$values,
    trait_cols = set$traits,
    id_col = "id",
    G = set$G,
    P = set$P,
    lower_is_better = set$lower,
    center_traits = FALSE,
    scale_traits = FALSE,
    n_select = 4L,
    main_trait = "yield"
  )
  list(
    Smith_Hazel = do.call(selection_index, c(
      common,
      list(method = "smith_hazel", economic_weights = set$economic)
    )),
    Base = do.call(selection_index, c(
      common,
      list(method = "base", economic_weights = set$economic)
    )),
    Pesek_Baker = do.call(selection_index, c(
      common,
      list(
        method = "pesek_baker",
        desired_gains = set$desired,
        aggregate_weights = set$economic
      )
    )),
    Yamada = do.call(selection_index, c(
      common,
      list(
        method = "yamada",
        desired_gains = set$desired,
        aggregate_weights = set$economic
      )
    )),
    Mulamba_Mock = do.call(selection_index, c(
      common,
      list(method = "mulamba_mock")
    )),
    Elston = do.call(selection_index, c(
      common,
      list(method = "elston", culling_thresholds = set$limits)
    )),
    Independent_culling = do.call(selection_index, c(
      common,
      list(
        method = "independent_culling",
        culling_thresholds = set$limits
      )
    )),
    Tandem = do.call(selection_index, c(
      common,
      list(method = "tandem", tandem_order = set$traits)
    ))
  )
}

test_that("all eight classical selection families execute their contracts", {
  set <- method_coverage_fixture()
  fits <- fit_all_classical_families(set)
  expected_methods <- c(
    Smith_Hazel = "smith_hazel",
    Base = "base",
    Pesek_Baker = "pesek_baker",
    Yamada = "yamada",
    Mulamba_Mock = "mulamba_mock",
    Elston = "elston",
    Independent_culling = "independent_culling",
    Tandem = "tandem"
  )

  expect_setequal(names(fits), names(expected_methods))
  for (label in names(fits)) {
    fit <- fits[[label]]
    expect_true(
      inherits(fit, "desiredgainr_index"),
      info = paste(label, "did not return a desiredgainr_index")
    )
    expect_identical(fit$method, unname(expected_methods[[label]]), info = label)
    expect_identical(fit$n_selected, 4L, info = label)
    expect_equal(nrow(fit$ranking), nrow(set$values), info = label)
    expect_equal(nrow(fit$selected), 4L, info = label)
    expect_identical(anyDuplicated(fit$selected$id), 0L, info = label)
    expect_true(all(is.finite(fit$observed_differential$Differential)),
      info = label
    )
  }

  expect_equal(
    fits$Pesek_Baker$coefficients,
    fits$Yamada$coefficients,
    tolerance = 1e-10
  )
  expect_null(fits$Mulamba_Mock$evaluation)
  expect_null(fits$Elston$evaluation)
  expect_null(fits$Independent_culling$evaluation)
  expect_null(fits$Tandem$evaluation)
  expect_true(all(fits$Elston$culling_report$passed ==
    is.finite(fits$Elston$ranking$score[match(
      fits$Elston$culling_report$id,
      fits$Elston$ranking$id
    )])))
})

test_that("tandem selection applies the declared stage order", {
  values <- data.frame(
    id = LETTERS[1:8],
    yield = 8:1,
    disease = c(8, 1, 2, 3, 0, 0, 0, 0)
  )
  fit <- selection_index(
    values,
    c("yield", "disease"),
    id_col = "id",
    method = "tandem",
    lower_is_better = "disease",
    tandem_order = c("yield", "disease"),
    center_traits = FALSE,
    scale_traits = FALSE,
    n_select = 2L
  )
  expect_setequal(fit$selected$id, c("B", "C"))
  expect_error(
    selection_index(
      values, c("yield", "disease"), id_col = "id",
      method = "tandem", tandem_order = "yield",
      center_traits = FALSE, scale_traits = FALSE
    ),
    "n_select is required"
  )
  expect_error(
    selection_index(
      values, c("yield", "disease"), id_col = "id",
      method = "tandem", tandem_order = c("yield", "unknown"),
      center_traits = FALSE, scale_traits = FALSE, n_select = 2L
    ),
    "tandem_order"
  )
})

test_that("all five restricted index families execute their constraints", {
  set <- method_coverage_fixture()
  methods <- c(
    "kempthorne_nordskog", "restricted_smith_hazel",
    "tallis", "mallard", "harville"
  )
  fits <- lapply(methods, function(method) {
    proportional <- method %in% c("tallis", "mallard", "harville")
    restricted_index(
      set$values,
      set$traits,
      id_col = "id",
      method = method,
      G = set$G,
      P = set$P,
      economic_weights = set$economic,
      restricted_traits = if (proportional) NULL else "height",
      target_gains = if (proportional) set$desired else NULL,
      lower_is_better = set$lower,
      center_traits = FALSE,
      scale_traits = FALSE,
      n_select = 4L
    )
  })
  names(fits) <- methods

  for (method in methods) {
    fit <- fits[[method]]
    expect_true(
      inherits(fit, "desiredgainr_index"),
      info = paste(method, "did not return a desiredgainr_index")
    )
    expect_identical(fit$method, method, info = method)
    expect_true(
      fit$constraint$largest_violation < 1e-8,
      info = paste(method, "did not satisfy its declared restriction")
    )
    expect_identical(fit$n_selected, 4L, info = method)
  }
  for (method in c("tallis", "mallard", "harville")) {
    expect_true(
      inherits(
        fits[[method]]$constraint$satoh_response,
        "desiredgainr_restricted_response"
      ),
      info = paste(method, "did not return a Satoh response assessment")
    )
    expect_true(
      inherits(
        fits[[method]]$restricted_breeding_values,
        "desiredgainr_restricted_bv"
      ),
      info = paste(method, "did not return restricted breeding values")
    )
  }

  unrestricted <- selection_index(
    set$values,
    set$traits,
    id_col = "id",
    method = "smith_hazel",
    G = set$G,
    P = set$P,
    economic_weights = set$economic,
    lower_is_better = set$lower,
    center_traits = FALSE,
    scale_traits = FALSE,
    n_select = 4L
  )
  comparison <- compare_selection_methods(
    c(list(Unrestricted = unrestricted), fits),
    target_gains = set$desired
  )
  expect_equal(nrow(comparison$summary), 6L)
  expect_setequal(
    comparison$summary$Family,
    c("smith_hazel", methods)
  )
})

test_that("general economic and desired-gain indices are both covered", {
  set <- method_coverage_fixture()
  information <- selection_information(
    set$values[set$traits],
    P = set$P,
    C = set$G,
    G = set$G
  )
  economic <- generalized_index(
    information,
    set$economic,
    method = "economic",
    n_select = 4L,
    lower_is_better = set$lower
  )
  desired <- generalized_index(
    information,
    set$desired,
    method = "desired_gain",
    n_select = 4L,
    lower_is_better = set$lower
  )
  expect_s3_class(economic, "desiredgainr_generalized_index")
  expect_s3_class(desired, "desiredgainr_generalized_index")
  expect_true(all(is.finite(economic$expected_response)))
  expect_true(all(is.finite(desired$expected_response)))
  expect_equal(
    desired$expected_response / desired$objective,
    setNames(
      rep(desired$expected_response[[1L]] / desired$objective[[1L]], 3L),
      names(desired$objective)
    ),
    tolerance = 1e-9
  )
})

test_that("comparison covers every classical family and decision output", {
  set <- method_coverage_fixture()
  fits <- fit_all_classical_families(set)
  comparison <- compare_selection_methods(fits, target_gains = set$desired)

  expect_s3_class(comparison, "desiredgainr_method_comparison")
  expect_identical(comparison$models, names(fits))
  expect_equal(nrow(comparison$summary), 8L)
  expect_equal(nrow(comparison$responses), 8L * length(set$traits))
  expect_identical(dim(comparison$rank_correlation), c(8L, 8L))
  expect_identical(dim(comparison$selected_jaccard), c(8L, 8L))
  expect_identical(dim(comparison$selected_overlap), c(8L, 8L))
  expect_equal(unname(diag(comparison$rank_correlation)), rep(1, 8))
  expect_equal(unname(diag(comparison$selected_jaccard)), rep(1, 8))
  expect_equal(unname(diag(comparison$selected_overlap)), rep(4, 8))
  expect_equal(
    comparison$rank_correlation,
    t(comparison$rank_correlation),
    tolerance = 1e-12
  )
  expect_equal(
    comparison$selected_jaccard,
    t(comparison$selected_jaccard),
    tolerance = 1e-12
  )
  expect_true(all(comparison$fairness$Satisfied))
  expect_true(all(is.finite(comparison$responses$Observed_differential)))

  linear <- comparison$summary$Family %in%
    c("smith_hazel", "base", "pesek_baker", "yamada")
  expect_true(all(is.finite(
    comparison$summary$Worst_expected_attainment[linear]
  )))
  expect_true(all(is.na(
    comparison$summary$Worst_expected_attainment[!linear]
  )))
  expect_true(all(is.finite(
    comparison$summary$Mahalanobis_alignment[linear]
  )))
  expect_match(
    paste(capture.output(print(comparison)), collapse = "\n"),
    "Methods:"
  )

  without_target <- compare_selection_methods(fits)
  expect_true(all(is.na(without_target$responses$Target)))
  expect_error(compare_selection_methods(fits[1L]), "at least two")
  expect_error(
    compare_selection_methods(unname(fits)),
    "unique, non-empty names"
  )
  expect_error(
    compare_selection_methods(fits, target_gains = c(
      yield = 1, disease = -0.5, height = 0.3
    )),
    "non-negative"
  )

  altered <- fits
  altered$Base$transformation$scale[[1L]] <- 2
  expect_error(
    compare_selection_methods(altered, target_gains = set$desired),
    "same direction, centring, and scaling"
  )
  reordered <- fits
  reordered$Base$candidate_id <- rev(reordered$Base$candidate_id)
  reordered$Base$ranking <- reordered$Base$ranking[
    match(reordered$Base$candidate_id, reordered$Base$ranking$id)
  ]
  expect_s3_class(
    compare_selection_methods(reordered, target_gains = set$desired),
    "desiredgainr_method_comparison"
  )

  expect_error(
    compare_selection_methods(c(fits[1L], Invalid = list(list()))),
    "generalized_index"
  )
})

test_that("comparison accepts low-assumption methods fitted without G", {
  set <- method_coverage_fixture()
  smith <- selection_index(
    set$values, set$traits,
    id_col = "id", method = "smith_hazel",
    G = set$G, P = set$P, economic_weights = set$economic,
    lower_is_better = set$lower,
    center_traits = FALSE, scale_traits = FALSE, n_select = 4L
  )
  rank_sum <- selection_index(
    set$values, set$traits,
    id_col = "id", method = "mulamba_mock",
    lower_is_better = set$lower,
    center_traits = FALSE, scale_traits = FALSE, n_select = 4L
  )

  expect_silent(comparison <- compare_selection_methods(
    list(Smith_Hazel = smith, Rank_sum = rank_sum)
  ))
  expect_s3_class(comparison, "desiredgainr_method_comparison")
  expect_false(comparison$fairness$Satisfied[
    comparison$fairness$Condition ==
      "One common genetic covariance for response geometry"
  ])
  expect_true(all(is.finite(comparison$rank_correlation)))
})

test_that("one objective compares classical, general, DGSI, and QGSI fits", {
  set <- method_coverage_fixture()
  values <- set$values
  names(values)[names(values) == "id"] <- "GenoID"
  common <- list(
    values = values,
    trait_cols = set$traits,
    id_col = "GenoID",
    G = set$G,
    P = set$P,
    lower_is_better = set$lower,
    center_traits = FALSE,
    scale_traits = FALSE,
    n_select = 4L
  )
  smith <- do.call(selection_index, c(common, list(
    method = "smith_hazel", economic_weights = set$economic
  )))
  desired <- do.call(selection_index, c(common, list(
    method = "yamada", desired_gains = set$desired,
    aggregate_weights = set$economic
  )))
  restricted <- restricted_index(
    values, set$traits,
    id_col = "GenoID", method = "restricted_smith_hazel",
    G = set$G, P = set$P, economic_weights = set$economic,
    restricted_traits = "height", lower_is_better = set$lower,
    center_traits = FALSE, scale_traits = FALSE, n_select = 4L
  )
  information <- selection_information(
    values[set$traits], P = set$P, C = set$G, G = set$G
  )
  information$candidate_id <- values$GenoID
  general <- generalized_index(
    information, set$economic,
    method = "economic", n_select = 4L,
    lower_is_better = set$lower
  )
  candidate_sd <- vapply(values[set$traits], stats::sd, numeric(1L))
  dgsi_target <- set$desired * sqrt(diag(set$G)) / candidate_sd
  dgsi <- run_dgsi(
    init_data = values["GenoID"], cand_data = values,
    trait_cols = set$traits, dg = dgsi_target,
    G = set$G, P = set$P, lower_is_better = set$lower,
    n_select = 4L, n_iter = 20L, n_rep = 2L, seed = 91L
  )
  zero_W <- matrix(0, length(set$traits), length(set$traits),
    dimnames = list(set$traits, set$traits)
  )
  qgsi_W <- zero_W
  diag(qgsi_W) <- c(-0.05, -0.03, -0.02)
  qgsi <- run_qgsi(
    init_data = values["GenoID"], gebv_data = values,
    trait_cols = set$traits, linear_weights = set$economic,
    W = qgsi_W, Gamma = set$G, lower_is_better = set$lower,
    center_traits = FALSE, scale_traits = FALSE, n_select = 4L
  )
  objective <- comparison_objective(
    desired_gains = set$desired,
    aggregate_weights = set$economic,
    W = zero_W,
    G = set$G,
    gain_units = "genetic_sd"
  )
  comparison <- compare_selection_methods(
    list(
      Smith_Hazel = smith,
      Desired_gain = desired,
      Restricted = restricted,
      General = general,
      DGSI = dgsi,
      QGSI = qgsi
    ),
    objective = objective,
    validation_data = values
  )

  expect_s3_class(comparison, "desiredgainr_method_comparison")
  expect_equal(nrow(comparison$summary), 6L)
  expect_equal(nrow(comparison$responses), 6L * length(set$traits))
  expect_equal(nrow(comparison$validation_responses), 6L * length(set$traits))
  expect_equal(nrow(comparison$validation_utility), 6L)
  expect_true(all(is.finite(comparison$summary$Validation_utility_response)))
  expect_true(all(is.finite(comparison$summary$Mahalanobis_alignment)))
  expect_true(all(comparison$fairness$Satisfied))
  expect_equal(
    comparison$responses[
      comparison$responses$Method == "Smith_Hazel",
      "Expected_response"
    ],
    comparison$responses[
      comparison$responses$Method == "General",
      "Expected_response"
    ],
    tolerance = 1e-10
  )
  expect_match(
    comparison$summary$Expected_response_basis[
      comparison$summary$Method == "QGSI"
    ],
    "approximation"
  )

  curved_W <- zero_W
  diag(curved_W) <- c(-0.08, -0.04, -0.03)
  curved_objective <- comparison_objective(
    desired_gains = set$desired,
    aggregate_weights = set$economic,
    W = curved_W,
    G = set$G,
    gain_units = "genetic_sd"
  )
  curved <- compare_selection_methods(
    list(Smith_Hazel = smith, QGSI = qgsi),
    objective = curved_objective,
    validation_data = values
  )
  direction <- c(yield = 1, disease = -1, height = -1)
  validation <- sweep(
    as.matrix(values[set$traits]), 2L,
    direction / sqrt(diag(set$G)), "*"
  )
  utility <- as.numeric(validation %*% set$economic) +
    rowSums((validation %*% curved_W) * validation)
  selected <- values$GenoID %in% smith$selected$id
  manual_response <- mean(utility[selected]) - mean(utility)
  expect_equal(
    curved$summary$Validation_utility_response[
      curved$summary$Method == "Smith_Hazel"
    ],
    manual_response,
    tolerance = 1e-12
  )
  expect_true(all(is.finite(curved$summary$Common_merit_response)))
  smith_expected <- curved$responses$Expected_response[
    curved$responses$Method == "Smith_Hazel"
  ]
  expect_equal(
    curved$summary$Common_merit_response[
      curved$summary$Method == "Smith_Hazel"
    ],
    sum(set$economic * smith_expected),
    tolerance = 1e-12
  )
})

test_that("comparison objectives enforce one named scale and covariance", {
  set <- method_coverage_fixture()
  objective <- comparison_objective(
    desired_gains = set$desired,
    G = set$G,
    gain_units = "genetic_sd"
  )
  expect_s3_class(objective, "desiredgainr_comparison_objective")
  expect_error(
    comparison_objective(desired_gains = set$desired,
      gain_units = "genetic_sd"
    ),
    "requires G"
  )
  expect_error(
    comparison_objective(
      desired_gains = set$desired,
      aggregate_weights = c(yield = 1, disease = 1)
    ),
    "same traits"
  )
})

test_that("cross-family response comparison is invariant to fitted scaling", {
  set <- method_coverage_fixture()
  raw <- selection_index(
    set$values, set$traits,
    id_col = "id", method = "smith_hazel",
    G = set$G, P = set$P, economic_weights = set$economic,
    lower_is_better = set$lower,
    center_traits = FALSE, scale_traits = FALSE, n_select = 4L
  )
  sample_sd <- vapply(set$values[set$traits], stats::sd, numeric(1L))
  scaled <- selection_index(
    set$values, set$traits,
    id_col = "id", method = "smith_hazel",
    G = set$G, P = set$P,
    economic_weights = set$economic * sample_sd,
    lower_is_better = set$lower,
    center_traits = TRUE, scale_traits = TRUE, n_select = 4L
  )
  objective <- comparison_objective(
    aggregate_weights = set$economic,
    G = set$G,
    gain_units = "genetic_sd"
  )
  comparison <- compare_selection_methods(
    list(Raw = raw, Scaled = scaled),
    objective = objective,
    validation_data = set$values
  )
  raw_response <- comparison$responses$Expected_response[
    comparison$responses$Method == "Raw"
  ]
  scaled_response <- comparison$responses$Expected_response[
    comparison$responses$Method == "Scaled"
  ]
  expect_equal(raw_response, scaled_response, tolerance = 1e-10)
  expect_equal(comparison$rank_correlation["Raw", "Scaled"], 1,
    tolerance = 1e-12
  )
  expect_true(all(comparison$fairness$Satisfied))
})

test_that("comparison reports undefined geometry and constant rankings", {
  set <- method_coverage_fixture()
  first <- selection_index(
    set$values, set$traits,
    id_col = "id", method = "smith_hazel",
    G = set$G, P = set$P, economic_weights = set$economic,
    lower_is_better = set$lower,
    center_traits = FALSE, scale_traits = FALSE, n_select = 4L
  )
  constant <- first
  constant$ranking$score <- 1
  constant$ranking$rank <- 1

  singular_G <- set$G
  singular_G[3L, ] <- singular_G[2L, ]
  singular_G[, 3L] <- singular_G[, 2L]
  dimnames(singular_G) <- list(set$traits, set$traits)
  objective <- comparison_objective(
    desired_gains = set$desired,
    G = singular_G,
    gain_units = "trait"
  )

  expect_silent(comparison <- compare_selection_methods(
    list(Regular = first, Constant = constant),
    objective = objective,
    validation_data = set$values
  ))
  expect_true(all(is.na(comparison$summary$Mahalanobis_alignment)))
  expect_false(comparison$fairness$Satisfied[
    comparison$fairness$Condition ==
      "Common genetic covariance is invertible"
  ])
  expect_true(is.na(comparison$rank_correlation["Regular", "Constant"]))
  expect_true(is.na(comparison$rank_correlation["Constant", "Constant"]))
})

public_api_test_registry <- c(
  bend_covariance = "test-audit-2026-07-31.R",
  breeding_scenario = "test-literature-extensions.R",
  candidate_score_uncertainty = "test-literature-extensions.R",
  compare_dg_and_qgsi = "test-basic-workflows.R",
  compare_selection_methods = "test-public-api-coverage.R",
  define_desired_gain_intervals = "test-objective-tools.R",
  dosage_diagnostics = "test-simulation.R",
  draw_covariance_pairs = "test-multi-environment.R",
  effective_weights = "test-objective-tools.R",
  estimate_genetic_covariance = "test-genetic-covariance.R",
  evaluate_index = "test-index-families.R",
  evaluate_restricted_response = "test-literature-extensions.R",
  expand_environments = "test-multi-environment.R",
  founder_haplotypes = "test-simulation.R",
  founder_population = "test-simulation.R",
  gain_feasibility = "test-objective-tools.R",
  generalized_index = "test-public-api-coverage.R",
  haplotypes_from_inbred_dosage = "test-simulation.R",
  implied_desired_gains = "test-objective-tools.R",
  implied_economic_weights = "test-objective-tools.R",
  import_covariance = "test-gate-b.R",
  index_information_efficiency = "test-literature-extensions.R",
  index_uncertainty = "test-audit-2026-07-31.R",
  matrix_diagnostics = "test-index-families.R",
  open_desiredgain_guide = "test-breeder-guide.R",
  optimize_desired_gains = "test-optimisation.R",
  propagate_covariance_uncertainty = "test-optimisation.R",
  rahimi_debnath_2023 = "test-published-reproduction.R",
  restricted_breeding_values = "test-literature-extensions.R",
  restricted_index = "test-public-api-coverage.R",
  retrospective_weights = "test-objective-tools.R",
  run_dgsi = "test-basic-workflows.R",
  run_dgsi_qgsi_pipeline = "test-basic-workflows.R",
  run_qgsi = "test-basic-workflows.R",
  run_qgsi_desired_gain = "test-edge-cases.R",
  selection_index = "test-public-api-coverage.R",
  selection_information = "test-literature-extensions.R",
  simulate_selection_cycles = "test-simulation.R",
  simulation_calibration = "test-literature-extensions.R",
  stress_test_desired_gains = "test-literature-extensions.R",
  suggest_desired_gains = "test-population-suggestions.R",
  weight_sensitivity = "test-objective-tools.R",
  widen_environments = "test-multi-environment.R"
)
if ("comparison_objective" %in% getNamespaceExports("DesiredGainR")) {
  public_api_test_registry <- c(
    public_api_test_registry,
    comparison_objective = "test-public-api-coverage.R"
  )
}

test_that("the public API coverage registry is complete", {
  expect_setequal(
    names(public_api_test_registry),
    getNamespaceExports("DesiredGainR")
  )
  for (function_name in names(public_api_test_registry)) {
    test_file <- testthat::test_path(
      public_api_test_registry[[function_name]]
    )
    expect_true(file.exists(test_file), info = function_name)
    source <- readLines(test_file, warn = FALSE)
    expect_true(
      any(grepl(paste0(function_name, "("), source, fixed = TRUE)),
      info = paste(function_name, "is absent from", basename(test_file))
    )
  }
})

test_that("every exported function has a generated manual alias", {
  root <- testthat::test_path("..", "..")
  man_dir <- file.path(root, "man")
  testthat::skip_if_not(
    dir.exists(man_dir),
    "Generated manual sources are absent from the installed test context"
  )
  rd_files <- list.files(man_dir, pattern = "[.]Rd$", full.names = TRUE)
  manual_source <- paste(vapply(
    rd_files,
    function(path) paste(readLines(path, warn = FALSE), collapse = "\n"),
    character(1L)
  ), collapse = "\n")
  for (function_name in getNamespaceExports("DesiredGainR")) {
    expect_true(
      grepl(
        paste0("\\alias{", function_name, "}"),
        manual_source,
        fixed = TRUE
      ),
      info = paste("Missing manual alias for", function_name)
    )
  }
})

test_that("the complete breeder guide exposes every exported function", {
  root <- testthat::test_path("..", "..")
  guide <- file.path(
    root, "inst", "guide", "DesiredGainR_Breeder_Guide.Rmd"
  )
  testthat::skip_if_not(
    file.exists(guide),
    "Breeder Guide source is absent from the installed test context"
  )
  guide_source <- paste(readLines(guide, warn = FALSE), collapse = "\n")
  for (function_name in getNamespaceExports("DesiredGainR")) {
    expect_true(
      grepl(paste0(function_name, "()"), guide_source, fixed = TRUE),
      info = paste("Breeder Guide omits", function_name)
    )
  }
})

test_that("all index families are named in tests and breeder documentation", {
  classical <- c(
    "smith_hazel", "base", "pesek_baker", "yamada", "mulamba_mock",
    "elston", "independent_culling", "tandem"
  )
  restricted <- c(
    "kempthorne_nordskog", "restricted_smith_hazel",
    "tallis", "mallard", "harville"
  )
  expect_setequal(eval(formals(selection_index)$method), classical)
  expect_setequal(eval(formals(restricted_index)$method), restricted)
  expect_setequal(
    eval(formals(generalized_index)$method),
    c("economic", "desired_gain")
  )

  root <- testthat::test_path("..", "..")
  vignette <- file.path(
    root, "vignettes", "DesiredGainR-index-families.Rmd"
  )
  testthat::skip_if_not(
    file.exists(vignette),
    "Vignette sources are absent from the installed test context"
  )
  documentation <- paste(readLines(vignette, warn = FALSE), collapse = "\n")
  for (method in c(classical, restricted)) {
    expect_true(grepl(method, documentation, fixed = TRUE), info = method)
  }
  expect_match(documentation, "run_dgsi()", fixed = TRUE)
  expect_match(documentation, "run_qgsi()", fixed = TRUE)
})
