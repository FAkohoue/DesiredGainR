# Calibrated simulation scenarios and robustness analysis.

#' Describe a breeding simulation scenario
#'
#' A scenario records the biological assumptions that define a simulation.
#' Keeping these assumptions in one object reduces long function calls. It also
#' makes comparisons auditable.
#'
#' The programme can represent self-pollinated, outcrossing, clonal, or
#' testcross evaluation. Testcross describes the source of general combining
#' ability information. Parent allocation and cross design remain with
#' HapBlockR.
#'
#' @param label Short scenario name.
#' @param programme One of self, outcross, clonal, or testcross.
#' @param architecture Named list. Supported fields include
#'   qtl_per_chromosome, markers_per_chromosome, effect_distribution,
#'   qtl_shape, dominance_degree, and dominance_variance.
#' @param evaluation Named list describing phenotype, GEBV, or GCA information.
#' @param environments Named list describing the target environments and their
#'   genetic correlation structure.
#'
#' @return An object of class desiredgainr_scenario.
#'
#' @references
#' Beavis WD, Mahama AA, Suza W (2023). Simulation Modeling. In
#' *Quantitative Genetics for Plant Breeding*. Iowa State University Digital
#' Press.
#'
#' @export
breeding_scenario <- function(
  label,
  programme = c("self", "outcross", "clonal", "testcross"),
  architecture = list(),
  evaluation = list(),
  environments = list()
) {
  if (!is.character(label) || length(label) != 1L || !nzchar(label)) {
    stop("label must contain one short scenario name.", call. = FALSE)
  }
  programme <- match.arg(programme)
  if (!is.list(architecture) || !is.list(evaluation) ||
    !is.list(environments)) {
    stop("architecture, evaluation, and environments must be named lists.",
      call. = FALSE
    )
  }
  architecture <- .dgr_merge_named_list(
    list(
      qtl_per_chromosome = 100L,
      markers_per_chromosome = NULL,
      effect_distribution = "normal",
      qtl_shape = 1,
      dominance_degree = NULL,
      dominance_variance = NULL
    ),
    architecture, "architecture"
  )
  architecture$qtl_per_chromosome <- .dgr_positive_integer(
    architecture$qtl_per_chromosome, "architecture$qtl_per_chromosome"
  )
  if (!is.null(architecture$markers_per_chromosome)) {
    architecture$markers_per_chromosome <- .dgr_positive_integer(
      architecture$markers_per_chromosome,
      "architecture$markers_per_chromosome"
    )
  }
  architecture$effect_distribution <- match.arg(
    architecture$effect_distribution, c("normal", "gamma")
  )
  if (!is.numeric(architecture$qtl_shape) ||
    length(architecture$qtl_shape) != 1L ||
    !is.finite(architecture$qtl_shape) ||
    architecture$qtl_shape <= 0) {
    stop("architecture$qtl_shape must be positive.", call. = FALSE)
  }
  if (identical(programme, "clonal") &&
    is.null(architecture$dominance_degree)) {
    stop("A clonal scenario requires architecture$dominance_degree.",
      call. = FALSE
    )
  }

  evaluation <- .dgr_merge_named_list(
    list(type = if (identical(programme, "testcross")) "gca" else "phenotype"),
    evaluation, "evaluation"
  )
  evaluation$type <- match.arg(evaluation$type, c("phenotype", "gebv", "gca"))
  if (identical(programme, "testcross") && !identical(evaluation$type, "gca")) {
    stop("A testcross scenario uses evaluation$type equal to gca.",
      call. = FALSE
    )
  }
  environments <- .dgr_merge_named_list(
    list(names = "target", genetic_correlation = NULL),
    environments, "environments"
  )
  if (!is.character(environments$names) || !length(environments$names) ||
    anyDuplicated(environments$names)) {
    stop("environments$names must contain unique environment names.",
      call. = FALSE
    )
  }

  result <- list(
    label = label,
    programme = programme,
    architecture = architecture,
    evaluation = evaluation,
    environments = environments,
    scope = if (identical(programme, "testcross")) {
      paste(
        "The scenario represents GCA or testcross information.",
        "HapBlockR remains responsible for parent and cross decisions."
      )
    } else {
      "The scenario compares desired-gain objectives under truncation selection."
    }
  )
  class(result) <- c("desiredgainr_scenario", "list")
  result
}

#' Audit the calibration of an AlphaSimR setup
#'
#' The report compares the target and realised covariance. It also checks
#' heritability, marker coverage, and marker overlap with quantitative trait
#' loci. Each check receives a clear status.
#'
#' @param setup Object returned by [founder_population()].
#' @param tolerances Named list with variance, correlation, and heritability
#'   limits.
#'
#' @return An object of class desiredgainr_calibration.
#'
#' @export
simulation_calibration <- function(
  setup,
  tolerances = list(variance = 0.05, correlation = 0.05, heritability = 0.05)
) {
  if (!inherits(setup, "desiredgainr_sim_setup")) {
    stop("setup must come from founder_population().", call. = FALSE)
  }
  tolerances <- .dgr_merge_named_list(
    list(variance = 0.05, correlation = 0.05, heritability = 0.05),
    tolerances, "tolerances"
  )
  if (any(!is.finite(unlist(tolerances))) || any(unlist(tolerances) <= 0)) {
    stop("Every calibration tolerance must be positive.", call. = FALSE)
  }
  target <- setup$G_target
  realised <- setup$G_realised
  variance_error <- abs(diag(realised) - diag(target)) / diag(target)
  correlation_error <- abs(stats::cov2cor(realised) - stats::cov2cor(target))
  max_correlation_error <- if (nrow(correlation_error) > 1L) {
    max(correlation_error[upper.tri(correlation_error)])
  } else {
    0
  }

  SP <- setup$SP
  realised_h2 <- tryCatch(
    {
      genetic <- if (identical(setup$heritability_type, "broad")) {
        as.numeric(SP$varG)
      } else {
        as.numeric(SP$varA)
      }
      residual <- if (is.matrix(SP$varE)) diag(SP$varE) else as.numeric(SP$varE)
      genetic / (as.numeric(SP$varG) + residual)
    },
    error = function(e) rep(NA_real_, length(setup$h2))
  )
  names(realised_h2) <- names(setup$h2)
  heritability_error <- abs(realised_h2 - setup$h2)

  checks <- data.table::rbindlist(list(
    data.table::data.table(
      Check = paste0("variance_", names(variance_error)),
      Estimate = as.numeric(variance_error),
      Limit = tolerances$variance,
      Status = ifelse(variance_error <= tolerances$variance, "pass", "review")
    ),
    data.table::data.table(
      Check = "largest_correlation_error",
      Estimate = max_correlation_error,
      Limit = tolerances$correlation,
      Status = ifelse(
        max_correlation_error <= tolerances$correlation, "pass", "review"
      )
    ),
    data.table::data.table(
      Check = paste0("heritability_", names(heritability_error)),
      Estimate = as.numeric(heritability_error),
      Limit = tolerances$heritability,
      Status = ifelse(
        is.na(heritability_error), "unavailable",
        ifelse(heritability_error <= tolerances$heritability, "pass", "review")
      )
    ),
    data.table::data.table(
      Check = "neutral_marker_panel",
      Estimate = if (is.null(setup$marker_panel)) 0 else length(setup$marker_panel),
      Limit = 1,
      Status = if (is.null(setup$marker_panel)) "review" else "pass"
    ),
    data.table::data.table(
      Check = "marker_qtl_disjoint",
      Estimate = if (is.null(setup$marker_panel)) {
        NA_real_
      } else {
        as.numeric(!isTRUE(setup$marker_qtl_overlap))
      },
      Limit = 1,
      Status = if (is.null(setup$marker_panel)) {
        "unavailable"
      } else if (isTRUE(setup$marker_qtl_overlap)) {
        "review"
      } else {
        "pass"
      }
    )
  ), use.names = TRUE)

  result <- list(
    checks = checks,
    passed = all(checks$Status %in% c("pass", "unavailable")),
    target_G = target,
    realised_G = realised,
    realised_h2 = realised_h2,
    largest_variance_error = max(variance_error),
    largest_correlation_error = max_correlation_error,
    scenario = setup$scenario,
    interpretation = paste(
      "A review status identifies an assumption that needs inspection before",
      "the setup supports a breeding recommendation."
    )
  )
  class(result) <- c("desiredgainr_calibration", "list")
  result
}

#' Stress-test a desired-gain vector across simulation scenarios
#'
#' The function applies matched replicate seeds across calibrated setups. It
#' reports mean gain, Monte Carlo standard error, tail risk, target attainment,
#' and regret. Replication can continue until the utility standard error reaches
#' a stated precision.
#'
#' @param setups Named list of objects returned by [founder_population()].
#' @param desired_gains Named desired-gain direction.
#' @param options Named list. Fields include simulation, minimum_gains,
#'   min_replicates, max_replicates, batch_size, utility_mcse, and seed.
#'
#' @return An object of class desiredgainr_stress_test.
#'
#' @export
stress_test_desired_gains <- function(setups, desired_gains, options = list()) {
  if (inherits(setups, "desiredgainr_sim_setup")) setups <- list(base = setups)
  if (!is.list(setups) || !length(setups) ||
    any(!vapply(setups, inherits, logical(1L), "desiredgainr_sim_setup"))) {
    stop("setups must contain founder_population() results.", call. = FALSE)
  }
  if (is.null(names(setups)) || any(!nzchar(names(setups)))) {
    names(setups) <- paste0("scenario_", seq_along(setups))
  }
  options <- .dgr_merge_named_list(
    list(
      simulation = list(n_cycles = 5L),
      minimum_gains = NULL,
      min_replicates = 20L,
      max_replicates = 100L,
      batch_size = 10L,
      utility_mcse = NULL,
      seed = 42L
    ),
    options, "options"
  )
  options$min_replicates <- .dgr_positive_integer(
    options$min_replicates, "options$min_replicates"
  )
  options$max_replicates <- .dgr_positive_integer(
    options$max_replicates, "options$max_replicates"
  )
  options$batch_size <- .dgr_positive_integer(
    options$batch_size, "options$batch_size"
  )
  if (options$min_replicates < 2L) {
    stop("min_replicates must be at least two for Monte Carlo uncertainty.",
      call. = FALSE
    )
  }
  if (options$max_replicates < options$min_replicates) {
    stop("max_replicates must reach min_replicates.", call. = FALSE)
  }
  if (!is.null(options$utility_mcse) &&
    (!is.numeric(options$utility_mcse) ||
      length(options$utility_mcse) != 1L ||
      !is.finite(options$utility_mcse) ||
      options$utility_mcse <= 0)) {
    stop("options$utility_mcse must be positive.", call. = FALSE)
  }
  reserved <- intersect(
    names(options$simulation), c("setup", "desired_gains", "seed")
  )
  if (length(reserved)) {
    stop(
      "options$simulation must leave setup, desired_gains, and seed to the ",
      "stress-test controller.",
      call. = FALSE
    )
  }

  traits <- setups[[1L]]$trait_cols
  if (any(vapply(
    setups, function(x) !identical(x$trait_cols, traits),
    logical(1L)
  ))) {
    stop("Every setup must use the same ordered traits.", call. = FALSE)
  }
  desired_gains <- .dgr_named_vector(desired_gains, traits, "desired_gains")
  minimum_gains <- if (is.null(options$minimum_gains)) {
    NULL
  } else {
    .dgr_named_vector(options$minimum_gains, traits, "minimum_gains")
  }
  lower <- options$simulation$lower_is_better
  signs <- stats::setNames(rep(1, length(traits)), traits)
  if (length(lower)) signs[lower] <- -1
  direction <- abs(desired_gains)
  direction <- direction / sqrt(sum(direction^2))

  outcomes <- lapply(setups, function(x) {
    list(gain = matrix(numeric(0),
      nrow = 0, ncol = length(traits),
      dimnames = list(NULL, traits)
    ), utility = numeric(0))
  })
  total <- 0L
  repeat {
    additional <- min(options$batch_size, options$max_replicates - total)
    if (additional <= 0L) break
    replicate_ids <- total + seq_len(additional)
    for (scenario in seq_along(setups)) {
      setup <- setups[[scenario]]
      genetic_sd <- sqrt(diag(setup$G_target))
      for (replicate_id in replicate_ids) {
        arguments <- c(
          list(
            setup = setup,
            desired_gains = desired_gains,
            seed = options$seed + replicate_id - 1L
          ),
          options$simulation
        )
        run <- do.call(simulate_selection_cycles, arguments)
        gain <- run$cumulative_gain[traits]
        oriented_sd <- gain * signs / genetic_sd
        outcomes[[scenario]]$gain <- rbind(outcomes[[scenario]]$gain, gain)
        outcomes[[scenario]]$utility <- c(
          outcomes[[scenario]]$utility,
          sum(oriented_sd * direction)
        )
      }
    }
    total <- total + additional
    if (total < options$min_replicates) next
    if (is.null(options$utility_mcse)) break
    utility_error <- vapply(outcomes, function(x) {
      stats::sd(x$utility) / sqrt(length(x$utility))
    }, numeric(1L))
    if (all(utility_error <= options$utility_mcse)) break
  }

  details <- data.table::rbindlist(lapply(seq_along(setups), function(i) {
    gain <- outcomes[[i]]$gain
    data.table::data.table(
      Scenario = names(setups)[i],
      Trait = traits,
      Mean_gain = colMeans(gain),
      SD = apply(gain, 2L, stats::sd),
      MCSE = apply(gain, 2L, stats::sd) / sqrt(nrow(gain)),
      Lower_05 = apply(gain, 2L, stats::quantile, probs = 0.05),
      Upper_95 = apply(gain, 2L, stats::quantile, probs = 0.95),
      Unfavourable_probability = colMeans(sweep(gain, 2L, signs, "*") < 0)
    )
  }))
  summary <- data.table::rbindlist(lapply(seq_along(setups), function(i) {
    gain <- outcomes[[i]]$gain
    oriented <- sweep(gain, 2L, signs, "*")
    joint <- if (is.null(minimum_gains)) {
      NA_real_
    } else {
      mean(apply(
        sweep(oriented, 2L, abs(minimum_gains), FUN = ">="),
        1L, all
      ))
    }
    utility <- outcomes[[i]]$utility
    data.table::data.table(
      Scenario = names(setups)[i],
      Replicates = length(utility),
      Mean_utility = mean(utility),
      Utility_MCSE = stats::sd(utility) / sqrt(length(utility)),
      Joint_minimum_probability = joint,
      Any_unfavourable_probability = mean(apply(oriented < 0, 1L, any))
    )
  }))
  data.table::set(
    summary,
    j = "Regret",
    value = max(summary$Mean_utility) - summary$Mean_utility
  )
  result <- list(
    summary = summary,
    trait_results = details,
    replicate_outcomes = outcomes,
    desired_gains = desired_gains,
    minimum_gains = minimum_gains,
    options = options,
    matched_replicate_seeds = TRUE,
    random_number_alignment = paste(
      "Replicates use matched seeds across scenarios.",
      "Exact common-random-number alignment requires the same simulator call path."
    ),
    precision_reached = is.null(options$utility_mcse) ||
      all(summary$Utility_MCSE <= options$utility_mcse),
    interpretation = paste(
      "Regret measures the loss relative to the strongest tested scenario.",
      "Tail risks and target probabilities retain the stochastic uncertainty."
    )
  )
  class(result) <- c("desiredgainr_stress_test", "list")
  result
}

#' @export
print.desiredgainr_scenario <- function(x, ...) {
  cat("<desiredgainr_scenario>\n")
  cat("  Label:", x$label, "\n")
  cat("  Programme:", x$programme, "\n")
  cat("  Evaluation:", x$evaluation$type, "\n")
  cat("  QTL per chromosome:", x$architecture$qtl_per_chromosome, "\n")
  cat("  QTL effects:", x$architecture$effect_distribution, "\n")
  invisible(x)
}

#' @export
print.desiredgainr_calibration <- function(x, ...) {
  cat("<desiredgainr_calibration>\n")
  cat("  Overall status:", if (x$passed) "pass" else "review", "\n")
  print(x$checks)
  invisible(x)
}

#' @export
print.desiredgainr_stress_test <- function(x, ...) {
  cat("<desiredgainr_stress_test>\n")
  cat("  Matched replicate seeds: yes\n")
  cat(
    "  Precision target reached:",
    if (x$precision_reached) "yes" else "review", "\n"
  )
  print(x$summary)
  invisible(x)
}
