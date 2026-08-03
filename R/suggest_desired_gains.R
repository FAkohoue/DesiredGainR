# Population-driven desired-gain recommendations.

.dgr_population_genetic_covariance <- function(setup) {
  .dgr_covariance(setup$G_target, setup$trait_cols, "setup$G_target")
}

# G_target is the biological covariance estimate on which breeder thresholds
# and the exact selection-index geometry are defined. G_realised is a finite-
# QTL simulation calibration result; substituting it would make the scientific
# estimand depend on an arbitrary simulator seed.
.dgr_population_model_diagnostic <- function(
  setup, variance_tolerance = 0.15, correlation_tolerance = 0.10
) {
  target <- .dgr_population_genetic_covariance(setup)
  realised <- setup$G_realised
  if (is.null(realised)) {
    return(list(
      status = "not_assessed", calibrated = NA,
      reason = "setup$G_realised is unavailable."
    ))
  }
  realised <- .dgr_covariance(
    realised, setup$trait_cols, "setup$G_realised"
  )
  variance_error <- abs(diag(realised) - diag(target)) / diag(target)
  correlation_error <- suppressWarnings(abs(
    stats::cov2cor(realised) - stats::cov2cor(target)
  ))
  diag(correlation_error) <- 0
  maximum_variance_error <- if (all(is.finite(variance_error))) {
    max(variance_error)
  } else {
    Inf
  }
  maximum_correlation_error <- if (all(is.finite(correlation_error))) {
    max(correlation_error)
  } else {
    Inf
  }
  calibrated <- maximum_variance_error <= variance_tolerance &&
    maximum_correlation_error <= correlation_tolerance
  list(
    status = if (calibrated) "adequate" else "inadequate",
    calibrated = calibrated,
    maximum_relative_variance_error = maximum_variance_error,
    maximum_absolute_correlation_error = maximum_correlation_error,
    variance_error_by_trait = variance_error,
    variance_tolerance = variance_tolerance,
    correlation_tolerance = correlation_tolerance,
    reason = if (calibrated) {
      "The simulated finite-QTL architecture reproduces G_target within tolerance."
    } else {
      paste(
        "The simulated finite-QTL architecture does not reproduce G_target",
        "within the declared variance/correlation tolerances."
      )
    }
  )
}

# Construct the phenotypic covariance implied by a simulation setup.
.dgr_setup_phenotypic_covariance <- function(setup, G) {
  if (!is.null(setup$residual_covariance)) {
    return(G + setup$residual_covariance)
  }
  if (isTRUE(setup$dominance) && identical(setup$heritability_type, "narrow")) {
    stop(
      "The current phenotypic covariance cannot be recovered exactly from a ",
      "dominance setup using narrow-sense h2 alone. Rebuild the setup with ",
      "heritability = 'broad' or an explicit residual_covariance before ",
      "requesting a population-driven mathematical recommendation.",
      call. = FALSE
    )
  }
  P <- G
  diag(P) <- diag(G) / as.numeric(setup$h2)
  P
}

# Solve min_x x' B x subject to x >= lower.
#
# B is positive definite, so the KKT conditions are necessary and sufficient.
# The optimiser enumerates active bound sets. DesiredGainR's direction search
# supports at most 15 traits, making this exact active-set solution practical
# while avoiding convergence-tolerance dependence in a recommendation.
.dgr_quadratic_lower_bound <- function(B, lower, tolerance = 1e-9) {
  B <- as.matrix(B)
  lower <- as.numeric(lower)
  p <- length(lower)
  if (!identical(dim(B), c(p, p)) || any(!is.finite(B)) ||
    any(!is.finite(lower))) {
    stop("Invalid quadratic lower-bound problem.", call. = FALSE)
  }
  if (p > 15L) {
    stop("The exact active-set solver supports at most 15 traits.",
      call. = FALSE
    )
  }
  B <- (B + t(B)) / 2
  eigenvalues <- eigen(B, symmetric = TRUE, only.values = TRUE)$values
  if (min(eigenvalues) <= 0) {
    stop("The standardised response ellipsoid must be positive definite.",
      call. = FALSE
    )
  }

  scale <- max(1, max(abs(B)), max(abs(lower)))
  tol <- tolerance * scale
  solved <- .dgr_box_qp(
    H = 2 * B, f = numeric(p), lower = lower,
    upper = rep(Inf, p), tolerance = tolerance
  )
  best <- solved$solution
  best_value <- solved$objective
  gradient <- as.numeric(B %*% best)
  active <- abs(best - lower) <= 50 * tol
  kkt_residual <- max(
    max(pmax(lower - best, 0)),
    if (any(!active)) max(abs(gradient[!active])) else 0,
    if (any(active)) max(pmax(-gradient[active], 0)) else 0
  )
  list(
    solution = best,
    objective = best_value,
    active = active,
    kkt_residual = kkt_residual,
    eigenvalues = eigenvalues,
    globally_optimal = identical(solved$status, "solved") &&
      kkt_residual <= 100 * tol,
    solver_status = solved$status
  )
}

# Exact one-cycle population limits in favourable genetic-SD units.
.dgr_population_feasibility <- function(
  setup, minimum_gains, lower_is_better, n_parents
) {
  traits <- setup$trait_cols
  direction <- .dgr_direction(traits, lower_is_better)
  S <- diag(direction, nrow = length(traits))
  dimnames(S) <- list(traits, traits)
  G_current <- .dgr_population_genetic_covariance(setup)
  G <- S %*% G_current %*% S
  P_raw <- .dgr_setup_phenotypic_covariance(setup, G_current)
  P <- S %*% P_raw %*% S
  dimnames(G) <- dimnames(P) <- list(traits, traits)
  G_inverse <- .dgr_inverse(G, "G")$inverse
  A <- G_inverse %*% P %*% G_inverse
  D <- diag(sqrt(diag(G)), nrow = length(traits))
  dimnames(D) <- list(traits, traits)
  B <- D %*% A %*% D
  B <- (B + t(B)) / 2

  n_available <- AlphaSimR::nInd(setup$founder_pop)
  if (n_parents > n_available) {
    stop("programme$n_parents exceeds the founder candidates available.",
      call. = FALSE
    )
  }
  proportion <- n_parents / n_available
  intensity <- .dgr_intensity(proportion)

  threshold <- .dgr_quadratic_lower_bound(B, minimum_gains)
  common <- .dgr_quadratic_lower_bound(B, rep(1, length(traits)))
  threshold_scale <- if (threshold$objective > 0) {
    intensity / sqrt(threshold$objective)
  } else {
    Inf
  }
  maximum_common_gain <- if (common$objective > 0) {
    intensity / sqrt(common$objective)
  } else {
    Inf
  }
  names(threshold$solution) <- names(common$solution) <- traits
  threshold_response <- threshold$solution * threshold_scale
  common_response <- common$solution * maximum_common_gain
  limiting <- traits[common$active]

  list(
    theorem = paste(
      "For a linear index, R = i G b / sqrt(b' P b), hence",
      "R' G^-1 P G^-1 R = i^2. In favourable genetic-SD coordinates z,",
      "the package solves min z' B z subject to z >= the requested bounds.",
      "Positive-definiteness makes the problem strictly convex; the reported",
      "active-set solution satisfies the necessary and sufficient KKT",
      "conditions and is therefore the unique global optimum."
    ),
    genetic_covariance_source = paste(
      "setup$G_target (the biological covariance estimate; G_realised is",
      "used only as a simulation-calibration diagnostic)"
    ),
    first_cycle_selection_proportion = proportion,
    first_cycle_selection_intensity = intensity,
    minimum_gains = stats::setNames(minimum_gains, traits),
    minimum_required_intensity = sqrt(threshold$objective),
    minimum_feasible_in_one_cycle =
      intensity + 1e-10 >= sqrt(threshold$objective),
    minimum_cost_direction = threshold$solution,
    minimum_direction_attainable_response = threshold_response,
    minimum_kkt_residual = threshold$kkt_residual,
    maximum_common_gain = maximum_common_gain,
    maximum_common_direction = common$solution,
    maximum_common_response = common_response,
    limiting_traits = limiting,
    maximum_common_kkt_residual = common$kkt_residual,
    qualification = paste(
      "This is an exact one-cycle result under the supplied G and P and the",
      "normal-truncation selection intensity. Multi-cycle segregation, drift,",
      "recombination, finite-family sampling and variance erosion are assessed",
      "separately by simulation."
    )
  )
}

.dgr_jeffreys_binomial <- function(successes, trials, level) {
  alpha <- (1 - level) / 2
  shapes <- c(successes + 0.5, trials - successes + 0.5)
  c(
    mean = shapes[1L] / sum(shapes),
    lower = stats::qbeta(alpha, shapes[1L], shapes[2L]),
    upper = stats::qbeta(1 - alpha, shapes[1L], shapes[2L])
  )
}

# Bonferroni-adjusted exact one-sided Clopper-Pearson lower bound. The union
# bound does not require independence among directions, so it remains valid
# under common random numbers for a fixed, predeclared collection of vectors.
.dgr_simultaneous_binomial_lower <- function(successes, trials, alpha, tests) {
  if (successes <= 0) {
    return(0)
  }
  stats::qbeta(alpha / tests, successes, trials - successes + 1)
}

.dgr_simultaneous_binomial_upper <- function(successes, trials, alpha, tests) {
  if (successes >= trials) {
    return(1)
  }
  stats::qbeta(1 - alpha / tests, successes + 1, trials - successes)
}

# Distribution-free simultaneous lower confidence bound for a median. For the
# population median theta, P(X_(k) <= theta) is at least
# P{Binomial(n, 1/2) >= k}. We select the largest k whose non-coverage is no
# greater than alpha/tests.
.dgr_simultaneous_median_lower <- function(x, alpha, tests) {
  x <- sort(as.numeric(x))
  n <- length(x)
  admissible <- which(vapply(seq_len(n), function(k) {
    stats::pbinom(k - 1L, n, 0.5) <= alpha / tests
  }, logical(1L)))
  if (!length(admissible)) {
    return(-Inf)
  }
  x[max(admissible)]
}

.dgr_paired_best_probability <- function(values, seed, draws = 2000L) {
  n_candidates <- length(values)
  n <- min(vapply(values, length, integer(1L)))
  if (n < 2L) {
    return(rep(NA_real_, n_candidates))
  }
  if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    stats::runif(1L)
  }
  entry_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit(assign(".Random.seed", entry_seed, envir = globalenv()), add = TRUE)
  set.seed(.dgr_seed(seed))
  wins <- numeric(n_candidates)
  for (draw in seq_len(draws)) {
    sampled <- sample.int(n, n, replace = TRUE)
    scores <- vapply(values, function(x) mean(x[sampled]), numeric(1L))
    best <- which(scores >= max(scores) - 1e-12)
    wins[best] <- wins[best] + 1 / length(best)
  }
  wins / draws
}

.dgr_population_candidate_summary <- function(
  directions, replicates, traits, minimum_gains, scale, direction,
  level, seed
) {
  n_candidates <- length(replicates)
  alpha <- 1 - level
  joint_values <- worst_values <- vector("list", n_candidates)
  result <- data.table::data.table(candidate = seq_len(n_candidates))
  for (trait in traits) {
    data.table::set(result,
      j = paste0("d_", trait),
      value = directions[, trait]
    )
  }
  joint_mean <- joint_lower <- joint_upper <- joint_exact_lower <-
    joint_exact_upper <-
    mean_worst <- median_worst <- median_exact_lower <- numeric(n_candidates)
  per_trait <- per_trait_lower <- per_trait_upper <- matrix(
    NA_real_, n_candidates, length(traits),
    dimnames = list(NULL, traits)
  )
  mean_gain <- median_gain <- per_trait
  for (i in seq_len(n_candidates)) {
    raw <- as.matrix(replicates[[i]])[traits, , drop = FALSE]
    favourable <- sweep(raw, 1L, direction, "*")
    z <- sweep(favourable, 1L, scale, "/")
    achieved <- sweep(z, 1L, minimum_gains, ">=")
    joint <- as.numeric(colSums(achieved) == length(traits))
    worst <- apply(z, 2L, min)
    posterior <- .dgr_jeffreys_binomial(sum(joint), length(joint), level)
    joint_values[[i]] <- joint
    worst_values[[i]] <- worst
    joint_mean[i] <- posterior["mean"]
    joint_lower[i] <- posterior["lower"]
    joint_upper[i] <- posterior["upper"]
    joint_exact_lower[i] <- .dgr_simultaneous_binomial_lower(
      sum(joint), length(joint), alpha, n_candidates
    )
    joint_exact_upper[i] <- .dgr_simultaneous_binomial_upper(
      sum(joint), length(joint), alpha, n_candidates
    )
    mean_worst[i] <- mean(worst)
    median_worst[i] <- stats::median(worst)
    median_exact_lower[i] <- .dgr_simultaneous_median_lower(
      worst, alpha, n_candidates
    )
    mean_gain[i, ] <- rowMeans(z)
    median_gain[i, ] <- apply(z, 1L, stats::median)
    for (trait in seq_along(traits)) {
      trait_posterior <- .dgr_jeffreys_binomial(
        sum(achieved[trait, ]), ncol(achieved), level
      )
      per_trait[i, trait] <- trait_posterior["mean"]
      per_trait_lower[i, trait] <- trait_posterior["lower"]
      per_trait_upper[i, trait] <- trait_posterior["upper"]
    }
  }
  result[, `:=`(
    joint_probability = joint_mean,
    joint_probability_lower = joint_lower,
    joint_probability_upper = joint_upper,
    joint_exact_simultaneous_lower = joint_exact_lower,
    joint_exact_simultaneous_upper = joint_exact_upper,
    mean_worst_trait_gain = mean_worst,
    median_worst_trait_gain = median_worst,
    median_worst_exact_simultaneous_lower = median_exact_lower,
    bootstrap_selection_frequency_minimum = .dgr_paired_best_probability(
      joint_values, .dgr_shift_seed(seed, 1L)
    ),
    bootstrap_selection_frequency_balanced = .dgr_paired_best_probability(
      worst_values, .dgr_shift_seed(seed, 2L)
    )
  )]
  for (trait in traits) {
    data.table::set(result,
      j = paste0("minimum_probability_", trait),
      value = per_trait[, trait]
    )
    data.table::set(
      result,
      j = paste0("minimum_probability_lower_", trait),
      value = per_trait_lower[, trait]
    )
    data.table::set(
      result,
      j = paste0("minimum_probability_upper_", trait),
      value = per_trait_upper[, trait]
    )
    data.table::set(result,
      j = paste0("mean_gain_", trait),
      value = mean_gain[, trait]
    )
    data.table::set(result,
      j = paste0("median_gain_", trait),
      value = median_gain[, trait]
    )
  }
  list(table = result, joint = joint_values, worst = worst_values)
}

.dgr_merge_named_list <- function(defaults, supplied, name) {
  if (!is.list(supplied) || is.null(names(supplied)) && length(supplied)) {
    stop(name, " must be a named list.", call. = FALSE)
  }
  unknown <- setdiff(names(supplied), names(defaults))
  if (length(unknown)) {
    stop("Unknown ", name, " field(s): ", paste(unknown, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  utils::modifyList(defaults, supplied)
}

.dgr_direction_angle <- function(x, y) {
  cosine <- sum(x * y) / sqrt(sum(x^2) * sum(y^2))
  acos(pmin(1, pmax(-1, cosine))) * 180 / pi
}

.dgr_shift_seed <- function(seed, offset) {
  as.integer((abs(as.double(seed)) + as.double(offset)) %%
    (.Machine$integer.max - 1) + 1)
}

.dgr_support_status <- function(row, required_probability) {
  if (row$joint_exact_simultaneous_lower >= required_probability) {
    "supported"
  } else if (row$joint_exact_simultaneous_upper < required_probability) {
    "not_supported"
  } else {
    "uncertain"
  }
}

.dgr_rebuild_population_setup <- function(setup, G, E = NULL, seed) {
  suppressWarnings(founder_population(
    founders = setup$founders,
    G = G,
    h2 = setup$h2,
    residual_covariance = E,
    n_qtl_per_chromosome = setup$n_qtl_per_chromosome,
    n_markers_per_chromosome = setup$n_markers_per_chromosome,
    dominance_degree = setup$dominance_degree,
    dominance_variance = setup$dominance_degree_variance,
    heritability = setup$heritability_type,
    seed = seed
  ))
}

# Independent outer validation across newly sampled finite-QTL architectures
# and, when degrees of freedom are supplied, covariance sampling distributions.
.dgr_population_outer_validation <- function(
  setup, directions, traits, minimum_gains, scale, direction, programme,
  n_cycles, uncertainty, level, seed
) {
  if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    stats::runif(1L)
  }
  entry_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit(assign(".Random.seed", entry_seed, envir = globalenv()), add = TRUE)
  set.seed(.dgr_shift_seed(seed, 290000L))
  validate_setups <- function(setups, label, seed_offset) {
    if (!length(setups)) {
      return(list(status = "not_assessed", reason = paste(label, "draws = 0")))
    }
    usable <- Filter(Negate(is.null), setups)
    if (!length(usable)) {
      return(list(status = "failed", reason = paste(
        "No", label, "draw could be realised as a simulation architecture."
      )))
    }
    outcomes <- lapply(seq_len(nrow(directions)), function(candidate) {
      values <- vapply(seq_along(usable), function(draw) {
        .dgr_evaluate_direction(
          directions[candidate, ], usable[[draw]],
          c(
            list(n_cycles = n_cycles), programme,
            list(lower_is_better = names(direction)[direction < 0])
          ),
          n_replicates = 1L,
          base_seed = .dgr_shift_seed(seed, seed_offset + draw * 1009L),
          include_diversity = FALSE,
          direction_scale = scale
        )$replicates[, 1L]
      }, numeric(length(traits)))
      rownames(values) <- traits
      values
    })
    summary <- .dgr_population_candidate_summary(
      directions, outcomes, traits, minimum_gains, scale, direction,
      level, .dgr_shift_seed(seed, seed_offset)
    )$table
    list(
      status = "assessed", label = label, requested_draws = length(setups),
      usable_draws = length(usable), failed_draws = length(setups) - length(usable),
      results = summary, outcomes = outcomes
    )
  }

  architecture_n <- uncertainty$architecture_draws
  architecture_setups <- lapply(seq_len(architecture_n), function(draw) {
    tryCatch(
      .dgr_rebuild_population_setup(
        setup, setup$G_target,
        E = setup$residual_covariance,
        seed = .dgr_shift_seed(seed, 300000L + draw)
      ),
      error = function(e) NULL
    )
  })
  architecture <- validate_setups(
    architecture_setups, "genetic-architecture", 310000L
  )

  covariance <- list(
    status = "not_assessed",
    reason = "Covariance degrees of freedom and draws were not supplied."
  )
  if (uncertainty$covariance_draws > 0L) {
    if (is.null(uncertainty$genetic_df)) {
      stop("control$uncertainty$genetic_df is required when covariance_draws > 0.",
        call. = FALSE
      )
    }
    P <- .dgr_setup_phenotypic_covariance(setup, setup$G_target)
    draws <- draw_covariance_pairs(
      setup$G_target,
      P = P,
      genetic_df = uncertainty$genetic_df,
      residual_df = uncertainty$residual_df,
      n_draws = uncertainty$covariance_draws
    )
    covariance_setups <- lapply(seq_along(draws), function(draw) {
      tryCatch(
        .dgr_rebuild_population_setup(
          setup, draws[[draw]]$G, draws[[draw]]$E,
          seed = .dgr_shift_seed(seed, 400000L + draw)
        ),
        error = function(e) NULL
      )
    })
    covariance <- validate_setups(
      covariance_setups, "covariance", 410000L
    )
  }
  list(architecture = architecture, covariance = covariance)
}

#' Suggest desired-gain directions from the current population
#'
#' Gives a breeder who cannot state a complete desired-gain vector two
#' population-driven recommendations without inventing economic weights:
#'
#' 1. a **minimum-attainment direction**, chosen to maximise robust evidence
#'    that every trait reaches its own breeder-specified minimum; and
#' 2. a **maximum-balanced direction**, chosen to maximise the typical gain of
#'    the worst-responding trait in genetic standard deviations.
#'
#' @details
#' `minimum_gains` is always stated in favourable genetic standard deviations
#' calculated from `setup$G_target`, the biological genetic-covariance estimate
#' supplied by the breeder. `setup$G_realised` is used only to diagnose whether
#' the finite-QTL simulation reproduces that target; using it for the units
#' would make a breeder's threshold depend on the simulator seed. Thus
#' `c(yield = 1, disease = 0.5)` asks for at least one estimated
#' genetic SD of yield improvement and half an SD of disease reduction when
#' `lower_is_better = "disease"`. Values may differ by trait and zero means
#' non-decline, not omission. At least one value must be positive.
#'
#' The public call is deliberately short. `programme` contains arguments of
#' [simulate_selection_cycles()] that describe the breeding programme;
#' `control` contains computational settings. Unknown fields are rejected so a
#' misspelling cannot silently change a recommendation.
#'
#' # Exact feasibility mathematics
#'
#' For selection intensity \eqn{i}, genetic covariance \eqn{G}, phenotypic
#' covariance \eqn{P}, and index coefficients \eqn{b}, the response is
#' \deqn{R = iGb / \sqrt{b^\mathsf{T}Pb}}
#' and all attainable one-cycle responses satisfy
#' \deqn{R^\mathsf{T}G^{-1}PG^{-1}R=i^2.}
#' After orienting traits and dividing each response by its estimated genetic SD,
#' DesiredGainR solves the strictly convex programme
#' \deqn{\min_z z^\mathsf{T}Bz \quad \hbox{subject to}\quad z_j\ge m_j.}
#' The exact active-set solver checks the KKT conditions, which are necessary
#' and sufficient because \eqn{B} is positive definite. It also solves the same
#' problem with every \eqn{m_j=1}; scaling its solution to the planned selection
#' intensity gives the mathematically greatest common lower bound attainable
#' across traits in one cycle.
#'
#' # Multi-cycle recommendation
#'
#' The analytical directions and a surrogate-assisted Pareto search are
#' evaluated through forward simulation. For each direction, joint success is
#' the event that every favourable cumulative response reaches its own
#' `minimum_gains` threshold at `control$n_cycles`. The table reports a
#' Jeffreys beta-binomial posterior probability and an exact, one-sided
#' Clopper-Pearson lower confidence bound with Bonferroni family-wise coverage
#' across all evaluated directions. No independence among directions is
#' required for that bound.
#'
#' "Highest" has no unique meaning for several traits without preferences. The
#' least preferential definition used here is maximin: maximise the
#' worst-responding trait after all traits are expressed in genetic SD. The
#' recommendation is ranked by a simultaneous, distribution-free lower
#' confidence bound for the median worst-trait gain, obtained from binomial
#' order-statistic theory. It does not assume normal simulation outcomes.
#'
#' Search and inference are separated in three stages. The adaptive search
#' discovers directions. Once that candidate set is locked, every direction is
#' rerun with new common random numbers; all confidence bounds and both choices
#' use this independent screening sample. Finally, the two locked
#' recommendations are rerun with a third seed stream. `confirmation` is
#' therefore an independent assessment of the reported vectors rather than
#' reuse of either optimisation or screening noise.
#'
#' @param setup A simulation setup from [founder_population()].
#' @param minimum_gains Named non-negative vector of trait-specific minimum
#'   favourable gains, in estimated `G_target` genetic standard deviations.
#' @param lower_is_better Traits for which a reduction is favourable.
#' @param programme Named list of breeding-programme arguments accepted by
#'   [simulate_selection_cycles()]. Common entries are `mating_system`,
#'   `n_parents`, `n_crosses`, and `n_progeny_per_cross`.
#' @param control Named list of computational settings: `n_cycles` (5),
#'   `budget` (60), `n_initial` (`NULL`), `n_replicates` (50),
#'   `n_candidates` (2000), `screening_replicates` (50),
#'   `confirmation_replicates` (200), `confirmation_finalists` (5),
#'   `required_success_probability` (0.8), `search_starts` (3),
#'   `search_angle_tolerance` (10 degrees), `probability_level` (0.95),
#'   `uncertainty` (a named list), `model` (a named list), `seed` (42),
#'   `checkpoint` (`NULL`), and `verbose` (`TRUE`). The uncertainty list accepts
#'   `architecture_draws` (30), `covariance_draws` (0), `genetic_df`, and
#'   `residual_df`. A supported result remains explicitly conditional on the
#'   point covariance unless covariance draws and defensible degrees of freedom
#'   are supplied. Smaller simulation settings are useful for examples, not
#'   decisions.
#'
#' @return An object of class `desiredgainr_gain_suggestion`. The two
#'   recommendations contain a unit-norm desired-gain direction in favourable
#'   genetic-SD space and a signed vector in original trait units. The object
#'   also contains exact one-cycle feasibility, all simulated candidate
#'   summaries, independent confirmation, and full search provenance.
#'
#' @references
#' Brown LD, Cai TT, DasGupta A (2001). Interval estimation for a binomial
#' proportion. *Statistical Science* 16:101--133.
#' \doi{10.1214/ss/1009213286}
#'
#' Casella G, Berger RL (2002). *Statistical Inference*, 2nd ed. Duxbury.
#'
#' Pesek J, Baker RJ (1969). Desired improvement in relation to selection
#' indices. *Canadian Journal of Plant Science* 49:803--804.
#' \doi{10.4141/cjps69-137}
#'
#' Yang W-N, Nelson BL (1991). Using common random numbers and control variates
#' in multiple-comparison procedures. *Operations Research* 39:583--591.
#' \doi{10.1287/opre.39.4.583}
#'
#' @seealso [define_desired_gain_intervals()], [optimize_desired_gains()],
#'   [gain_feasibility()]
#' @export
suggest_desired_gains <- function(
  setup,
  minimum_gains,
  lower_is_better = NULL,
  programme = list(),
  control = list()
) {
  .dgr_require_alphasimr("suggest_desired_gains()")
  if (!inherits(setup, "desiredgainr_sim_setup")) {
    stop("setup must be created by founder_population().", call. = FALSE)
  }
  traits <- setup$trait_cols
  minimum_gains <- .dgr_named_vector(
    minimum_gains, traits, "minimum_gains"
  )
  if (any(!is.finite(minimum_gains)) || any(minimum_gains < 0) ||
    !any(minimum_gains > 0)) {
    stop("minimum_gains must be finite and non-negative, with at least one ",
      "positive trait-specific minimum.",
      call. = FALSE
    )
  }
  direction <- .dgr_direction(traits, lower_is_better)

  programme_defaults <- list(
    mating_system = "self", n_parents = 20L, n_crosses = 50L,
    n_progeny_per_cross = 10L, n_selfing_generations = 3L,
    use_doubled_haploids = FALSE, reestimate_index = TRUE,
    n_clonal_replicates = 1L, n_threads = 1L,
    prediction = list(method = "phenotype")
  )
  programme <- .dgr_merge_named_list(
    programme_defaults, programme, "programme"
  )
  control_defaults <- list(
    n_cycles = 5L, budget = 60L, n_initial = NULL, n_replicates = 50L,
    n_candidates = 2000L, screening_replicates = 50L,
    confirmation_replicates = 200L, confirmation_finalists = 5L,
    required_success_probability = 0.8,
    search_starts = 3L, search_angle_tolerance = 10,
    probability_level = 0.95,
    uncertainty = list(
      architecture_draws = 30L, covariance_draws = 0L,
      genetic_df = NULL, residual_df = NULL
    ),
    model = list(variance_tolerance = 0.15, correlation_tolerance = 0.10),
    seed = 42L, checkpoint = NULL, verbose = TRUE
  )
  control <- .dgr_merge_named_list(control_defaults, control, "control")
  control$uncertainty <- .dgr_merge_named_list(
    control_defaults$uncertainty, control$uncertainty,
    "control$uncertainty"
  )
  control$model <- .dgr_merge_named_list(
    control_defaults$model, control$model, "control$model"
  )
  control$n_replicates <- .dgr_positive_integer(
    control$n_replicates, "control$n_replicates"
  )
  control$confirmation_replicates <- .dgr_positive_integer(
    control$confirmation_replicates, "control$confirmation_replicates"
  )
  control$screening_replicates <- .dgr_positive_integer(
    control$screening_replicates, "control$screening_replicates"
  )
  control$confirmation_finalists <- .dgr_positive_integer(
    control$confirmation_finalists, "control$confirmation_finalists"
  )
  control$search_starts <- .dgr_positive_integer(
    control$search_starts, "control$search_starts"
  )
  control$uncertainty$architecture_draws <- .dgr_non_negative_integer(
    control$uncertainty$architecture_draws,
    "control$uncertainty$architecture_draws"
  )
  control$uncertainty$covariance_draws <- .dgr_non_negative_integer(
    control$uncertainty$covariance_draws,
    "control$uncertainty$covariance_draws"
  )
  if (control$n_replicates < 20L) {
    warning("Fewer than 20 search replicates give weak recommendation ",
      "stability; use the default or more for a decision.",
      call. = FALSE
    )
  }
  if (control$confirmation_replicates < 50L) {
    warning("Fewer than 50 independent confirmation replicates give wide ",
      "uncertainty; use the default or more for a decision.",
      call. = FALSE
    )
  }
  if (control$screening_replicates < 30L) {
    warning("Fewer than 30 independent screening replicates give unstable ",
      "simultaneous rankings; use the default or more for a decision.",
      call. = FALSE
    )
  }
  required_probability <- control$required_success_probability
  if (!is.numeric(required_probability) || length(required_probability) != 1L ||
    !is.finite(required_probability) || required_probability <= 0 ||
    required_probability >= 1) {
    stop("control$required_success_probability must lie strictly between zero and one.",
      call. = FALSE
    )
  }
  if (!is.numeric(control$search_angle_tolerance) ||
    length(control$search_angle_tolerance) != 1L ||
    !is.finite(control$search_angle_tolerance) ||
    control$search_angle_tolerance <= 0 ||
    control$search_angle_tolerance >= 90) {
    stop("control$search_angle_tolerance must be between 0 and 90 degrees.",
      call. = FALSE
    )
  }
  for (field in c("variance_tolerance", "correlation_tolerance")) {
    value <- control$model[[field]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value < 0) {
      stop("control$model$", field,
        " must be one finite non-negative number.",
        call. = FALSE
      )
    }
  }
  level <- control$probability_level
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
    level <= 0 || level >= 1) {
    stop("control$probability_level must lie strictly between zero and one.",
      call. = FALSE
    )
  }
  seed <- .dgr_seed(control$seed)
  programme$n_parents <- .dgr_positive_integer(
    programme$n_parents, "programme$n_parents"
  )

  feasibility <- .dgr_population_feasibility(
    setup, minimum_gains, lower_is_better, programme$n_parents
  )
  model_diagnostic <- .dgr_population_model_diagnostic(
    setup, control$model$variance_tolerance,
    control$model$correlation_tolerance
  )
  population_setup <- setup
  searches <- lapply(seq_len(control$search_starts), function(start) {
    checkpoint <- control$checkpoint
    if (!is.null(checkpoint) && control$search_starts > 1L) {
      checkpoint <- paste0(checkpoint, ".start-", start, ".rds")
    }
    search_call <- c(
      list(
        setup = population_setup, n_cycles = control$n_cycles, mode = "pareto",
        budget = control$budget, n_initial = control$n_initial,
        n_replicates = control$n_replicates, include_diversity = FALSE,
        non_negative = TRUE, desired_gain_units = "genetic_sd",
        probability_level = level, n_candidates = control$n_candidates,
        checkpoint = checkpoint,
        seed = .dgr_shift_seed(seed, (start - 1L) * 100003L),
        verbose = control$verbose, lower_is_better = lower_is_better
      ),
      programme
    )
    do.call(optimize_desired_gains, search_call)
  })
  search <- if (length(searches) == 1L) searches[[1L]] else searches

  scale <- sqrt(diag(population_setup$G_target))
  names(scale) <- traits
  direction_blocks <- lapply(searches, `[[`, "directions")
  raw_directions <- do.call(rbind, direction_blocks)
  direction_keys <- apply(round(raw_directions, 10L), 1L, paste, collapse = ":")
  keep_unique <- !duplicated(direction_keys)
  directions <- raw_directions[keep_unique, , drop = FALSE]
  source_start <- rep(seq_along(direction_blocks), vapply(
    direction_blocks, nrow, integer(1L)
  ))[keep_unique]
  source_evaluation <- sequence(vapply(
    direction_blocks, nrow, integer(1L)
  ))[keep_unique]
  unique_keys <- direction_keys[keep_unique]
  start_candidates <- lapply(direction_blocks, function(block) {
    keys <- apply(round(block, 10L), 1L, paste, collapse = ":")
    unique(match(keys, unique_keys))
  })
  analytical <- rbind(
    feasibility$minimum_cost_direction,
    feasibility$maximum_common_direction
  )
  analytical <- sweep(analytical, 1L, sqrt(rowSums(analytical^2)), "/")
  colnames(analytical) <- traits
  for (i in seq_len(nrow(analytical))) {
    distances <- sqrt(rowSums(sweep(directions, 2L, analytical[i, ], "-")^2))
    if (min(distances) > 1e-8) {
      directions <- rbind(directions, analytical[i, ])
      source_start <- c(source_start, 0L)
      source_evaluation <- c(source_evaluation, i)
    }
  }
  colnames(directions) <- traits

  # The adaptive search has seen its own response streams. Lock all directions
  # before generating a fresh screening stream, otherwise even an exact
  # binomial interval would inherit optimiser selection bias.
  seed_modulus <- .Machine$integer.max - max(
    control$screening_replicates,
    control$confirmation_replicates
  ) - 10L
  screening_seed <- as.integer(
    (abs(as.double(seed)) + 10000001) %% seed_modulus
  )
  screening_replicates <- lapply(seq_len(nrow(directions)), function(index) {
    .dgr_evaluate_direction(
      directions[index, ], population_setup,
      c(
        list(n_cycles = control$n_cycles), programme,
        list(lower_is_better = lower_is_better)
      ),
      control$screening_replicates,
      base_seed = screening_seed,
      include_diversity = FALSE, direction_scale = scale
    )$replicates
  })
  summary <- .dgr_population_candidate_summary(
    directions, screening_replicates, traits, minimum_gains, scale, direction,
    level, screening_seed
  )
  table <- summary$table
  data.table::set(table, j = "search_start", value = source_start)
  data.table::set(table, j = "search_evaluation", value = source_evaluation)
  minimum_order <- order(
    -table$joint_exact_simultaneous_lower,
    -table$joint_probability,
    -table$median_worst_trait_gain
  )
  balanced_order <- order(
    -table$median_worst_exact_simultaneous_lower,
    -table$median_worst_trait_gain,
    -table$mean_worst_trait_gain
  )

  # Confirm several predeclared finalists, not merely the screening winner.
  # The final choice is made only from the new confirmation stream.
  finalist_count <- min(control$confirmation_finalists, nrow(table))
  locked <- unique(c(
    utils::head(minimum_order, finalist_count),
    utils::head(balanced_order, finalist_count)
  ))

  search_diagnostics <- lapply(seq_len(control$search_starts), function(start) {
    available <- start_candidates[[start]]
    half <- unique(match(
      apply(
        round(direction_blocks[[start]][
          seq_len(ceiling(nrow(direction_blocks[[start]]) / 2)), ,
          drop = FALSE
        ], 10L),
        1L, paste,
        collapse = ":"
      ),
      unique_keys
    ))
    full_min <- available[order(
      -table$joint_exact_simultaneous_lower[available],
      -table$joint_probability[available]
    )[1L]]
    full_bal <- available[order(
      -table$median_worst_exact_simultaneous_lower[available],
      -table$median_worst_trait_gain[available]
    )[1L]]
    half_min <- half[order(
      -table$joint_exact_simultaneous_lower[half],
      -table$joint_probability[half]
    )[1L]]
    half_bal <- half[order(
      -table$median_worst_exact_simultaneous_lower[half],
      -table$median_worst_trait_gain[half]
    )[1L]]
    data.frame(
      start = start,
      minimum_candidate = full_min,
      balanced_candidate = full_bal,
      minimum_half_budget_angle = .dgr_direction_angle(
        directions[full_min, ], directions[half_min, ]
      ),
      balanced_half_budget_angle = .dgr_direction_angle(
        directions[full_bal, ], directions[half_bal, ]
      )
    )
  })
  search_diagnostics <- do.call(rbind, search_diagnostics)
  pairwise_max <- function(indices) {
    if (length(indices) < 2L) {
      return(NA_real_)
    }
    max(utils::combn(indices, 2L, function(pair) {
      .dgr_direction_angle(directions[pair[1L], ], directions[pair[2L], ])
    }))
  }
  between_angles <- c(
    pairwise_max(search_diagnostics$minimum_candidate),
    pairwise_max(search_diagnostics$balanced_candidate)
  )
  maximum_between_start_angle <- if (all(is.na(between_angles))) {
    NA_real_
  } else {
    max(between_angles, na.rm = TRUE)
  }
  maximum_half_budget_angle <- max(
    search_diagnostics$minimum_half_budget_angle,
    search_diagnostics$balanced_half_budget_angle,
    na.rm = TRUE
  )
  search_status <- if (control$search_starts < 2L) {
    "not_assessed"
  } else if (maximum_between_start_angle <= control$search_angle_tolerance &&
    maximum_half_budget_angle <= control$search_angle_tolerance) {
    "resolved"
  } else {
    "unresolved"
  }
  search_stability <- list(
    status = search_status,
    angle_tolerance_degrees = control$search_angle_tolerance,
    maximum_between_start_angle = maximum_between_start_angle,
    maximum_half_vs_full_budget_angle = maximum_half_budget_angle,
    by_start = search_diagnostics
  )

  confirmation_seed <- as.integer(
    (abs(as.double(seed)) + 20000003) %% seed_modulus
  )
  confirmation_replicates <- lapply(locked, function(index) {
    .dgr_evaluate_direction(
      directions[index, ], population_setup,
      c(
        list(n_cycles = control$n_cycles), programme,
        list(lower_is_better = lower_is_better)
      ),
      control$confirmation_replicates,
      base_seed = confirmation_seed,
      include_diversity = FALSE, direction_scale = scale
    )$replicates
  })
  confirmation_summary <- .dgr_population_candidate_summary(
    directions[locked, , drop = FALSE], confirmation_replicates,
    traits, minimum_gains, scale, direction, level,
    confirmation_seed
  )$table
  data.table::set(confirmation_summary, j = "candidate", value = locked)
  minimum_row <- order(
    -confirmation_summary$joint_exact_simultaneous_lower,
    -confirmation_summary$joint_probability,
    -confirmation_summary$median_worst_trait_gain
  )[1L]
  balanced_row <- order(
    -confirmation_summary$median_worst_exact_simultaneous_lower,
    -confirmation_summary$median_worst_trait_gain,
    -confirmation_summary$mean_worst_trait_gain
  )[1L]
  minimum_index <- confirmation_summary$candidate[minimum_row]
  balanced_index <- confirmation_summary$candidate[balanced_row]

  outer_validation <- .dgr_population_outer_validation(
    setup, directions[locked, , drop = FALSE], traits, minimum_gains,
    scale, direction, programme, control$n_cycles, control$uncertainty,
    level, seed
  )
  for (component in c("architecture", "covariance")) {
    if (identical(outer_validation[[component]]$status, "assessed")) {
      data.table::set(
        outer_validation[[component]]$results,
        j = "candidate", value = locked
      )
    }
  }

  validation_for <- function(component, index) {
    validation <- outer_validation[[component]]
    if (!identical(validation$status, "assessed")) {
      return(validation)
    }
    row <- validation$results[validation$results$candidate == index]
    list(
      status = if (validation$failed_draws > 0L) {
        "failed"
      } else {
        .dgr_support_status(row, required_probability)
      },
      requested_draws = validation$requested_draws,
      usable_draws = validation$usable_draws,
      failed_draws = validation$failed_draws,
      summary = as.list(row)
    )
  }
  make_recommendation <- function(index) {
    favourable <- directions[index, ]
    names(favourable) <- traits
    confirmation_row <- confirmation_summary[
      confirmation_summary[["candidate"]] == index
    ]
    decision_status <- .dgr_support_status(
      confirmation_row, required_probability
    )
    architecture <- validation_for("architecture", index)
    covariance <- validation_for("covariance", index)
    operational_status <- if (!identical(model_diagnostic$status, "adequate")) {
      paste0("model_", model_diagnostic$status)
    } else if (!identical(search_status, "resolved")) {
      paste0("search_", search_status)
    } else if (!identical(decision_status, "supported")) {
      decision_status
    } else if (identical(architecture$status, "not_supported")) {
      "architecture_not_supported"
    } else if (identical(architecture$status, "uncertain")) {
      "architecture_uncertain"
    } else if (!identical(architecture$status, "supported")) {
      paste0("architecture_", architecture$status)
    } else if (identical(covariance$status, "not_supported")) {
      "covariance_not_supported"
    } else if (identical(covariance$status, "uncertain")) {
      "covariance_uncertain"
    } else if (identical(covariance$status, "not_assessed")) {
      "supported_conditional_on_covariance"
    } else if (!identical(covariance$status, "supported")) {
      paste0("covariance_", covariance$status)
    } else {
      "supported"
    }
    list(
      candidate = index,
      desired_gain_direction = favourable,
      desired_gains_trait_units = favourable * scale * direction,
      screening = as.list(table[index]),
      confirmation = as.list(confirmation_row),
      decision_status = decision_status,
      operational_status = operational_status,
      recommended = operational_status %in% c(
        "supported", "supported_conditional_on_covariance"
      ),
      architecture_validation = architecture,
      covariance_validation = covariance
    )
  }
  minimum_recommendation <- make_recommendation(minimum_index)
  balanced_recommendation <- make_recommendation(balanced_index)

  result <- list(
    minimum_gains = minimum_gains,
    gain_units = "genetic_sd",
    lower_is_better = lower_is_better,
    minimum_recommendation = minimum_recommendation,
    maximum_balanced_recommendation = balanced_recommendation,
    analytical_feasibility = feasibility,
    model_diagnostic = model_diagnostic,
    candidate_results = table,
    screening_replicate_outcomes = screening_replicates,
    confirmation = confirmation_summary,
    outer_validation = outer_validation,
    search_stability = search_stability,
    search = search,
    programme = programme,
    control = control,
    interpretation = paste(
      "The desired-gain vectors are directions, not promised responses.",
      "The minimum recommendation targets the breeder's trait-specific",
      "thresholds. The maximum-balanced recommendation uses a symmetric",
      "maximin rule in genetic-SD units because no economic or preference",
      "weights were supplied. A vector is recommended only when its exact",
      "simultaneous confirmation bound, model diagnostic, search stability",
      "and requested outer validations support it. Covariance support remains",
      "conditional unless degrees of freedom were supplied."
    )
  )
  class(result) <- c("desiredgainr_gain_suggestion", "list")
  result
}

#' @export
print.desiredgainr_gain_suggestion <- function(x, ...) {
  cat("<desiredgainr_gain_suggestion>\n")
  cat("  Trait-specific minima (favourable genetic SD):\n")
  print(round(x$minimum_gains, 4L))
  cat(sprintf(
    "  Exact one-cycle feasibility of all minima: %s\n",
    if (isTRUE(x$analytical_feasibility$minimum_feasible_in_one_cycle)) {
      "yes"
    } else {
      "no"
    }
  ))
  cat(sprintf(
    "  Exact maximum common one-cycle gain: %.3f genetic SD\n",
    x$analytical_feasibility$maximum_common_gain
  ))
  cat("  Simulation-model calibration:", x$model_diagnostic$status, "\n")
  cat("  Multi-start search stability:", x$search_stability$status, "\n")
  if (length(x$analytical_feasibility$limiting_traits)) {
    cat(
      "  Limiting trait(s):",
      paste(x$analytical_feasibility$limiting_traits, collapse = ", "), "\n"
    )
  }
  show_recommendation <- function(label, recommendation) {
    cat("  ", label, " candidate [", recommendation$operational_status,
      "]:\n",
      sep = ""
    )
    print(round(recommendation$desired_gain_direction, 4L))
    confirmation <- recommendation$confirmation
    cat(sprintf(
      "    Independent joint probability of all minima: %.1f%%",
      100 * confirmation$joint_probability
    ))
    cat(sprintf(
      " (%.1f%% to %.1f%%)\n",
      100 * confirmation$joint_probability_lower,
      100 * confirmation$joint_probability_upper
    ))
    cat(sprintf(
      "    Exact simultaneous probability bounds: %.1f%% to %.1f%%\n",
      100 * confirmation$joint_exact_simultaneous_lower,
      100 * confirmation$joint_exact_simultaneous_upper
    ))
    cat(sprintf(
      "    Independent median worst-trait gain: %.3f genetic SD\n",
      confirmation$median_worst_trait_gain
    ))
  }
  show_recommendation("Minimum-attainment", x$minimum_recommendation)
  show_recommendation(
    "Maximum-balanced", x$maximum_balanced_recommendation
  )
  cat("  A candidate labelled uncertain/unresolved is not a recommendation.\n")
  cat("  Probabilities are conditional on the declared programme and model.\n")
  invisible(x)
}
