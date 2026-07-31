# Surrogate-assisted search for the desired-gain direction that performs best
# over several breeding cycles.
#
# The simulation is the objective function and is never replaced by a cheaper
# approximation. The Gaussian process only decides where to spend the next
# simulation, and it retains an exploration term everywhere, so no region of
# the search space can be permanently excluded by the surrogate.

#' Halton quasi-random sequence
#'
#' @param n Number of points.
#' @param d Number of dimensions.
#' @param skip Number of leading points discarded, which improves the spread of
#'   the first few points.
#'
#' @return An `n` by `d` matrix of points in the unit hypercube.
#' @noRd
.dgr_halton <- function(n, d, skip = 20L) {
  primes <- c(2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47)
  if (d > length(primes)) {
    stop("Quasi-random initialisation supports at most ", length(primes),
         " dimensions.", call. = FALSE)
  }
  van_der_corput <- function(index, base) {
    result <- numeric(length(index))
    fraction <- 1
    remaining <- index
    while (any(remaining > 0)) {
      fraction <- fraction / base
      result <- result + fraction * (remaining %% base)
      remaining <- remaining %/% base
    }
    result
  }
  index <- seq_len(n) + skip
  vapply(seq_len(d), function(j) van_der_corput(index, primes[j]), numeric(n))
}

#' Map unit-hypercube points onto directions
#'
#' Only the direction of a desired-gain vector affects the index, because the
#' attainable magnitude is fixed by selection intensity. The search space is
#' therefore the unit sphere rather than the whole of trait space.
#'
#' @param U Matrix of points in the unit hypercube.
#' @param non_negative Whether to restrict to the non-negative orthant.
#'
#' @return A matrix of unit-norm rows.
#' @noRd
.dgr_directions_from_cube <- function(U, non_negative = TRUE) {
  Z <- stats::qnorm(pmin(pmax(U, 1e-6), 1 - 1e-6))
  if (isTRUE(non_negative)) Z <- abs(Z)
  norms <- sqrt(rowSums(Z^2))
  norms[norms == 0] <- 1
  sweep(Z, 1L, norms, "/")
}

#' Identify the non-dominated rows of an objective matrix
#'
#' @param objectives Matrix with one row per solution and one column per
#'   objective, all to be maximised.
#'
#' @return A logical vector marking the non-dominated rows.
#' @noRd
.dgr_non_dominated <- function(objectives) {
  n <- nrow(objectives)
  keep <- rep(TRUE, n)
  for (i in seq_len(n)) {
    if (!keep[i]) next
    dominated <- apply(objectives, 1L, function(other) {
      all(other >= objectives[i, ]) && any(other > objectives[i, ])
    })
    if (any(dominated)) keep[i] <- FALSE
  }
  keep
}

#' Evaluate one desired-gain direction by simulation
#'
#' @param direction Unit-norm desired-gain direction.
#' @param setup Simulation setup.
#' @param arguments List of arguments forwarded to
#'   [simulate_selection_cycles()].
#' @param n_replicates Independent simulation replicates.
#' @param base_seed Seed of the first replicate. The same seed sequence is used
#'   for every direction, so that comparisons between directions share their
#'   stochasticity.
#' @param include_diversity Whether to append the negated mean relationship as
#'   an additional objective.
#'
#' @return A list with the mean objective vector and its standard error.
#' @noRd
.dgr_evaluate_direction <- function(
    direction, setup, arguments, n_replicates, base_seed, include_diversity
) {
  trait_cols <- setup$trait_cols
  named_direction <- stats::setNames(as.numeric(direction), trait_cols)
  replicate_values <- vapply(seq_len(n_replicates), function(replicate) {
    simulation <- do.call(
      simulate_selection_cycles,
      c(
        list(
          setup = setup,
          desired_gains = named_direction,
          seed = base_seed + replicate - 1L
        ),
        arguments
      )
    )
    outcome <- simulation$cumulative_gain[trait_cols]
    if (isTRUE(include_diversity)) {
      outcome <- c(outcome, diversity = -simulation$final_relationship)
    }
    as.numeric(outcome)
  }, numeric(length(trait_cols) + as.integer(include_diversity)))

  if (is.null(dim(replicate_values))) {
    replicate_values <- matrix(replicate_values, ncol = n_replicates)
  }
  list(
    mean = rowMeans(replicate_values),
    sd = if (n_replicates > 1L) {
      apply(replicate_values, 1L, stats::sd) / sqrt(n_replicates)
    } else {
      rep(NA_real_, nrow(replicate_values))
    }
  )
}

#' Search for the desired-gain direction giving the best multi-cycle outcome
#'
#' A single-cycle response calculation cannot distinguish desired-gain
#' directions beyond what the achievable-response ellipsoid already states, and
#' [gain_feasibility()] provides that answer exactly and without simulation.
#' Over several cycles the ranking can change, however, because truncation
#' selection erodes the genetic variance that response depends on, drift
#' accumulates in a finite population, and an antagonistic correlation can
#' drive a secondary trait past an unacceptable level. This function therefore
#' searches for the direction that performs best once those processes are
#' represented.
#'
#' @details
#' # The simulation is the objective, never the surrogate
#'
#' Each candidate direction is evaluated by running
#' [simulate_selection_cycles()]. A Gaussian process is fitted to the
#' accumulated results, but it is used only to decide where the next simulation
#' should be spent. It never filters, screens, or replaces an evaluation.
#' Because the acquisition function retains an exploration term across the
#' whole search space, no region can be permanently excluded by the surrogate,
#' and a direction that the surrogate expects to be poor can still be visited.
#'
#' # Search space
#'
#' Only the direction of a desired-gain vector affects the index, since the
#' attainable magnitude is fixed by selection intensity. The domain is
#' therefore the unit sphere, with `p - 1` free parameters for `p` traits,
#' which is small enough to cover densely for the trait numbers breeders use.
#' By default the search is restricted to the non-negative orthant, meaning
#' improvement is sought in every trait; set `non_negative = FALSE` to admit
#' directions that deliberately concede ground on a trait.
#'
#' # Ranking modes
#'
#' \describe{
#'   \item{`"pareto"`}{Returns the non-dominated set of multi-cycle outcomes
#'     and the direction generating each. This is the default, because it is
#'     the only mode that does not require the economic weights that breeders
#'     find hardest to state; choosing a point on a frontier is an easier
#'     judgement, and it is weight elicitation by revealed preference. The
#'     search uses randomised augmented-Chebyshev scalarisation, which
#'     converges on the frontier while needing only single-objective
#'     improvement.}
#'   \item{`"economic"`}{Maximises `sum(economic_weights * cumulative_gain)`.}
#'   \item{`"target"`}{Minimises the distance between the cumulative gain and a
#'     stated absolute target, reporting the per-trait shortfall.}
#'   \item{`"constrained"`}{Maximises the gain in `focal_trait` subject to
#'     `gain_floors` on the remaining traits, weighting improvement by the
#'     probability that every floor is met.}
#' }
#'
#' # Noise and the frontier
#'
#' Simulation output is stochastic. Every direction is evaluated with the same
#' sequence of seeds, so that comparisons between directions share their
#' stochasticity, which removes far more comparison variance than increasing
#' the replicate count. The reported frontier is computed from Gaussian-process
#' posterior means rather than from raw replicate averages, because a frontier
#' built from raw draws is populated by fortunate runs.
#'
#' # Interpreting the result
#'
#' No single best direction is returned. The optimum is conditional on the
#' supplied genetic covariance matrix, the founder germplasm and the programme
#' parameters, so the result reports a stability region: every direction whose
#' posterior outcome lies within `stability_tolerance` of the best. Choose
#' within that region on other grounds.
#'
#' @param setup An object from [founder_population()].
#' @param n_cycles Number of selection cycles per evaluation.
#' @param mode Ranking mode. See Details.
#' @param budget Total number of simulated directions, including the initial
#'   design.
#' @param n_initial Size of the quasi-random initial design. Defaults to ten
#'   times the number of free parameters.
#' @param n_replicates Simulation replicates per direction.
#' @param economic_weights Named weights required by `mode = "economic"`.
#' @param target_gains Named absolute targets required by `mode = "target"`.
#' @param focal_trait,gain_floors Objective and constraints required by
#'   `mode = "constrained"`.
#' @param include_diversity Whether to treat diversity, measured as the negated
#'   mean relationship among selected parents, as an additional Pareto
#'   objective. Family balancing and coancestry control can cost a large share
#'   of nominal gain, so the trade-off is reported rather than assumed away.
#' @param non_negative Whether to restrict the search to directions seeking
#'   improvement in every trait.
#' @param stability_tolerance Relative tolerance defining the stability region.
#' @param n_candidates Size of the quasi-random pool over which the acquisition
#'   function is maximised at each iteration.
#' @param checkpoint Optional file path. When supplied, the accumulated
#'   evaluations are written after every simulation and reloaded automatically
#'   on a subsequent call, so that an interrupted search is not lost.
#' @param seed Random seed. The caller's random number generator state is
#'   restored on exit.
#' @param verbose Whether to report progress.
#' @param ... Further arguments passed to [simulate_selection_cycles()], such
#'   as `mating_system`, `n_parents` and `lower_is_better`.
#'
#' @return An object of class `desiredgainr_optimisation` containing the
#'   evaluated directions and outcomes, the posterior-mean Pareto set, the
#'   recommended region, and the search diagnostics.
#'
#' @seealso [simulate_selection_cycles()], [gain_feasibility()]
#' @export
optimize_desired_gains <- function(
    setup,
    n_cycles = 5L,
    mode = c("pareto", "economic", "target", "constrained"),
    budget = 60L,
    n_initial = NULL,
    n_replicates = 3L,
    economic_weights = NULL,
    target_gains = NULL,
    focal_trait = NULL,
    gain_floors = NULL,
    include_diversity = TRUE,
    non_negative = TRUE,
    stability_tolerance = 0.05,
    n_candidates = 2000L,
    checkpoint = NULL,
    seed = 42L,
    verbose = TRUE,
    ...
) {
  mode <- match.arg(mode)
  .dgr_require_alphasimr("optimize_desired_gains()")
  if (!inherits(setup, "desiredgainr_sim_setup")) {
    stop("setup must be created by founder_population().", call. = FALSE)
  }
  trait_cols <- setup$trait_cols
  p <- length(trait_cols)
  if (p < 2L) {
    stop("At least two traits are required to search a desired-gain direction.",
         call. = FALSE)
  }
  budget <- .dgr_positive_integer(budget, "budget")
  n_replicates <- .dgr_positive_integer(n_replicates, "n_replicates")
  n_candidates <- .dgr_positive_integer(n_candidates, "n_candidates")
  if (is.null(n_initial)) n_initial <- min(budget, 10L * (p - 1L))
  n_initial <- .dgr_positive_integer(n_initial, "n_initial")
  if (n_initial > budget) {
    stop("n_initial cannot exceed budget.", call. = FALSE)
  }

  objective_names <- trait_cols
  if (isTRUE(include_diversity) && mode == "pareto") {
    objective_names <- c(objective_names, "diversity")
  } else {
    include_diversity <- FALSE
  }

  if (mode == "economic") {
    economic_weights <- .dgr_named_vector(
      economic_weights, trait_cols, "economic_weights"
    )
  } else if (mode == "target") {
    target_gains <- .dgr_named_vector(
      target_gains, trait_cols, "target_gains"
    )
  } else if (mode == "constrained") {
    if (is.null(focal_trait) || !focal_trait %in% trait_cols) {
      stop("focal_trait must name one of the traits.", call. = FALSE)
    }
    other_traits <- setdiff(trait_cols, focal_trait)
    gain_floors <- .dgr_named_vector(
      gain_floors, other_traits, "gain_floors"
    )
  }

  if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    stats::runif(1L)
  }
  entry_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit(assign(".Random.seed", entry_seed, envir = globalenv()), add = TRUE)
  set.seed(as.integer(seed))

  simulation_arguments <- c(list(n_cycles = n_cycles), list(...))
  candidate_pool <- .dgr_directions_from_cube(
    .dgr_halton(n_candidates, p, skip = 500L), non_negative
  )
  colnames(candidate_pool) <- trait_cols

  evaluated_directions <- NULL
  evaluated_objectives <- NULL
  if (!is.null(checkpoint) && file.exists(checkpoint)) {
    saved <- readRDS(checkpoint)
    if (identical(saved$trait_cols, trait_cols) &&
        identical(saved$objective_names, objective_names)) {
      evaluated_directions <- saved$directions
      evaluated_objectives <- saved$objectives
      if (isTRUE(verbose)) {
        message(sprintf(
          "Resuming from %d checkpointed evaluation(s).",
          nrow(evaluated_directions)
        ))
      }
    } else {
      warning("The checkpoint does not match this problem and was ignored.",
              call. = FALSE)
    }
  }

  evaluate_and_store <- function(direction) {
    outcome <- .dgr_evaluate_direction(
      direction, setup, simulation_arguments, n_replicates,
      base_seed = as.integer(seed) + 1000L, include_diversity
    )
    evaluated_directions <<- rbind(evaluated_directions, direction)
    evaluated_objectives <<- rbind(evaluated_objectives, outcome$mean)
    if (!is.null(checkpoint)) {
      saveRDS(
        list(
          trait_cols = trait_cols, objective_names = objective_names,
          directions = evaluated_directions, objectives = evaluated_objectives
        ),
        checkpoint
      )
    }
    invisible(NULL)
  }

  already <- if (is.null(evaluated_directions)) 0L else nrow(evaluated_directions)
  if (already < n_initial) {
    initial <- .dgr_directions_from_cube(
      .dgr_halton(n_initial, p), non_negative
    )
    for (i in seq(already + 1L, n_initial)) {
      if (isTRUE(verbose)) {
        message(sprintf("Initial design %d of %d.", i, n_initial))
      }
      evaluate_and_store(initial[i, , drop = FALSE])
    }
  }

  scalarise <- function(objectives, weights = NULL) {
    switch(
      mode,
      economic = as.numeric(objectives[, trait_cols, drop = FALSE] %*%
                              economic_weights),
      target = {
        scale <- pmax(abs(target_gains), 1e-8)
        deviation <- sweep(
          objectives[, trait_cols, drop = FALSE], 2L, target_gains, "-"
        )
        -sqrt(rowSums(sweep(deviation, 2L, scale, "/")^2))
      },
      constrained = objectives[, focal_trait],
      pareto = {
        # Augmented Chebyshev scalarisation on objectives normalised to the
        # unit interval, drawn afresh each iteration so that the search sweeps
        # the frontier rather than one corner of it.
        ranges <- apply(objectives, 2L, range)
        spread <- pmax(ranges[2L, ] - ranges[1L, ], 1e-12)
        normalised <- sweep(
          sweep(objectives, 2L, ranges[1L, ], "-"), 2L, spread, "/"
        )
        weighted <- sweep(normalised, 2L, weights, "*")
        apply(weighted, 1L, min) + 0.05 * rowSums(weighted)
      }
    )
  }

  iteration_log <- list()
  while (nrow(evaluated_directions) < budget) {
    colnames(evaluated_objectives) <- objective_names
    chebyshev_weights <- if (mode == "pareto") {
      draw <- stats::rexp(length(objective_names))
      draw / sum(draw)
    } else {
      NULL
    }
    response <- scalarise(evaluated_objectives, chebyshev_weights)

    model <- .dgr_gp_fit(evaluated_directions, response)
    prediction <- .dgr_gp_predict(model, candidate_pool)
    acquisition <- .dgr_expected_improvement(
      prediction$mean, prediction$sd, max(response)
    )

    if (mode == "constrained") {
      for (trait in names(gain_floors)) {
        constraint_model <- .dgr_gp_fit(
          evaluated_directions, evaluated_objectives[, trait]
        )
        constraint_prediction <- .dgr_gp_predict(
          constraint_model, candidate_pool
        )
        acquisition <- acquisition * .dgr_probability_feasible(
          constraint_prediction$mean, constraint_prediction$sd,
          gain_floors[[trait]]
        )
      }
    }

    chosen <- which.max(acquisition)
    iteration_log[[length(iteration_log) + 1L]] <- data.table::data.table(
      evaluation = nrow(evaluated_directions) + 1L,
      acquisition = acquisition[chosen],
      best_response = max(response),
      gp_nlml = model$nlml
    )
    if (isTRUE(verbose)) {
      message(sprintf(
        "Evaluation %d of %d, expected improvement %.4g.",
        nrow(evaluated_directions) + 1L, budget, acquisition[chosen]
      ))
    }
    evaluate_and_store(candidate_pool[chosen, , drop = FALSE])
  }

  colnames(evaluated_directions) <- trait_cols
  colnames(evaluated_objectives) <- objective_names

  # The frontier and the recommendation are read from posterior means rather
  # than raw replicate averages, so that fortunate runs do not populate them.
  smoothed <- vapply(objective_names, function(objective) {
    model <- .dgr_gp_fit(evaluated_directions, evaluated_objectives[, objective])
    .dgr_gp_predict(model, evaluated_directions)$mean
  }, numeric(nrow(evaluated_directions)))
  if (is.null(dim(smoothed))) {
    smoothed <- matrix(smoothed, ncol = length(objective_names))
  }
  colnames(smoothed) <- objective_names

  results <- data.table::data.table(
    evaluation = seq_len(nrow(evaluated_directions))
  )
  for (trait in trait_cols) {
    data.table::set(results, j = paste0("d_", trait),
                    value = evaluated_directions[, trait])
  }
  for (objective in objective_names) {
    data.table::set(results, j = paste0("observed_", objective),
                    value = evaluated_objectives[, objective])
    data.table::set(results, j = paste0("posterior_", objective),
                    value = smoothed[, objective])
  }

  pareto_flag <- .dgr_non_dominated(smoothed)
  data.table::set(results, j = "pareto_optimal", value = pareto_flag)

  final_scalar <- if (mode == "pareto") {
    rep(NA_real_, nrow(smoothed))
  } else {
    scalarise(smoothed)
  }
  stability <- rep(NA, nrow(smoothed))
  recommended <- NULL
  if (mode != "pareto") {
    best_value <- max(final_scalar)
    threshold <- best_value - stability_tolerance * abs(best_value)
    stability <- final_scalar >= threshold
    recommended <- evaluated_directions[which.max(final_scalar), ]
    names(recommended) <- trait_cols
    data.table::set(results, j = "objective", value = final_scalar)
    data.table::set(results, j = "within_stability_region", value = stability)
  }

  result <- list(
    mode = mode,
    trait_cols = trait_cols,
    objective_names = objective_names,
    n_cycles = n_cycles,
    n_replicates = n_replicates,
    budget = budget,
    n_initial = n_initial,
    results = results,
    directions = evaluated_directions,
    observed_objectives = evaluated_objectives,
    posterior_objectives = smoothed,
    pareto_set = results[pareto_optimal == TRUE],
    recommended_direction = recommended,
    stability_proportion = if (mode == "pareto") NA_real_ else mean(stability),
    stability_tolerance = stability_tolerance,
    iteration_log = data.table::rbindlist(iteration_log),
    simulation_arguments = simulation_arguments,
    seed = as.integer(seed),
    interpretation = paste(
      "Only the direction of a desired-gain vector is optimised, because the",
      "attainable magnitude is fixed by selection intensity. The reported",
      "optimum is conditional on the supplied genetic covariance, founder",
      "germplasm and programme parameters, so the stability region rather",
      "than the single best direction should guide the decision."
    )
  )
  class(result) <- c("desiredgainr_optimisation", "list")
  result
}

#' @export
print.desiredgainr_optimisation <- function(x, ...) {
  cat("<desiredgainr_optimisation>\n")
  cat(sprintf(
    "  Mode: %s   Cycles: %d   Evaluations: %d (%d replicates each)\n",
    x$mode, x$n_cycles, nrow(x$directions), x$n_replicates
  ))
  cat("  Objectives:", paste(x$objective_names, collapse = ", "), "\n")
  if (identical(x$mode, "pareto")) {
    cat(sprintf(
      "  Non-dominated directions: %d of %d\n",
      nrow(x$pareto_set), nrow(x$directions)
    ))
    cat("  No single direction is recommended; choose a point on the",
        "frontier.\n")
  } else {
    cat("  Recommended direction:\n")
    print(round(x$recommended_direction, 4L))
    cat(sprintf(
      "  Directions within %.0f%% of the optimum: %.1f%%\n",
      100 * x$stability_tolerance, 100 * x$stability_proportion
    ))
  }
  invisible(x)
}
