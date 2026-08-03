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
      " dimensions.",
      call. = FALSE
    )
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

# Map hypercube points into a breeder-defined box, then retain only direction.
.dgr_directions_from_intervals <- function(U, intervals) {
  lower <- intervals[, "lower"]
  upper <- intervals[, "upper"]
  points <- sweep(U, 2L, upper - lower, "*")
  points <- sweep(points, 2L, lower, "+")
  norms <- sqrt(rowSums(points^2))
  if (any(norms <= .Machine$double.eps)) {
    stop("desired_gain_intervals generated an all-zero desired-gain vector. ",
      "At least one trait must be bounded away from zero.",
      call. = FALSE
    )
  }
  sweep(points, 1L, norms, "/")
}

.dgr_validate_gain_intervals <- function(x, trait_cols, non_negative) {
  if (is.null(x)) {
    return(NULL)
  }
  if (inherits(x, "desiredgainr_gain_intervals")) {
    intervals <- as.matrix(x[, c("lower", "upper"), drop = FALSE])
    rownames(intervals) <- x$trait
  } else {
    intervals <- as.matrix(x)
    if (is.null(rownames(intervals)) ||
      !all(c("lower", "upper") %in% colnames(intervals))) {
      stop("desired_gain_intervals must come from ",
        "define_desired_gain_intervals(), or be a trait-named matrix/data ",
        "frame with lower and upper columns in the favourable direction.",
        call. = FALSE
      )
    }
    intervals <- intervals[, c("lower", "upper"), drop = FALSE]
  }
  if (!setequal(rownames(intervals), trait_cols)) {
    stop("desired_gain_intervals must contain exactly these traits: ",
      paste(trait_cols, collapse = ", "), ".",
      call. = FALSE
    )
  }
  intervals <- intervals[trait_cols, , drop = FALSE]
  storage.mode(intervals) <- "double"
  if (any(!is.finite(intervals)) || any(intervals[, "lower"] >
    intervals[, "upper"])) {
    stop("desired_gain_intervals must have finite bounds with lower <= upper.",
      call. = FALSE
    )
  }
  if (isTRUE(non_negative) && any(intervals < 0)) {
    stop("Negative favourable-direction interval bounds require ",
      "non_negative = FALSE.",
      call. = FALSE
    )
  }
  if (all(intervals[, "upper"] == 0)) {
    stop("desired_gain_intervals cannot restrict every trait to zero.",
      call. = FALSE
    )
  }
  if (!any(intervals[, "upper"] > intervals[, "lower"])) {
    stop("desired_gain_intervals fixes a single direction; use ",
      "simulate_selection_cycles() for that vector instead of optimisation.",
      call. = FALSE
    )
  }
  intervals
}

.dgr_representative_gains <- function(directions, intervals) {
  directions <- as.matrix(directions)
  if (is.null(intervals)) {
    return(directions)
  }
  representatives <- t(apply(directions, 1L, function(direction) {
    lower_scale <- 0
    upper_scale <- Inf
    for (trait in seq_along(direction)) {
      if (direction[trait] > .Machine$double.eps) {
        lower_scale <- max(
          lower_scale, intervals[trait, "lower"] / direction[trait]
        )
        upper_scale <- min(
          upper_scale, intervals[trait, "upper"] / direction[trait]
        )
      } else if (direction[trait] < -.Machine$double.eps) {
        lower_scale <- max(
          lower_scale, intervals[trait, "upper"] / direction[trait]
        )
        upper_scale <- min(
          upper_scale, intervals[trait, "lower"] / direction[trait]
        )
      } else if (intervals[trait, "lower"] > 0 ||
        intervals[trait, "upper"] < 0) {
        stop("A searched direction does not intersect desired_gain_intervals.",
          call. = FALSE
        )
      }
    }
    if (lower_scale > upper_scale + 1e-10 || !is.finite(lower_scale)) {
      stop("A searched direction does not intersect desired_gain_intervals.",
        call. = FALSE
      )
    }
    scale <- if (is.finite(upper_scale)) {
      (lower_scale + upper_scale) / 2
    } else {
      max(lower_scale, 1)
    }
    direction * scale
  }))
  dimnames(representatives) <- dimnames(directions)
  representatives
}

.dgr_desired_gain_scale <- function(setup, units) {
  G <- setup$G_target
  scale <- switch(units,
    trait = rep(1, nrow(G)),
    genetic_sd = sqrt(diag(G)),
    phenotypic_sd = {
      if (!is.null(setup$residual_covariance)) {
        sqrt(diag(G + setup$residual_covariance))
      } else {
        sqrt(diag(G) / as.numeric(setup$h2))
      }
    }
  )
  names(scale) <- setup$trait_cols
  if (any(!is.finite(scale)) || any(scale <= 0)) {
    stop("The requested desired-gain units cannot be derived from this setup.",
      call. = FALSE
    )
  }
  scale
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

# Fingerprint every input that determines an optimiser result.
#
# A checkpoint that matches only on trait names can resume evaluations produced
# under different founders, a different genetic covariance, a different mating
# system or a different cycle count, and present them as this run's frontier.
# The fields below are everything a resumed evaluation must share for the
# pooled results to describe one problem.
#
# Written as plain comments, not roxygen: a roxygen block placed immediately
# after another block's tag merges with it, and @noRd followed by text is an
# error. Internal helpers need no manual page.
.dgr_optimiser_fingerprint <- function(
  setup, trait_cols, objective_names, simulation_arguments, n_replicates,
  seed, include_diversity, non_negative, mode, objective_parameters,
  search_parameters
) {
  digest_of <- function(x) {
    # Hash the complete serialization. Truncating a large founder object made
    # changes at late loci invisible and allowed incompatible checkpoints to
    # collide.
    digest::digest(x, algo = "sha256", serialize = TRUE)
  }
  list(
    trait_cols = trait_cols,
    objective_names = objective_names,
    genetic_target = digest_of(round(unname(as.matrix(setup$G_target)), 10L)),
    genetic_realised = digest_of(
      if (is.null(setup$G_realised)) {
        NULL
      } else {
        round(unname(as.matrix(setup$G_realised)), 10L)
      }
    ),
    heritability = digest_of(round(unname(setup$h2), 10L)),
    heritability_type = setup$heritability_type,
    founder_haplotypes = digest_of(setup$founders$haplotypes),
    genetic_map = digest_of(setup$founders$gen_map),
    n_qtl_per_chromosome = setup$n_qtl_per_chromosome,
    marker_panel = digest_of(sort(as.character(setup$marker_panel))),
    dominance = setup$dominance,
    setup_seed = setup$seed,
    simulation_arguments = digest_of(
      simulation_arguments[order(names(simulation_arguments))]
    ),
    n_replicates = n_replicates,
    seed = seed,
    include_diversity = isTRUE(include_diversity),
    non_negative = isTRUE(non_negative),
    mode = mode,
    objective_parameters = digest_of(objective_parameters),
    search_parameters = digest_of(search_parameters),
    desiredgainr_version = as.character(
      utils::packageVersion("DesiredGainR")
    ),
    alphasimr_version = as.character(utils::packageVersion("AlphaSimR")),
    r_version = paste(R.version$major, R.version$minor, sep = ".")
  )
}

#' Compute Monte Carlo error after applying the actual scalarisation
#'
#' The surrogate is fitted to a scalar summary of the objective vector, so the
#' replicate-level objective vectors are scalarised first, preserving their
#' covariance and the non-linearity of target-distance and Chebyshev modes.
#'
#' @noRd
.dgr_scalarised_error <- function(
  replicates, scalarise, weights = NULL, ranges = NULL
) {
  if (is.null(replicates) || !length(replicates)) {
    return(NULL)
  }
  vapply(replicates, function(values) {
    values <- as.matrix(values)
    if (ncol(values) < 2L) {
      return(NA_real_)
    }
    scalar_replicates <- scalarise(
      t(values),
      weights = weights, ranges = ranges
    )
    stats::sd(scalar_replicates) / sqrt(length(scalar_replicates))
  }, numeric(1L))
}

# Convert replicate outcomes into interval-attainment probabilities.
#
# A success is the joint event that every favourable, standardised cumulative
# response reaches its breeder-defined lower bound. With s successes in n
# independent simulation replicates, Jeffreys' invariant binomial prior gives
# Beta(s + 1/2, n - s + 1/2). The resulting interval stays finite at s = 0 and
# s = n, unlike the plug-in normal interval.
.dgr_interval_attainment <- function(
  replicates, trait_cols, intervals, outcome_scale, favourable_direction,
  level = 0.95, augmentation = 0.05, bootstrap_seed = 1L
) {
  alpha_tail <- (1 - level) / 2
  n_directions <- length(replicates)
  per_trait_probability <- matrix(
    NA_real_,
    nrow = n_directions, ncol = length(trait_cols),
    dimnames = list(NULL, trait_cols)
  )
  per_trait_lower <- per_trait_upper <- per_trait_probability
  joint_mean <- joint_lower <- joint_upper <- expected_attainment <-
    numeric(n_directions)
  replicate_scores <- vector("list", n_directions)

  beta_summary <- function(successes, trials) {
    shape_one <- successes + 0.5
    shape_two <- trials - successes + 0.5
    c(
      mean = shape_one / (shape_one + shape_two),
      lower = stats::qbeta(alpha_tail, shape_one, shape_two),
      upper = stats::qbeta(1 - alpha_tail, shape_one, shape_two)
    )
  }

  width <- intervals[, "upper"] - intervals[, "lower"]
  width[width <= .Machine$double.eps] <- pmax(
    abs(intervals[width <= .Machine$double.eps, "upper"]), 1
  )
  for (i in seq_len(n_directions)) {
    values <- as.matrix(replicates[[i]])[trait_cols, , drop = FALSE]
    favourable <- sweep(values, 1L, favourable_direction[trait_cols], "*")
    standardised <- sweep(favourable, 1L, outcome_scale[trait_cols], "/")
    trials <- ncol(standardised)
    achieved <- sweep(
      standardised, 1L, intervals[, "lower"], ">="
    )
    joint <- colSums(achieved) == length(trait_cols)
    joint_summary <- beta_summary(sum(joint), trials)
    joint_mean[i] <- joint_summary["mean"]
    joint_lower[i] <- joint_summary["lower"]
    joint_upper[i] <- joint_summary["upper"]
    for (trait in seq_along(trait_cols)) {
      trait_summary <- beta_summary(
        sum(achieved[trait, ]), trials
      )
      per_trait_probability[i, trait] <- trait_summary["mean"]
      per_trait_lower[i, trait] <- trait_summary["lower"]
      per_trait_upper[i, trait] <- trait_summary["upper"]
    }
    attainment <- sweep(
      sweep(standardised, 1L, intervals[, "lower"], "-"),
      1L, width, "/"
    )
    # The worst-satisfied trait is primary. The small bounded average term
    # breaks ties without allowing exceptional gain in one trait to conceal a
    # failure in another.
    scores <- apply(attainment, 2L, min) +
      augmentation * colMeans(pmin(attainment, 1))
    replicate_scores[[i]] <- scores
    expected_attainment[i] <- mean(scores)
  }

  # Shared replicate k uses the same random-number stream for every direction.
  # A paired bootstrap therefore estimates which vector is best without
  # discarding the induced cross-direction covariance.
  minimum_replicates <- min(vapply(replicate_scores, length, integer(1L)))
  bootstrap_selection_frequency <- rep(NA_real_, n_directions)
  if (minimum_replicates > 1L) {
    if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      stats::runif(1L)
    }
    entry_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", entry_seed, envir = globalenv()), add = TRUE)
    set.seed(.dgr_seed(bootstrap_seed))
    wins <- numeric(n_directions)
    for (draw in seq_len(1000L)) {
      sampled <- sample.int(
        minimum_replicates, minimum_replicates,
        replace = TRUE
      )
      draw_scores <- vapply(replicate_scores, function(x) {
        mean(x[seq_len(minimum_replicates)][sampled])
      }, numeric(1L))
      best <- which(draw_scores >= max(draw_scores) - 1e-12)
      wins[best] <- wins[best] + 1 / length(best)
    }
    bootstrap_selection_frequency <- wins / 1000
  }

  list(
    joint_probability = joint_mean,
    joint_lower = joint_lower,
    joint_upper = joint_upper,
    per_trait_probability = per_trait_probability,
    per_trait_lower = per_trait_lower,
    per_trait_upper = per_trait_upper,
    expected_attainment = expected_attainment,
    bootstrap_selection_frequency = bootstrap_selection_frequency,
    level = level,
    definition = paste(
      "Joint success means that every favourable cumulative response reaches",
      "the lower bound of its desired-gain interval at the declared cycle",
      "horizon. Probabilities use a Jeffreys Beta(1/2, 1/2) binomial posterior",
      "over independent simulation replicates. They are conditional on the",
      "founder population, covariance matrices and simulation assumptions.",
      "Bootstrap selection frequency is a resampling stability diagnostic,",
      "not the posterior probability that a direction is truly best."
    )
  )
}

# Evaluate one desired-gain direction by simulation.
#
#   direction          Unit-norm desired-gain direction.
#   setup              Simulation setup from founder_population().
#   arguments          Forwarded to simulate_selection_cycles().
#   n_replicates       Independent simulation replicates.
#   base_seed          Seed of the first replicate. The same seed sequence is
#                      used for every direction, so comparisons between
#                      directions share their stochasticity.
#   include_diversity  Whether to append the negated mean relationship.
#
# Returns the mean objective vector, its Monte Carlo standard error, the
# replicate-level outcomes and their covariance.
.dgr_evaluate_direction <- function(
  direction, setup, arguments, n_replicates, base_seed, include_diversity,
  direction_scale = NULL
) {
  trait_cols <- setup$trait_cols
  if (is.null(direction_scale)) {
    direction_scale <- stats::setNames(rep(1, length(trait_cols)), trait_cols)
  }
  named_direction <- stats::setNames(
    as.numeric(direction) * direction_scale[trait_cols], trait_cols
  )
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
  rownames(replicate_values) <- c(
    trait_cols, if (isTRUE(include_diversity)) "diversity"
  )
  # The Monte Carlo error is known per direction and differs between them, so
  # it is retained rather than recomputed or discarded. A surrogate that fits
  # one homogeneous nugget instead treats a precisely evaluated direction and a
  # noisy one as equally trustworthy.
  standard_error <- if (n_replicates > 1L) {
    apply(replicate_values, 1L, stats::sd) / sqrt(n_replicates)
  } else {
    rep(NA_real_, nrow(replicate_values))
  }
  list(
    mean = rowMeans(replicate_values),
    sd = standard_error,
    replicates = replicate_values,
    covariance = if (n_replicates > 1L) {
      stats::cov(t(replicate_values)) / n_replicates
    } else {
      NULL
    },
    n_replicates = n_replicates
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
#' When `desired_gain_intervals` is supplied, the search is restricted to rays
#' that intersect the breeder's acceptable interval box. This lets the breeder
#' say, for example, that yield improvement should lie between 0.5 and 1.5
#' genetic standard deviations while disease reduction should lie between
#' 0.25 and 1.0, without claiming that one exact ratio is known. The intervals
#' constrain the requested response pattern; they do not guarantee that the
#' corresponding absolute gains are attainable.
#'
#' # Ranking modes
#'
#' \describe{
#'   \item{`"pareto"`}{Returns the non-dominated set of multi-cycle outcomes
#'     and the direction generating each. This is the default because it makes
#'     no scalar preference assumption; choosing a point on a frontier can be
#'     easier than stating weights. The search uses randomised
#'     augmented-Chebyshev scalarisation, which converges on the frontier while
#'     needing only single-objective improvement.}
#'   \item{`"economic"`}{Maximises `sum(economic_weights * cumulative_gain)`.}
#'   \item{`"target"`}{Minimises the distance between the cumulative gain and a
#'     stated absolute target, reporting the per-trait shortfall.}
#'   \item{`"constrained"`}{Maximises the gain in `focal_trait` subject to
#'     `gain_floors` on the remaining traits, weighting improvement by the
#'     probability that every floor is met.}
#'   \item{`"interval"`}{Requires `desired_gain_intervals`. A high-gain
#'     replicate is one in which every favourable cumulative response reaches
#'     its interval's lower bound. Each evaluated vector receives a joint and
#'     per-trait posterior success probability. The recommendation maximises
#'     the lower credible bound of the joint probability, then expected
#'     balanced attainment, so a spectacular response in one trait cannot hide
#'     failure in another. No economic weights are used.}
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
#' # Interval-mode probability statement
#'
#' Let `R_j` be the favourable cumulative response for trait `j`, expressed in
#' the interval's units, and let `L_j` be its lower bound. A simulation
#' replicate achieves high genetic gain when
#' \deqn{R_j \ge L_j \quad\hbox{for every trait }j.}
#' If this joint event occurs in `s` of `n` independent replicates, interval
#' mode reports the Jeffreys binomial posterior
#' \deqn{p \mid s,n \sim \mathrm{Beta}(s+1/2,n-s+1/2).}
#' The recommended vector maximises the lower credible bound for `p`, with the
#' posterior mean and an augmented maximin attainment score as tie-breakers.
#' Consequently a direction supported by a few fortunate runs cannot outrank a
#' more precisely supported direction merely through its point estimate.
#'
#' The probability is conditional, not universal: it integrates the stochastic
#' segregation, mating and selection represented by the simulation while
#' holding the supplied founders, covariance matrices and programme parameters
#' fixed. Use [propagate_covariance_uncertainty()] to assess sensitivity to
#' estimated covariance matrices.
#'
#' # Interpreting the result
#'
#' Pareto mode returns no single best direction. Scalar modes return a
#' recommendation plus a stability region rather than presenting the numerical
#' maximiser as certain. Every recommendation is conditional on the supplied
#' genetic covariance matrix, founder germplasm and programme parameters.
#'
#' @param setup An object from [founder_population()].
#' @param n_cycles Number of selection cycles per evaluation.
#' @param mode Ranking mode. See Details.
#' @param budget Total number of simulated directions, including the initial
#'   design.
#' @param n_initial Size of the quasi-random initial design. Defaults to ten
#'   times the number of free parameters.
#' @param n_replicates Simulation replicates per direction. When omitted,
#'   interval mode uses 50 and other modes use 3. Interval mode warns below 20;
#'   at least 50 are recommended when leading vectors have similar success
#'   probabilities.
#' @param economic_weights Named weights required by `mode = "economic"`.
#' @param target_gains Named absolute targets required by `mode = "target"`.
#' @param focal_trait,gain_floors Objective and constraints required by
#'   `mode = "constrained"`.
#' @param include_diversity Whether to treat diversity, measured as the negated
#'   mean relationship among selected parents, as an additional Pareto
#'   objective. Family balancing and coancestry control can cost a large share
#'   of nominal gain, so the trade-off is reported rather than assumed away.
#'   Defaults to `FALSE` because a default founder setup has no verified marker
#'   panel.
#' @param allow_unverified_diversity Experimental override allowing a marker
#'   panel whose disjointness from QTL could not be verified. The override is
#'   included in checkpoint provenance; overlap known to be present is never
#'   allowed.
#' @param non_negative Whether to restrict the search to directions seeking
#'   improvement in every trait.
#' @param desired_gain_intervals Optional acceptable bounds created by
#'   [define_desired_gain_intervals()], or a trait-named matrix/data frame with
#'   `lower` and `upper` columns expressed in the favourable direction. The
#'   search considers only desired-gain directions whose ray intersects this
#'   box.
#' @param desired_gain_units Units in which searched directions and plain
#'   `desired_gain_intervals` are expressed. The default, genetic standard
#'   deviations, makes relative gains comparable across traits. An object from
#'   [define_desired_gain_intervals()] carries its own units and overrides this
#'   argument.
#' @param probability_level Credible level for interval-mode success
#'   probabilities. These probabilities describe simulation variation
#'   conditional on the supplied population and covariance parameters.
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
#' @references
#' Brown LD, Cai TT, DasGupta A (2001). Interval estimation for a binomial
#' proportion. *Statistical Science* 16:101--133.
#' \doi{10.1214/ss/1009213286}
#'
#' Wierzbicki AP (1982). A mathematical basis for satisficing decision making.
#' *Mathematical Modelling* 3:391--405.
#' \doi{10.1016/0270-0255(82)90038-0}
#'
#' Yang W-N, Nelson BL (1991). Using common random numbers and control variates
#' in multiple-comparison procedures. *Operations Research* 39:583--591.
#' \doi{10.1287/opre.39.4.583}
#'
#' @seealso [simulate_selection_cycles()], [gain_feasibility()]
#' @export
optimize_desired_gains <- function(
  setup,
  n_cycles = 5L,
  mode = c("pareto", "interval", "economic", "target", "constrained"),
  budget = 60L,
  n_initial = NULL,
  n_replicates = 3L,
  economic_weights = NULL,
  target_gains = NULL,
  focal_trait = NULL,
  gain_floors = NULL,
  include_diversity = FALSE,
  allow_unverified_diversity = FALSE,
  non_negative = TRUE,
  desired_gain_intervals = NULL,
  desired_gain_units = c("genetic_sd", "phenotypic_sd", "trait"),
  probability_level = 0.95,
  stability_tolerance = 0.05,
  n_candidates = 2000L,
  checkpoint = NULL,
  seed = 42L,
  verbose = TRUE,
  ...
) {
  mode <- match.arg(mode)
  if (mode == "interval" && missing(n_replicates)) {
    n_replicates <- 50L
  }
  desired_gain_units <- match.arg(desired_gain_units)
  .dgr_require_alphasimr("optimize_desired_gains()")
  if (!inherits(setup, "desiredgainr_sim_setup")) {
    stop("setup must be created by founder_population().", call. = FALSE)
  }
  trait_cols <- setup$trait_cols
  p <- length(trait_cols)
  if (p < 2L) {
    stop("At least two traits are required to search a desired-gain direction.",
      call. = FALSE
    )
  }
  if (inherits(desired_gain_intervals, "desiredgainr_gain_intervals")) {
    desired_gain_units <- attr(desired_gain_intervals, "gain_units")
    desired_gain_units <- match.arg(
      desired_gain_units, c("genetic_sd", "phenotypic_sd", "trait")
    )
  }
  intervals <- .dgr_validate_gain_intervals(
    desired_gain_intervals, trait_cols, non_negative
  )
  if (mode == "interval" && is.null(intervals)) {
    stop("mode = 'interval' requires desired_gain_intervals.",
      call. = FALSE
    )
  }
  if (!is.numeric(probability_level) || length(probability_level) != 1L ||
    !is.finite(probability_level) || probability_level <= 0 ||
    probability_level >= 1) {
    stop("probability_level must be one number strictly between zero and one.",
      call. = FALSE
    )
  }
  interval_horizon <- if (inherits(
    desired_gain_intervals, "desiredgainr_gain_intervals"
  )) {
    attr(desired_gain_intervals, "horizon_cycles")
  } else {
    NULL
  }
  if (!is.null(interval_horizon) && interval_horizon != n_cycles) {
    stop("desired_gain_intervals were defined for ", interval_horizon,
      " cycles but n_cycles = ", n_cycles, ".",
      call. = FALSE
    )
  }
  budget <- .dgr_positive_integer(budget, "budget")
  n_replicates <- .dgr_positive_integer(n_replicates, "n_replicates")
  if (mode == "interval" && n_replicates < 20L) {
    warning("Interval-attainment probabilities based on fewer than 20 ",
      "replicates will have wide credible intervals. Use at least 20 for a ",
      "decision and 50 or more when candidate vectors are close.",
      call. = FALSE
    )
  }
  n_cycles <- .dgr_positive_integer(n_cycles, "n_cycles")
  seed <- .dgr_seed(seed)
  # Replicate seeds are derived by adding offsets, so the base must leave room
  # for them inside the integer range rather than overflowing silently.
  if (abs(seed) > .Machine$integer.max - 1000L - budget * n_replicates) {
    stop("seed is too close to the integer limit for the replicate seeds ",
      "derived from it to remain distinct. Use a smaller seed.",
      call. = FALSE
    )
  }
  n_candidates <- .dgr_positive_integer(n_candidates, "n_candidates")
  if (is.null(n_initial)) n_initial <- min(budget, 10L * (p - 1L))
  n_initial <- .dgr_positive_integer(n_initial, "n_initial")
  if (n_initial > budget) {
    stop("n_initial cannot exceed budget.", call. = FALSE)
  }
  if (budget > n_initial + n_candidates) {
    stop("budget cannot exceed n_initial + n_candidates for the fixed ",
      "quasi-random search pool. Increase n_candidates or reduce budget.",
      call. = FALSE
    )
  }

  objective_names <- trait_cols
  if (isTRUE(include_diversity) && mode == "pareto") {
    # Diversity is only a defensible objective when it is measured on loci
    # that are not under selection. Without a verified marker panel it falls
    # back to all segregating sites, which include the quantitative trait
    # loci, and the diversity axis then partly measures genetic gain itself:
    # a direction gets penalised for working. Refuse rather than produce a
    # frontier whose second axis cannot be interpreted.
    if (is.null(setup$marker_panel) || !length(setup$marker_panel)) {
      stop(
        "include_diversity = TRUE requires a marker panel that excludes the ",
        "quantitative trait loci, and this setup has none.\n",
        "  Rebuild it with founder_population(n_markers_per_chromosome = ...), ",
        "or set include_diversity = FALSE.\n",
        "  Without a panel, diversity would be measured on all segregating ",
        "sites including the QTL, so the diversity axis would partly measure ",
        "the gain it is meant to be traded against.",
        call. = FALSE
      )
    }
    if (identical(setup$marker_qtl_overlap, TRUE)) {
      stop("The marker panel overlaps quantitative trait loci and cannot be ",
        "used as a neutral diversity objective.",
        call. = FALSE
      )
    }
    if (!identical(setup$marker_qtl_overlap, FALSE) &&
      !isTRUE(allow_unverified_diversity)) {
      stop(
        "The marker panel could not be verified as disjoint from the ",
        "quantitative trait loci. Set allow_unverified_diversity = TRUE only ",
        "for an explicitly experimental analysis; the resulting diversity ",
        "axis may partly measure response at QTL rather than coancestry.",
        call. = FALSE
      )
    }
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
  direction_scale <- .dgr_desired_gain_scale(setup, desired_gain_units)

  if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    stats::runif(1L)
  }
  entry_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit(assign(".Random.seed", entry_seed, envir = globalenv()), add = TRUE)
  set.seed(seed)

  simulation_arguments <- c(list(n_cycles = n_cycles), list(...))
  favourable_direction <- .dgr_direction(
    trait_cols,
    if (is.null(simulation_arguments$lower_is_better)) {
      NULL
    } else {
      simulation_arguments$lower_is_better
    }
  )
  make_directions <- function(U) {
    if (is.null(intervals)) {
      .dgr_directions_from_cube(U, non_negative)
    } else {
      .dgr_directions_from_intervals(U, intervals)
    }
  }
  candidate_pool <- make_directions(.dgr_halton(
    n_candidates, p,
    skip = 500L
  ))
  colnames(candidate_pool) <- trait_cols

  # Everything that determines a result goes into the fingerprint. Matching on
  # trait names alone let a checkpoint resume evaluations from a different
  # scientific problem -- different founders, different G, a different mating
  # system or cycle count -- and present them as this run's.
  fingerprint <- .dgr_optimiser_fingerprint(
    setup = setup, trait_cols = trait_cols, objective_names = objective_names,
    simulation_arguments = simulation_arguments, n_replicates = n_replicates,
    seed = seed, include_diversity = include_diversity,
    non_negative = non_negative, mode = mode,
    objective_parameters = list(
      economic_weights = economic_weights,
      target_gains = target_gains,
      focal_trait = focal_trait,
      gain_floors = gain_floors,
      desired_gain_intervals = intervals,
      desired_gain_units = desired_gain_units,
      direction_scale = direction_scale,
      probability_level = probability_level,
      interval_horizon = interval_horizon
    ),
    search_parameters = list(
      n_initial = n_initial, n_candidates = n_candidates,
      stability_tolerance = stability_tolerance,
      allow_unverified_diversity = isTRUE(allow_unverified_diversity)
    )
  )

  evaluated_directions <- NULL
  evaluated_objectives <- NULL
  evaluated_errors <- NULL
  evaluated_counts <- NULL
  evaluated_replicates <- NULL
  if (!is.null(checkpoint) && file.exists(checkpoint)) {
    saved <- readRDS(checkpoint)
    if (identical(saved$fingerprint, fingerprint)) {
      evaluated_directions <- saved$directions
      evaluated_objectives <- saved$objectives
      evaluated_errors <- saved$errors
      evaluated_counts <- saved$counts
      evaluated_replicates <- saved$replicates
      if (is.null(evaluated_replicates) ||
        length(evaluated_replicates) != nrow(evaluated_directions)) {
        stop("The checkpoint predates replicate-level uncertainty storage and ",
          "cannot be resumed safely. Start a new checkpoint.",
          call. = FALSE
        )
      }
      if (isTRUE(verbose)) {
        message(sprintf(
          "Resuming from %d checkpointed evaluation(s).",
          nrow(evaluated_directions)
        ))
      }
    } else {
      differing <- names(fingerprint)[
        !vapply(names(fingerprint), function(field) {
          identical(saved$fingerprint[[field]], fingerprint[[field]])
        }, logical(1L))
      ]
      stop(
        "The checkpoint at '", checkpoint, "' was produced under different ",
        "conditions and cannot be resumed.\n",
        "  Differing: ",
        paste(if (length(differing)) differing else "unknown", collapse = ", "),
        ".\n  Resuming would mix evaluations from two different problems in ",
        "one frontier. Delete the file to start again, or point checkpoint ",
        "at a new path.",
        call. = FALSE
      )
    }
  }

  # Directions are matched to a tolerance so that a repeated proposal is
  # treated as additional replication of the same point rather than as a new
  # one. Expected improvement at an already-observed point can stay positive
  # under a fitted nugget, so without this a finite budget can be spent
  # re-evaluating one direction while the record shows several.
  match_evaluated <- function(direction) {
    if (is.null(evaluated_directions) || !nrow(evaluated_directions)) {
      return(NA_integer_)
    }
    distances <- sqrt(colSums(
      (t(evaluated_directions) - as.numeric(direction))^2
    ))
    closest <- which.min(distances)
    if (distances[closest] < 1e-8) as.integer(closest) else NA_integer_
  }

  evaluate_and_store <- function(direction) {
    existing <- match_evaluated(direction)
    seed_offset <- if (is.na(existing)) 0L else evaluated_counts[existing]
    outcome <- .dgr_evaluate_direction(
      direction, setup, simulation_arguments, n_replicates,
      base_seed = seed + 1000L + seed_offset, include_diversity,
      direction_scale = direction_scale
    )
    if (is.na(existing)) {
      evaluated_directions <<- rbind(evaluated_directions, direction)
      evaluated_objectives <<- rbind(evaluated_objectives, outcome$mean)
      evaluated_errors <<- rbind(evaluated_errors, outcome$sd)
      evaluated_counts <<- c(evaluated_counts, n_replicates)
      evaluated_replicates[[length(evaluated_replicates) + 1L]] <<-
        outcome$replicates
    } else {
      # A genuine repeat receives new seeds. Pool its replicate outcomes and
      # recompute uncertainty from those observations; deterministic replays
      # of the same seeds must never masquerade as extra information.
      combined <- cbind(evaluated_replicates[[existing]], outcome$replicates)
      total_n <- ncol(combined)
      evaluated_replicates[[existing]] <<- combined
      evaluated_objectives[existing, ] <<- rowMeans(combined)
      evaluated_errors[existing, ] <<- apply(combined, 1L, stats::sd) /
        sqrt(total_n)
      evaluated_counts[existing] <<- total_n
      if (isTRUE(verbose)) {
        message(sprintf(
          "Direction already evaluated; pooled as %d replicates.", total_n
        ))
      }
    }
    if (!is.null(checkpoint)) {
      saveRDS(
        list(
          fingerprint = fingerprint,
          trait_cols = trait_cols, objective_names = objective_names,
          directions = evaluated_directions,
          objectives = evaluated_objectives,
          errors = evaluated_errors,
          counts = evaluated_counts,
          replicates = evaluated_replicates
        ),
        checkpoint
      )
    }
    invisible(NULL)
  }

  already <- if (is.null(evaluated_directions)) 0L else nrow(evaluated_directions)
  if (already < n_initial) {
    initial <- make_directions(.dgr_halton(n_initial, p))
    for (i in seq(already + 1L, n_initial)) {
      if (isTRUE(verbose)) {
        message(sprintf("Initial design %d of %d.", i, n_initial))
      }
      evaluate_and_store(initial[i, , drop = FALSE])
    }
  }

  scalarise <- function(objectives, weights = NULL, ranges = NULL) {
    objectives <- as.matrix(objectives)
    switch(mode,
      economic = as.numeric(objectives[, trait_cols, drop = FALSE] %*%
        economic_weights),
      interval = {
        favourable <- sweep(
          objectives[, trait_cols, drop = FALSE], 2L,
          favourable_direction, "*"
        )
        standardised <- sweep(
          favourable, 2L, direction_scale, "/"
        )
        width <- intervals[, "upper"] - intervals[, "lower"]
        width[width <= .Machine$double.eps] <- pmax(
          abs(intervals[width <= .Machine$double.eps, "upper"]), 1
        )
        attainment <- sweep(
          sweep(standardised, 2L, intervals[, "lower"], "-"),
          2L, width, "/"
        )
        apply(attainment, 1L, min) +
          0.05 * rowMeans(pmin(attainment, 1))
      },
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
        favourable_objectives <- objectives
        favourable_objectives[, trait_cols] <- sweep(
          favourable_objectives[, trait_cols, drop = FALSE], 2L,
          favourable_direction, "*"
        )
        if (is.null(ranges)) {
          ranges <- apply(favourable_objectives, 2L, range)
        }
        spread <- pmax(ranges[2L, ] - ranges[1L, ], 1e-12)
        normalised <- sweep(
          sweep(favourable_objectives, 2L, ranges[1L, ], "-"),
          2L, spread, "/"
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
    pareto_ranges <- if (mode == "pareto") {
      favourable_objectives <- evaluated_objectives
      favourable_objectives[, trait_cols] <- sweep(
        favourable_objectives[, trait_cols, drop = FALSE], 2L,
        favourable_direction, "*"
      )
      apply(favourable_objectives, 2L, range)
    } else {
      NULL
    }
    response <- scalarise(
      evaluated_objectives, chebyshev_weights,
      ranges = pareto_ranges
    )

    # The Monte Carlo error of each evaluation is known, and differs between
    # directions because a direction that erodes variance quickly gives noisier
    # replicates. Supplying it as a per-point nugget stops the surrogate from
    # treating a precisely evaluated direction and a noisy one as equally
    # trustworthy, which is what a single fitted homogeneous nugget does.
    response_error <- .dgr_scalarised_error(
      evaluated_replicates, scalarise, chebyshev_weights,
      ranges = pareto_ranges
    )
    response_error[!is.finite(response_error)] <- 0
    model <- .dgr_gp_fit(
      evaluated_directions, response,
      noise_variance = response_error^2
    )
    prediction <- .dgr_gp_predict(model, candidate_pool)
    acquisition <- .dgr_expected_improvement(
      prediction$mean, prediction$sd, max(response)
    )
    # A direction already evaluated cannot be improved upon by evaluating it
    # again under this acquisition, so it is excluded from the pool rather than
    # competing with unvisited points.
    if (!is.null(evaluated_directions) && nrow(evaluated_directions)) {
      distances <- apply(candidate_pool, 1L, function(point) {
        min(sqrt(colSums((t(evaluated_directions) - point)^2)))
      })
      acquisition[distances < 1e-8] <- -Inf
    }

    if (mode == "constrained") {
      for (trait in names(gain_floors)) {
        constraint_model <- .dgr_gp_fit(
          evaluated_directions, evaluated_objectives[, trait],
          noise_variance = if (is.null(evaluated_errors)) {
            NULL
          } else {
            evaluated_errors[, trait]^2
          }
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

    if (!any(is.finite(acquisition))) {
      stop(
        "No unevaluated candidate direction remains in the acquisition pool. ",
        "Increase n_candidates or resume with a larger, new checkpoint; the ",
        "optimizer stopped instead of repeatedly evaluating one direction.",
        call. = FALSE
      )
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
  representative_gains <- .dgr_representative_gains(
    evaluated_directions, intervals
  )
  colnames(evaluated_objectives) <- objective_names
  if (!is.null(evaluated_errors)) {
    colnames(evaluated_errors) <- objective_names
  }

  # The frontier and the recommendation are read from posterior means rather
  # than raw replicate averages, so that fortunate runs do not populate them.
  smoothed <- vapply(objective_names, function(objective) {
    model <- .dgr_gp_fit(
      evaluated_directions, evaluated_objectives[, objective],
      noise_variance = if (is.null(evaluated_errors)) {
        NULL
      } else {
        evaluated_errors[, objective]^2
      }
    )
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
    data.table::set(results,
      j = paste0("d_", trait),
      value = evaluated_directions[, trait]
    )
    data.table::set(results,
      j = paste0("requested_", trait),
      value = representative_gains[, trait]
    )
  }
  for (objective in objective_names) {
    data.table::set(results,
      j = paste0("observed_", objective),
      value = evaluated_objectives[, objective]
    )
    data.table::set(results,
      j = paste0("posterior_", objective),
      value = smoothed[, objective]
    )
  }

  dominance_objectives <- smoothed
  dominance_objectives[, trait_cols] <- sweep(
    dominance_objectives[, trait_cols, drop = FALSE], 2L,
    favourable_direction, "*"
  )
  pareto_flag <- .dgr_non_dominated(dominance_objectives)
  data.table::set(results, j = "pareto_optimal", value = pareto_flag)

  # A frontier drawn as a set of points implies the membership is certain. It
  # is not: each coordinate carries Monte Carlo error, and a direction sitting
  # just inside the frontier may be there by luck. Resampling the evaluations
  # within their own standard errors gives the probability that each direction
  # is non-dominated, which is what a breeder choosing from the frontier needs.
  frontier_probability <- rep(NA_real_, nrow(smoothed))
  minimum_replicates <- if (length(evaluated_replicates)) {
    min(vapply(evaluated_replicates, ncol, integer(1L)))
  } else {
    0L
  }
  if (minimum_replicates > 1L) {
    set.seed(as.integer(
      (abs(as.double(seed)) + 600001) %% .Machine$integer.max
    ))
    n_frontier_draws <- 400L
    membership <- matrix(
      0,
      nrow = nrow(smoothed), ncol = n_frontier_draws
    )
    for (draw in seq_len(n_frontier_draws)) {
      # Replicate k uses the same simulation seed for every direction. A paired
      # bootstrap therefore preserves both within-direction objective
      # covariance and common-random-number covariance across directions.
      sampled <- sample.int(minimum_replicates, minimum_replicates,
        replace = TRUE
      )
      perturbation <- t(vapply(
        evaluated_replicates,
        function(values) {
          values <- values[, seq_len(minimum_replicates), drop = FALSE]
          rowMeans(values[, sampled, drop = FALSE]) - rowMeans(values)
        },
        numeric(length(objective_names))
      ))
      perturbed <- smoothed + perturbation
      perturbed[, trait_cols] <- sweep(
        perturbed[, trait_cols, drop = FALSE], 2L,
        favourable_direction, "*"
      )
      membership[, draw] <- .dgr_non_dominated(perturbed)
    }
    frontier_probability <- rowMeans(membership)
  }
  data.table::set(
    results,
    j = "pareto_probability", value = frontier_probability
  )
  for (objective in objective_names) {
    data.table::set(
      results,
      j = paste0("se_", objective),
      value = if (is.null(evaluated_errors)) {
        NA_real_
      } else {
        evaluated_errors[, objective]
      }
    )
  }
  data.table::set(
    results,
    j = "n_replicates",
    value = if (is.null(evaluated_counts)) NA_integer_ else evaluated_counts
  )

  interval_attainment <- NULL
  if (mode == "interval") {
    interval_attainment <- .dgr_interval_attainment(
      evaluated_replicates, trait_cols, intervals, direction_scale,
      favourable_direction,
      level = probability_level,
      bootstrap_seed = as.integer(
        (abs(as.double(seed)) + 700001) %% .Machine$integer.max
      )
    )
    data.table::set(
      results,
      j = "joint_high_gain_probability",
      value = interval_attainment$joint_probability
    )
    data.table::set(
      results,
      j = "joint_probability_lower",
      value = interval_attainment$joint_lower
    )
    data.table::set(
      results,
      j = "joint_probability_upper",
      value = interval_attainment$joint_upper
    )
    data.table::set(
      results,
      j = "expected_balanced_attainment",
      value = interval_attainment$expected_attainment
    )
    data.table::set(
      results,
      j = "bootstrap_selection_frequency",
      value = interval_attainment$bootstrap_selection_frequency
    )
    for (trait in trait_cols) {
      data.table::set(
        results,
        j = paste0("high_gain_probability_", trait),
        value = interval_attainment$per_trait_probability[, trait]
      )
      data.table::set(
        results,
        j = paste0("high_gain_probability_lower_", trait),
        value = interval_attainment$per_trait_lower[, trait]
      )
      data.table::set(
        results,
        j = paste0("high_gain_probability_upper_", trait),
        value = interval_attainment$per_trait_upper[, trait]
      )
    }
  }

  final_scalar <- if (mode == "pareto") {
    rep(NA_real_, nrow(smoothed))
  } else if (mode == "interval") {
    interval_attainment$joint_lower
  } else {
    scalarise(smoothed)
  }
  stability <- rep(NA, nrow(smoothed))
  recommended <- NULL
  recommended_gains <- NULL
  if (mode != "pareto") {
    best_index <- if (mode == "interval") {
      order(
        -interval_attainment$joint_lower,
        -interval_attainment$joint_probability,
        -interval_attainment$expected_attainment
      )[1L]
    } else {
      which.max(final_scalar)
    }
    best_value <- final_scalar[best_index]
    threshold <- best_value - stability_tolerance * abs(best_value)
    stability <- final_scalar >= threshold
    recommended <- evaluated_directions[best_index, ]
    names(recommended) <- trait_cols
    recommended_gains <- representative_gains[best_index, ]
    names(recommended_gains) <- trait_cols
    data.table::set(results, j = "objective", value = final_scalar)
    data.table::set(results, j = "within_stability_region", value = stability)
  }
  recommendation_order <- if (mode == "interval") {
    order(
      -results[["joint_probability_lower"]],
      -results[["joint_high_gain_probability"]],
      -results[["expected_balanced_attainment"]]
    )
  } else {
    NULL
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
    desired_gains = representative_gains,
    simulation_directions = sweep(
      evaluated_directions, 2L, direction_scale, "*"
    ),
    desired_gain_units = desired_gain_units,
    desired_gain_intervals = intervals,
    observed_objectives = evaluated_objectives,
    posterior_objectives = smoothed,
    pareto_set = results[pareto_optimal == TRUE],
    recommended_direction = recommended,
    recommended_desired_gains = recommended_gains,
    recommended_probability = if (mode == "interval") {
      c(
        mean = interval_attainment$joint_probability[best_index],
        lower = interval_attainment$joint_lower[best_index],
        upper = interval_attainment$joint_upper[best_index]
      )
    } else {
      NULL
    },
    recommendations = if (mode == "interval") {
      results[recommendation_order]
    } else {
      NULL
    },
    interval_attainment = interval_attainment,
    stability_proportion = if (mode == "pareto") NA_real_ else mean(stability),
    stability_tolerance = stability_tolerance,
    iteration_log = data.table::rbindlist(iteration_log),
    simulation_arguments = simulation_arguments,
    seed = seed,
    monte_carlo_error = evaluated_errors,
    replicate_counts = evaluated_counts,
    replicate_outcomes = evaluated_replicates,
    fingerprint = fingerprint,
    uncertainty = list(
      frontier_probability_available = !all(is.na(frontier_probability)),
      note = paste(
        "pareto_probability is the proportion of resampled frontiers on which",
        "a direction remained non-dominated under a paired bootstrap of the",
        "shared-seed replicate outcomes. It preserves covariance among",
        "objectives and directions and represents simulation error only.",
        "It does NOT represent uncertainty in G and P, which are treated as",
        "fixed throughout this search; use index_uncertainty() to assess that",
        "separately, and treat a frontier from a single covariance estimate as",
        "conditional on it."
      )
    ),
    interpretation = paste(
      "Only the direction of a desired-gain vector is optimised, because the",
      "attainable magnitude is fixed by selection intensity. The reported",
      "desired_gains are representative points on the searched rays in",
      desired_gain_units, "units; they are requests, not guaranteed realised",
      "responses. The ranking mode defines what 'best' means. The reported",
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
  cat("  Desired-gain directions expressed in:", x$desired_gain_units, "\n")
  if (!is.null(x$desired_gain_intervals)) {
    cat("  Search restricted to breeder-defined desired-gain intervals.\n")
  }
  if (identical(x$mode, "pareto")) {
    cat(sprintf(
      "  Non-dominated directions: %d of %d\n",
      nrow(x$pareto_set), nrow(x$directions)
    ))
    cat(
      "  No single direction is recommended; choose a point on the",
      "frontier.\n"
    )
  } else {
    if (!is.null(x$recommended_desired_gains)) {
      cat("  Recommended desired-gain request (", x$desired_gain_units,
        "):\n",
        sep = ""
      )
      print(round(x$recommended_desired_gains, 4L))
    } else {
      cat("  Recommended direction:\n")
      print(round(x$recommended_direction, 4L))
    }
    if (identical(x$mode, "interval")) {
      probability <- x$recommended_probability
      cat(sprintf(
        "  Joint probability of reaching every lower bound: %.1f%%\n",
        100 * probability["mean"]
      ))
      cat(sprintf(
        "  %.0f%% credible interval: %.1f%% to %.1f%%\n",
        100 * x$interval_attainment$level,
        100 * probability["lower"], 100 * probability["upper"]
      ))
    }
    cat(sprintf(
      "  Directions within %.0f%% of the optimum: %.1f%%\n",
      100 * x$stability_tolerance, 100 * x$stability_proportion
    ))
  }
  invisible(x)
}
