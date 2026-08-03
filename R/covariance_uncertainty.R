# Propagating covariance-estimation uncertainty through the simulation and the
# direction optimiser.
#
# index_uncertainty() answers the question for a single-cycle index. This file
# answers it for the multi-cycle recommendation, which is harder because the
# genetic covariance is not a matrix the simulation multiplies by: it is the
# architecture the founder population is built from. Resampling it therefore
# means rebuilding the founders, which is why this is an outer loop rather than
# a perturbation.

#' Draw jointly admissible genetic and residual covariance matrices
#'
#' A pair \eqn{(\mathbf{G}^{*}, \mathbf{P}^{*})} is only admissible when
#' \eqn{\mathbf{P}^{*} - \mathbf{G}^{*}} is itself a covariance matrix.
#' Resampling the two independently produces pairs that violate this, and an
#' index built on such a pair reports a heritability above one. Drawing the
#' genetic and residual matrices separately and reassembling
#' \eqn{\mathbf{P}^{*} = \mathbf{G}^{*} + \mathbf{E}^{*}} makes every draw
#' admissible by construction.
#'
#' @param G Genetic covariance matrix.
#' @param P Phenotypic covariance matrix. When `NULL`, only `G` is drawn and
#'   the result carries no residual.
#' @param genetic_df,residual_df Degrees of freedom. See [index_uncertainty()]
#'   for what these count; they are the number of independent genetic units,
#'   not the number of plots.
#' @param n_draws Number of draws.
#'
#' @return A list of draws, each with `G`, `E` and `P`.
#' @export
#'
#' @examples
#' traits <- c("yield", "protein")
#' G <- matrix(c(1.0, 0.2, 0.2, 0.5), 2, dimnames = list(traits, traits))
#' P <- matrix(c(2.5, 0.4, 0.4, 1.2), 2, dimnames = list(traits, traits))
#' draws <- draw_covariance_pairs(G, P,
#'   genetic_df = 40, residual_df = 200,
#'   n_draws = 5
#' )
#' length(draws)
#' # Every draw is admissible by construction.
#' all(vapply(draws, function(d) {
#'   min(eigen(d$P - d$G, symmetric = TRUE, only.values = TRUE)$values) > 0
#' }, logical(1)))
#'
#' @seealso [index_uncertainty()], [optimize_desired_gains()]
draw_covariance_pairs <- function(
  G, P = NULL, genetic_df, residual_df = NULL, n_draws = 20L
) {
  trait_cols <- colnames(G)
  if (is.null(trait_cols)) {
    stop("G must carry trait names.", call. = FALSE)
  }
  p <- ncol(G)
  n_draws <- .dgr_positive_integer(n_draws, "n_draws")
  genetic_df <- .dgr_wishart_df(genetic_df, p, "genetic_df")

  G_draws <- .dgr_wishart_draws(G, genetic_df, n_draws, "G")
  E_draws <- NULL
  if (!is.null(P)) {
    .dgr_check_compatible(G, P, "draw_covariance_pairs()")
    residual_df <- .dgr_wishart_df(
      if (is.null(residual_df)) genetic_df else residual_df, p, "residual_df"
    )
    E <- P - G
    dimnames(E) <- dimnames(G)
    E_draws <- .dgr_wishart_draws(E, residual_df, n_draws, "P - G")
  }

  lapply(seq_len(n_draws), function(draw) {
    G_star <- G_draws[, , draw]
    dimnames(G_star) <- list(trait_cols, trait_cols)
    if (is.null(E_draws)) {
      return(list(G = G_star, E = NULL, P = NULL, draw = draw))
    }
    E_star <- E_draws[, , draw]
    dimnames(E_star) <- list(trait_cols, trait_cols)
    P_star <- G_star + E_star
    dimnames(P_star) <- list(trait_cols, trait_cols)
    list(G = G_star, E = E_star, P = P_star, draw = draw)
  })
}

#' Propagate covariance uncertainty into a desired-gain recommendation
#'
#' [optimize_desired_gains()] conditions on one genetic covariance estimate.
#' Its frontier therefore carries Monte Carlo error but treats the covariance
#' as known, so a direction that is only best because of how `G` happened to be
#' estimated is indistinguishable from one that is robustly best.
#'
#' This function repeats the evaluation of a fixed set of directions across
#' draws from the sampling distribution of the covariance, rebuilding the
#' founder population for each draw, and reports how far the recommendation
#' moves.
#'
#' @details
#' # Why the founders are rebuilt
#'
#' The genetic covariance is not a matrix the simulation multiplies by. It is
#' the architecture the trait is built from, so a different `G` means different
#' quantitative trait locus effects and a different founder population. Each
#' draw therefore calls [founder_population()] again. That is what makes this
#' expensive, and why the directions are supplied rather than searched: running
#' a full Bayesian optimisation inside every covariance draw would multiply an
#' already slow search by `n_covariance_draws`.
#'
#' The intended workflow is to run [optimize_desired_gains()] once and pass the
#' complete optimisation object here. The trait-unit directions actually used
#' by the simulator are then selected automatically; this avoids mistaking
#' genetic-SD search coordinates for original trait units.
#'
#' # The two sources of uncertainty are separated
#'
#' Within one covariance draw, replicate simulations differ by Monte Carlo
#' error alone. Across draws, they differ by both. The variance decomposition
#' reported as `variance_components` splits them, because they have different
#' remedies: Monte Carlo error falls with more replicates, covariance
#' uncertainty does not and can only be reduced by estimating `G` better.
#'
#' If the covariance component dominates, running more simulation replicates is
#' wasted effort.
#'
#' # What to read
#'
#' When `rank_weights` is supplied, `rank_churn` is the Spearman correlation
#' between the declared scalar ranking under each draw and under the point
#' estimate. Raw objectives are never averaged across incompatible units.
#' `frontier_membership` is the
#' proportion of draws in which each direction was non-dominated. A direction
#' that is on the frontier under the point estimate but in only a third of the
#' draws is not a robust recommendation.
#'
#' @param setup A setup from [founder_population()], used for the founder
#'   genomes, heritabilities and marker panel. Its `G_target` is the point
#'   estimate around which the draws are taken.
#' @param directions An object returned by [optimize_desired_gains()], or a
#'   matrix of desired-gain directions in original trait units, one per row,
#'   with columns named by trait. Passing the optimisation object is safer when
#'   its search used genetic- or phenotypic-SD units because the stored
#'   trait-unit simulation directions are used automatically.
#' @param genetic_df Degrees of freedom for the genetic covariance. This
#'   governs the width of everything reported and counts independent genetic
#'   units, not plots.
#' @param P Optional phenotypic covariance. When supplied, the residual is
#'   resampled too and each draw is admissible by construction.
#' @param residual_df Degrees of freedom for the residual covariance.
#' @param n_covariance_draws Number of covariance draws.
#' @param n_replicates Simulation replicates within each draw.
#' @param n_cycles Selection cycles per evaluation.
#' @param include_diversity Whether to append diversity as an objective.
#' @param allow_unverified_diversity Experimental override for a marker panel
#'   whose QTL disjointness could not be verified. Known overlap is rejected.
#' @param rank_weights Optional named scalarisation weights for every objective.
#'   Required to compute rank churn; when absent, rank churn is reported as not
#'   requested rather than averaging objectives in their raw units.
#' @param verbose Whether to report progress.
#' @param seed Random seed. The caller's stream is restored on exit.
#' @param ... Further arguments passed to [simulate_selection_cycles()].
#'
#' @return An object of class `desiredgainr_covariance_uncertainty`.
#'
#' @examples
#' \donttest{
#' # See the "Uncertainty" section of the simulation vignette; this example
#' # needs AlphaSimR and takes several minutes.
#' }
#'
#' @seealso [optimize_desired_gains()], [index_uncertainty()],
#'   [draw_covariance_pairs()]
#' @export
propagate_covariance_uncertainty <- function(
  setup,
  directions,
  genetic_df,
  P = NULL,
  residual_df = NULL,
  n_covariance_draws = 20L,
  n_replicates = 3L,
  n_cycles = 5L,
  include_diversity = FALSE,
  allow_unverified_diversity = FALSE,
  rank_weights = NULL,
  verbose = FALSE,
  seed = 42L,
  ...
) {
  .dgr_require_alphasimr("propagate_covariance_uncertainty()")
  if (!inherits(setup, "desiredgainr_sim_setup")) {
    stop("setup must be created by founder_population().", call. = FALSE)
  }
  trait_cols <- setup$trait_cols
  if (inherits(directions, "desiredgainr_optimisation")) {
    if (is.null(directions$simulation_directions)) {
      stop("This optimisation object does not contain trait-unit simulation ",
        "directions. Re-run optimize_desired_gains() with the current package.",
        call. = FALSE
      )
    }
    directions <- directions$simulation_directions
  }
  directions <- as.matrix(directions)
  if (is.null(colnames(directions))) {
    if (ncol(directions) != length(trait_cols)) {
      stop("directions must have one column per trait.", call. = FALSE)
    }
    colnames(directions) <- trait_cols
  }
  absent <- setdiff(trait_cols, colnames(directions))
  if (length(absent)) {
    stop("directions is missing traits: ", paste(absent, collapse = ", "),
      call. = FALSE
    )
  }
  directions <- directions[, trait_cols, drop = FALSE]
  if (!nrow(directions)) {
    stop("At least one direction must be supplied.", call. = FALSE)
  }

  n_covariance_draws <- .dgr_positive_integer(
    n_covariance_draws, "n_covariance_draws"
  )
  n_replicates <- .dgr_positive_integer(n_replicates, "n_replicates")
  if (n_replicates < 2L) {
    stop("n_replicates must be at least 2 to separate Monte Carlo variance ",
      "from covariance-estimation variance.",
      call. = FALSE
    )
  }
  n_cycles <- .dgr_positive_integer(n_cycles, "n_cycles")
  seed <- .dgr_seed(seed)

  if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    stats::runif(1L)
  }
  entry_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit(assign(".Random.seed", entry_seed, envir = globalenv()), add = TRUE)
  set.seed(seed)

  objective_names <- trait_cols
  if (isTRUE(include_diversity)) {
    if (is.null(setup$marker_panel) || !length(setup$marker_panel)) {
      stop("include_diversity = TRUE requires a marker panel disjoint from ",
        "the quantitative trait loci; rebuild the setup with ",
        "n_markers_per_chromosome.",
        call. = FALSE
      )
    }
    if (identical(setup$marker_qtl_overlap, TRUE)) {
      stop("The marker panel overlaps QTL and cannot be used as a neutral ",
        "diversity objective.",
        call. = FALSE
      )
    }
    if (!identical(setup$marker_qtl_overlap, FALSE) &&
      !isTRUE(allow_unverified_diversity)) {
      stop("Marker-QTL disjointness was not verified. Set ",
        "allow_unverified_diversity = TRUE only for an explicitly ",
        "experimental analysis.",
        call. = FALSE
      )
    }
    objective_names <- c(objective_names, "diversity")
  }

  if (!is.null(rank_weights)) {
    rank_weights <- .dgr_named_vector(
      rank_weights, objective_names, "rank_weights"
    )
  }

  draws <- draw_covariance_pairs(
    G = setup$G_target, P = P, genetic_df = genetic_df,
    residual_df = residual_df, n_draws = n_covariance_draws
  )
  simulation_arguments <- c(list(n_cycles = n_cycles), list(...))

  # outcomes[direction, objective, draw]
  outcomes <- array(
    NA_real_,
    dim = c(nrow(directions), length(objective_names), n_covariance_draws),
    dimnames = list(NULL, objective_names, NULL)
  )
  within_draw_variance <- array(
    NA_real_,
    dim = c(nrow(directions), length(objective_names), n_covariance_draws)
  )
  rebuild_failures <- 0L

  for (draw_index in seq_along(draws)) {
    if (isTRUE(verbose)) {
      message(sprintf(
        "Covariance draw %d of %d.", draw_index, n_covariance_draws
      ))
    }
    drawn <- draws[[draw_index]]
    # The founder population is rebuilt because G defines the trait
    # architecture rather than multiplying into it. A draw that cannot be
    # realised as a trait architecture is recorded and skipped rather than
    # silently replaced by the point estimate.
    drawn_setup <- tryCatch(
      suppressWarnings(founder_population(
        founders = setup$founders,
        G = drawn$G,
        h2 = setup$h2,
        residual_covariance = drawn$E,
        n_qtl_per_chromosome = setup$n_qtl_per_chromosome,
        n_markers_per_chromosome = setup$n_markers_per_chromosome,
        dominance_degree = setup$dominance_degree,
        dominance_variance = setup$dominance_degree_variance,
        heritability = setup$heritability_type,
        seed = setup$seed
      )),
      error = function(e) NULL
    )
    if (is.null(drawn_setup)) {
      rebuild_failures <- rebuild_failures + 1L
      next
    }

    for (row in seq_len(nrow(directions))) {
      evaluation <- .dgr_evaluate_direction(
        directions[row, , drop = FALSE], drawn_setup, simulation_arguments,
        n_replicates,
        base_seed = seed + 1000L * draw_index,
        include_diversity
      )
      outcomes[row, , draw_index] <- evaluation$mean
      within_draw_variance[row, , draw_index] <- evaluation$sd^2
    }
  }

  usable <- which(apply(outcomes, 3L, function(slice) all(is.finite(slice))))
  if (length(usable) < 2L) {
    stop("Fewer than two covariance draws produced usable results, so ",
      "nothing can be said about the spread. ", rebuild_failures,
      " draw(s) could not be realised as a trait architecture.",
      call. = FALSE
    )
  }
  outcomes <- outcomes[, , usable, drop = FALSE]
  within_draw_variance <- within_draw_variance[, , usable, drop = FALSE]

  covariance_draw_mean <- apply(outcomes, c(1L, 2L), mean)
  colnames(covariance_draw_mean) <- objective_names

  # Evaluate the original covariance setup explicitly. The average of sampled
  # draws is an uncertainty summary, not the point estimate. When P is given,
  # rebuild once with its point residual so the reference and the draws use the
  # same phenotype model.
  point_setup <- setup
  if (!is.null(P)) {
    point_setup <- suppressWarnings(founder_population(
      founders = setup$founders,
      G = setup$G_target,
      h2 = setup$h2,
      residual_covariance = P - setup$G_target,
      n_qtl_per_chromosome = setup$n_qtl_per_chromosome,
      n_markers_per_chromosome = setup$n_markers_per_chromosome,
      dominance_degree = setup$dominance_degree,
      dominance_variance = setup$dominance_degree_variance,
      heritability = setup$heritability_type,
      seed = setup$seed
    ))
  }
  point_estimate <- t(vapply(seq_len(nrow(directions)), function(row) {
    .dgr_evaluate_direction(
      directions[row, , drop = FALSE], point_setup, simulation_arguments,
      n_replicates,
      base_seed = seed + 500000L, include_diversity
    )$mean
  }, numeric(length(objective_names))))
  colnames(point_estimate) <- objective_names

  # Variance decomposition. The total spread of a direction's outcome across
  # draws contains both the covariance uncertainty and the Monte Carlo error
  # that every evaluation carries; subtracting the mean within-draw variance
  # leaves the part attributable to the covariance estimate alone.
  total_variance <- apply(outcomes, c(1L, 2L), stats::var)
  monte_carlo_variance <- apply(within_draw_variance, c(1L, 2L), mean,
    na.rm = TRUE
  )
  covariance_variance <- pmax(total_variance - monte_carlo_variance, 0)

  # Frontier membership and rank churn under the drawn covariances.
  draw_matrix <- function(draw_index) {
    matrix(
      outcomes[, , draw_index],
      nrow = nrow(directions),
      ncol = length(objective_names),
      dimnames = list(NULL, objective_names)
    )
  }
  membership <- vapply(seq_along(usable), function(draw_index) {
    .dgr_non_dominated(draw_matrix(draw_index))
  }, logical(nrow(directions)))
  if (is.null(dim(membership))) {
    membership <- matrix(membership, nrow = nrow(directions))
  }
  rank_churn <- if (is.null(rank_weights) || nrow(directions) < 2L) {
    rep(NA_real_, length(usable))
  } else {
    reference_rank <- rank(-as.numeric(point_estimate %*% rank_weights))
    vapply(seq_along(usable), function(draw_index) {
      drawn_rank <- rank(-as.numeric(draw_matrix(draw_index) %*% rank_weights))
      suppressWarnings(
        stats::cor(reference_rank, drawn_rank, method = "spearman")
      )
    }, numeric(1L))
  }

  summary_table <- data.table::data.table(
    direction = seq_len(nrow(directions))
  )
  for (trait in trait_cols) {
    data.table::set(summary_table,
      j = paste0("d_", trait),
      value = directions[, trait]
    )
  }
  for (objective in objective_names) {
    index <- match(objective, objective_names)
    data.table::set(summary_table,
      j = paste0("mean_", objective),
      value = covariance_draw_mean[, index]
    )
    data.table::set(summary_table,
      j = paste0("point_", objective),
      value = point_estimate[, index]
    )
    data.table::set(summary_table,
      j = paste0("sd_covariance_", objective),
      value = sqrt(covariance_variance[, index])
    )
    data.table::set(summary_table,
      j = paste0("sd_montecarlo_", objective),
      value = sqrt(monte_carlo_variance[, index])
    )
  }
  data.table::set(summary_table,
    j = "frontier_membership",
    value = rowMeans(membership)
  )

  dominant <- if (mean(covariance_variance, na.rm = TRUE) >
    mean(monte_carlo_variance, na.rm = TRUE)) {
    "covariance estimation"
  } else {
    "Monte Carlo simulation"
  }

  result <- list(
    trait_cols = trait_cols,
    objective_names = objective_names,
    directions = directions,
    point_estimate = point_estimate,
    covariance_draw_mean = covariance_draw_mean,
    genetic_df = genetic_df,
    residual_df = residual_df,
    n_covariance_draws = length(usable),
    n_requested_draws = n_covariance_draws,
    n_rebuild_failures = rebuild_failures,
    n_replicates = n_replicates,
    n_cycles = n_cycles,
    outcomes = outcomes,
    summary = summary_table,
    frontier_membership = rowMeans(membership),
    rank_churn = list(
      spearman = rank_churn,
      mean = if (all(is.na(rank_churn))) {
        NA_real_
      } else {
        mean(rank_churn, na.rm = TRUE)
      },
      minimum = if (all(is.na(rank_churn))) {
        NA_real_
      } else {
        min(rank_churn, na.rm = TRUE)
      },
      weights = rank_weights,
      note = if (is.null(rank_weights)) {
        "Not computed: supply rank_weights to declare a unit-aware scalarisation."
      } else if (nrow(directions) < 2L) {
        "Not computed: rank correlation requires at least two directions."
      } else {
        "Computed from the user-declared scalarisation weights."
      }
    ),
    variance_components = list(
      covariance = covariance_variance,
      monte_carlo = monte_carlo_variance,
      total = total_variance,
      dominant_source = dominant,
      note = paste(
        "The dominant source of uncertainty here is", dominant,
        "error. Monte Carlo error falls as the square root of the replicate",
        "count; covariance-estimation error does not fall with more",
        "simulation and can only be reduced by estimating G from more",
        "independent genetic units. If the covariance component dominates,",
        "running more replicates is wasted effort."
      )
    ),
    seed = seed,
    interpretation = paste(
      "Each covariance draw rebuilds the founder population, because G is the",
      "architecture the traits are built from rather than a matrix the",
      "simulation multiplies by. A direction on the frontier under the point",
      "estimate but in only a minority of draws is not a robust",
      "recommendation, however clean the point-estimate frontier looks."
    )
  )
  class(result) <- c("desiredgainr_covariance_uncertainty", "list")
  result
}

#' @export
print.desiredgainr_covariance_uncertainty <- function(x, ...) {
  cat("<desiredgainr_covariance_uncertainty>\n")
  cat(sprintf(
    "  %d directions, %d usable covariance draws of %d, %d replicates each\n",
    nrow(x$directions), x$n_covariance_draws, x$n_requested_draws,
    x$n_replicates
  ))
  if (x$n_rebuild_failures > 0L) {
    cat(sprintf(
      "  %d draw(s) could not be realised as a trait architecture.\n",
      x$n_rebuild_failures
    ))
  }
  cat(sprintf(
    "  Genetic df: %d   Rank correlation with the point estimate: %.3f (min %.3f)\n",
    x$genetic_df, x$rank_churn$mean, x$rank_churn$minimum
  ))
  cat("  Dominant uncertainty:", x$variance_components$dominant_source, "\n")
  cat("\n  Frontier membership across covariance draws:\n")
  print(round(x$frontier_membership, 3L))
  cat("\n  ", x$variance_components$note, "\n")
  invisible(x)
}
