#' Propagate sampling error in the covariance matrices into an index
#'
#' Every quantity a selection index reports treats \eqn{\mathbf{G}} and
#' \eqn{\mathbf{P}} as known constants. They are not. They are estimates, often
#' from a trial with far fewer independent families than the covariance matrix
#' has free parameters, and the index coefficients are a non-linear function of
#' them. An index whose coefficients cannot be distinguished from zero is not a
#' different index from one whose coefficients are sharply determined, yet the
#' point estimate presents them identically.
#'
#' This function resamples the covariance matrices from their approximate
#' sampling distribution, refits the index on each draw, and reports how much
#' the answer moves.
#'
#' @details
#' # What is resampled
#'
#' A covariance matrix estimated from \eqn{\nu} independent multivariate normal
#' vectors has a Wishart sampling distribution, so draws are taken as
#'
#' \deqn{\mathbf{G}^{*} \sim \mathcal{W}_p(\nu_g, \mathbf{G}/\nu_g),}{
#' G* ~ Wishart_p(nu_g, G / nu_g),}
#'
#' which has expectation \eqn{\mathbf{G}}. When `residual_df` is supplied the
#' residual matrix \eqn{\mathbf{E} = \mathbf{P} - \mathbf{G}} is drawn the same
#' way and the phenotypic matrix reassembled as
#' \eqn{\mathbf{P}^{*} = \mathbf{G}^{*} + \mathbf{E}^{*}}. That structure
#' matters: \eqn{\mathbf{G}} and \eqn{\mathbf{P}} are not independently
#' estimated quantities, because both are built from the same records, and
#' drawing them independently would break a correlation that damps the movement
#' in \eqn{\mathbf{P}^{-1}\mathbf{G}}.
#'
#' # Choosing the degrees of freedom
#'
#' This is the input that governs the width of every interval below, and it is
#' the one users get wrong.
#'
#' `genetic_df` is the number of **independent genetic units** that informed the
#' genetic covariance, not the number of plots and not the number of
#' observations. In a half-sib trial it is on the order of the number of
#' families; in a diallel, the number of parents; in a clonally replicated
#' trial, the number of distinct clones. Replication within a genotype sharpens
#' the residual matrix, not the genetic one.
#'
#' Supplying the number of plots where the number of families was meant will
#' understate every interval by roughly the square root of the replication.
#' When the design does not yield a clean answer, run the function twice at the
#' plausible extremes and report both.
#'
#' # What the Wishart assumption buys and what it costs
#'
#' It is an approximation. A restricted maximum likelihood estimate from an
#' unbalanced trial with a relationship matrix is not a sample covariance from
#' independent normal vectors, and its true sampling distribution has no closed
#' form. The Wishart matches the mean exactly, matches the variance of a
#' balanced design, and enforces positive semidefiniteness on every draw, which
#' a naive normal perturbation of the elements does not.
#'
#' Where the fitting software reports an asymptotic variance-covariance matrix
#' of the variance components, that is a better basis than this. Where it does
#' not, which is the common case, this is a defensible substitute provided the
#' result is described as what it is.
#'
#' # What is reported
#'
#' The coefficient intervals answer how well determined the index is. The
#' second block answers the question that actually matters, which is what the
#' estimation error costs.
#'
#' For `smith_hazel` and `base`, that cost is the classical relative
#' efficiency: under each drawn truth, the accuracy of the fitted index divided
#' by the accuracy of the index that would have been optimal for that truth.
#'
#' For `pesek_baker` and `yamada` it is sharper. The Pesek-Baker coefficients
#' \eqn{\mathbf{b} = \mathbf{G}^{-1}\mathbf{d}} produce a correlated response
#' proportional to \eqn{\mathbf{G}\mathbf{b} = \mathbf{d}} exactly, so an index
#' fitted on the true \eqn{\mathbf{G}} delivers the desired gains in the
#' requested proportions by construction. Any departure under a drawn truth is
#' therefore attributable to estimation error alone, and `alignment`, the
#' cosine between the achieved response \eqn{\mathbf{G}^{*}\hat{\mathbf{b}}}
#' and the desired gains \eqn{\mathbf{d}}, measures it directly. An alignment
#' interval whose lower bound sits well below one means the desired gains were
#' specified more precisely than the data can deliver them.
#'
#' @param index A fitted object from [selection_index()] using one of
#'   `"smith_hazel"`, `"base"`, `"pesek_baker"` or `"yamada"`.
#' @param genetic_df Degrees of freedom for the genetic covariance matrix. See
#'   the guidance above; this is the number of independent genetic units.
#' @param residual_df Degrees of freedom for the residual covariance matrix
#'   \eqn{\mathbf{P} - \mathbf{G}}. When `NULL`, `P` is held fixed and the
#'   reported uncertainty is attributable to `G` alone.
#' @param n_draws Number of resampling draws.
#' @param level Width of the reported intervals.
#' @param seed Optional seed. The random number generator state is restored on
#'   exit.
#'
#' @return An object of class `desiredgainr_uncertainty` holding the
#'   coefficient summary, the estimation-cost summary, the rank stability, and
#'   the draws themselves.
#'
#' @examples
#' set.seed(1)
#' traits <- c("yield", "protein")
#' G <- matrix(c(1.0, 0.2, 0.2, 0.5), 2, dimnames = list(traits, traits))
#' P <- matrix(c(2.5, 0.4, 0.4, 1.2), 2, dimnames = list(traits, traits))
#' values <- as.data.frame(matrix(
#'   stats::rnorm(80),
#'   ncol = 2, dimnames = list(paste0("G", 1:40), traits)
#' ))
#' fit <- selection_index(
#'   values, traits,
#'   method = "pesek_baker", G = G, P = P,
#'   desired_gains = c(yield = 1, protein = 0.5), n_select = 8
#' )
#' index_uncertainty(fit, genetic_df = 30, residual_df = 120, n_draws = 200)
#'
#' @references
#' Hayes, J.F. and Hill, W.G. (1981) Modification of estimates of parameters in
#' the construction of genetic selection indices. *Biometrics* 37, 483-493.
#'
#' Sales, J. and Hill, W.G. (1976) Effect of sampling errors on efficiency of
#' selection indices. *Animal Production* 22, 1-17.
#'
#' @seealso [selection_index()], [bend_covariance()], [matrix_diagnostics()],
#'   [predict.desiredgainr_index()]
#' @export
index_uncertainty <- function(
  index,
  genetic_df,
  residual_df = NULL,
  n_draws = 500L,
  level = 0.95,
  seed = NULL
) {
  if (!inherits(index, "desiredgainr_index")) {
    stop("index must be an object returned by selection_index().",
      call. = FALSE
    )
  }
  if (!identical(index$strategy, "index")) {
    stop("Resampling applies to the covariance-based families only; ",
      "method = '", index$method, "' derives no coefficients from G.",
      call. = FALSE
    )
  }
  if (is.null(index$G)) {
    stop("index was fitted without G, so there is no covariance estimate ",
      "whose sampling error could be propagated.",
      call. = FALSE
    )
  }
  if (is.null(index$scaled_values) || is.null(index$objective)) {
    stop("index was produced by an earlier version of selection_index() that ",
      "did not retain the values needed to refit. Refit the index.",
      call. = FALSE
    )
  }
  n_draws <- .dgr_positive_integer(n_draws, "n_draws")
  if (!is.numeric(level) || length(level) != 1L ||
    level <= 0 || level >= 1) {
    stop("level must lie strictly between 0 and 1.", call. = FALSE)
  }

  traits <- index$trait_cols
  p <- length(traits)
  G <- index$G
  P <- index$P
  X <- index$scaled_values
  objective <- index$objective
  method <- index$method
  b_hat <- index$coefficients
  gain_method <- method %in% c("pesek_baker", "yamada")
  if (is.null(P)) {
    # `base` and `pesek_baker` can be fitted without P. The efficiency and
    # response calculations still need one, so fall back to the candidates'
    # covariance and say so, rather than silently reporting nothing.
    P <- stats::cov(X)
    dimnames(P) <- list(traits, traits)
    P_note <- paste(
      "The index was fitted without P, so the candidates' own covariance was",
      "used for the response calculations. It is not a population phenotypic",
      "covariance, and the reported cost of estimation error inherits that."
    )
  } else {
    P_note <- NULL
  }

  genetic_df <- .dgr_wishart_df(genetic_df, p, "genetic_df")
  if (!is.null(residual_df)) {
    if (is.null(index$P)) {
      stop("residual_df was supplied but the index was fitted without P, so ",
        "there is no residual covariance to resample.",
        call. = FALSE
      )
    }
    residual_df <- .dgr_wishart_df(residual_df, p, "residual_df")
  } else if (!is.null(index$P)) {
    # With P held fixed, a drawn G can exceed it in some direction, which is
    # not a possible pair of population values. Such draws are counted and
    # excluded rather than silently contributing impossible coefficients, and
    # the mode is described as local sensitivity analysis in the result.
    .dgr_check_compatible(G, P, "index_uncertainty()")
  }

  if (!is.null(seed)) {
    if (!exists(".Random.seed", envir = globalenv())) stats::runif(1L)
    saved <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", saved, envir = globalenv()), add = TRUE)
    set.seed(seed)
  }

  G_draws <- .dgr_wishart_draws(G, genetic_df, n_draws, "G")
  E_draws <- NULL
  residual_note <- NULL
  if (!is.null(residual_df)) {
    E <- P - G
    dimnames(E) <- dimnames(G)
    residual_diagnostics <- matrix_diagnostics(E, "P - G")
    if (residual_diagnostics$minimum_eigenvalue <=
      1e-10 * max(1, residual_diagnostics$maximum_eigenvalue)) {
      residual_note <- paste(
        "P - G is not positive definite, which means the supplied G and P are",
        "mutually inconsistent: some linear combination of traits has a",
        "genetic variance at least as large as its phenotypic variance.",
        "The residual matrix was bent before resampling; check the estimates."
      )
      warning(residual_note, call. = FALSE)
      E <- bend_covariance(E, method = "asrgenomics")$G
    }
    E_draws <- .dgr_wishart_draws(E, residual_df, n_draws, "P - G")
  }

  coefficient_draws <- matrix(
    NA_real_,
    nrow = n_draws, ncol = p, dimnames = list(NULL, traits)
  )
  metric <- rep(NA_real_, n_draws)
  rank_correlation <- rep(NA_real_, n_draws)
  overlap <- rep(NA_real_, n_draws)
  scores_hat <- as.numeric(X %*% b_hat)
  n_select <- index$n_select
  selected_hat <- if (!is.null(n_select)) {
    order(-scores_hat, index$candidate_id)[seq_len(n_select)]
  } else {
    NULL
  }
  d_direction <- if (gain_method) objective / sqrt(sum(objective^2)) else NULL

  n_inadmissible <- 0L
  for (draw in seq_len(n_draws)) {
    G_star <- G_draws[, , draw]
    dimnames(G_star) <- list(traits, traits)
    P_star <- if (is.null(E_draws)) {
      # P fixed: the drawn G may exceed it in some direction, giving a pair
      # that could not both be true. Those draws are discarded.
      residual_values <- eigen(P - G_star,
        symmetric = TRUE,
        only.values = TRUE
      )$values
      if (min(residual_values) < -1e-8 * max(1, max(abs(diag(P))))) {
        n_inadmissible <- n_inadmissible + 1L
        next
      }
      P
    } else {
      # P* = G* + E* is admissible by construction, since E* is a covariance
      # matrix and P* - G* = E*.
      star <- G_star + E_draws[, , draw]
      dimnames(star) <- list(traits, traits)
      star
    }
    refit <- tryCatch(
      suppressWarnings(
        .dgr_index_coefficients(method, G_star, P_star, objective, X, traits)
      ),
      error = function(e) NULL
    )
    if (is.null(refit)) next
    b_star <- refit$coefficients
    coefficient_draws[draw, ] <- b_star

    # The cost of estimation error: what the fitted index achieves when this
    # draw is taken as the truth, against what was achievable under it.
    metric[draw] <- if (gain_method) {
      response <- as.numeric(G_star %*% b_hat)
      denominator <- sqrt(sum(response^2))
      if (denominator > 0) sum(response * d_direction) / denominator else NA_real_
    } else {
      .dgr_accuracy_ratio(b_hat, b_star, G_star, P_star, objective)
    }

    scores_star <- as.numeric(X %*% b_star)
    rank_correlation[draw] <- suppressWarnings(
      stats::cor(scores_hat, scores_star, method = "spearman")
    )
    if (!is.null(selected_hat)) {
      selected_star <- order(-scores_star, index$candidate_id)[seq_len(n_select)]
      overlap[draw] <- length(intersect(selected_hat, selected_star)) / n_select
    }
  }

  converged <- stats::complete.cases(coefficient_draws)
  if (sum(converged) < 2L) {
    stop("Fewer than two draws produced usable coefficients, which means G ",
      "is too ill-conditioned for resampling at genetic_df = ",
      genetic_df, ". Inspect matrix_diagnostics(G) first.",
      call. = FALSE
    )
  }
  probabilities <- c((1 - level) / 2, 1 - (1 - level) / 2)

  coefficients <- data.table::data.table(
    Trait = traits,
    Estimate = as.numeric(b_hat),
    SD = apply(coefficient_draws, 2L, stats::sd, na.rm = TRUE),
    Lower = apply(coefficient_draws, 2L, stats::quantile,
      probs = probabilities[1L], na.rm = TRUE
    ),
    Upper = apply(coefficient_draws, 2L, stats::quantile,
      probs = probabilities[2L], na.rm = TRUE
    ),
    Sign_stability = vapply(seq_len(p), function(j) {
      column <- coefficient_draws[, j]
      mean(sign(column) == sign(b_hat[j]), na.rm = TRUE)
    }, numeric(1L))
  )
  data.table::set(
    coefficients,
    j = "Relative_SD",
    value = ifelse(
      abs(coefficients$Estimate) > 0,
      coefficients$SD / abs(coefficients$Estimate),
      NA_real_
    )
  )

  cost <- list(
    quantity = if (gain_method) "alignment" else "relative_efficiency",
    mean = mean(metric, na.rm = TRUE),
    lower = unname(stats::quantile(metric, probabilities[1L], na.rm = TRUE)),
    upper = unname(stats::quantile(metric, probabilities[2L], na.rm = TRUE)),
    draws = metric,
    interpretation = if (gain_method) {
      paste(
        "Cosine between the response the fitted index achieves under each",
        "drawn covariance matrix and the requested desired gains. An index",
        "fitted on the true G attains 1 exactly, so any shortfall is",
        "estimation error."
      )
    } else {
      paste(
        "Accuracy of the fitted index under each drawn covariance matrix,",
        "divided by the accuracy of the index that would have been optimal",
        "for it. Bounded above by 1."
      )
    }
  )

  result <- list(
    method = method,
    trait_cols = traits,
    genetic_df = genetic_df,
    residual_df = residual_df,
    n_draws = n_draws,
    n_usable = sum(converged),
    n_inadmissible = n_inadmissible,
    mode = if (is.null(residual_df)) {
      "local sensitivity analysis (P held fixed)"
    } else {
      "joint resampling of G and the residual"
    },
    mode_note = if (is.null(residual_df)) {
      paste(
        "P was held fixed while G was resampled, so this is a local",
        "sensitivity analysis rather than a sample from the joint sampling",
        "distribution. Some drawn G exceed the fixed P in one or more",
        "directions, giving a pair that could not both be true; those draws",
        "were discarded and are counted in n_inadmissible. Supply residual_df",
        "to resample G and the residual jointly, which is admissible by",
        "construction."
      )
    } else {
      NULL
    },
    level = level,
    coefficients = coefficients,
    estimation_cost = cost,
    rank_stability = list(
      spearman_mean = mean(rank_correlation, na.rm = TRUE),
      spearman_lower = unname(stats::quantile(
        rank_correlation, probabilities[1L],
        na.rm = TRUE
      )),
      selection_overlap_mean = if (is.null(selected_hat)) {
        NA_real_
      } else {
        mean(overlap, na.rm = TRUE)
      },
      selection_overlap_lower = if (is.null(selected_hat)) {
        NA_real_
      } else {
        unname(stats::quantile(overlap, probabilities[1L], na.rm = TRUE))
      },
      n_select = n_select
    ),
    coefficient_draws = coefficient_draws,
    residual_note = residual_note,
    phenotypic_note = P_note,
    provenance = paste(
      "Wishart resampling of G", if (is.null(residual_df)) {
        "with P held fixed"
      } else {
        "and of P - G, recombined as P* = G* + E*"
      },
      "at", genetic_df, "genetic degrees of freedom. The interval widths are",
      "governed by that number; verify it describes independent genetic units",
      "rather than plots."
    )
  )
  class(result) <- c("desiredgainr_uncertainty", "list")
  result
}

#' @export
print.desiredgainr_uncertainty <- function(x, ...) {
  cat("<desiredgainr_uncertainty>\n")
  cat("  Method:", x$method, "\n")
  cat(sprintf(
    "  Draws: %d usable of %d  Genetic df: %d%s\n",
    x$n_usable, x$n_draws, x$genetic_df,
    if (is.null(x$residual_df)) {
      "  (P held fixed)"
    } else {
      paste0("  Residual df: ", x$residual_df)
    }
  ))
  cat(sprintf("  %.0f%% intervals\n", 100 * x$level))
  cat("  Mode:", x$mode, "\n")
  if (x$n_inadmissible > 0L) {
    cat(sprintf(
      "  Discarded %d draw(s) where P - G was not a covariance matrix.\n",
      x$n_inadmissible
    ))
  }
  cat("\n  Coefficients:\n")
  print(x$coefficients)
  cat(sprintf(
    "\n  %s: %.3f [%.3f, %.3f]\n",
    x$estimation_cost$quantity, x$estimation_cost$mean,
    x$estimation_cost$lower, x$estimation_cost$upper
  ))
  cat(sprintf(
    "  Rank correlation with the point estimate: %.3f (lower %.3f)\n",
    x$rank_stability$spearman_mean, x$rank_stability$spearman_lower
  ))
  if (!is.na(x$rank_stability$selection_overlap_mean)) {
    cat(sprintf(
      "  Selected set retained: %.1f%% of %d (lower %.1f%%)\n",
      100 * x$rank_stability$selection_overlap_mean,
      x$rank_stability$n_select,
      100 * x$rank_stability$selection_overlap_lower
    ))
  }
  for (note in c(x$mode_note, x$phenotypic_note, x$residual_note)) {
    cat("\n  Note:", note, "\n")
  }
  invisible(x)
}

#' Validate degrees of freedom for a Wishart draw
#' @noRd
.dgr_wishart_df <- function(df, p, name) {
  if (!is.numeric(df) || length(df) != 1L || !is.finite(df)) {
    stop(name, " must be a single finite number.", call. = FALSE)
  }
  df <- as.integer(round(df))
  if (df < p) {
    stop(
      name, " is ", df, " but the covariance matrix has ", p, " traits. ",
      "A Wishart draw on fewer degrees of freedom than traits is singular, ",
      "and a covariance matrix estimated from fewer independent units than ",
      "traits carries no usable information about the correlations.",
      call. = FALSE
    )
  }
  if (df < 3 * p) {
    warning(
      name, " = ", df, " for ", p, " traits. The covariance matrix has ",
      p * (p + 1) / 2, " free parameters, so the resulting intervals will be ",
      "wide. That width is a property of the estimate, not of this function.",
      call. = FALSE
    )
  }
  df
}

#' Draw Wishart matrices with the supplied matrix as the expectation
#' @noRd
.dgr_wishart_draws <- function(M, df, n_draws, name) {
  diagnostics <- matrix_diagnostics(M, name)
  if (!diagnostics$positive_definite) {
    stop(name, " must be positive definite to be resampled; its smallest ",
      "eigenvalue is ", format(diagnostics$minimum_eigenvalue, digits = 3),
      ". Repair it with bend_covariance() first.",
      call. = FALSE
    )
  }
  # Scaling by the degrees of freedom makes the expectation of each draw equal
  # to M, since E[W_p(df, S)] = df * S.
  stats::rWishart(n_draws, df = df, Sigma = M / df)
}

#' Accuracy of one coefficient vector relative to another under a given truth
#' @noRd
.dgr_accuracy_ratio <- function(b_hat, b_star, G, P, a) {
  accuracy <- function(b) {
    numerator <- as.numeric(t(b) %*% G %*% a)
    denominator <- sqrt(as.numeric(t(b) %*% P %*% b)) *
      sqrt(as.numeric(t(a) %*% G %*% a))
    if (!is.finite(denominator) || denominator <= 0) {
      return(NA_real_)
    }
    numerator / denominator
  }
  best <- accuracy(b_star)
  if (!is.finite(best) || best <= 0) {
    return(NA_real_)
  }
  accuracy(b_hat) / best
}
