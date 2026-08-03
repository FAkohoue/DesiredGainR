# Standard criteria for evaluating and comparing selection indices.
#
# These are the four criteria used by Rahimi and Debnath (2023) and reported by
# the established selection-index software (Genes, MIX, RIndSel, SelAction),
# together with the per-trait expected responses.

#' Evaluate a selection index against the standard criteria
#'
#' A selection index cannot be judged from its coefficients. Rahimi and Debnath
#' (2023) evaluate every index on four criteria, and this function reproduces
#' them so that alternative index families and alternative objectives can be
#' compared on identical data.
#'
#' @details
#' Let \eqn{b} be the index coefficients, \eqn{a} the aggregate weights
#' defining net merit \eqn{H = a^\mathsf{T}g}, and \eqn{k} the standardised
#' selection intensity. The criteria are:
#'
#' \describe{
#'   \item{\eqn{R_{HI}}}{The correlation between the index and net merit,
#'     \eqn{b^\mathsf{T}Ga / \sqrt{b^\mathsf{T}Pb \cdot a^\mathsf{T}Ga}}.
#'     Maximising it maximises the response in aggregate merit.}
#'   \item{\eqn{\Delta H}}{The expected gain in aggregate merit,
#'     \eqn{k R_{HI} \sqrt{a^\mathsf{T}Ga}}.}
#'   \item{\eqn{\Delta_j}}{The expected response for each trait,
#'     \eqn{k (Gb)_j / \sqrt{b^\mathsf{T}Pb}}.}
#'   \item{RE}{Efficiency relative to direct selection on a single main trait.}
#'   \item{\eqn{CV_I}}{The coefficient of variation of the index scores.}
#'   \item{\eqn{h^2_I}}{The heritability of the index treated as a trait,
#'     \eqn{b^\mathsf{T}Gb / b^\mathsf{T}Pb}, whose square root is the
#'     correlation between the index score and the genetic value it estimates.
#'     This is the only accuracy available for the desired-gain families
#'     without first supplying economic weights.}
#' }
#'
#' # Interpreting relative efficiency
#'
#' Relative efficiency below 1 is not a defect. Every value reported by Rahimi
#' and Debnath was below 1, because a multi-trait index deliberately trades
#' response in the main trait for response in the remainder of the objective.
#' Hence RE should be read alongside the per-trait responses that direct
#' selection would sacrifice, and not as a pass-or-fail test.
#'
#' @param coefficients Named numeric vector of index coefficients, in the same
#'   trait space as `G` and `P`.
#' @param G Genetic variance-covariance matrix, named by trait.
#' @param P Phenotypic variance-covariance matrix, named by trait.
#' @param aggregate_weights Named weights defining net merit. When `NULL`, the
#'   criteria that depend on net merit are returned as `NA`.
#' @param selection_intensity Standardised selection intensity.
#' @param main_trait Trait against which relative efficiency is computed.
#' @param scores Optional vector of realised index scores, used for
#'   \eqn{CV_I}.
#'
#' @return An object of class `desiredgainr_evaluation`.
#'
#' @references
#' Rahimi M, Debnath S (2023). *Scientific Reports* 13:18977.
#' \doi{10.1038/s41598-023-46368-6}
#'
#' @seealso [selection_index()]
#' @export
evaluate_index <- function(
  coefficients,
  G,
  P,
  aggregate_weights = NULL,
  selection_intensity = NA_real_,
  main_trait = names(coefficients)[1L],
  scores = NULL
) {
  trait_cols <- names(coefficients)
  if (is.null(trait_cols) || anyDuplicated(trait_cols)) {
    stop("coefficients must be a numeric vector with unique trait names.",
      call. = FALSE
    )
  }
  G <- .dgr_covariance(G, trait_cols, "G")
  P <- .dgr_covariance(P, trait_cols, "P")
  b <- .dgr_named_vector(coefficients, trait_cols, "coefficients")
  if (!main_trait %in% trait_cols) {
    stop("main_trait must name one of the traits.", call. = FALSE)
  }
  # Every criterion below is a ratio of a genetic quantity to a phenotypic one.
  # If P - G is not a covariance matrix, those ratios can exceed their natural
  # bounds and the result looks valid while being impossible.
  .dgr_check_compatible(G, P, "evaluate_index()")

  index_variance <- as.numeric(crossprod(b, P %*% b))
  if (!is.finite(index_variance) || index_variance <= 0) {
    stop("The index has non-positive variance.", call. = FALSE)
  }
  index_sd <- sqrt(index_variance)

  # Per-trait expected response, i k Gb / sd(I).
  response <- as.numeric(selection_intensity * (G %*% b) / index_sd)
  names(response) <- trait_cols

  # Direct selection on the main trait alone gives k G[m, m] / sqrt(P[m, m]).
  direct <- selection_intensity * G[main_trait, main_trait] /
    sqrt(P[main_trait, main_trait])
  relative_efficiency <- if (is.finite(direct) && direct != 0) {
    response[[main_trait]] / direct
  } else {
    NA_real_
  }

  # Heritability of the index itself. The index is a linear combination of
  # phenotypes and so has a genetic variance b'Gb and a phenotypic variance
  # b'Pb of its own. Their ratio is the heritability of the index as a trait,
  # and its square root is the correlation between the index score and the
  # genetic value the index estimates, which is the accuracy of selection on
  # the index irrespective of how net merit is defined. Reported because it is
  # the one accuracy measure available for the desired-gain families without
  # first supplying economic weights.
  index_genetic_variance <- as.numeric(crossprod(b, G %*% b))
  index_heritability <- if (index_genetic_variance >= 0) {
    # Compatibility has been validated, so a ratio above one can only be
    # rounding error and is clamped. Anything larger raises an error there.
    .dgr_clamp_correlation(
      index_genetic_variance / index_variance, "The index heritability"
    )
  } else {
    NA_real_
  }
  index_accuracy <- if (is.finite(index_heritability) &&
    index_heritability >= 0) {
    sqrt(index_heritability)
  } else {
    NA_real_
  }

  merit_correlation <- NA_real_
  merit_gain <- NA_real_
  merit_sd <- NA_real_
  if (!is.null(aggregate_weights)) {
    a <- .dgr_named_vector(aggregate_weights, trait_cols, "aggregate_weights")
    merit_variance <- as.numeric(crossprod(a, G %*% a))
    if (merit_variance > 0) {
      merit_sd <- sqrt(merit_variance)
      merit_correlation <- .dgr_clamp_correlation(
        as.numeric(crossprod(b, G %*% a)) / (index_sd * merit_sd),
        "R_HI, the correlation between the index and net merit"
      )
      merit_gain <- selection_intensity * merit_correlation * merit_sd
    }
  }

  # The coefficient of variation of an index is only interpretable when the
  # index has a meaningful non-zero location. Standardising or centring the
  # traits puts the index mean at zero, so the ratio diverges and must be
  # withheld rather than reported as an enormous number.
  coefficient_of_variation <- NA_real_
  cv_note <- "Not requested."
  if (!is.null(scores) && length(scores) > 1L) {
    mean_score <- mean(scores)
    sd_score <- stats::sd(scores)
    if (is.finite(sd_score) && sd_score > 0 &&
      abs(mean_score) > 1e-6 * sd_score) {
      coefficient_of_variation <- 100 * sd_score / abs(mean_score)
      cv_note <- "Reported on the supplied index scores."
    } else {
      cv_note <- paste(
        "Withheld: the index mean is negligible relative to its standard",
        "deviation, which happens whenever the traits were centred or",
        "standardised. A coefficient of variation is undefined there."
      )
    }
  }

  result <- list(
    R_HI = merit_correlation,
    delta_H = merit_gain,
    expected_response = response,
    RE = relative_efficiency,
    CV_I = coefficient_of_variation,
    CV_I_note = cv_note,
    h2_index = index_heritability,
    accuracy_index = index_accuracy,
    index_variance = index_variance,
    index_genetic_variance = index_genetic_variance,
    index_sd = index_sd,
    merit_sd = merit_sd,
    selection_intensity = selection_intensity,
    main_trait = main_trait,
    response_table = data.table::data.table(
      Trait = trait_cols,
      Coefficient = as.numeric(b),
      Expected_response = as.numeric(response)
    ),
    interpretation = paste(
      "Relative efficiency below 1 means the index yields less response in",
      main_trait, "than direct selection on it, which is the intended",
      "trade-off of multi-trait selection rather than a defect."
    )
  )
  class(result) <- c("desiredgainr_evaluation", "list")
  result
}

#' @export
print.desiredgainr_evaluation <- function(x, ...) {
  cat("<desiredgainr_evaluation>\n")
  cat(sprintf("  R_HI (index vs net merit): %s\n", format(x$R_HI, digits = 5)))
  cat(sprintf(
    "  Delta H (aggregate gain):  %s\n",
    format(x$delta_H, digits = 5)
  ))
  cat(sprintf(
    "  RE vs direct selection on %s: %s\n",
    x$main_trait, format(x$RE, digits = 5)
  ))
  cat(sprintf(
    "  Index heritability: %s  Accuracy: %s\n",
    format(x$h2_index, digits = 4),
    format(x$accuracy_index, digits = 4)
  ))
  cat(sprintf("  CV of index: %s%%\n", format(x$CV_I, digits = 4)))
  cat("  Expected response per trait:\n")
  print(round(x$expected_response, 4L))
  invisible(x)
}

#' Bend a covariance matrix to positive definiteness
#'
#' Multi-trait restricted maximum likelihood frequently returns a genetic
#' covariance matrix that is not positive definite, because the number of
#' variance components grows quadratically in the number of traits while the
#' information available to estimate them does not. Such a matrix cannot be
#' inverted, so every index built on it fails.
#'
#' Rejecting it outright, which is what DesiredGainR previously did, sends the
#' user to an ad hoc repair outside the package and outside the audit trail.
#' Bending performs the repair explicitly and records what was changed.
#'
#' @details
#' # The algorithm
#'
#' The default `"correlation"` method repairs the correlation matrix by
#' eigenvalue clipping and then restores the original trait standard
#' deviations. It therefore preserves every marginal variance, which is the
#' appropriate constraint for a trait covariance. `"asrgenomics"` retains the
#' earlier relationship-matrix-oriented `G.tuneup()` route for reproducibility.
#'
#' `G.tuneup()` bends by calling `Matrix::nearPD()`, which implements Higham's
#' (2002) alternating-projections algorithm. That iteration alternately
#' projects onto the set of positive semidefinite matrices and onto the set
#' satisfying the remaining constraints, converging on the nearest admissible
#' matrix in the Frobenius norm. It is a well-tested implementation used across
#' animal and plant breeding, which is the reason to delegate to it rather than
#' write a local substitute.
#'
#' # What G.tuneup fixes, and what this function therefore cannot expose
#'
#' `G.tuneup()` hard-wires two of `nearPD()`'s arguments. Both are consequences
#' for a trait covariance rather than options, so they are stated here rather
#' than passed through.
#'
#' `keepDiag = FALSE`. **The trait variances are not preserved.** Bending will
#' move them, so the heritabilities implied by the bent matrix are not the ones
#' supplied. Compare `diag()` before and after, which `adjustment` reports as
#' `max_variance_change`, and if that number is material to the argument being
#' made, say so when reporting the result.
#'
#' `posd.tol = 1e-2`. Eigenvalues are floored at a hundredth of the largest, so
#' the bent matrix has a condition number near 100 or below. For the
#' \eqn{n \times n} genomic relationship matrix `G.tuneup()` was written for,
#' near-singular by construction, that is the right target. For a
#' \eqn{p \times p} trait covariance it is aggressive: condition numbers of
#' \eqn{10^3} to \eqn{10^4} are ordinary and genuine there rather than
#' artefacts, so bending will smooth away correlation structure the data
#' support. Run [matrix_diagnostics()] on both matrices and check how far the
#' conditioning was moved.
#'
#' The `eig_tol` argument is passed through to `G.tuneup()`, but note that it
#' governs `nearPD()`'s `eig.tol` and not `posd.tol`, so it does not lift the
#' condition-number cap.
#'
#' # Interpreting the result
#'
#' Report the returned `adjustment` alongside any result computed from a bent
#' matrix. A large adjustment means the estimate contained little information
#' about the correlation structure, and conclusions that depend on that
#' structure should be treated accordingly. Bending is a repair of an estimate,
#' not an estimate.
#'
#' @param M Symmetric covariance matrix. Trait names are required on both
#'   dimensions and must match, because `G.tuneup()` validates them.
#' @param method `"correlation"` (default) to preserve trait variances, or
#'   `"asrgenomics"` to reproduce the earlier aggressive relationship-matrix
#'   repair.
#' @param eig_tol Relative eigenvalue tolerance passed to `G.tuneup()`. See the
#'   note above on why this does not control the condition-number cap.
#' @param digits Rounding applied by `G.tuneup()` to the bent matrix.
#'
#' @return A list of class `desiredgainr_bent_covariance` holding the bent
#'   matrix, the diagnostics before and after, the reciprocal condition numbers
#'   reported by `G.tuneup()`, and the size of the adjustment.
#'
#' @references
#' Higham, N.J. (2002) Computing the nearest correlation matrix: a problem from
#' finance. *IMA Journal of Numerical Analysis* 22, 329-343.
#' \doi{10.1093/imanum/22.3.329}
#'
#' Nazarian, A. and Gezan, S.A. (2016) GenoMatrix: a software package for
#' pedigree-based and genomic prediction analyses on complex traits.
#' *Journal of Heredity* 107, 372-379.
#'
#' @examples
#' \donttest{
#' traits <- c("t1", "t2", "t3")
#' # A matrix implying a correlation structure that cannot exist: t1 and t2 are
#' # strongly positively correlated, t2 and t3 likewise, yet t1 and t3 are
#' # strongly negatively correlated.
#' M <- matrix(
#'   c(
#'     1.0, 0.9, 0.9,
#'     0.9, 1.0, 0.9,
#'     0.9, 0.9, 1.0
#'   ), 3,
#'   dimnames = list(traits, traits)
#' )
#' M["t1", "t3"] <- M["t3", "t1"] <- -0.9
#' matrix_diagnostics(M)$positive_definite
#'
#' if (requireNamespace("ASRgenomics", quietly = TRUE)) {
#'   bent <- bend_covariance(M)
#'   bent
#'   # The variances moved, because G.tuneup() does not hold them fixed.
#'   diag(bent$G)
#' }
#' }
#'
#' @seealso [matrix_diagnostics()], [estimate_genetic_covariance()]
#' @export
bend_covariance <- function(
  M, method = c("correlation", "asrgenomics"), eig_tol = 1e-6, digits = 8
) {
  method <- match.arg(method)
  if (identical(method, "asrgenomics") &&
    !requireNamespace("ASRgenomics", quietly = TRUE)) {
    stop(
      "method = 'asrgenomics' requires the ASRgenomics package, which ",
      "provides G.tuneup().\n",
      "  Install it with install.packages('ASRgenomics').",
      call. = FALSE
    )
  }
  if (!is.matrix(M) || nrow(M) != ncol(M)) {
    stop("M must be a square matrix.", call. = FALSE)
  }
  if (any(!is.finite(M))) {
    stop("M must contain only finite values.", call. = FALSE)
  }
  if (!is.numeric(eig_tol) || length(eig_tol) != 1L ||
    !is.finite(eig_tol) || eig_tol <= 0) {
    stop("eig_tol must be a single positive number.", call. = FALSE)
  }
  if (is.null(rownames(M)) || is.null(colnames(M))) {
    stop("M must carry trait names on both dimensions, because G.tuneup() ",
      "requires them.",
      call. = FALSE
    )
  }
  if (!identical(rownames(M), colnames(M))) {
    stop("The row and column names of M must match.", call. = FALSE)
  }

  original <- (M + t(M)) / 2
  dimnames(original) <- dimnames(M)
  before <- matrix_diagnostics(original, "supplied matrix")
  if (before$maximum_eigenvalue <= 0) {
    stop("M has no positive eigenvalue and cannot be bent into a covariance ",
      "matrix.",
      call. = FALSE
    )
  }

  if (identical(method, "correlation")) {
    if (any(diag(original) <= 0)) {
      stop("method = 'correlation' requires positive trait variances.",
        call. = FALSE
      )
    }
    trait_sd <- sqrt(diag(original))
    correlation <- stats::cov2cor(original)
    decomposition <- eigen(correlation, symmetric = TRUE)
    floor_value <- eig_tol * max(decomposition$values)
    clipped <- pmax(decomposition$values, floor_value)
    repaired_correlation <- decomposition$vectors %*%
      (clipped * t(decomposition$vectors))
    repaired_correlation <- stats::cov2cor(repaired_correlation)
    bent <- outer(trait_sd, trait_sd) * repaired_correlation
    dimnames(bent) <- dimnames(original)
    fit_rcn <- c(
      before = if (is.finite(before$condition_number)) {
        1 / before$condition_number
      } else {
        0
      },
      after = 1 / matrix_diagnostics(bent, "bent matrix")$condition_number
    )
    determinant_before <- det(original)
    raised_threshold <- floor_value
    provenance <- paste(
      "The trait correlation matrix was repaired by clipping eigenvalues",
      "below", format(floor_value, digits = 3), "and renormalising its",
      "diagonal to one. The original marginal trait variances were restored",
      "exactly. This is a transparent repair of an estimate, not an estimate."
    )
  } else {
    # G.tuneup() insists that exactly one of blend, bend and align is requested,
    # and messages its diagnostics; those are captured in the result instead.
    fit <- ASRgenomics::G.tuneup(
      G = original, bend = TRUE, eig.tol = eig_tol, rcn = TRUE,
      digits = digits, sparseform = FALSE, determinant = TRUE,
      message = FALSE
    )
    bent <- as.matrix(fit$Gb)
    bent <- (bent + t(bent)) / 2
    dimnames(bent) <- dimnames(original)
    fit_rcn <- c(before = fit$rcn0, after = fit$rcnb)
    determinant_before <- fit$det0
    raised_threshold <- 1e-2 * before$maximum_eigenvalue
    provenance <- paste(
      "Bent by ASRgenomics::G.tuneup(bend = TRUE), which calls",
      "Matrix::nearPD(). G.tuneup() fixes keepDiag = FALSE and posd.tol =",
      "1e-2, so this relationship-matrix repair changes trait variances and",
      "caps the condition number near 100."
    )
  }
  after <- matrix_diagnostics(bent, "bent matrix")

  if (!after$positive_definite) {
    warning(
      "Bending did not produce a positive definite matrix. The supplied ",
      "estimate is too far from admissible for G.tuneup() to repair; ",
      "re-estimate the covariance on fewer traits, or supply the matrix from ",
      "a model that constrains it to be admissible.",
      call. = FALSE
    )
  }

  variance_change <- max(abs(diag(bent) - diag(original)))
  result <- list(
    G = bent,
    original = original,
    before = before,
    after = after,
    method = method,
    reciprocal_condition = fit_rcn,
    determinant_before = determinant_before,
    n_eigenvalues_raised = sum(
      before$eigenvalues < raised_threshold
    ),
    adjustment = list(
      max_absolute_change = max(abs(bent - original)),
      # Undefined when the supplied matrix has a non-positive diagonal, which
      # happens when a residual matrix is bent rather than a covariance.
      max_correlation_change = tryCatch(
        suppressWarnings(max(abs(
          stats::cov2cor(bent) - stats::cov2cor(original)
        ))),
        error = function(e) NA_real_
      ),
      max_variance_change = variance_change,
      max_relative_variance_change = if (all(diag(original) > 0)) {
        max(abs(diag(bent) - diag(original)) / diag(original))
      } else {
        NA_real_
      },
      frobenius_relative = sqrt(sum((bent - original)^2)) /
        sqrt(sum(original^2))
    ),
    provenance = provenance
  )
  class(result) <- c("desiredgainr_bent_covariance", "list")
  result
}

#' @export
print.desiredgainr_bent_covariance <- function(x, ...) {
  cat("<desiredgainr_bent_covariance>\n")
  cat("  Method:", x$method, "\n")
  cat(sprintf(
    "  Smallest eigenvalue: %.3g -> %.3g\n",
    x$before$minimum_eigenvalue, x$after$minimum_eigenvalue
  ))
  cat(sprintf(
    "  Condition number: %.4g -> %.4g\n",
    x$before$condition_number, x$after$condition_number
  ))
  cat(sprintf(
    "  Reciprocal condition number: %.4g -> %.4g\n",
    x$reciprocal_condition[["before"]], x$reciprocal_condition[["after"]]
  ))
  cat(sprintf(
    "  Largest correlation change: %s\n",
    format(x$adjustment$max_correlation_change, digits = 4)
  ))
  cat(sprintf(
    "  Largest variance change: %s (%s of the supplied variance)\n",
    format(x$adjustment$max_variance_change, digits = 4),
    format(x$adjustment$max_relative_variance_change, digits = 3)
  ))
  cat(sprintf(
    "  Relative Frobenius distance moved: %s\n",
    format(x$adjustment$frobenius_relative, digits = 4)
  ))
  cat(
    "  Variances preserved:",
    if (identical(x$method, "correlation")) "yes" else "no", "\n"
  )
  invisible(x)
}

#' Report the numerical conditioning of a covariance matrix
#'
#' Selection-index coefficients are obtained by inverting covariance matrices,
#' and an ill-conditioned matrix produces coefficients that are numerically
#' meaningless while remaining superficially plausible. Rahimi and Debnath
#' (2023) report \eqn{R_{HI} = 0.0018} for a desired-gain index on
#' unstandardised maize traits spanning plant height in centimetres, grain
#' counts, and yield; that collapse is consistent with an inversion performed
#' on a severely ill-conditioned genetic covariance matrix.
#'
#' Therefore inspect conditioning before trusting any index built on original
#' trait scales, and standardise the traits when the reciprocal condition
#' number is small.
#'
#' @param M Symmetric covariance matrix.
#' @param name Object name used in messages.
#'
#' @return A list giving the eigenvalues, the reciprocal condition number, the
#'   numerical rank, and a logical flag for positive definiteness.
#'
#' @export
matrix_diagnostics <- function(M, name = "matrix") {
  if (!is.matrix(M) || nrow(M) != ncol(M)) {
    stop(name, " must be a square matrix.", call. = FALSE)
  }
  M <- (M + t(M)) / 2
  values <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
  largest <- max(abs(values))
  reciprocal_condition <- if (largest > 0) min(abs(values)) / largest else 0
  list(
    name = name,
    dimension = nrow(M),
    eigenvalues = values,
    minimum_eigenvalue = min(values),
    maximum_eigenvalue = max(values),
    reciprocal_condition = reciprocal_condition,
    condition_number = if (reciprocal_condition > 0) {
      1 / reciprocal_condition
    } else {
      Inf
    },
    numerical_rank = sum(values > 1e-10 * max(1, largest)),
    positive_definite = min(values) > 1e-10 * max(1, largest)
  )
}
