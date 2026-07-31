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
         call. = FALSE)
  }
  G <- .dgr_covariance(G, trait_cols, "G")
  P <- .dgr_covariance(P, trait_cols, "P")
  b <- .dgr_named_vector(coefficients, trait_cols, "coefficients")
  if (!main_trait %in% trait_cols) {
    stop("main_trait must name one of the traits.", call. = FALSE)
  }

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

  merit_correlation <- NA_real_
  merit_gain <- NA_real_
  merit_sd <- NA_real_
  if (!is.null(aggregate_weights)) {
    a <- .dgr_named_vector(aggregate_weights, trait_cols, "aggregate_weights")
    merit_variance <- as.numeric(crossprod(a, G %*% a))
    if (merit_variance > 0) {
      merit_sd <- sqrt(merit_variance)
      merit_correlation <- as.numeric(crossprod(b, G %*% a)) /
        (index_sd * merit_sd)
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
    index_variance = index_variance,
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
  cat(sprintf("  Delta H (aggregate gain):  %s\n",
              format(x$delta_H, digits = 5)))
  cat(sprintf("  RE vs direct selection on %s: %s\n",
              x$main_trait, format(x$RE, digits = 5)))
  cat(sprintf("  CV of index: %s%%\n", format(x$CV_I, digits = 4)))
  cat("  Expected response per trait:\n")
  print(round(x$expected_response, 4L))
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
