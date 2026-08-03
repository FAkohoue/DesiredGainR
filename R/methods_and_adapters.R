# Consistent S3 methods for fitted objects, plus validation and provenance for
# covariance matrices extracted from upstream mixed-model software.

#' @export
coef.desiredgainr_index <- function(object, ...) {
  object$coefficients
}

#' @export
coef.desired_gain_index <- function(object, ...) {
  object$coefficients
}

#' Summarise a fitted selection index
#'
#' Returns one tidy table per question a reader asks of an index: what the
#' coefficients are, what response is expected, and what the index is worth.
#' The tables are data frames rather than printed text so that they can be
#' written into a report without reformatting.
#'
#' @param object A fitted index from [selection_index()] or
#'   [restricted_index()].
#' @param ... Unused.
#'
#' @return An object of class `desiredgainr_index_summary`.
#'
#' @examples
#' set.seed(1)
#' traits <- c("yield", "protein")
#' G <- matrix(c(1.0, 0.2, 0.2, 0.5), 2, dimnames = list(traits, traits))
#' P <- matrix(c(2.5, 0.4, 0.4, 1.2), 2, dimnames = list(traits, traits))
#' values <- as.data.frame(matrix(
#'   stats::rnorm(40),
#'   ncol = 2, dimnames = list(paste0("g", 1:20), traits)
#' ))
#' fit <- selection_index(
#'   values, traits,
#'   method = "smith_hazel", G = G, P = P,
#'   economic_weights = c(yield = 2, protein = 1), n_select = 5
#' )
#' summary(fit)
#'
#' @export
summary.desiredgainr_index <- function(object, ...) {
  traits <- object$trait_cols
  evaluation <- object$evaluation

  coefficients <- data.table::data.table(
    Trait = traits,
    Coefficient = as.numeric(object$coefficients),
    Effective_weight = if (is.null(object$effective_weights)) {
      NA_real_
    } else {
      as.numeric(object$effective_weights$weights)[seq_along(traits)]
    },
    Aggregate_weight = if (is.null(object$aggregate_weights)) {
      NA_real_
    } else {
      as.numeric(object$aggregate_weights)
    }
  )

  response <- if (is.null(evaluation)) {
    NULL
  } else {
    data.table::data.table(
      Trait = traits,
      Expected_response = as.numeric(evaluation$expected_response),
      Observed_differential = if (is.null(object$observed_differential)) {
        NA_real_
      } else {
        object$observed_differential$Differential
      }
    )
  }

  criteria <- if (is.null(evaluation)) {
    NULL
  } else {
    data.table::data.table(
      Criterion = c(
        "R_HI", "Delta_H", "RE", "h2_index", "accuracy_index", "CV_I"
      ),
      Value = c(
        evaluation$R_HI, evaluation$delta_H, evaluation$RE,
        evaluation$h2_index, evaluation$accuracy_index, evaluation$CV_I
      ),
      Meaning = c(
        "Correlation between the index and net merit",
        "Expected gain in aggregate merit",
        paste(
          "Efficiency relative to direct selection on",
          evaluation$main_trait
        ),
        "Heritability of the index treated as a trait",
        "Correlation between the index and the value it estimates",
        "Coefficient of variation of the index (percent)"
      )
    )
  }

  result <- list(
    method = object$method,
    n_candidates = object$n_candidates,
    n_selected = object$n_selected,
    n_requested = object$n_select,
    selection_intensity = object$selection_intensity,
    transformation = object$transformation,
    coefficients = coefficients,
    response = response,
    criteria = criteria,
    constraint = object$constraint,
    culling_report = object$culling_report
  )
  class(result) <- c("desiredgainr_index_summary", "list")
  result
}

#' @export
print.desiredgainr_index_summary <- function(x, ...) {
  cat("<desiredgainr_index_summary>\n")
  cat("  Method:", x$method, "\n")
  cat(sprintf(
    "  Candidates: %d   Selected: %s   Intensity: %s\n",
    x$n_candidates,
    if (is.null(x$n_selected)) "none" else x$n_selected,
    format(x$selection_intensity, digits = 4)
  ))
  if (!is.null(x$n_requested) && !identical(x$n_requested, x$n_selected)) {
    cat(sprintf(
      "  NOTE: %d were requested but only %d could be selected.\n",
      x$n_requested, x$n_selected
    ))
  }
  cat(
    "  Traits centred:", isTRUE(x$transformation$centred),
    " scaled:", isTRUE(x$transformation$scaled),
    paste0("(", x$transformation$scale_by, ")"), "\n"
  )
  cat("\n  Coefficients:\n")
  print(x$coefficients)
  if (!is.null(x$response)) {
    cat("\n  Response:\n")
    print(x$response)
  }
  if (!is.null(x$criteria)) {
    cat("\n  Criteria:\n")
    print(x$criteria)
  }
  if (!is.null(x$constraint)) {
    cat("\n  Constraint:", x$constraint$method, "\n")
    cat(
      "  Largest violation:",
      format(x$constraint$largest_violation, digits = 3), "\n"
    )
  }
  invisible(x)
}

#' @export
summary.desired_gain_index <- function(object, ...) {
  traits <- object$trait_cols
  required_named <- c(
    "coefficients", "desired_gain", "optimised_d",
    "realised_response"
  )
  malformed <- required_named[!vapply(required_named, function(component) {
    value <- object[[component]]
    is.numeric(value) && all(traits %in% names(value))
  }, logical(1L))]
  if (length(malformed)) {
    stop("The fitted DGSI object has missing or unnamed trait components: ",
      paste(malformed, collapse = ", "), ". Refit it with the current ",
      "DesiredGainR version rather than silently returning empty columns.",
      call. = FALSE
    )
  }
  theoretical <- object$theoretical_response
  result <- list(
    coefficients = data.table::data.table(
      Trait = traits,
      Coefficient = as.numeric(object$coefficients)
    ),
    response = data.table::data.table(
      Trait = traits,
      Requested_gain = as.numeric(object$desired_gain[traits]),
      Optimised_gain = as.numeric(object$optimised_d[traits]),
      Empirical_differential = as.numeric(object$realised_response[traits]),
      Transmitted_response = if (is.null(theoretical)) {
        NA_real_
      } else {
        as.numeric(theoretical$standardised[traits])
      }
    ),
    optimism = object$optimism,
    compatibility = object$covariance_provenance$compatibility,
    eligibility = object$eligibility
  )
  class(result) <- c("desiredgainr_dgsi_summary", "list")
  result
}

#' @export
print.desiredgainr_dgsi_summary <- function(x, ...) {
  cat("<desiredgainr_dgsi_summary>\n")
  cat("  All gains and responses are in candidate standard-deviation units.\n")
  cat("  Empirical_differential is the difference among selected candidates;\n")
  cat("  Transmitted_response is what the next generation inherits.\n\n")
  print(x$response)
  cat("\n  Coefficients:\n")
  print(x$coefficients)
  if (!is.null(x$optimism)) {
    cat("\n  Replicate chosen on:", x$optimism$selection_rule, "\n")
  }
  if (!is.null(x$compatibility)) {
    cat("  P - G status:", x$compatibility$status, "\n")
  }
  invisible(x)
}

# ---------------------------------------------------------------------------
# Upstream covariance validation
# ---------------------------------------------------------------------------

#' Validate and record an extracted genetic covariance matrix
#'
#' Mixed-model packages report multi-trait variance components in different
#' shapes. Except for simple named ASReml component vectors, this function does
#' not extract a fitted object: the caller supplies the intended covariance
#' term and this function validates it. The errors it prevents are not subtle:
#' a covariance in the wrong
#' trait order, a correlation matrix mistaken for a covariance, or a
#' genetic-by-environment block read as a genetic one. This function converts
#' the common formats and validates what it produces.
#'
#' @details
#' # What is checked
#'
#' The trait names are matched to `trait_cols` and reordered, not assumed to
#' already be in that order. Symmetry, finiteness and positive semidefiniteness
#' are verified. Where `P` is supplied as well, the pair is checked for
#' admissibility, because a genetic covariance exceeding the phenotypic one in
#' any direction is what produces heritabilities above one downstream.
#'
#' # Source labels
#'
#' \describe{
#'   \item{`"asreml"`}{A named vector or matrix of variance components from
#'     `summary(fit)$varcomp`, or an explicit covariance matrix.}
#'   \item{`"sommer"`}{A matrix such as `fit$sigma[[term]]`, extracted by the
#'     caller.}
#'   \item{`"breedR"`, `"bglr"`, `"stagewise"`}{A covariance matrix extracted by
#'     the caller; this function validates and orders it.}
#'   \item{`"matrix"`}{Any matrix, for the case where the extraction has
#'     already been done.}
#' }
#'
#' Because the upstream APIs change between versions, this function accepts the
#' extracted matrix rather than the fitted object wherever possible. That keeps
#' the contract stable: the caller is responsible for pulling the right term
#' out of their model, and this function is responsible for checking it is
#' usable.
#'
#' @param x A covariance matrix, or a named vector of variance components.
#' @param trait_cols Trait names, in the order the rest of the analysis uses.
#' @param source Which upstream package produced `x`; see Details.
#' @param estimand What `x` represents, recorded in the provenance.
#' @param P Optional phenotypic covariance, checked for admissibility with `x`.
#' @param is_correlation Whether `x` is a correlation matrix that must be
#'   rescaled by `variances`.
#' @param variances Named trait variances, required when `is_correlation` is
#'   `TRUE`.
#'
#' @return A list of class `desiredgainr_imported_covariance` holding the
#'   ordered matrix, its diagnostics and its provenance.
#'
#' @examples
#' traits <- c("yield", "protein")
#' # A correlation matrix and separate variances, which is how several
#' # packages report multi-trait results.
#' correlation <- matrix(c(1, 0.3, 0.3, 1), 2,
#'   dimnames = list(traits, traits)
#' )
#' import_covariance(
#'   correlation, traits,
#'   source = "matrix",
#'   is_correlation = TRUE, variances = c(yield = 1.4, protein = 0.6)
#' )
#'
#' @seealso [estimate_genetic_covariance()], [matrix_diagnostics()],
#'   [bend_covariance()]
#' @export
import_covariance <- function(
  x,
  trait_cols,
  source = c("matrix", "asreml", "sommer", "breedR", "bglr", "stagewise"),
  estimand = "additive genetic covariance",
  P = NULL,
  is_correlation = FALSE,
  variances = NULL
) {
  source <- match.arg(source)
  if (!is.character(trait_cols) || anyDuplicated(trait_cols)) {
    stop("trait_cols must contain unique trait names.", call. = FALSE)
  }
  p <- length(trait_cols)

  if (is.vector(x) && !is.matrix(x)) {
    # A named vector of components, as ASReml-R reports them. Real component
    # labels often wrap traits in at(), us() or model-term syntax, so match the
    # supplied trait names within the complete label instead of assuming one
    # delimiter convention.
    if (is.null(names(x))) {
      stop("A vector of variance components must be named so that its ",
        "entries can be matched to traits.",
        call. = FALSE
      )
    }
    matrix_form <- matrix(NA_real_, p, p, dimnames = list(trait_cols, trait_cols))
    for (entry in names(x)) {
      hits <- trait_cols[vapply(
        trait_cols, function(trait) grepl(trait, entry, fixed = TRUE),
        logical(1L)
      )]
      # Prefer longer names when one trait name is a substring of another.
      hits <- hits[!vapply(hits, function(hit) {
        any(nchar(hits) > nchar(hit) & grepl(hit, hits, fixed = TRUE))
      }, logical(1L))]
      if (length(hits) == 1L) {
        matrix_form[hits, hits] <- x[[entry]]
      } else if (length(hits) == 2L) {
        matrix_form[hits[1L], hits[2L]] <- x[[entry]]
        matrix_form[hits[2L], hits[1L]] <- x[[entry]]
      }
    }
    if (anyNA(matrix_form)) {
      missing_cells <- which(is.na(matrix_form), arr.ind = TRUE)
      stop(
        "The variance components do not cover every trait pair; ",
        nrow(missing_cells), " entries could not be filled, including ",
        trait_cols[missing_cells[1L, 1L]], " with ",
        trait_cols[missing_cells[1L, 2L]],
        ". Check that the component names contain the trait names.",
        call. = FALSE
      )
    }
    x <- matrix_form
  }

  x <- as.matrix(x)
  if (is.null(dimnames(x)) || is.null(rownames(x))) {
    if (nrow(x) != p) {
      stop("The matrix has no trait names and its dimension (", nrow(x),
        ") does not match the ", p, " traits supplied, so its rows cannot ",
        "be identified.",
        call. = FALSE
      )
    }
    warning("The matrix carried no trait names; the supplied order of ",
      "trait_cols was assumed. Verify it against the fitted model.",
      call. = FALSE
    )
    dimnames(x) <- list(trait_cols, trait_cols)
  }

  absent <- setdiff(trait_cols, rownames(x))
  if (length(absent)) {
    stop("The imported matrix is missing traits: ",
      paste(absent, collapse = ", "),
      call. = FALSE
    )
  }
  # Reordering rather than assuming the order is the single most valuable
  # thing this function does: a silently transposed trait order produces an
  # index that is wrong in a way nothing downstream can detect.
  x <- x[trait_cols, trait_cols, drop = FALSE]

  if (isTRUE(is_correlation)) {
    if (is.null(variances)) {
      stop("variances must be supplied when is_correlation = TRUE.",
        call. = FALSE
      )
    }
    variances <- .dgr_named_vector(variances, trait_cols, "variances")
    if (any(variances <= 0)) {
      stop("Every trait variance must be positive.", call. = FALSE)
    }
    if (max(abs(diag(x) - 1)) > 1e-6) {
      stop("is_correlation = TRUE but the diagonal is not 1; supply the ",
        "matrix as a covariance instead.",
        call. = FALSE
      )
    }
    scaling <- sqrt(variances)
    x <- sweep(sweep(x, 1L, scaling, "*"), 2L, scaling, "*")
    dimnames(x) <- list(trait_cols, trait_cols)
  } else if (max(abs(diag(x) - 1)) < 1e-8 && p > 1L) {
    warning(
      "Every diagonal element is 1. If this is a correlation matrix rather ",
      "than a covariance, set is_correlation = TRUE and supply variances; ",
      "otherwise the trait scales will be wrong throughout the analysis.",
      call. = FALSE
    )
  }

  x <- (x + t(x)) / 2
  if (any(!is.finite(x))) {
    stop("The imported covariance contains non-finite values.", call. = FALSE)
  }
  diagnostics <- matrix_diagnostics(x, "imported covariance")
  if (!diagnostics$positive_definite) {
    warning(
      "The imported covariance is not positive definite (smallest eigenvalue ",
      format(diagnostics$minimum_eigenvalue, digits = 3), "). Repair it with ",
      "bend_covariance() and report the adjustment.",
      call. = FALSE
    )
  }
  if (!is.null(P)) {
    P <- .dgr_covariance(P, trait_cols, "P")
    .dgr_check_compatible(x, P, "import_covariance()")
  }

  result <- list(
    covariance = x,
    trait_cols = trait_cols,
    diagnostics = diagnostics,
    provenance = list(
      source = source,
      estimand = estimand,
      rescaled_from_correlation = isTRUE(is_correlation),
      reordered_to = trait_cols,
      checked_against_P = !is.null(P)
    )
  )
  class(result) <- c("desiredgainr_imported_covariance", "list")
  result
}

#' @export
print.desiredgainr_imported_covariance <- function(x, ...) {
  cat("<desiredgainr_imported_covariance>\n")
  cat("  Source:", x$provenance$source, "\n")
  cat("  Estimand:", x$provenance$estimand, "\n")
  cat("  Traits:", paste(x$trait_cols, collapse = ", "), "\n")
  cat(sprintf(
    "  Condition number: %.4g   Positive definite: %s\n",
    x$diagnostics$condition_number,
    if (x$diagnostics$positive_definite) "yes" else "no"
  ))
  if (isTRUE(x$provenance$rescaled_from_correlation)) {
    cat("  Rescaled from a correlation matrix using the supplied variances.\n")
  }
  invisible(x)
}
