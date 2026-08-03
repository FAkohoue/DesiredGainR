#' Estimate a working genetic covariance matrix
#'
#' Estimates a covariance matrix for use in a desired-gain selection index
#' when a covariance matrix from a fitted multivariate genetic model is not
#' available. The estimand depends on `method`; the function never relabels the
#' covariance of predictions as total genetic covariance.
#'
#' @param values A data frame or numeric matrix containing one row per
#'   genotype and one column per trait.
#' @param trait_cols Character vector naming and ordering the trait columns.
#'   For a matrix, `colnames(values)` are used when `trait_cols` is omitted.
#' @param method Estimation method. See Details.
#' @param prediction_error_covariance Prediction-error covariance (PEV)
#'   information for `method = "pev_corrected"`. Supply either one common
#'   trait-by-trait matrix or a trait-by-trait-by-genotype array.
#' @param prediction_se Prediction standard errors for
#'   `method = "se_diagonal_corrected"`. Supply a named vector with one value
#'   per trait or a genotype-by-trait matrix or data frame.
#' @param relationship_matrix Genotype relationship matrix for
#'   `method = "relationship_adjusted"`.
#' @param ids Genotype identifiers in the row order of `values`. Required when
#'   `relationship_matrix` is supplied; the matrix is matched by its dimnames.
#' @param eigen_tolerance Relative eigenvalue tolerance used for positive
#'   semidefiniteness and numerical-rank checks.
#' @param symmetry_tolerance Relative tolerance used to check matrix symmetry.
#'
#' @details
#' `method = "prediction_covariance"` returns the sample covariance of the
#' supplied predictions. For BLUPs or GEBVs this is a covariance of predicted
#' genetic values and is generally shrunken relative to total genetic
#' covariance.
#'
#' `method = "pev_corrected"` uses
#' \deqn{\widehat{G} = \operatorname{Cov}(\widehat{u}) +
#' \overline{\operatorname{PEV}}.}
#' This follows the BLUP identity
#' \eqn{\operatorname{PEV} = \operatorname{Var}(u) -
#' \operatorname{Var}(\widehat{u})}. It approximates total genetic covariance
#' when the supplied values are compatible multivariate BLUPs, their PEV
#' matrices are on the same scale, and the candidate sample and prediction
#' model support that identity. Full cross-trait PEV matrices are required.
#'
#' `method = "se_diagonal_corrected"` adds mean squared prediction standard
#' errors to the diagonal of the prediction covariance. Because standard
#' errors contain no cross-trait prediction-error covariance, off-diagonal
#' entries remain covariances of predictions. The result is therefore a
#' partial working approximation, not a complete estimate of total genetic
#' covariance.
#'
#' `method = "adjusted_means_surrogate"` returns the covariance of
#' across-environment adjusted means and labels it as a working surrogate for
#' genetic covariance. Covarrubias-Pazaran (2021) recommends this for
#' operational use when no genetic covariance matrix is available, on the
#' grounds that a mixed-model adjustment removes most environmental and design
#' variation. However, the adjusted means still carry estimation error, so the
#' result is not an estimate of additive genetic covariance and will differ
#' from the covariance governing response in the next generation. Hence it
#' should be reported as a surrogate and revisited once a fitted multi-trait
#' model becomes available.
#'
#' `method = "relationship_adjusted"` estimates the covariance of predictions
#' after accounting for a supplied relationship matrix:
#' \deqn{\widehat{G}_{pred} =
#' X_c^\mathsf{T} K^+ X_c / \operatorname{rank}(K).}
#' It does not undo BLUP or GEBV shrinkage. Estimating total genetic covariance
#' for related candidates requires prediction-error covariance from the fitted
#' multivariate model, including the relevant relationship structure.
#'
#' Raw phenotypes contain residual and environmental covariance. Their
#' covariance cannot identify genetic covariance without a genetic model,
#' replication, or relationship information, and should not be passed to this
#' function as genetic predictions. Best linear unbiased estimates (BLUEs) and
#' adjusted means occupy an intermediate position: they are not genetic
#' covariance estimates either, but `method = "adjusted_means_surrogate"`
#' makes their use explicit and auditable rather than leaving breeders to
#' substitute them silently.
#'
#' @references
#' Covarrubias-Pazaran G (2021). *Practical implementation of selection
#' indices.* CGIAR Excellence in Breeding.
#'
#' @return An object of class `desiredgainr_covariance_estimate` containing
#'   `G`, the estimated matrix; `estimand`; `provenance`; `assumptions`;
#'   component matrices; and numerical diagnostics.
#' @export
estimate_genetic_covariance <- function(
  values,
  trait_cols = NULL,
  method = c(
    "prediction_covariance",
    "pev_corrected",
    "se_diagonal_corrected",
    "relationship_adjusted",
    "adjusted_means_surrogate"
  ),
  prediction_error_covariance = NULL,
  prediction_se = NULL,
  relationship_matrix = NULL,
  ids = NULL,
  eigen_tolerance = 1e-8,
  symmetry_tolerance = 1e-8
) {
  method <- match.arg(method)
  if (!is.numeric(eigen_tolerance) || length(eigen_tolerance) != 1L ||
    !is.finite(eigen_tolerance) || eigen_tolerance <= 0) {
    stop("eigen_tolerance must be a positive finite scalar.", call. = FALSE)
  }
  if (!is.numeric(symmetry_tolerance) ||
    length(symmetry_tolerance) != 1L ||
    !is.finite(symmetry_tolerance) || symmetry_tolerance <= 0) {
    stop("symmetry_tolerance must be a positive finite scalar.", call. = FALSE)
  }

  prepared <- .dgr_genetic_covariance_values(values, trait_cols)
  X <- prepared$values
  trait_cols <- prepared$trait_cols
  n <- nrow(X)
  p <- ncol(X)
  prediction_covariance <- stats::cov(X)
  dimnames(prediction_covariance) <- list(trait_cols, trait_cols)
  correction <- matrix(
    0, p, p,
    dimnames = list(trait_cols, trait_cols)
  )
  estimand <- "covariance of supplied genetic predictions"
  source <- "sample covariance of supplied genetic predictions"
  assumptions <- c(
    "Rows represent the target population of candidate genotypes.",
    "Trait columns are genetic predictions on a common, interpretable scale."
  )
  relationship_rank <- NA_integer_

  if (method == "pev_corrected") {
    if (is.null(prediction_error_covariance)) {
      stop(
        paste(
          "prediction_error_covariance is required for",
          "method = 'pev_corrected'."
        ),
        call. = FALSE
      )
    }
    correction <- .dgr_average_pev(
      prediction_error_covariance,
      trait_cols = trait_cols,
      n = n,
      eigen_tolerance = eigen_tolerance,
      symmetry_tolerance = symmetry_tolerance
    )
    estimand <- "approximate total genetic covariance"
    source <- paste(
      "sample covariance of supplied genetic predictions plus",
      "mean prediction-error covariance"
    )
    assumptions <- c(
      assumptions,
      paste(
        "Values are compatible multivariate BLUPs and PEV matrices are from",
        "the same fitted model and scale."
      ),
      paste(
        "The BLUP covariance identity is applicable to the target candidate",
        "sample."
      )
    )
  } else if (method == "se_diagonal_corrected") {
    if (is.null(prediction_se)) {
      stop(
        "prediction_se is required for method = 'se_diagonal_corrected'.",
        call. = FALSE
      )
    }
    se <- .dgr_prediction_se(prediction_se, trait_cols, n)
    diag(correction) <- colMeans(se^2)
    estimand <- paste(
      "partially corrected working covariance;",
      "off-diagonals are covariances of predictions"
    )
    source <- paste(
      "sample covariance of supplied genetic predictions plus mean squared",
      "prediction standard errors on the diagonal"
    )
    assumptions <- c(
      assumptions,
      paste(
        "Prediction standard errors are from the same fitted model and scale",
        "as the supplied predictions."
      ),
      paste(
        "Cross-trait prediction-error covariances are unavailable;",
        "off-diagonal genetic covariances are not corrected."
      )
    )
  } else if (method == "adjusted_means_surrogate") {
    estimand <- paste(
      "working surrogate for genetic covariance;",
      "covariance of across-environment adjusted means"
    )
    source <- paste(
      "covariance of across-environment adjusted means used as a practical",
      "surrogate for G, following Covarrubias-Pazaran (2021)"
    )
    assumptions <- c(
      assumptions,
      paste(
        "Values are across-environment adjusted means from a mixed model, so",
        "that most environmental and design variation has been removed."
      ),
      paste(
        "The residual estimation error remaining in the adjusted means is",
        "small relative to the genetic covariance among them."
      ),
      paste(
        "This surrogate is recommended for operational use by the CGIAR",
        "Excellence in Breeding guideline when no genetic covariance matrix",
        "is available. It is not an estimate of additive genetic covariance,",
        "and it will differ from the covariance governing response in the",
        "next generation."
      )
    )
  } else if (method == "relationship_adjusted") {
    relationship <- .dgr_relationship_inverse(
      relationship_matrix = relationship_matrix,
      ids = ids,
      n = n,
      eigen_tolerance = eigen_tolerance,
      symmetry_tolerance = symmetry_tolerance
    )
    Xc <- sweep(X, 2L, colMeans(X), "-")
    prediction_covariance <- crossprod(
      Xc, relationship$inverse %*% Xc
    ) / relationship$rank
    prediction_covariance <- (
      prediction_covariance + t(prediction_covariance)
    ) / 2
    dimnames(prediction_covariance) <- list(trait_cols, trait_cols)
    relationship_rank <- relationship$rank
    source <- paste(
      "relationship-adjusted covariance of supplied genetic predictions",
      "using a spectral Moore-Penrose inverse"
    )
    assumptions <- c(
      assumptions,
      paste(
        "The supplied relationship matrix describes covariance among rows",
        "and is on the intended relationship scale."
      ),
      "This adjustment does not recover variance removed by prediction shrinkage."
    )
  }

  G <- prediction_covariance + correction
  G <- (G + t(G)) / 2
  dimnames(G) <- list(trait_cols, trait_cols)
  .dgr_check_psd(G, "estimated G", tolerance = eigen_tolerance)
  eigenvalues <- eigen(G, symmetric = TRUE, only.values = TRUE)$values

  result <- list(
    G = G,
    method = method,
    estimand = estimand,
    provenance = source,
    assumptions = assumptions,
    prediction_covariance = prediction_covariance,
    correction = correction,
    diagnostics = list(
      n_genotypes = n,
      n_traits = p,
      minimum_eigenvalue = min(eigenvalues),
      maximum_eigenvalue = max(eigenvalues),
      numerical_rank = sum(
        eigenvalues > eigen_tolerance * max(1, max(abs(eigenvalues)))
      ),
      relationship_rank = relationship_rank,
      full_cross_trait_pev = identical(method, "pev_corrected")
    )
  )
  class(result) <- c("desiredgainr_covariance_estimate", "list")
  result
}

#' @export
print.desiredgainr_covariance_estimate <- function(x, ...) {
  cat("DesiredGainR covariance estimate\n")
  cat("  Method:", x$method, "\n")
  cat("  Estimand:", x$estimand, "\n")
  cat(
    "  Genotypes:", x$diagnostics$n_genotypes,
    " Traits:", x$diagnostics$n_traits, "\n"
  )
  invisible(x)
}

.dgr_genetic_covariance_values <- function(values, trait_cols) {
  if (is.data.frame(values)) {
    if (is.null(trait_cols) || !length(trait_cols) ||
      anyDuplicated(trait_cols)) {
      stop(
        "trait_cols must contain unique trait names for data-frame values.",
        call. = FALSE
      )
    }
    absent <- setdiff(trait_cols, names(values))
    if (length(absent)) {
      stop(
        "values is missing trait columns: ",
        paste(absent, collapse = ", "),
        call. = FALSE
      )
    }
    X <- as.matrix(as.data.frame(values)[, trait_cols, drop = FALSE])
  } else if (is.matrix(values)) {
    X <- values
    if (is.null(trait_cols)) {
      trait_cols <- colnames(X)
    }
    if (is.null(trait_cols) || length(trait_cols) != ncol(X) ||
      anyDuplicated(trait_cols)) {
      stop(
        paste(
          "A matrix must have unique column names or be accompanied by",
          "trait_cols."
        ),
        call. = FALSE
      )
    }
    if (!is.null(colnames(X))) {
      absent <- setdiff(trait_cols, colnames(X))
      if (length(absent)) {
        stop(
          "values dimnames do not contain: ",
          paste(absent, collapse = ", "),
          call. = FALSE
        )
      }
      X <- X[, trait_cols, drop = FALSE]
    }
  } else {
    stop("values must be a data frame or numeric matrix.", call. = FALSE)
  }
  storage.mode(X) <- "double"
  if (nrow(X) < 2L) {
    stop("values must contain at least two genotypes.", call. = FALSE)
  }
  if (any(!is.finite(X))) {
    stop("values must contain only finite numeric trait values.", call. = FALSE)
  }
  colnames(X) <- trait_cols
  list(values = X, trait_cols = trait_cols)
}

.dgr_average_pev <- function(
  value,
  trait_cols,
  n,
  eigen_tolerance,
  symmetry_tolerance
) {
  p <- length(trait_cols)
  if (is.matrix(value)) {
    matrices <- list(value)
  } else if (is.array(value) && length(dim(value)) == 3L &&
    all(dim(value)[1:2] == p) && dim(value)[3] == n) {
    matrices <- lapply(seq_len(n), function(i) value[, , i, drop = TRUE])
  } else {
    stop(
      paste(
        "prediction_error_covariance must be a trait-by-trait matrix or a",
        "trait-by-trait-by-genotype array."
      ),
      call. = FALSE
    )
  }
  matrices <- lapply(seq_along(matrices), function(i) {
    matrix_i <- .dgr_named_symmetric_matrix(
      matrices[[i]],
      names = trait_cols,
      label = paste0("prediction_error_covariance[, , ", i, "]"),
      symmetry_tolerance = symmetry_tolerance
    )
    .dgr_check_psd(
      matrix_i,
      paste0("prediction_error_covariance[, , ", i, "]"),
      tolerance = eigen_tolerance
    )
    matrix_i
  })
  Reduce("+", matrices) / length(matrices)
}

.dgr_prediction_se <- function(value, trait_cols, n) {
  p <- length(trait_cols)
  if (is.numeric(value) && is.null(dim(value))) {
    if (is.null(names(value)) || !all(trait_cols %in% names(value))) {
      stop(
        "A prediction_se vector must be named for every trait.",
        call. = FALSE
      )
    }
    value <- matrix(value[trait_cols], nrow = 1L)
  } else if (is.data.frame(value)) {
    absent <- setdiff(trait_cols, names(value))
    if (length(absent)) {
      stop(
        "prediction_se is missing: ",
        paste(absent, collapse = ", "),
        call. = FALSE
      )
    }
    value <- as.matrix(
      as.data.frame(value)[, trait_cols, drop = FALSE]
    )
  } else if (is.matrix(value)) {
    if (!is.null(colnames(value))) {
      absent <- setdiff(trait_cols, colnames(value))
      if (length(absent)) {
        stop(
          "prediction_se dimnames do not contain: ",
          paste(absent, collapse = ", "),
          call. = FALSE
        )
      }
      value <- value[, trait_cols, drop = FALSE]
    }
  } else {
    stop(
      "prediction_se must be a named vector, data frame, or matrix.",
      call. = FALSE
    )
  }
  if (!nrow(value) %in% c(1L, n) || ncol(value) != p) {
    stop(
      "prediction_se must have one row or one row per genotype and one column per trait.",
      call. = FALSE
    )
  }
  storage.mode(value) <- "double"
  if (any(!is.finite(value)) || any(value < 0)) {
    stop("prediction_se must contain finite non-negative values.",
      call. = FALSE
    )
  }
  colnames(value) <- trait_cols
  value
}

.dgr_relationship_inverse <- function(
  relationship_matrix,
  ids,
  n,
  eigen_tolerance,
  symmetry_tolerance
) {
  if (is.null(relationship_matrix)) {
    stop(
      "relationship_matrix is required for method = 'relationship_adjusted'.",
      call. = FALSE
    )
  }
  if (is.null(ids) || length(ids) != n || anyNA(ids) ||
    anyDuplicated(ids)) {
    stop(
      "ids must contain one unique, non-missing identifier per genotype.",
      call. = FALSE
    )
  }
  if (!is.matrix(relationship_matrix) ||
    any(dim(relationship_matrix) != n) ||
    is.null(rownames(relationship_matrix)) ||
    is.null(colnames(relationship_matrix))) {
    stop(
      paste(
        "relationship_matrix must be a named square matrix with one row and",
        "column per genotype."
      ),
      call. = FALSE
    )
  }
  ids <- as.character(ids)
  if (!all(ids %in% rownames(relationship_matrix)) ||
    !all(ids %in% colnames(relationship_matrix))) {
    stop(
      "relationship_matrix dimnames must contain every genotype id.",
      call. = FALSE
    )
  }
  K <- relationship_matrix[ids, ids, drop = FALSE]
  K <- .dgr_named_symmetric_matrix(
    K, ids, "relationship_matrix", symmetry_tolerance
  )
  decomposition <- eigen(K, symmetric = TRUE)
  largest <- max(1, max(abs(decomposition$values)))
  if (min(decomposition$values) < -eigen_tolerance * largest) {
    stop("relationship_matrix must be positive semidefinite.",
      call. = FALSE
    )
  }
  retained <- decomposition$values > eigen_tolerance * largest
  if (!any(retained)) {
    stop("relationship_matrix has zero numerical rank.", call. = FALSE)
  }
  vectors <- decomposition$vectors[, retained, drop = FALSE]
  inverse <- vectors %*%
    (t(vectors) / decomposition$values[retained])
  list(inverse = inverse, rank = sum(retained))
}

.dgr_named_symmetric_matrix <- function(
  value,
  names,
  label,
  symmetry_tolerance
) {
  p <- length(names)
  if (!is.matrix(value) || any(dim(value) != p)) {
    stop(label, " must be a ", p, " x ", p, " matrix.", call. = FALSE)
  }
  if (xor(is.null(rownames(value)), is.null(colnames(value)))) {
    stop(label, " must have both row and column names, or neither.",
      call. = FALSE
    )
  }
  if (!is.null(rownames(value))) {
    if (!all(names %in% rownames(value)) ||
      !all(names %in% colnames(value))) {
      stop(label, " dimnames do not contain all required names.",
        call. = FALSE
      )
    }
    value <- value[names, names, drop = FALSE]
  }
  storage.mode(value) <- "double"
  if (any(!is.finite(value))) {
    stop(label, " must contain only finite values.", call. = FALSE)
  }
  scale <- max(1, max(abs(value)))
  if (max(abs(value - t(value))) > symmetry_tolerance * scale) {
    stop(label, " must be symmetric.", call. = FALSE)
  }
  value <- (value + t(value)) / 2
  dimnames(value) <- list(names, names)
  value
}
