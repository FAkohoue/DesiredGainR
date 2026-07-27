#' Fit and evaluate a quadratic genomic selection index
#'
#' `run_qgsi()` implements the quadratic genomic selection index (QGSI)
#' \deqn{I_{qg,i} = w^\mathsf{T}\hat{\gamma}_i +
#' \hat{\gamma}_i^\mathsf{T}W\hat{\gamma}_i,}
#' where \eqn{\hat{\gamma}_i} is the vector of trait genomic estimated
#' breeding values (GEBVs), \eqn{w} contains the linear economic weights, and
#' \eqn{W} is the symmetric matrix of economic weights for squared and
#' cross-product genetic merit. Desired gains are not economic weights and are
#' not used to construct either \eqn{w} or \eqn{W}.
#'
#' The genomic covariance matrix \eqn{\Gamma} may be supplied, estimated by
#' \deqn{\hat{\Gamma} = g^{-1}\hat{\gamma}^{\mathsf{T}}\hat{\gamma},}
#' or estimated with a supplied genomic relationship matrix \eqn{\Phi} as
#' \deqn{\hat{\Gamma} =
#' g^{-1}\hat{\gamma}^{\mathsf{T}}\Phi^{-1}\hat{\gamma}.}
#' A Moore-Penrose inverse is used, and reported, when \eqn{\Phi} is
#' positive-semidefinite but rank deficient.
#'
#' @param init_data Candidate identifiers and metadata.
#' @param gebv_data Candidate identifiers and one GEBV column per trait.
#' @param trait_cols Unique names of the trait GEBV columns.
#' @param linear_weights Named linear economic-weight vector \eqn{w}.
#' @param W Symmetric matrix of economic weights for squared and cross-product
#'   terms. For traits `i` and `j`, the coefficient multiplying
#'   \eqn{\gamma_i\gamma_j} is `2 * W[i, j]`.
#' @param id_col Candidate identifier column.
#' @param reference_gebv_data Optional reference-population GEBVs used to
#'   estimate \eqn{\Gamma}. The candidate GEBVs are used when this is `NULL`.
#' @param Gamma Optional genomic covariance matrix in the original trait
#'   coordinates. When supplied, it is transformed consistently with the GEBVs.
#' @param relationship_matrix Optional named genomic relationship matrix
#'   \eqn{\Phi} for the reference genotypes. It is used only when `Gamma` is
#'   `NULL`; row and column names must match reference genotype identifiers.
#' @param true_G Optional true genetic covariance matrix. This is normally
#'   available only in simulation. Supplying it enables MSPE and squared
#'   index-merit correlation calculations; these quantities are not fabricated
#'   from empirical data when `true_G` is absent.
#' @param lower_is_better Traits for which smaller original values are
#'   favourable. Their GEBVs are multiplied by -1. `linear_weights` and `W`
#'   must describe this favourable-direction trait space.
#' @param center_traits Whether to centre GEBVs using reference means. The QGSI
#'   theory assumes zero-mean GEBVs, so `TRUE` is the default.
#' @param scale_traits Whether to divide GEBVs by reference standard deviations.
#'   Supplied weights must refer to the resulting scale.
#' @param missing_policy Explicit missing-value policy.
#' @param n_select Optional exact number of top-ranked candidates to select.
#' @param selection_proportion Optional selected proportion in `(0, 1]`.
#'   Specify at most one of `n_select` and `selection_proportion`.
#' @param return_contributions Whether to return candidate-specific linear and
#'   quadratic contribution tables.
#' @param relationship_tolerance Relative eigenvalue tolerance used for the
#'   Moore-Penrose inverse of `relationship_matrix`.
#' @param symmetry_tolerance Relative tolerance used when validating symmetric
#'   matrices.
#' @param debug Whether to print progress messages.
#'
#' @return An object of class `quadratic_genomic_index`. It contains rankings,
#'   optional selection decisions, score contributions, the genomic covariance
#'   estimate and provenance, and QGSI theoretical parameters. Expected
#'   per-trait gains are the linear-regression gains in Supplementary Equation
#'   16 of Ceron-Rojas et al. (2026); `observed_selection_differential` instead
#'   reports the selected candidates' GEBV shift and is not labelled realised
#'   genetic gain.
#'
#' @references
#' Ceron-Rojas JJ, Montesinos-Lopez OA, Montesinos-Lopez A, et al. (2026).
#' Nonlinear genomic selection index accelerates multi-trait crop improvement.
#' \emph{Nature Communications}, 17, 1991.
#' \doi{10.1038/s41467-026-69890-3}
#'
#' @examples
#' traits <- c("yield", "height")
#' gebv <- data.frame(
#'   GenoID = paste0("G", 1:8),
#'   yield = c(-1.2, -0.7, -0.2, 0.1, 0.3, 0.7, 1.0, 1.4),
#'   height = c(0.8, 0.3, -0.1, -0.4, 0.5, -0.7, -0.2, -0.6)
#' )
#' W <- matrix(
#'   c(0.10, -0.02, -0.02, -0.05),
#'   2, dimnames = list(traits, traits)
#' )
#' fit <- run_qgsi(
#'   init_data = gebv["GenoID"],
#'   gebv_data = gebv,
#'   trait_cols = traits,
#'   linear_weights = c(yield = 1, height = 0.2),
#'   W = W,
#'   n_select = 2
#' )
#' fit$theoretical_parameters
#' fit$expected_gain_per_trait
#' fit$observed_selection_differential
#'
#' @export
run_qgsi <- function(
    init_data,
    gebv_data,
    trait_cols,
    linear_weights,
    W,
    id_col = "GenoID",
    reference_gebv_data = NULL,
    Gamma = NULL,
    relationship_matrix = NULL,
    true_G = NULL,
    lower_is_better = NULL,
    center_traits = TRUE,
    scale_traits = FALSE,
    missing_policy = c("error", "complete_cases", "mean_impute"),
    n_select = NULL,
    selection_proportion = NULL,
    return_contributions = TRUE,
    relationship_tolerance = 1e-8,
    symmetry_tolerance = sqrt(.Machine$double.eps),
    debug = FALSE
) {
  missing_policy <- match.arg(missing_policy)
  .dgr_qg_traits(trait_cols)
  .dgr_qg_scalar(
    relationship_tolerance, "relationship_tolerance",
    lower = 0, lower_open = TRUE
  )
  .dgr_qg_scalar(
    symmetry_tolerance, "symmetry_tolerance",
    lower = 0, lower_open = FALSE
  )
  if (missing(W)) {
    stop(
      paste(
        "W is required. It must contain scientifically or economically",
        "justified squared and cross-product weights."
      ),
      call. = FALSE
    )
  }
  if (!is.null(Gamma) && !is.null(relationship_matrix)) {
    stop(
      "Supply either Gamma or relationship_matrix, not both.",
      call. = FALSE
    )
  }

  prepared <- .dgr_prepare_values(
    init_data = init_data,
    cand_data = gebv_data,
    ref_data = reference_gebv_data,
    validation_data = NULL,
    trait_cols = trait_cols,
    id_col = id_col,
    lower_is_better = lower_is_better,
    scale_traits = scale_traits,
    centre_traits = center_traits,
    missing_policy = missing_policy
  )
  X <- prepared$candidate_matrix
  Xref <- prepared$reference_matrix
  n_candidates <- nrow(X)

  weights <- .dgr_named_vector(
    linear_weights, trait_cols, "linear_weights"
  )
  W <- .dgr_qg_symmetric_matrix(
    W, trait_cols, "W", symmetry_tolerance
  )

  transform <- diag(
    prepared$direction / prepared$scale,
    nrow = length(trait_cols)
  )
  dimnames(transform) <- list(trait_cols, trait_cols)

  if (is.null(Gamma)) {
    gamma_fit <- .dgr_qg_estimate_gamma(
      X = Xref,
      ids = prepared$reference_ids,
      relationship_matrix = relationship_matrix,
      trait_cols = trait_cols,
      tolerance = relationship_tolerance,
      symmetry_tolerance = symmetry_tolerance
    )
    Gamma <- gamma_fit$Gamma
    gamma_provenance <- gamma_fit$provenance
  } else {
    Gamma <- .dgr_qg_symmetric_matrix(
      Gamma, trait_cols, "Gamma", symmetry_tolerance
    )
    Gamma <- transform %*% Gamma %*% transform
    dimnames(Gamma) <- list(trait_cols, trait_cols)
    gamma_provenance <- list(
      source = "user supplied",
      equation = "not estimated by DesiredGainR",
      reference_n = NA_integer_,
      relationship_rank = NA_integer_,
      relationship_dimension = NA_integer_
    )
  }
  .dgr_check_psd(Gamma, "Gamma")

  G_analysis <- NULL
  if (!is.null(true_G)) {
    G_analysis <- .dgr_qg_symmetric_matrix(
      true_G, trait_cols, "true_G", symmetry_tolerance
    )
    G_analysis <- transform %*% G_analysis %*% transform
    dimnames(G_analysis) <- list(trait_cols, trait_cols)
    .dgr_check_psd(G_analysis, "true_G")
  }

  linear_by_trait <- sweep(X, 2L, weights, "*")
  linear_score <- rowSums(linear_by_trait)
  quadratic_score <- rowSums((X %*% W) * X)
  total_score <- linear_score + quadratic_score

  selection <- .dgr_qg_select(
    score = total_score,
    ids = prepared$init_data[[id_col]],
    n_select = n_select,
    selection_proportion = selection_proportion
  )

  ranked <- data.table::copy(prepared$init_data)
  ranked[, LinearPart := linear_score]
  ranked[, QuadraticPart := quadratic_score]
  ranked[, QGSI := total_score]
  ranked[, Rank := data.table::frank(-QGSI, ties.method = "average")]
  ranked[, Selected := selection$selected]
  data.table::setorderv(
    ranked,
    cols = c("QGSI", id_col),
    order = c(-1L, 1L)
  )

  linear_contributions <- NULL
  quadratic_contributions <- NULL
  if (isTRUE(return_contributions)) {
    linear_contributions <- data.table::as.data.table(linear_by_trait)
    data.table::setnames(
      linear_contributions,
      paste0("Linear_", trait_cols)
    )
    linear_contributions[, (id_col) := prepared$init_data[[id_col]]]
    data.table::setcolorder(linear_contributions, id_col)

    pairs <- which(upper.tri(W, diag = TRUE), arr.ind = TRUE)
    pair_values <- vapply(seq_len(nrow(pairs)), function(k) {
      i <- pairs[k, 1L]
      j <- pairs[k, 2L]
      multiplier <- if (i == j) 1 else 2
      multiplier * W[i, j] * X[, i] * X[, j]
    }, numeric(n_candidates))
    if (is.null(dim(pair_values))) {
      pair_values <- matrix(pair_values, ncol = 1L)
    }
    colnames(pair_values) <- vapply(seq_len(nrow(pairs)), function(k) {
      paste0(
        "Quadratic_", trait_cols[pairs[k, 1L]], "_x_",
        trait_cols[pairs[k, 2L]]
      )
    }, character(1))
    quadratic_contributions <- data.table::as.data.table(pair_values)
    quadratic_contributions[, (id_col) := prepared$init_data[[id_col]]]
    data.table::setcolorder(quadratic_contributions, id_col)
  }

  theory <- .dgr_qg_theory(
    weights = weights,
    W = W,
    Gamma = Gamma,
    true_G = G_analysis,
    selection_intensity = selection$normal_selection_intensity,
    trait_cols = trait_cols
  )

  observed_differential <- NULL
  if (selection$selected_n > 0L) {
    selected_X <- X[selection$selected, , drop = FALSE]
    differential <- colMeans(selected_X) - colMeans(X)
    observed_differential <- data.table::data.table(
      Trait = trait_cols,
      Mean_all = colMeans(X),
      Mean_selected = colMeans(selected_X),
      Observed_GEBV_differential = as.numeric(differential)
    )
  }

  result <- list(
    method = "QGSI",
    call = match.call(),
    trait_cols = trait_cols,
    linear_weights = weights,
    W = W,
    Gamma = Gamma,
    covariance_provenance = gamma_provenance,
    true_G = G_analysis,
    transformation = prepared$transformation,
    missing_data = prepared$missing_data,
    selection = selection[c(
      "requested", "selected_n", "selection_proportion",
      "normal_selection_intensity", "tie_break"
    )],
    ranked_geno = ranked,
    selected_geno = ranked[Selected == TRUE],
    component_summary = data.table::data.table(
      Component = c("LinearPart", "QuadraticPart", "QGSI"),
      Mean = c(mean(linear_score), mean(quadratic_score), mean(total_score)),
      SD = c(
        stats::sd(linear_score), stats::sd(quadratic_score),
        stats::sd(total_score)
      ),
      Min = c(min(linear_score), min(quadratic_score), min(total_score)),
      Max = c(max(linear_score), max(quadratic_score), max(total_score))
    ),
    theoretical_parameters = theory$parameters,
    expected_gain_per_trait = theory$expected_gain_per_trait,
    observed_selection_differential = observed_differential,
    linear_contributions = linear_contributions,
    quadratic_contributions = quadratic_contributions,
    weight_diagnostics = list(
      eigenvalues = stats::setNames(
        eigen(W, symmetric = TRUE, only.values = TRUE)$values,
        paste0("eigen_", seq_along(trait_cols))
      ),
      off_diagonal_coefficient_convention =
        "The coefficient of gamma_i * gamma_j is 2 * W[i, j].",
      interpretation = paste(
        "Negative diagonal curvature favours intermediate values;",
        "positive diagonal curvature favours extremes, conditional on",
        "the linear and cross-product terms."
      )
    ),
    contribution_scope = paste(
      "QGSI contributions are candidate-specific. They are not a single",
      "global additive marker-effect vector."
    )
  )
  class(result) <- c("quadratic_genomic_index", "list")
  .desiredgainr_dbg(
    debug,
    "Scored %d candidates with QGSI; Gamma source: %s.",
    n_candidates,
    gamma_provenance$source
  )
  result
}

.dgr_qg_traits <- function(trait_cols) {
  if (!is.character(trait_cols) || !length(trait_cols) ||
      anyNA(trait_cols) || any(!nzchar(trait_cols)) ||
      anyDuplicated(trait_cols)) {
    stop("trait_cols must contain unique, non-missing trait names.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.dgr_qg_scalar <- function(x, name, lower, lower_open) {
  valid <- is.numeric(x) && length(x) == 1L && is.finite(x)
  if (valid) {
    valid <- if (lower_open) x > lower else x >= lower
  }
  if (!valid) {
    relation <- if (lower_open) "greater than" else "at least"
    stop(name, " must be a finite scalar ", relation, " ", lower, ".",
         call. = FALSE)
  }
  invisible(TRUE)
}

.dgr_qg_symmetric_matrix <- function(
    x, trait_cols, name, tolerance
) {
  p <- length(trait_cols)
  if (!is.matrix(x) || any(dim(x) != p)) {
    stop(name, " must be a ", p, " x ", p, " matrix.", call. = FALSE)
  }
  if (xor(is.null(rownames(x)), is.null(colnames(x)))) {
    stop(name, " must have both row and column names, or neither.",
         call. = FALSE)
  }
  if (!is.null(rownames(x))) {
    absent_rows <- setdiff(trait_cols, rownames(x))
    absent_cols <- setdiff(trait_cols, colnames(x))
    if (length(absent_rows) || length(absent_cols)) {
      stop(name, " dimnames must contain every trait.", call. = FALSE)
    }
    x <- x[trait_cols, trait_cols, drop = FALSE]
  }
  storage.mode(x) <- "double"
  if (any(!is.finite(x))) {
    stop(name, " must contain finite values.", call. = FALSE)
  }
  asymmetry <- max(abs(x - t(x)))
  matrix_scale <- max(1, max(abs(x)))
  if (asymmetry > tolerance * matrix_scale) {
    stop(name, " must be symmetric within symmetry_tolerance.",
         call. = FALSE)
  }
  x <- (x + t(x)) / 2
  dimnames(x) <- list(trait_cols, trait_cols)
  x
}

.dgr_qg_estimate_gamma <- function(
    X,
    ids,
    relationship_matrix,
    trait_cols,
    tolerance,
    symmetry_tolerance
) {
  X_cov <- sweep(X, 2L, colMeans(X), "-")
  g <- nrow(X_cov)
  if (is.null(relationship_matrix)) {
    Gamma <- crossprod(X_cov) / g
    dimnames(Gamma) <- list(trait_cols, trait_cols)
    return(list(
      Gamma = Gamma,
      provenance = list(
        source = "estimated from reference GEBVs",
        equation = "Ceron-Rojas et al. (2026), Supplementary Equation 19.2",
        reference_n = g,
        divisor = g,
        covariance_centring = "reference GEBV column means removed",
        relationship_rank = NA_integer_,
        relationship_dimension = NA_integer_
      )
    ))
  }

  if (is.null(ids) || anyNA(ids) || anyDuplicated(ids)) {
    stop(
      paste(
        "reference_gebv_data must contain unique, non-missing id_col values",
        "when relationship_matrix is supplied."
      ),
      call. = FALSE
    )
  }
  K <- relationship_matrix
  if (!is.matrix(K) || any(dim(K) != g)) {
    stop(
      "relationship_matrix must be square with one row per reference genotype.",
      call. = FALSE
    )
  }
  if (is.null(rownames(K)) || is.null(colnames(K))) {
    stop(
      "relationship_matrix must have genotype row and column names.",
      call. = FALSE
    )
  }
  id_names <- as.character(ids)
  if (!all(id_names %in% rownames(K)) ||
      !all(id_names %in% colnames(K))) {
    stop(
      "relationship_matrix dimnames must contain every reference genotype.",
      call. = FALSE
    )
  }
  K <- K[id_names, id_names, drop = FALSE]
  K <- .dgr_qg_symmetric_matrix(
    K, id_names, "relationship_matrix", symmetry_tolerance
  )
  decomposition <- eigen(K, symmetric = TRUE)
  largest <- max(abs(decomposition$values))
  negative_limit <- tolerance * max(1, largest)
  if (min(decomposition$values) < -negative_limit) {
    stop("relationship_matrix must be positive semidefinite.",
         call. = FALSE)
  }
  retained <- decomposition$values > tolerance * max(1, largest)
  if (!any(retained)) {
    stop("relationship_matrix has zero numerical rank.", call. = FALSE)
  }
  vectors <- decomposition$vectors[, retained, drop = FALSE]
  projected <- crossprod(vectors, X_cov)
  inverse_projected <- sweep(
    projected, 1L, decomposition$values[retained], "/"
  )
  K_inverse_X <- vectors %*% inverse_projected
  Gamma <- crossprod(X_cov, K_inverse_X) / g
  Gamma <- (Gamma + t(Gamma)) / 2
  dimnames(Gamma) <- list(trait_cols, trait_cols)
  list(
    Gamma = Gamma,
    provenance = list(
      source = "relationship-adjusted estimate from reference GEBVs",
      equation = "Ceron-Rojas et al. (2026), Supplementary Equation 19.1",
      reference_n = g,
      divisor = g,
      covariance_centring = "reference GEBV column means removed",
      relationship_rank = sum(retained),
      relationship_dimension = g,
      inverse = if (all(retained)) {
        "spectral inverse"
      } else {
        "spectral Moore-Penrose inverse"
      },
      eigenvalue_tolerance = tolerance
    )
  )
}

.dgr_qg_select <- function(
    score, ids, n_select, selection_proportion
) {
  n <- length(score)
  if (!is.null(n_select) && !is.null(selection_proportion)) {
    stop(
      "Specify at most one of n_select and selection_proportion.",
      call. = FALSE
    )
  }
  requested <- !is.null(n_select) || !is.null(selection_proportion)
  if (!requested) {
    return(list(
      requested = FALSE,
      selected = rep(FALSE, n),
      selected_n = 0L,
      selection_proportion = NA_real_,
      normal_selection_intensity = NA_real_,
      tie_break = "No selection requested."
    ))
  }
  if (!is.null(n_select)) {
    if (!is.numeric(n_select) || length(n_select) != 1L ||
        !is.finite(n_select) || n_select != as.integer(n_select) ||
        n_select < 1L || n_select > n) {
      stop("n_select must be an integer between 1 and the candidate count.",
           call. = FALSE)
    }
    selected_n <- as.integer(n_select)
  } else {
    if (!is.numeric(selection_proportion) ||
        length(selection_proportion) != 1L ||
        !is.finite(selection_proportion) ||
        selection_proportion <= 0 || selection_proportion > 1) {
      stop("selection_proportion must be in (0, 1].", call. = FALSE)
    }
    selected_n <- min(n, max(1L, ceiling(selection_proportion * n)))
  }
  id_key <- as.character(ids)
  ordering <- order(-score, id_key, na.last = NA)
  selected <- rep(FALSE, n)
  selected[ordering[seq_len(selected_n)]] <- TRUE
  selected_fraction <- selected_n / n
  intensity <- if (selected_fraction >= 1) {
    0
  } else {
    cutoff <- stats::qnorm(1 - selected_fraction)
    stats::dnorm(cutoff) / selected_fraction
  }
  list(
    requested = TRUE,
    selected = selected,
    selected_n = selected_n,
    selection_proportion = selected_fraction,
    normal_selection_intensity = intensity,
    tie_break = paste(
      "Exact selection count; QGSI score descending, then",
      "candidate identifier ascending."
    )
  )
}

.dgr_qg_trace <- function(x) {
  sum(diag(x))
}

.dgr_qg_theory <- function(
    weights, W, Gamma, true_G, selection_intensity, trait_cols
) {
  linear_variance <- as.numeric(
    crossprod(weights, Gamma %*% weights)
  )
  quadratic_variance <- 2 * .dgr_qg_trace(
    W %*% Gamma %*% W %*% Gamma
  )
  numerical_scale <- max(1, abs(linear_variance), abs(quadratic_variance))
  if (linear_variance < -1e-10 * numerical_scale ||
      quadratic_variance < -1e-10 * numerical_scale) {
    stop(
      "The supplied weights and Gamma produced an invalid negative variance.",
      call. = FALSE
    )
  }
  linear_variance <- max(0, linear_variance)
  quadratic_variance <- max(0, quadratic_variance)
  index_variance <- linear_variance + quadratic_variance
  index_mean <- .dgr_qg_trace(W %*% Gamma)

  expected_response <- if (is.finite(selection_intensity)) {
    selection_intensity * sqrt(index_variance)
  } else {
    NA_real_
  }
  expected_gain <- rep(NA_real_, length(trait_cols))
  if (is.finite(selection_intensity) && linear_variance > 0) {
    expected_gain <- as.numeric(
      selection_intensity * Gamma %*% weights /
        sqrt(linear_variance)
    )
  }

  merit_variance <- NA_real_
  merit_index_covariance <- NA_real_
  accuracy_squared <- NA_real_
  mspe <- NA_real_
  accuracy_available <- !is.null(true_G)
  if (accuracy_available) {
    merit_variance <- as.numeric(
      crossprod(weights, true_G %*% weights)
    ) + 2 * .dgr_qg_trace(W %*% true_G %*% W %*% true_G)
    merit_index_covariance <- linear_variance +
      2 * .dgr_qg_trace(W %*% Gamma %*% W %*% true_G)
    mspe <- merit_variance + index_variance -
      2 * merit_index_covariance
    if (abs(mspe) <= 1e-10 * max(1, merit_variance, index_variance)) {
      mspe <- 0
    }
    if (merit_variance > 0) {
      accuracy_squared <- index_variance / merit_variance
    }
  }

  list(
    parameters = list(
      model_expected_index = index_mean,
      linear_index_variance = linear_variance,
      quadratic_index_variance = quadratic_variance,
      total_index_variance = index_variance,
      index_standard_deviation = sqrt(index_variance),
      selection_intensity = selection_intensity,
      expected_net_merit_response = expected_response,
      true_merit_variance = merit_variance,
      merit_index_covariance = merit_index_covariance,
      squared_index_merit_correlation = accuracy_squared,
      mean_squared_prediction_error = mspe,
      accuracy_and_mspe_available = accuracy_available,
      assumptions = c(
        "GEBVs represent the stated transformed trait space.",
        "Gamma is the genomic covariance of those GEBVs.",
        paste(
          "Normal-selection response uses the reported selection intensity;",
          "it is a model-based expectation."
        ),
        if (accuracy_available) {
          "Accuracy and MSPE use the supplied true_G."
        } else {
          paste(
            "Accuracy and MSPE are not reported because true_G was not",
            "supplied; Gamma is not relabelled as true genetic covariance."
          )
        }
      )
    ),
    expected_gain_per_trait = data.table::data.table(
      Trait = trait_cols,
      Expected_Genetic_Gain = expected_gain,
      Basis = paste(
        "Ceron-Rojas et al. (2026), Supplementary Equation 16;",
        "linear regression of breeding value on the index"
      )
    )
  )
}

#' Print a quadratic genomic selection index
#'
#' @param x A `quadratic_genomic_index` object.
#' @param ... Unused.
#'
#' @return `x`, invisibly.
#' @keywords internal
#' @export
print.quadratic_genomic_index <- function(x, ...) {
  cat("<quadratic_genomic_index>\n")
  cat("  Candidates:", nrow(x$ranked_geno), "\n")
  cat("  Traits:", length(x$trait_cols), "\n")
  cat("  Gamma:", x$covariance_provenance$source, "\n")
  if (isTRUE(x$selection$requested)) {
    cat(
      "  Selected:", x$selection$selected_n,
      sprintf("(%.1f%%)\n", 100 * x$selection$selection_proportion)
    )
  } else {
    cat("  Selection: ranking only\n")
  }
  cat(
    "  Model index SD:",
    format(x$theoretical_parameters$index_standard_deviation, digits = 5),
    "\n"
  )
  invisible(x)
}
