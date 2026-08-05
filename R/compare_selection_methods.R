# Comparative evaluation of multi-trait selection methods.

#' Define one objective for a cross-family comparison
#'
#' A comparison objective fixes the yardstick before fitted methods are
#' inspected. Desired gains define a response target. `aggregate_weights` and
#' `W` define one common utility
#' \deqn{U(g)=a^\mathsf{T}g+g^\mathsf{T}Wg.}
#' Every method is evaluated against this same objective.
#'
#' @param desired_gains Optional named favourable desired-gain vector.
#' @param aggregate_weights Optional named linear aggregate-weight vector.
#' @param W Optional symmetric matrix of squared and cross-product utility
#'   weights. Its names must match the objective traits.
#' @param G Optional genetic covariance matrix in original trait coordinates.
#'   It is required when `gain_units = "genetic_sd"`. It also supplies one
#'   common geometry for target alignment and common-merit evaluation.
#' @param gain_units Coordinate system for the objective. Use either original
#'   trait units or genetic standard deviations. `desired_gains` are expressed
#'   in these coordinates. `aggregate_weights` and `W` must act on gains in
#'   the same coordinates.
#'
#' @return An object of class `desiredgainr_comparison_objective`.
#'
#' @examples
#' traits <- c("yield", "disease")
#' G <- matrix(c(1, -0.2, -0.2, 0.8), 2,
#'   dimnames = list(traits, traits)
#' )
#' comparison_objective(
#'   desired_gains = c(yield = 1, disease = 0.5),
#'   G = G,
#'   gain_units = "genetic_sd"
#' )
#'
#' @seealso [compare_selection_methods()]
#' @export
comparison_objective <- function(
  desired_gains = NULL,
  aggregate_weights = NULL,
  W = NULL,
  G = NULL,
  gain_units = c("trait", "genetic_sd")
) {
  gain_units <- match.arg(gain_units)
  supplied <- list(desired_gains, aggregate_weights)
  named <- supplied[!vapply(supplied, is.null, logical(1L))]
  trait_sets <- lapply(named, function(x) {
    if (!is.numeric(x) || !length(x) || is.null(names(x)) ||
      any(!is.finite(x)) || any(!nzchar(names(x))) || anyDuplicated(names(x))) {
      stop(
        "Objective vectors must be finite numeric vectors with unique names.",
        call. = FALSE
      )
    }
    names(x)
  })
  if (!is.null(W)) {
    if (!is.matrix(W) || !is.numeric(W) || nrow(W) != ncol(W) ||
      is.null(rownames(W)) || is.null(colnames(W)) ||
      !identical(rownames(W), colnames(W)) || any(!is.finite(W))) {
      stop("W must be a finite square matrix with matching names.",
        call. = FALSE
      )
    }
    if (!isTRUE(all.equal(W, t(W), tolerance = sqrt(.Machine$double.eps)))) {
      stop("W must be symmetric.", call. = FALSE)
    }
    trait_sets[[length(trait_sets) + 1L]] <- rownames(W)
  }
  if (!is.null(G)) {
    if (!is.matrix(G) || !is.numeric(G) || nrow(G) != ncol(G) ||
      is.null(rownames(G)) || is.null(colnames(G)) ||
      !identical(rownames(G), colnames(G)) || any(!is.finite(G))) {
      stop("G must be a finite square covariance matrix with matching names.",
        call. = FALSE
      )
    }
    trait_sets[[length(trait_sets) + 1L]] <- rownames(G)
  }
  if (!length(trait_sets)) {
    stop("Supply desired_gains, aggregate_weights, W, or G.", call. = FALSE)
  }
  traits <- trait_sets[[1L]]
  if (!all(vapply(trait_sets, setequal, logical(1L), y = traits))) {
    stop("Every comparison-objective component must name the same traits.",
      call. = FALSE
    )
  }
  reorder_vector <- function(x) if (is.null(x)) NULL else x[traits]
  desired_gains <- reorder_vector(desired_gains)
  aggregate_weights <- reorder_vector(aggregate_weights)
  if (!is.null(desired_gains) &&
    (any(desired_gains < 0) || !any(desired_gains > 0))) {
    stop(
      "desired_gains must contain favourable non-negative gains and at least one positive value.",
      call. = FALSE
    )
  }
  if (!is.null(W)) W <- W[traits, traits, drop = FALSE]
  if (!is.null(G)) {
    G <- G[traits, traits, drop = FALSE]
    .dgr_check_psd(G, "comparison objective G")
  }
  if (identical(gain_units, "genetic_sd") && is.null(G)) {
    stop("gain_units = 'genetic_sd' requires G.", call. = FALSE)
  }
  result <- list(
    traits = traits,
    desired_gains = desired_gains,
    aggregate_weights = aggregate_weights,
    W = W,
    G = G,
    gain_units = gain_units
  )
  class(result) <- c("desiredgainr_comparison_objective", "list")
  result
}

.dgr_comparison_ids <- function(model, table, preferred = NULL) {
  if (!is.null(model$candidate_id)) {
    return(as.character(model$candidate_id))
  }
  if (!is.null(preferred) && preferred %in% names(table)) {
    return(as.character(table[[preferred]]))
  }
  candidates <- setdiff(
    names(table),
    c(
      "SelectionIndex", "QGSI", "LinearPart", "QuadraticPart", "Rank",
      "Selected", "Eligible"
    )
  )
  if (!length(candidates)) {
    stop("The fitted object does not expose candidate identifiers.",
      call. = FALSE
    )
  }
  as.character(table[[candidates[[1L]]]])
}

.dgr_comparison_adapter <- function(model) {
  if (inherits(model, "desiredgainr_index")) {
    traits <- model$trait_cols
    scale <- model$transformation$scale[traits]
    direction <- model$transformation$direction[traits]
    ids <- as.character(model$candidate_id)
    ranks <- stats::setNames(model$ranking$rank, model$ranking$id)[ids]
    scores <- stats::setNames(model$ranking$score, model$ranking$id)[ids]
    selected <- as.character(model$selected$id)
    expected <- if (is.null(model$evaluation)) {
      NULL
    } else {
      model$evaluation$expected_response[traits] * scale
    }
    observed <- if (is.null(model$observed_differential)) {
      NULL
    } else {
      stats::setNames(
        model$observed_differential$Differential[
          match(traits, model$observed_differential$Trait)
        ],
        traits
      ) * scale
    }
    G <- if (is.null(model$G)) {
      NULL
    } else {
      S <- diag(scale, nrow = length(traits))
      dimnames(S) <- list(traits, traits)
      S %*% model$G[traits, traits, drop = FALSE] %*% S
    }
    native <- model$evaluation
    return(list(
      family = model$method, strategy = model$strategy,
      traits = traits, direction = direction, ids = ids,
      score = as.numeric(scores), rank = as.numeric(ranks),
      selected = selected, expected = expected, observed = observed,
      G = G, G_source = "fitted genetic covariance",
      scale_to_trait = scale,
      selection_intensity = model$selection_intensity,
      expected_basis = if (is.null(expected)) {
        "unavailable"
      } else {
        "exact linear-index response"
      },
      native_R_HI = if (is.null(native)) NA_real_ else native$R_HI,
      native_Delta_H = if (is.null(native)) NA_real_ else native$delta_H,
      native_RE = if (is.null(native)) NA_real_ else native$RE,
      native_h2 = if (is.null(native)) NA_real_ else native$h2_index,
      native_accuracy = if (is.null(native)) NA_real_ else native$accuracy_index,
      native_MSPE = NA_real_
    ))
  }

  if (inherits(model, "desiredgainr_generalized_index")) {
    traits <- model$model$objective_names
    direction <- stats::setNames(rep(1, length(traits)), traits)
    direction[intersect(model$lower_is_better, traits)] <- -1
    ids <- as.character(model$candidate_id)
    ranks <- stats::setNames(model$ranking$rank, model$ranking$id)[ids]
    scores <- stats::setNames(model$ranking$score, model$ranking$id)[ids]
    D <- diag(direction, nrow = length(traits))
    dimnames(D) <- list(traits, traits)
    return(list(
      family = paste0("generalized_", model$method), strategy = "index",
      traits = traits, direction = direction, ids = ids,
      score = as.numeric(scores), rank = as.numeric(ranks),
      selected = as.character(model$selected$id),
      expected = model$expected_response[traits] * direction,
      observed = NULL,
      G = D %*% model$model$G[traits, traits, drop = FALSE] %*% D,
      G_source = "general information-model genetic covariance",
      scale_to_trait = stats::setNames(rep(1, length(traits)), traits),
      selection_intensity = model$selection_intensity,
      expected_basis = "exact general linear-index response",
      native_R_HI = model$accuracy,
      native_Delta_H = model$delta_merit,
      native_RE = NA_real_, native_h2 = NA_real_,
      native_accuracy = model$accuracy, native_MSPE = NA_real_
    ))
  }

  if (inherits(model, "desiredgainr_dgsi") ||
    inherits(model, "desired_gain_index")) {
    traits <- model$trait_cols
    scale <- model$transformation$scale[traits]
    direction <- model$transformation$direction[traits]
    table <- model$ranked_geno
    preferred <- if (!is.null(model$id_col)) model$id_col else "GenoID"
    row_ids <- if (preferred %in% names(table)) {
      as.character(table[[preferred]])
    } else {
      .dgr_comparison_ids(structure(list(), class = "list"), table)
    }
    ids <- if (is.null(model$candidate_id)) {
      row_ids
    } else {
      as.character(model$candidate_id)
    }
    rank_by_id <- stats::setNames(table$Rank, row_ids)
    score_by_id <- stats::setNames(table$SelectionIndex, row_ids)
    selected <- row_ids[table$Selected]
    expected <- if (is.null(model$theoretical_response)) {
      NULL
    } else {
      model$theoretical_response$analysis_units[traits] * scale
    }
    gain_scale <- model$candidate_sd_analysis
    observed <- if (is.null(gain_scale)) {
      NULL
    } else {
      model$realised_response[traits] * gain_scale[traits] * scale
    }
    S <- diag(scale, nrow = length(traits))
    dimnames(S) <- list(traits, traits)
    return(list(
      family = "dgsi_iterative_search", strategy = "index",
      traits = traits, direction = direction, ids = ids,
      score = as.numeric(score_by_id[ids]), rank = as.numeric(rank_by_id[ids]),
      selected = selected, expected = expected, observed = observed,
      G = S %*% model$G[traits, traits, drop = FALSE] %*% S,
      G_source = "DGSI fitted genetic covariance",
      scale_to_trait = scale,
      selection_intensity = if (is.null(model$theoretical_response)) {
        NA_real_
      } else {
        model$theoretical_response$selection_intensity
      },
      expected_basis = "exact response of the final linear index",
      native_R_HI = NA_real_, native_Delta_H = NA_real_, native_RE = NA_real_,
      native_h2 = NA_real_, native_accuracy = NA_real_, native_MSPE = NA_real_
    ))
  }

  if (inherits(model, "quadratic_genomic_index")) {
    traits <- model$trait_cols
    scale <- model$transformation$scale[traits]
    direction <- model$transformation$direction[traits]
    table <- model$ranked_geno
    preferred <- if (!is.null(model$id_col)) model$id_col else "GenoID"
    row_ids <- if (preferred %in% names(table)) {
      as.character(table[[preferred]])
    } else {
      .dgr_comparison_ids(structure(list(), class = "list"), table)
    }
    ids <- if (is.null(model$candidate_id)) {
      row_ids
    } else {
      as.character(model$candidate_id)
    }
    rank_by_id <- stats::setNames(table$Rank, row_ids)
    score_by_id <- stats::setNames(table$QGSI, row_ids)
    selected <- row_ids[table$Selected]
    gains <- model$expected_gain_per_trait
    expected <- stats::setNames(
      gains$Expected_Genetic_Gain[match(traits, gains$Trait)], traits
    ) * scale
    observed <- if (is.null(model$observed_selection_differential)) {
      NULL
    } else {
      stats::setNames(
        model$observed_selection_differential$Observed_GEBV_differential[
          match(traits, model$observed_selection_differential$Trait)
        ],
        traits
      ) * scale
    }
    G_analysis <- if (is.null(model$true_G)) model$Gamma else model$true_G
    S <- diag(scale, nrow = length(traits))
    dimnames(S) <- list(traits, traits)
    theory <- model$theoretical_parameters
    squared_accuracy <- theory$squared_index_merit_correlation
    q_accuracy <- if (length(squared_accuracy) == 1L &&
      is.finite(squared_accuracy)) {
      sqrt(squared_accuracy)
    } else {
      NA_real_
    }
    has_curvature <- any(abs(model$W) > sqrt(.Machine$double.eps))
    return(list(
      family = if (has_curvature) "qgsi" else "linear_genomic_index",
      strategy = if (has_curvature) "quadratic_index" else "index",
      traits = traits, direction = direction, ids = ids,
      score = as.numeric(score_by_id[ids]), rank = as.numeric(rank_by_id[ids]),
      selected = selected, expected = expected, observed = observed,
      G = S %*% G_analysis[traits, traits, drop = FALSE] %*% S,
      G_source = if (is.null(model$true_G)) {
        "genomic prediction covariance Gamma"
      } else {
        "supplied true_G"
      },
      scale_to_trait = scale,
      selection_intensity = model$selection$normal_selection_intensity,
      expected_basis = if (has_curvature) {
        "linear-regression approximation; the quadratic score is non-normal"
      } else {
        "exact linear genomic-index response (W = 0)"
      },
      native_R_HI = NA_real_,
      native_Delta_H = theory$expected_net_merit_response,
      native_RE = NA_real_, native_h2 = NA_real_,
      native_accuracy = q_accuracy,
      native_MSPE = theory$mean_squared_prediction_error
    ))
  }

  stop(
    paste(
      "Every model must come from selection_index(), restricted_index(),",
      "generalized_index(), run_dgsi(), or run_qgsi()."
    ),
    call. = FALSE
  )
}

.dgr_comparison_validation <- function(data, traits, ids, direction, unit_sd) {
  frame <- as.data.frame(data)
  absent <- setdiff(traits, names(frame))
  if (length(absent)) {
    stop("validation_data is missing objective traits: ",
      paste(absent, collapse = ", "),
      call. = FALSE
    )
  }
  id_candidates <- setdiff(names(frame), traits)
  matched <- id_candidates[vapply(id_candidates, function(name) {
    values <- as.character(frame[[name]])
    length(values) == length(ids) && !anyNA(values) && !anyDuplicated(values) &&
      setequal(values, ids)
  }, logical(1L))]
  if (length(matched) == 1L) {
    validation_ids <- as.character(frame[[matched]])
  } else if (!is.null(rownames(frame)) && setequal(rownames(frame), ids)) {
    validation_ids <- rownames(frame)
  } else {
    stop(
      "validation_data needs one identifier column, or row names, matching every candidate.",
      call. = FALSE
    )
  }
  Z <- as.matrix(frame[match(ids, validation_ids), traits, drop = FALSE])
  storage.mode(Z) <- "double"
  if (any(!is.finite(Z))) {
    stop("validation_data objective values must be finite.", call. = FALSE)
  }
  Z <- sweep(Z, 2L, direction[traits], "*")
  sweep(Z, 2L, unit_sd[traits], "/")
}

#' Compare selection methods on one common decision problem
#'
#' This function accepts results from [selection_index()], [restricted_index()],
#' [generalized_index()], [run_dgsi()], and [run_qgsi()]. Internal adapters put
#' response vectors in favourable original trait units. A
#' [comparison_objective()] can then place them in trait or genetic-standard-
#' deviation units and evaluate every method against one fixed target and
#' utility.
#'
#' @details
#' Coefficients are never compared across families. Candidate ranks, selected
#' sets, per-trait responses, target alignment, and common utility are compared.
#' QGSI expected gains retain their documented linear-regression approximation
#' when `W` contains curvature. The `W = 0` special case has the exact linear
#' response. A quadratic common utility is evaluated exactly only from
#' `validation_data`. A singular common `G` still permits response and decision
#' comparisons. Mahalanobis alignment and Satoh projection remain missing
#' because those criteria require an inverse.
#'
#' `target_gains` is retained for backward compatibility. It uses the analysis
#' scale of the first fitted `desiredgainr_index` and requires every model to
#' have the same transformation. New cross-family analyses should use
#' `objective`.
#'
#' @param models Named list containing at least two supported fitted objects.
#' @param target_gains Optional legacy target on the common fitted analysis
#'   scale.
#' @param objective Optional object returned by [comparison_objective()].
#' @param validation_data Optional data frame containing one common set of
#'   objective-trait values for every candidate. One non-trait identifier
#'   column, or the row names, must match the fitted candidate identifiers.
#'
#' @return An object of class `desiredgainr_method_comparison`. It contains
#'   method summaries, model-based responses, optional validation responses,
#'   fairness checks, rank correlations, and selected-set agreement.
#'
#' @seealso [comparison_objective()], [compare_dg_and_qgsi()],
#'   [evaluate_index()], [candidate_score_uncertainty()]
#' @export
compare_selection_methods <- function(
  models,
  target_gains = NULL,
  objective = NULL,
  validation_data = NULL
) {
  if (!is.list(models) || length(models) < 2L) {
    stop("models must contain at least two fitted indices.", call. = FALSE)
  }
  if (is.null(names(models)) || any(names(models) == "") ||
    anyDuplicated(names(models))) {
    stop("models must have unique, non-empty names.", call. = FALSE)
  }
  if (!is.null(target_gains) && !is.null(objective)) {
    stop("Supply either target_gains or objective, not both.", call. = FALSE)
  }
  if (!is.null(objective) &&
    !inherits(objective, "desiredgainr_comparison_objective")) {
    stop("objective must come from comparison_objective().", call. = FALSE)
  }

  adapted <- lapply(models, .dgr_comparison_adapter)
  traits <- adapted[[1L]]$traits
  if (!all(vapply(adapted, function(x) setequal(x$traits, traits), logical(1L)))) {
    stop("Every model must use the same objective traits.", call. = FALSE)
  }
  adapted <- lapply(adapted, function(x) {
    order <- match(traits, x$traits)
    x$traits <- traits
    x$direction <- x$direction[traits]
    x$expected <- if (is.null(x[["expected"]])) {
      NULL
    } else {
      x[["expected"]][traits]
    }
    x$observed <- if (is.null(x[["observed"]])) {
      NULL
    } else {
      x[["observed"]][traits]
    }
    fitted_G <- x[["G"]]
    x["G"] <- list(if (is.null(fitted_G)) {
      NULL
    } else {
      fitted_G[traits, traits, drop = FALSE]
    })
    x$scale_to_trait <- x$scale_to_trait[traits]
    x
  })
  direction <- adapted[[1L]]$direction
  if (!all(vapply(
    adapted, function(x) identical(x$direction, direction),
    logical(1L)
  ))) {
    stop("Every model must use the same favourable trait directions.",
      call. = FALSE
    )
  }
  ids <- adapted[[1L]]$ids
  if (!all(vapply(adapted, function(x) setequal(x$ids, ids), logical(1L)))) {
    stop("Every model must use the same candidates.", call. = FALSE)
  }
  adapted <- lapply(adapted, function(x) {
    order <- match(ids, x$ids)
    x$ids <- ids
    x$score <- x$score[order]
    x$rank <- x$rank[order]
    x
  })

  legacy <- !is.null(target_gains)
  if (legacy) {
    if (!all(vapply(models, inherits, logical(1L), "desiredgainr_index"))) {
      stop("Legacy target_gains supports selection_index() and",
      "restricted_index() only. Use objective for cross-family comparisons.",
        call. = FALSE
      )
    }
    reference_scale <- adapted[[1L]]$scale_to_trait
    same_scale <- all(vapply(adapted, function(x) {
      isTRUE(all.equal(
        x$scale_to_trait, reference_scale,
        check.attributes = FALSE
      ))
    }, logical(1L)))
    if (!same_scale) {
      stop("Every model must use the same direction,",
      "centring, and scaling for legacy target_gains.",
        call. = FALSE
      )
    }
    target_gains <- .dgr_named_vector(target_gains, traits, "target_gains")
    if (any(target_gains < 0) || !any(target_gains > 0)) {
      stop(
        "target_gains must contain favourable",
        "non-negative gains and at least one positive value.",
        call. = FALSE
      )
    }
    objective <- structure(list(
      traits = traits,
      desired_gains = target_gains * reference_scale,
      aggregate_weights = NULL, W = NULL, G = NULL, gain_units = "trait"
    ), class = c("desiredgainr_comparison_objective", "list"))
  }

  if (!is.null(objective) && !setequal(objective$traits, traits)) {
    stop("objective and fitted models must name the same traits.",
      call. = FALSE
    )
  }
  if (!is.null(objective)) {
    objective$desired_gains <- objective$desired_gains[traits]
    objective$aggregate_weights <- objective$aggregate_weights[traits]
    if (!is.null(objective$W)) objective$W <- objective$W[traits, traits]
    if (!is.null(objective$G)) objective$G <- objective$G[traits, traits]
  }

  unit_sd <- stats::setNames(rep(1, length(traits)), traits)
  common_G <- NULL
  common_G_source <- "unavailable"
  if (!is.null(objective) && !is.null(objective$G)) {
    D <- diag(direction, nrow = length(traits))
    dimnames(D) <- list(traits, traits)
    common_G <- D %*% objective$G %*% D
    common_G_source <- "comparison_objective()"
  } else {
    available <- which(vapply(
      adapted, function(x) !is.null(x[["G"]]), logical(1L)
    ))
    if (length(available)) {
      common_G <- adapted[[available[[1L]]]][["G"]]
      common_G_source <- adapted[[available[[1L]]]]$G_source
    }
  }
  gain_units <- if (is.null(objective)) "trait" else objective$gain_units
  if (identical(gain_units, "genetic_sd")) {
    unit_sd <- sqrt(diag(common_G))
    if (any(!is.finite(unit_sd) | unit_sd <= 0)) {
      stop("The common G has invalid genetic standard deviations.",
        call. = FALSE
      )
    }
  }
  if (!is.null(common_G)) {
    common_G <- sweep(sweep(common_G, 1L, unit_sd, "/"), 2L, unit_sd, "/")
  }
  common_G_invertible <- if (is.null(common_G)) {
    FALSE
  } else {
    eigenvalues <- eigen(common_G, symmetric = TRUE, only.values = TRUE)$values
    min(eigenvalues) > sqrt(.Machine$double.eps) * max(1, max(eigenvalues))
  }

  target <- if (is.null(objective)) NULL else objective$desired_gains
  weights <- if (is.null(objective)) NULL else objective$aggregate_weights
  W <- if (is.null(objective)) NULL else objective$W
  if (is.null(W) && !is.null(objective)) {
    W <- matrix(0, length(traits), length(traits),
      dimnames = list(traits, traits)
    )
  }

  finite_min <- function(x) {
    x <- x[is.finite(x)]
    if (length(x)) min(x) else NA_real_
  }
  finite_mean <- function(x) {
    x <- x[is.finite(x)]
    if (length(x)) mean(x) else NA_real_
  }
  safe_spearman <- function(x, y) {
    keep <- is.finite(x) & is.finite(y)
    if (sum(keep) < 2L || stats::sd(x[keep]) == 0 ||
      stats::sd(y[keep]) == 0) {
      return(NA_real_)
    }
    stats::cor(x[keep], y[keep], method = "spearman")
  }
  alignment <- function(response) {
    if (is.null(target) || !common_G_invertible || any(!is.finite(response))) {
      return(c(beta = NA_real_, alignment = NA_real_, residual = NA_real_))
    }
    value <- evaluate_restricted_response(response, target, common_G)
    c(
      beta = value$beta, alignment = value$mahalanobis_alignment,
      residual = value$mahalanobis_residual
    )
  }

  response_rows <- lapply(names(adapted), function(label) {
    x <- adapted[[label]]
    expected <- if (is.null(x[["expected"]])) {
      stats::setNames(rep(NA_real_, length(traits)), traits)
    } else {
      x[["expected"]] / unit_sd
    }
    observed <- if (is.null(x[["observed"]])) {
      stats::setNames(rep(NA_real_, length(traits)), traits)
    } else {
      x[["observed"]] / unit_sd
    }
    data.frame(
      Method = label, Trait = traits,
      Expected_response = as.numeric(expected),
      Observed_differential = as.numeric(observed),
      Target = if (is.null(target)) NA_real_ else as.numeric(target),
      Expected_attainment = if (is.null(target)) {
        NA_real_
      } else {
        ifelse(target > 0, expected / target, NA_real_)
      },
      Observed_attainment = if (is.null(target)) {
        NA_real_
      } else {
        ifelse(target > 0, observed / target, NA_real_)
      },
      Expected_response_basis = x$expected_basis,
      Units = gain_units,
      stringsAsFactors = FALSE
    )
  })
  responses <- data.table::rbindlist(response_rows)

  validation_responses <- NULL
  validation_utility <- NULL
  validation_matrix <- NULL
  if (!is.null(validation_data)) {
    validation_matrix <- .dgr_comparison_validation(
      validation_data, traits, ids, direction, unit_sd
    )
    has_utility <- !is.null(weights) || (!is.null(W) && any(W != 0))
    utility <- if (has_utility) rep(0, nrow(validation_matrix)) else NULL
    if (!is.null(weights)) {
      utility <- utility + as.numeric(validation_matrix %*% weights)
    }
    if (!is.null(W) && any(W != 0)) {
      utility <- utility + rowSums((validation_matrix %*% W) * validation_matrix)
    }
    validation_responses <- data.table::rbindlist(lapply(names(adapted), function(label) {
      selected <- ids %in% adapted[[label]]$selected
      delta <- if (any(selected)) {
        colMeans(validation_matrix[selected, , drop = FALSE]) -
          colMeans(validation_matrix)
      } else {
        rep(NA_real_, length(traits))
      }
      data.frame(
        Method = label, Trait = traits,
        Validation_response = as.numeric(delta), Units = gain_units
      )
    }))
    if (has_utility) {
      validation_utility <- data.table::rbindlist(lapply(names(adapted), function(label) {
        selected <- ids %in% adapted[[label]]$selected
        data.frame(
          Method = label,
          Validation_utility_response = if (any(selected)) {
            mean(utility[selected]) - mean(utility)
          } else {
            NA_real_
          },
          Validation_utility_rank_correlation = safe_spearman(
            adapted[[label]]$score, utility
          )
        )
      }))
    }
  }

  summary_rows <- lapply(names(adapted), function(label) {
    x <- adapted[[label]]
    method_responses <- responses[responses$Method == label]
    expected <- method_responses$Expected_response
    diagnostics <- alignment(stats::setNames(expected, traits))
    # This is the response in the linear component a'g. It remains meaningful
    # when the declared utility also contains g'Wg. The quadratic component
    # depends on second moments and cannot be recovered from mean trait gains.
    common_merit_response <- if (!is.null(weights) &&
      all(is.finite(expected))) {
      sum(weights * expected)
    } else {
      NA_real_
    }
    common_merit_correlation <- NA_real_
    if (is.finite(common_merit_response) && !is.null(common_G) &&
      grepl("^exact", x$expected_basis) &&
      is.finite(x$selection_intensity) && x$selection_intensity > 0) {
      merit_sd <- sqrt(as.numeric(crossprod(weights, common_G %*% weights)))
      if (merit_sd > 0) {
        common_merit_correlation <- .dgr_clamp_correlation(
          common_merit_response / (x$selection_intensity * merit_sd),
          "The common-merit correlation"
        )
      }
    }
    validation_row <- if (is.null(validation_utility)) {
      NULL
    } else {
      validation_utility[validation_utility$Method == label]
    }
    data.frame(
      Method = label, Family = x$family, Strategy = x$strategy,
      N_selected = length(x$selected),
      Selected_fraction = length(x$selected) / length(ids),
      Selection_intensity = x$selection_intensity,
      R_HI = x$native_R_HI, Delta_H = x$native_Delta_H,
      RE = x$native_RE, Index_heritability = x$native_h2,
      Index_accuracy = x$native_accuracy, MSPE = x$native_MSPE,
      Common_merit_response = common_merit_response,
      Common_merit_correlation = common_merit_correlation,
      Worst_expected_attainment = finite_min(method_responses$Expected_attainment),
      Mean_expected_attainment = finite_mean(method_responses$Expected_attainment),
      Worst_observed_attainment = finite_min(method_responses$Observed_attainment),
      Satoh_beta = diagnostics[["beta"]],
      Mahalanobis_alignment = diagnostics[["alignment"]],
      Mahalanobis_residual = diagnostics[["residual"]],
      Validation_utility_response = if (is.null(validation_row)) {
        NA_real_
      } else {
        validation_row$Validation_utility_response
      },
      Validation_utility_rank_correlation = if (is.null(validation_row)) {
        NA_real_
      } else {
        validation_row$Validation_utility_rank_correlation
      },
      Expected_response_basis = x$expected_basis,
      stringsAsFactors = FALSE
    )
  })
  summary <- data.table::rbindlist(summary_rows)

  rank_matrix <- vapply(adapted, `[[`, numeric(length(ids)), "rank")
  colnames(rank_matrix) <- names(adapted)
  rank_correlation <- outer(
    names(adapted), names(adapted),
    Vectorize(function(a, b) {
      safe_spearman(
        rank_matrix[, a], rank_matrix[, b]
      )
    })
  )
  dimnames(rank_correlation) <- list(names(adapted), names(adapted))
  selected_sets <- lapply(adapted, `[[`, "selected")
  selected_overlap <- outer(
    selected_sets, selected_sets,
    Vectorize(function(a, b) length(intersect(a, b)))
  )
  selected_jaccard <- outer(
    selected_sets, selected_sets,
    Vectorize(function(a, b) {
      total <- length(union(a, b))
      if (total) length(intersect(a, b)) / total else NA_real_
    })
  )
  dimnames(selected_overlap) <- dimnames(selected_jaccard) <-
    list(names(adapted), names(adapted))

  G_agreement <- if (is.null(common_G)) {
    FALSE
  } else {
    all(vapply(adapted, function(x) {
      fitted_G <- x[["G"]]
      if (is.null(fitted_G)) {
        return(FALSE)
      }
      candidate_G <- sweep(sweep(fitted_G, 1L, unit_sd, "/"), 2L, unit_sd, "/")
      isTRUE(all.equal(candidate_G, common_G,
        tolerance = 1e-7,
        check.attributes = FALSE
      ))
    }, logical(1L)))
  }
  counts <- vapply(adapted, function(x) length(x$selected), integer(1L))
  intensities <- vapply(adapted, `[[`, numeric(1L), "selection_intensity")
  fairness <- data.frame(
    Condition = c(
      "Same candidates and objective traits",
      "Same favourable trait directions",
      "Common response units",
      "Same number selected",
      "Same selection intensity",
      "One common genetic covariance for response geometry",
      "Common genetic covariance is invertible",
      "Validation identifiers aligned when supplied"
    ),
    Satisfied = c(
      TRUE, TRUE, TRUE,
      length(unique(counts)) == 1L,
      all(is.finite(intensities)) &&
        length(unique(signif(intensities, 12L))) == 1L,
      G_agreement,
      common_G_invertible,
      TRUE
    ),
    Interpretation = c(
      "Required and enforced by this function.",
      "Required and enforced by this function.",
      paste("Responses are reported in", gain_units, "units."),
      "Hard culling can retain fewer candidates than requested.",
      "Model-based responses require one common selection pressure.",
      paste("Alignment uses", common_G_source, "; disagreement is reported."),
      "Mahalanobis alignment and Satoh projection require an invertible covariance matrix.",
      if (is.null(validation_matrix)) {
        "No validation data were supplied. External validation evidence remains unavailable."
      } else {
        "Validation values provide the same external yardstick for every method."
      }
    )
  )

  result <- list(
    summary = summary, responses = responses,
    validation_responses = validation_responses,
    validation_utility = validation_utility,
    rank_correlation = rank_correlation,
    selected_jaccard = selected_jaccard,
    selected_overlap = selected_overlap,
    fairness = fairness,
    objective = objective,
    target_gains = target,
    response_units = gain_units,
    common_G = common_G,
    common_G_source = common_G_source,
    models = names(models),
    interpretation = paste(
      "Read per-trait response and target attainment first.",
      "Common-merit quantities use one fixed linear utility.",
      "Quadratic utility is compared on common validation values.",
      "QGSI model-based gains retain their stated approximation."
    )
  )
  class(result) <- c("desiredgainr_method_comparison", "list")
  result
}

#' @export
print.desiredgainr_method_comparison <- function(x, ...) {
  cat("<desiredgainr_method_comparison>\n")
  cat("  Methods:", paste(x$models, collapse = ", "), "\n")
  cat("  Response units:", x$response_units, "\n")
  failed <- x$fairness$Condition[!x$fairness$Satisfied]
  if (length(failed)) {
    cat("  Comparison conditions requiring attention:\n")
    cat(paste0("    * ", failed, collapse = "\n"), "\n")
  } else {
    cat("  All recorded comparison conditions are satisfied.\n")
  }
  columns <- c(
    "Method", "Family", "N_selected", "Common_merit_response",
    "Worst_expected_attainment", "Mahalanobis_alignment",
    "Validation_utility_response"
  )
  print(as.data.frame(x$summary)[, columns, drop = FALSE], row.names = FALSE)
  invisible(x)
}
