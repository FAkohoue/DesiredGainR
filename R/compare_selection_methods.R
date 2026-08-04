# Comparative evaluation of multi-trait selection methods.

#' Compare multi-trait selection methods on one decision problem
#'
#' A fair comparison uses the same candidates, trait orientation, measurement
#' scale, selection intensity, and definition of aggregate merit. This function
#' checks those conditions. It then places biological response, objective
#' attainment, ranking agreement, and selected-set agreement in one result.
#'
#' @details
#' The function compares fitted objects returned by [selection_index()]. It
#' keeps expected response separate from the observed differential among the
#' supplied candidates. Expected response follows selection-index theory.
#' Observed differential describes the present candidate set.
#'
#' When `target_gains` is supplied, it must use the analysis scale stored in the
#' fitted objects. For example, a model fitted to genetic-standard-deviation
#' values requires targets in genetic standard deviations. Positive values
#' mean improvement because trait direction has already been applied by
#' [selection_index()].
#'
#' Methods based on ranks or thresholds lack a closed-form expected response in
#' this implementation. Their expected-response fields remain missing. Their
#' observed differentials, ranks, and selected sets remain comparable.
#'
#' @param models Named list containing at least two fitted
#'   `desiredgainr_index` objects.
#' @param target_gains Optional named vector of favourable desired gains on the
#'   common fitted analysis scale.
#'
#' @return An object of class `desiredgainr_method_comparison`. The `summary`
#'   component contains one row per method. The `responses` component contains
#'   one row per method and trait. The three matrices are `rank_correlation`,
#'   `selected_jaccard`, and `selected_overlap`. The `fairness` component
#'   records the conditions required for interpretation.
#'
#' @references
#' Beavis W, Lamkey K, Mahama AA, Suza W (2023). Multiple Trait Selection.
#' In Suza WP and Lamkey KR, editors, *Quantitative Genetics for Plant
#' Breeding*. Iowa State University Digital Press.
#'
#' Cunningham EP (1969). The relative efficiencies of selection indexes.
#' *Acta Agriculturae Scandinavica* 19:45-48.
#'
#' @examples
#' traits <- c("yield", "disease")
#' values <- data.frame(
#'   id = paste0("L", 1:30),
#'   yield = seq(-1.5, 1.5, length.out = 30),
#'   disease = seq(1.5, -1.5, length.out = 30)
#' )
#' G <- matrix(c(1, -0.2, -0.2, 0.8), 2,
#'   dimnames = list(traits, traits)
#' )
#' P <- matrix(c(1.8, -0.2, -0.2, 1.4), 2,
#'   dimnames = list(traits, traits)
#' )
#' objective <- c(yield = 1, disease = 0.5)
#' fits <- list(
#'   Smith_Hazel = selection_index(
#'     values, traits, id_col = "id", method = "smith_hazel",
#'     G = G, P = P, economic_weights = objective,
#'     lower_is_better = "disease", scale_traits = FALSE, n_select = 6
#'   ),
#'   Base = selection_index(
#'     values, traits, id_col = "id", method = "base",
#'     G = G, P = P, economic_weights = objective,
#'     lower_is_better = "disease", scale_traits = FALSE, n_select = 6
#'   )
#' )
#' compare_selection_methods(fits)
#'
#' @seealso [selection_index()], [evaluate_index()],
#'   [candidate_score_uncertainty()], [index_information_efficiency()]
#' @export
compare_selection_methods <- function(models, target_gains = NULL) {
  if (!is.list(models) || length(models) < 2L) {
    stop("models must contain at least two fitted indices.", call. = FALSE)
  }
  if (is.null(names(models)) || any(names(models) == "") ||
    anyDuplicated(names(models))) {
    stop("models must have unique, non-empty names.", call. = FALSE)
  }
  valid <- vapply(
    models, inherits, logical(1L), what = "desiredgainr_index"
  )
  if (!all(valid)) {
    stop(
      "Every entry in models must come from selection_index(). Invalid: ",
      paste(names(models)[!valid], collapse = ", "),
      call. = FALSE
    )
  }

  reference <- models[[1L]]
  traits <- reference$trait_cols
  ids <- reference$candidate_id
  same_traits <- vapply(
    models, function(x) identical(x$trait_cols, traits), logical(1L)
  )
  if (!all(same_traits)) {
    stop("Every model must use the same trait order.", call. = FALSE)
  }
  same_candidates <- vapply(
    models, function(x) identical(x$candidate_id, ids), logical(1L)
  )
  if (!all(same_candidates)) {
    stop("Every model must use the same candidates in the same order.",
      call. = FALSE
    )
  }
  same_transformation <- vapply(models, function(x) {
    isTRUE(all.equal(
      x$transformation, reference$transformation,
      check.attributes = FALSE
    ))
  }, logical(1L))
  if (!all(same_transformation)) {
    stop(
      "Every model must use the same direction, centring, and scaling.",
      call. = FALSE
    )
  }

  target <- NULL
  if (!is.null(target_gains)) {
    target <- .dgr_named_vector(target_gains, traits, "target_gains")
    if (any(target < 0) || !any(target > 0)) {
      stop(
        "target_gains must contain favourable non-negative gains and at ",
        "least one positive value.",
        call. = FALSE
      )
    }
  }

  model_response <- function(model, expected = TRUE) {
    if (isTRUE(expected)) {
      if (is.null(model$evaluation)) {
        return(stats::setNames(rep(NA_real_, length(traits)), traits))
      }
      return(model$evaluation$expected_response[traits])
    }
    observed <- model$observed_differential
    if (is.null(observed)) {
      return(stats::setNames(rep(NA_real_, length(traits)), traits))
    }
    stats::setNames(observed$Differential[match(traits, observed$Trait)], traits)
  }

  response_rows <- lapply(names(models), function(label) {
    model <- models[[label]]
    expected <- model_response(model, TRUE)
    observed <- model_response(model, FALSE)
    data.frame(
      Method = label,
      Trait = traits,
      Expected_response = as.numeric(expected),
      Observed_differential = as.numeric(observed),
      Target = if (is.null(target)) NA_real_ else as.numeric(target),
      Expected_attainment = if (is.null(target)) {
        NA_real_
      } else {
        ifelse(target > 0, as.numeric(expected) / target, NA_real_)
      },
      Observed_attainment = if (is.null(target)) {
        NA_real_
      } else {
        ifelse(target > 0, as.numeric(observed) / target, NA_real_)
      },
      stringsAsFactors = FALSE
    )
  })
  responses <- data.table::rbindlist(response_rows)

  finite_min <- function(x) {
    x <- x[is.finite(x)]
    if (length(x)) min(x) else NA_real_
  }
  finite_mean <- function(x) {
    x <- x[is.finite(x)]
    if (length(x)) mean(x) else NA_real_
  }
  alignment <- function(model, response) {
    if (is.null(target) || any(!is.finite(response)) || is.null(model$G)) {
      return(c(beta = NA_real_, alignment = NA_real_, residual = NA_real_))
    }
    value <- evaluate_restricted_response(response, target, model$G)
    c(
      beta = value$beta,
      alignment = value$mahalanobis_alignment,
      residual = value$mahalanobis_residual
    )
  }

  summary_rows <- lapply(names(models), function(label) {
    model <- models[[label]]
    evaluation <- model$evaluation
    expected <- model_response(model, TRUE)
    diagnostics <- alignment(model, expected)
    method_responses <- responses[responses$Method == label]
    data.frame(
      Method = label,
      Family = model$method,
      Strategy = model$strategy,
      N_selected = model$n_selected,
      Selected_fraction = model$n_selected / model$n_candidates,
      Selection_intensity = model$selection_intensity,
      R_HI = if (is.null(evaluation)) NA_real_ else evaluation$R_HI,
      Delta_H = if (is.null(evaluation)) NA_real_ else evaluation$delta_H,
      RE = if (is.null(evaluation)) NA_real_ else evaluation$RE,
      Index_heritability = if (is.null(evaluation)) {
        NA_real_
      } else {
        evaluation$h2_index
      },
      Index_accuracy = if (is.null(evaluation)) {
        NA_real_
      } else {
        evaluation$accuracy_index
      },
      Worst_expected_attainment = finite_min(
        method_responses$Expected_attainment
      ),
      Mean_expected_attainment = finite_mean(
        method_responses$Expected_attainment
      ),
      Worst_observed_attainment = finite_min(
        method_responses$Observed_attainment
      ),
      Satoh_beta = diagnostics[["beta"]],
      Mahalanobis_alignment = diagnostics[["alignment"]],
      Mahalanobis_residual = diagnostics[["residual"]],
      stringsAsFactors = FALSE
    )
  })
  summary <- data.table::rbindlist(summary_rows)

  rank_matrix <- vapply(models, function(model) {
    candidate_rank <- stats::setNames(model$ranking$rank, model$ranking$id)
    as.numeric(candidate_rank[ids])
  }, numeric(length(ids)))
  colnames(rank_matrix) <- names(models)
  rank_correlation <- stats::cor(
    rank_matrix, method = "spearman", use = "pairwise.complete.obs"
  )

  selected_sets <- lapply(models, function(model) model$selected$id)
  selected_overlap <- outer(
    selected_sets, selected_sets,
    Vectorize(function(a, b) length(intersect(a, b)))
  )
  selected_jaccard <- outer(
    selected_sets, selected_sets,
    Vectorize(function(a, b) {
      union_size <- length(union(a, b))
      if (union_size) length(intersect(a, b)) / union_size else NA_real_
    })
  )
  dimnames(selected_overlap) <- dimnames(selected_jaccard) <-
    list(names(models), names(models))

  merit_models <- models[vapply(models, function(x) {
    !is.null(x$aggregate_weights) && !is.null(x$evaluation)
  }, logical(1L))]
  common_merit <- length(merit_models) >= 2L && all(vapply(
    merit_models[-1L],
    function(x) isTRUE(all.equal(
      x$aggregate_weights, merit_models[[1L]]$aggregate_weights,
      check.attributes = FALSE
    )),
    logical(1L)
  ))
  selected_counts <- vapply(models, `[[`, integer(1L), "n_selected")
  intensities <- vapply(models, function(x) x$selection_intensity, numeric(1L))
  fairness <- data.frame(
    Condition = c(
      "Same candidates and trait order",
      "Same direction, centring, and scaling",
      "Same number selected",
      "Same selection intensity",
      "Common aggregate merit among merit-based comparisons"
    ),
    Satisfied = c(
      TRUE,
      TRUE,
      length(unique(selected_counts)) == 1L,
      length(unique(signif(intensities, 12L))) == 1L,
      common_merit
    ),
    Interpretation = c(
      "Required and enforced by this function.",
      "Required and enforced by this function.",
      "A culling method may select fewer candidates when few pass every limit.",
      "Expected responses require one common selection pressure.",
      "R_HI and Delta_H answer one common question only under one merit vector."
    ),
    stringsAsFactors = FALSE
  )

  result <- list(
    summary = summary,
    responses = responses,
    rank_correlation = rank_correlation,
    selected_jaccard = selected_jaccard,
    selected_overlap = selected_overlap,
    fairness = fairness,
    target_gains = target,
    models = names(models),
    interpretation = paste(
      "Read per-trait response and target attainment first.",
      "Then examine common-merit criteria, ranking agreement, selected-set",
      "agreement, and uncertainty. A method is suitable when it represents",
      "the breeding objective and remains stable under credible perturbations."
    )
  )
  class(result) <- c("desiredgainr_method_comparison", "list")
  result
}

#' @export
print.desiredgainr_method_comparison <- function(x, ...) {
  cat("<desiredgainr_method_comparison>\n")
  cat("  Methods:", paste(x$models, collapse = ", "), "\n")
  failed <- x$fairness$Condition[!x$fairness$Satisfied]
  if (length(failed)) {
    cat("  Comparison conditions requiring attention:\n")
    cat(paste0("    * ", failed, collapse = "\n"), "\n")
  } else {
    cat("  All recorded comparison conditions are satisfied.\n")
  }
  columns <- c(
    "Method", "N_selected", "R_HI", "Delta_H", "RE",
    "Worst_expected_attainment", "Mahalanobis_alignment"
  )
  print(as.data.frame(x$summary)[, columns, drop = FALSE], row.names = FALSE)
  invisible(x)
}
