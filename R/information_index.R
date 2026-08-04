# Selection indices with distinct information and objective vectors.

.dgr_orient_objective <- function(
  objective, objective_names, lower_is_better, name
) {
  input <- .dgr_named_vector(objective, objective_names, name)
  signs <- stats::setNames(rep(1, length(objective_names)), objective_names)
  if (length(lower_is_better)) {
    unknown <- setdiff(lower_is_better, objective_names)
    if (length(unknown)) {
      stop("Unknown lower-is-better objective traits: ",
        paste(unknown, collapse = ", "), call. = FALSE
      )
    }
    if (any(input < 0)) {
      stop(
        name, " must contain favourable magnitudes when lower_is_better is used.",
        call. = FALSE
      )
    }
    signs[lower_is_better] <- -1
  }
  list(input = input, signed = input * signs, signs = signs)
}

#' Define the information used by a selection index
#'
#' Mixed-model and genomic analyses can use records, family means, genomic
#' estimated breeding values, and environment-specific predictions as
#' selection information. The breeding objective can contain a different set
#' of genetic quantities. Let \eqn{P=Var(x)}, \eqn{C=Cov(x,g)}, and
#' \eqn{G=Var(g)}. This function stores that model and checks its joint
#' covariance.
#'
#' @param values Candidate-by-information matrix or data frame.
#' @param P Covariance matrix of the information variables.
#' @param C Covariance between information variables and objective traits.
#' @param G Genetic covariance matrix of the objective traits.
#' @param id_col Optional candidate identifier column in values.
#'
#' @return An object of class desiredgainr_information.
#'
#' @references
#' Beavis WD, Lamkey K, Mahama AA, Suza W (2023). Multiple Trait Selection.
#' In *Quantitative Genetics for Plant Breeding*. Iowa State University Digital
#' Press.
#'
#' Henderson CR, Quaas RL (1976). Multiple trait evaluation using relatives'
#' records. *Journal of Animal Science* 43:1188-1197.
#'
#' @seealso [generalized_index()], [candidate_score_uncertainty()]
#' @export
selection_information <- function(values, P, C, G, id_col = NULL) {
  frame <- as.data.frame(values)
  information_names <- rownames(P)
  if (is.null(information_names) || !identical(information_names, colnames(P))) {
    stop("P must have matching information names on both dimensions.",
      call. = FALSE
    )
  }
  if (!length(information_names) || anyDuplicated(information_names)) {
    stop("P must contain unique information names.", call. = FALSE)
  }
  absent <- setdiff(information_names, names(frame))
  if (length(absent)) {
    stop("values is missing information columns: ",
      paste(absent, collapse = ", "), call. = FALSE
    )
  }
  X <- as.matrix(frame[, information_names, drop = FALSE])
  storage.mode(X) <- "double"
  if (!nrow(X) || any(!is.finite(X))) {
    stop("values must contain finite information values.", call. = FALSE)
  }

  P <- .dgr_covariance(P, information_names, "P")
  .dgr_check_psd(P, "P")
  objective_names <- rownames(G)
  if (is.null(objective_names) || !identical(objective_names, colnames(G)) ||
    anyDuplicated(objective_names)) {
    stop("G must have matching and unique objective names.", call. = FALSE)
  }
  G <- .dgr_covariance(G, objective_names, "G")
  .dgr_check_psd(G, "G")

  C <- as.matrix(C)
  storage.mode(C) <- "double"
  if (is.null(rownames(C)) || is.null(colnames(C))) {
    stop("C must carry information names by row and objective names by column.",
      call. = FALSE
    )
  }
  missing_rows <- setdiff(information_names, rownames(C))
  missing_columns <- setdiff(objective_names, colnames(C))
  if (length(missing_rows) || length(missing_columns)) {
    stop("C must cover every information variable and objective trait.",
      call. = FALSE
    )
  }
  C <- C[information_names, objective_names, drop = FALSE]
  if (any(!is.finite(C))) {
    stop("C must contain finite covariances.", call. = FALSE)
  }

  joint <- rbind(cbind(P, C), cbind(t(C), G))
  joint <- (joint + t(joint)) / 2
  joint_values <- eigen(joint, symmetric = TRUE, only.values = TRUE)$values
  tolerance <- 1e-8 * max(1, max(abs(diag(joint))))
  if (min(joint_values) < -tolerance) {
    stop(
      "P, C, and G do not form a valid joint covariance matrix. The smallest ",
      "joint eigenvalue is ", format(min(joint_values), digits = 4), ".",
      call. = FALSE
    )
  }

  if (is.null(id_col)) {
    candidate_id <- rownames(X)
    if (is.null(candidate_id)) candidate_id <- paste0("row_", seq_len(nrow(X)))
  } else {
    if (!id_col %in% names(frame)) {
      stop("id_col was not found in values.", call. = FALSE)
    }
    candidate_id <- as.character(frame[[id_col]])
  }
  if (anyNA(candidate_id) || anyDuplicated(candidate_id)) {
    stop("Candidate identifiers must be unique and complete.", call. = FALSE)
  }

  result <- list(
    values = X,
    candidate_id = candidate_id,
    information_names = information_names,
    objective_names = objective_names,
    P = P,
    C = C,
    G = G,
    joint_minimum_eigenvalue = min(joint_values),
    provenance = paste(
      length(information_names), "selection variables predict",
      length(objective_names), "objective traits. P = Var(x), C = Cov(x, g),",
      "and G = Var(g)."
    )
  )
  class(result) <- c("desiredgainr_information", "list")
  result
}

#' Fit an index from a general information model
#'
#' The economic solution is \eqn{b=P^{-1}Ca}. The desired-gain solution is
#' \eqn{b=P^{-1}C(C^\mathsf{T}P^{-1}C)^{-1}d}. The expected response is
#' \eqn{\Delta g=iC^\mathsf{T}b/\sqrt{b^\mathsf{T}Pb}}. Hence the records and
#' objective traits may differ in number and meaning.
#'
#' @param model Object returned by [selection_information()].
#' @param objective Named economic weights or desired gains. Desired gains must
#'   use the trait units of `model$G`. Convert genetic-standard-deviation gains
#'   by multiplying them by `sqrt(diag(model$G))`.
#' @param method Either economic or desired_gain.
#' @param n_select Optional number of candidates selected.
#' @param selection_intensity Optional standardised selection intensity.
#' @param lower_is_better Objective traits for which smaller original values
#'   are favourable. When supplied, objective entries are favourable
#'   magnitudes. The function gives these traits a negative sign in the
#'   original trait coordinates.
#'
#' @return An object of class desiredgainr_generalized_index.
#'
#' @seealso [selection_information()], [selection_index()]
#' @export
generalized_index <- function(
  model,
  objective,
  method = c("economic", "desired_gain"),
  n_select = NULL,
  selection_intensity = NULL,
  lower_is_better = NULL
) {
  if (!inherits(model, "desiredgainr_information")) {
    stop("model must come from selection_information().", call. = FALSE)
  }
  method <- match.arg(method)
  oriented_objective <- .dgr_orient_objective(
    objective, model$objective_names, lower_is_better, "objective"
  )
  objective_input <- oriented_objective$input
  objective <- oriented_objective$signed
  if (identical(method, "desired_gain") && all(objective == 0)) {
    stop("A desired-gain objective needs at least one non-zero entry.",
      call. = FALSE
    )
  }

  P_inverse <- .dgr_inverse(model$P, "P")$inverse
  if (identical(method, "economic")) {
    coefficients <- P_inverse %*% model$C %*% objective
  } else {
    estimable <- crossprod(model$C, P_inverse %*% model$C)
    estimable <- (estimable + t(estimable)) / 2
    dimnames(estimable) <- list(model$objective_names, model$objective_names)
    coefficients <- P_inverse %*% model$C %*%
      (.dgr_inverse(estimable, "C' P^-1 C")$inverse %*% objective)
  }
  coefficients <- as.numeric(coefficients)
  names(coefficients) <- model$information_names
  scores <- as.numeric(model$values %*% coefficients)

  n_candidates <- nrow(model$values)
  selected <- rep(FALSE, n_candidates)
  if (!is.null(n_select)) {
    n_select <- .dgr_positive_integer(n_select, "n_select")
    if (n_select > n_candidates) {
      stop("n_select exceeds the number of candidates.", call. = FALSE)
    }
    keep <- order(-scores, model$candidate_id)[seq_len(n_select)]
    selected[keep] <- TRUE
  }
  if (is.null(selection_intensity)) {
    selection_intensity <- if (is.null(n_select)) {
      NA_real_
    } else {
      .dgr_intensity(n_select / n_candidates)
    }
  }
  if (!is.numeric(selection_intensity) ||
    length(selection_intensity) != 1L ||
    (!is.na(selection_intensity) &&
      (!is.finite(selection_intensity) || selection_intensity < 0))) {
    stop("selection_intensity must be one non-negative finite value.",
      call. = FALSE
    )
  }

  index_variance <- as.numeric(crossprod(
    coefficients, model$P %*% coefficients
  ))
  if (!is.finite(index_variance) || index_variance <= 0) {
    stop("The fitted index has non-positive variance.", call. = FALSE)
  }
  index_sd <- sqrt(index_variance)
  response <- as.numeric(
    selection_intensity * crossprod(model$C, coefficients) / index_sd
  )
  names(response) <- model$objective_names

  accuracy <- NA_real_
  delta_merit <- NA_real_
  if (identical(method, "economic")) {
    merit_variance <- as.numeric(crossprod(
      objective, model$G %*% objective
    ))
    if (merit_variance > 0) {
      accuracy <- as.numeric(crossprod(
        coefficients, model$C %*% objective
      )) / (index_sd * sqrt(merit_variance))
      accuracy <- .dgr_clamp_correlation(
        accuracy, "The generalised index accuracy"
      )
      delta_merit <- selection_intensity * accuracy * sqrt(merit_variance)
    }
  }

  ranking <- data.table::data.table(
    id = model$candidate_id,
    score = scores,
    rank = data.table::frank(-scores, ties.method = "min"),
    selected = selected
  )
  data.table::setorderv(ranking, c("score", "id"), c(-1L, 1L))

  result <- list(
    method = method,
    model = model,
    objective = objective,
    objective_input = objective_input,
    lower_is_better = lower_is_better,
    coefficients = coefficients,
    scores = scores,
    candidate_id = model$candidate_id,
    ranking = ranking,
    selected = ranking[selected == TRUE],
    n_select = n_select,
    selection_intensity = selection_intensity,
    index_variance = index_variance,
    index_sd = index_sd,
    expected_response = response,
    accuracy = accuracy,
    delta_merit = delta_merit,
    interpretation = paste(
      "The coefficients weight selection information. Expected responses",
      "refer to the objective traits. Their names and dimensions can differ."
    )
  )
  class(result) <- c("desiredgainr_generalized_index", "list")
  result
}

#' Candidate-specific uncertainty in index scores and ranks
#'
#' For candidate \eqn{i}, the score standard error is
#' \eqn{\sqrt{b^\mathsf{T}PEC_i b}}. Monte Carlo draws then estimate selection
#' probabilities. Within-candidate correlations among variables are retained.
#' Prediction errors are treated as independent among candidates.
#'
#' @param index A fitted DesiredGainR coefficient-based index.
#' @param prediction_error_covariance One shared covariance matrix, one matrix
#'   per candidate in a list, or a three-dimensional array.
#' @param n_draws Number of Monte Carlo draws.
#' @param level Coverage of the score interval.
#' @param seed Optional random seed.
#'
#' @return An object of class desiredgainr_candidate_uncertainty.
#'
#' @export
candidate_score_uncertainty <- function(
  index,
  prediction_error_covariance,
  n_draws = 2000L,
  level = 0.95,
  seed = NULL
) {
  general <- inherits(index, "desiredgainr_generalized_index")
  conventional <- inherits(index, "desiredgainr_index")
  if (!general && !conventional) {
    stop("index must be a fitted DesiredGainR index.", call. = FALSE)
  }
  b <- index$coefficients
  if (is.null(b)) {
    stop("The fitted method has no coefficient vector.", call. = FALSE)
  }
  variable_names <- names(b)
  ids <- index$candidate_id
  scores <- if (general) {
    index$scores
  } else {
    raw <- numeric(length(ids))
    raw[match(index$ranking$id, ids)] <- index$ranking$score
    raw
  }
  n <- length(ids)
  n_draws <- .dgr_positive_integer(n_draws, "n_draws")
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
    level <= 0 || level >= 1) {
    stop("level must lie strictly between zero and one.", call. = FALSE)
  }

  dimensions <- dim(prediction_error_covariance)
  matrices <- if (is.list(prediction_error_covariance) &&
    !is.matrix(prediction_error_covariance)) {
    prediction_error_covariance
  } else if (length(dimensions) == 3L) {
    lapply(seq_len(dimensions[3L]), function(i) {
      prediction_error_covariance[, , i]
    })
  } else {
    rep(list(prediction_error_covariance), n)
  }
  if (length(matrices) != n) {
    stop("Supply one prediction error covariance for each candidate.",
      call. = FALSE
    )
  }
  score_variance <- vapply(seq_len(n), function(i) {
    M <- .dgr_covariance(matrices[[i]], variable_names, paste0("PEC[", i, "]"))
    .dgr_check_psd(M, paste0("PEC[", i, "]"))
    max(0, as.numeric(crossprod(b, M %*% b)))
  }, numeric(1L))
  score_se <- sqrt(score_variance)

  if (!is.null(seed)) {
    if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      stats::runif(1L)
    }
    saved_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", saved_seed, envir = globalenv()), add = TRUE)
    set.seed(seed)
  }
  draw_matrix <- matrix(stats::rnorm(n * n_draws), nrow = n)
  draw_matrix <- scores + score_se * draw_matrix
  n_select <- index$n_select
  selection_probability <- rep(NA_real_, n)
  if (!is.null(n_select)) {
    selected_draws <- matrix(FALSE, nrow = n, ncol = n_draws)
    for (draw in seq_len(n_draws)) {
      keep <- order(-draw_matrix[, draw], ids)[seq_len(n_select)]
      selected_draws[keep, draw] <- TRUE
    }
    selection_probability <- rowMeans(selected_draws)
  }
  alpha <- (1 - level) / 2
  summary <- data.table::data.table(
    id = ids,
    score = scores,
    score_se = score_se,
    lower = apply(draw_matrix, 1L, stats::quantile, probs = alpha),
    upper = apply(draw_matrix, 1L, stats::quantile, probs = 1 - alpha),
    selection_probability = selection_probability
  )
  data.table::setorderv(summary, c("score", "id"), c(-1L, 1L))
  result <- list(
    summary = summary,
    n_draws = n_draws,
    level = level,
    n_select = n_select,
    assumption = paste(
      "Prediction errors are correlated across variables within a candidate.",
      "Errors are independent among candidates."
    )
  )
  class(result) <- c("desiredgainr_candidate_uncertainty", "list")
  result
}

#' Cunningham information-deletion efficiency
#'
#' A value of 0.98 means that the index after deleting one information variable
#' retains 98 percent of the full optimum-index standard deviation.
#'
#' @param index A Smith-Hazel index or a general economic index.
#'
#' @return A data.table with one row per information variable.
#'
#' @references
#' Cunningham EP (1969). The relative efficiencies of selection indexes.
#' *Acta Agriculturae Scandinavica* 19:45-48.
#'
#' @export
index_information_efficiency <- function(index) {
  classical <- inherits(index, "desiredgainr_index") &&
    identical(index$method, "smith_hazel")
  general <- inherits(index, "desiredgainr_generalized_index") &&
    identical(index$method, "economic")
  if (!classical && !general) {
    stop(
      "Cunningham efficiency requires a Smith-Hazel or general economic index.",
      call. = FALSE
    )
  }
  P <- if (classical) index$P else index$model$P
  b <- index$coefficients
  inverse <- .dgr_inverse(P, "P")$inverse
  full <- as.numeric(crossprod(b, P %*% b))
  reduced <- pmax(0, full - b^2 / diag(inverse))
  efficiency <- sqrt(reduced / full)
  result <- data.table::data.table(
    Information = names(b),
    Full_index_efficiency = 1,
    Efficiency_after_deletion = as.numeric(efficiency),
    Proportional_loss = as.numeric(1 - efficiency),
    Response_variance_lost = as.numeric(1 - efficiency^2)
  )
  data.table::setorderv(result, "Efficiency_after_deletion")
  result
}

#' @export
print.desiredgainr_information <- function(x, ...) {
  cat("<desiredgainr_information>\n")
  cat("  Candidates:", nrow(x$values), "\n")
  cat("  Information variables:", length(x$information_names), "\n")
  cat("  Objective traits:", length(x$objective_names), "\n")
  cat("  Joint minimum eigenvalue:",
    format(x$joint_minimum_eigenvalue, digits = 4), "\n"
  )
  invisible(x)
}

#' @export
print.desiredgainr_generalized_index <- function(x, ...) {
  cat("<desiredgainr_generalized_index>\n")
  cat("  Method:", x$method, "\n")
  cat("  Information coefficients:\n")
  print(round(x$coefficients, 5L))
  cat("  Expected objective response:\n")
  print(round(x$expected_response, 5L))
  if (is.finite(x$accuracy)) {
    cat("  Accuracy for aggregate merit:", format(x$accuracy, digits = 4), "\n")
  }
  invisible(x)
}

#' @export
print.desiredgainr_candidate_uncertainty <- function(x, ...) {
  cat("<desiredgainr_candidate_uncertainty>\n")
  cat("  Monte Carlo draws:", x$n_draws, "\n")
  print(utils::head(x$summary, 10L))
  invisible(x)
}
