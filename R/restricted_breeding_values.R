# Restricted breeding values and the Satoh evaluation criterion.

#' Project breeding values into a restricted response space
#'
#' Satoh (2024) showed that restricted breeding values are linear projections
#' of ordinary breeding values. For constraints \eqn{C g=0}, the projection is
#'
#' \deqn{g_R=[I-GC^\mathsf{T}(CGC^\mathsf{T})^{-1}C]g.}{
#' g_R = [I - G C' (C G C')^-1 C] g.}
#'
#' A proportional desired-gain direction \eqn{d} defines a one-dimensional
#' space. Its projection is
#'
#' \deqn{g_R=d(d^\mathsf{T}G^{-1}d)^{-1}d^\mathsf{T}G^{-1}g
#' =\beta d.}{g_R = d (d' G^-1 d)^-1 d' G^-1 g = beta d.}
#'
#' The scalar \eqn{\beta} measures progress along the requested direction. It
#' combines response magnitude and proportional agreement in one quantity.
#'
#' @param breeding_values Candidate-by-trait matrix of breeding values or
#'   estimated breeding values.
#' @param G Genetic covariance matrix for the same traits.
#' @param direction Optional named proportional desired-gain direction.
#' @param constraint_matrix Optional restriction matrix with traits in columns.
#' @param id_col Optional identifier column in breeding_values.
#' @param lower_is_better Traits for which smaller original values are
#'   favourable. When supplied, direction entries are favourable magnitudes.
#'
#' @return An object of class desiredgainr_restricted_bv.
#'
#' @references
#' Satoh M (2024). Characteristics of restricted selection indices and
#' geometrical interpretation of restricted breeding values. *Journal of
#' Animal Breeding and Genetics* 141:353-363.
#' \doi{10.1111/jbg.12845}
#'
#' @seealso [restricted_index()], [evaluate_restricted_response()]
#' @export
restricted_breeding_values <- function(
  breeding_values,
  G,
  direction = NULL,
  constraint_matrix = NULL,
  id_col = NULL,
  lower_is_better = NULL
) {
  if (is.null(direction) == is.null(constraint_matrix)) {
    stop("Supply one of direction or constraint_matrix.", call. = FALSE)
  }
  frame <- as.data.frame(breeding_values)
  traits <- rownames(G)
  if (is.null(traits) || !identical(traits, colnames(G))) {
    stop("G must have matching trait names.", call. = FALSE)
  }
  absent <- setdiff(traits, names(frame))
  if (length(absent)) {
    stop("breeding_values is missing traits: ",
      paste(absent, collapse = ", "), call. = FALSE
    )
  }
  values <- as.matrix(frame[, traits, drop = FALSE])
  storage.mode(values) <- "double"
  if (!nrow(values) || any(!is.finite(values))) {
    stop("breeding_values must contain finite values.", call. = FALSE)
  }
  G <- .dgr_covariance(G, traits, "G")
  .dgr_check_psd(G, "G")

  ids <- if (is.null(id_col)) {
    value <- rownames(values)
    if (is.null(value)) paste0("row_", seq_len(nrow(values))) else value
  } else {
    if (!id_col %in% names(frame)) {
      stop("id_col was not found in breeding_values.", call. = FALSE)
    }
    as.character(frame[[id_col]])
  }
  if (anyNA(ids) || anyDuplicated(ids)) {
    stop("Candidate identifiers must be unique and complete.", call. = FALSE)
  }

  beta <- NULL
  if (!is.null(direction)) {
    oriented_direction <- .dgr_orient_objective(
      direction, traits, lower_is_better, "direction"
    )
    d <- oriented_direction$signed
    if (all(d == 0)) {
      stop("direction needs at least one non-zero value.", call. = FALSE)
    }
    G_inverse <- .dgr_inverse(G, "G")$inverse
    denominator <- as.numeric(crossprod(d, G_inverse %*% d))
    projection <- d %*% (t(d) %*% G_inverse) / denominator
    beta <- as.numeric(values %*% G_inverse %*% d / denominator)
    names(beta) <- ids
    restriction <- list(
      type = "proportional",
      direction = d,
      direction_input = oriented_direction$input,
      lower_is_better = lower_is_better
    )
  } else {
    C <- as.matrix(constraint_matrix)
    storage.mode(C) <- "double"
    if (is.null(colnames(C))) {
      if (ncol(C) != length(traits)) {
        stop("constraint_matrix must have one column per trait.", call. = FALSE)
      }
      colnames(C) <- traits
    } else {
      absent_C <- setdiff(traits, colnames(C))
      if (length(absent_C)) {
        stop("constraint_matrix is missing traits: ",
          paste(absent_C, collapse = ", "), call. = FALSE
        )
      }
      C <- C[, traits, drop = FALSE]
    }
    if (!nrow(C) || any(!is.finite(C))) {
      stop("constraint_matrix must contain finite restrictions.", call. = FALSE)
    }
    middle <- C %*% G %*% t(C)
    projection <- diag(length(traits)) - G %*% t(C) %*%
      .dgr_inverse(middle, "C G C'")$inverse %*% C
    restriction <- list(type = "linear", constraint_matrix = C)
  }
  dimnames(projection) <- list(traits, traits)
  projected <- values %*% t(projection)
  colnames(projected) <- traits
  rownames(projected) <- ids
  violation <- if (is.null(direction)) {
    apply(abs(projected %*% t(C)), 1L, max)
  } else {
    residual <- projected - beta * matrix(
      d, nrow = nrow(projected), ncol = length(d), byrow = TRUE
    )
    apply(abs(residual), 1L, max)
  }

  result <- list(
    ordinary = values,
    restricted = projected,
    beta = beta,
    projection = projection,
    restriction = restriction,
    largest_violation = max(violation),
    candidate_id = ids,
    interpretation = if (is.null(beta)) {
      paste(
        "Each breeding-value vector was projected into the linear restricted",
        "space. The projected values satisfy the supplied contrasts."
      )
    } else {
      paste(
        "Each restricted breeding-value vector equals beta times the desired",
        "direction. Larger beta means greater progress along that direction."
      )
    }
  )
  class(result) <- c("desiredgainr_restricted_bv", "list")
  result
}

#' Evaluate an achieved response against a desired-gain direction
#'
#' This function applies Satoh's one-dimensional restricted breeding-value
#' criterion to an achieved response. It reports progress along the direction,
#' proportional departure, and the residual Mahalanobis distance.
#'
#' @param response Named achieved response vector.
#' @param direction Named proportional desired-gain direction.
#' @param G Genetic covariance matrix.
#' @param lower_is_better Traits for which smaller original values are
#'   favourable. When supplied, direction entries are favourable magnitudes.
#'
#' @return An object of class desiredgainr_restricted_response.
#'
#' @export
evaluate_restricted_response <- function(
  response, direction, G, lower_is_better = NULL
) {
  traits <- rownames(G)
  if (is.null(traits) || !identical(traits, colnames(G))) {
    stop("G must have matching trait names.", call. = FALSE)
  }
  response <- .dgr_named_vector(response, traits, "response")
  direction <- .dgr_orient_objective(
    direction, traits, lower_is_better, "direction"
  )$signed
  if (all(direction == 0)) {
    stop("direction needs at least one non-zero value.", call. = FALSE)
  }
  G <- .dgr_covariance(G, traits, "G")
  G_inverse <- .dgr_inverse(G, "G")$inverse
  denominator <- as.numeric(crossprod(direction, G_inverse %*% direction))
  beta <- as.numeric(crossprod(direction, G_inverse %*% response)) /
    denominator
  projected <- beta * direction
  residual <- response - projected
  residual_distance <- sqrt(max(
    0, as.numeric(crossprod(residual, G_inverse %*% residual))
  ))
  response_distance <- sqrt(max(
    0, as.numeric(crossprod(response, G_inverse %*% response))
  ))
  alignment <- if (response_distance > 0) {
    as.numeric(crossprod(response, G_inverse %*% direction)) /
      (response_distance * sqrt(denominator))
  } else {
    NA_real_
  }
  result <- list(
    beta = beta,
    response = response,
    direction = direction,
    projected_response = projected,
    residual_response = residual,
    mahalanobis_residual = residual_distance,
    mahalanobis_alignment = alignment,
    interpretation = paste(
      "Beta measures progress along the desired-gain direction.",
      "The residual distance measures departure from the requested proportions."
    )
  )
  class(result) <- c("desiredgainr_restricted_response", "list")
  result
}

#' @export
print.desiredgainr_restricted_bv <- function(x, ...) {
  cat("<desiredgainr_restricted_bv>\n")
  cat("  Restriction:", x$restriction$type, "\n")
  cat("  Candidates:", nrow(x$restricted), "\n")
  cat("  Largest numerical violation:",
    format(x$largest_violation, digits = 4), "\n"
  )
  if (!is.null(x$beta)) {
    cat("  Beta range:",
      format(min(x$beta), digits = 4), "to",
      format(max(x$beta), digits = 4), "\n"
    )
  }
  invisible(x)
}

#' @export
print.desiredgainr_restricted_response <- function(x, ...) {
  cat("<desiredgainr_restricted_response>\n")
  cat("  Satoh beta:", format(x$beta, digits = 5), "\n")
  cat("  Mahalanobis alignment:",
    format(x$mahalanobis_alignment, digits = 5), "\n"
  )
  cat("  Mahalanobis residual:",
    format(x$mahalanobis_residual, digits = 5), "\n"
  )
  invisible(x)
}
