# Restricted and predetermined-proportional-gain selection indices.
#
# These are the closed-form relatives of what run_dgsi() approaches by
# stochastic search. Where a constraint can be written exactly, solving it
# exactly is preferable to searching for it: the answer is unique, reproducible,
# and carries no optimism from a best-of-n_rep choice.

#' Restricted and proportional-gain selection indices
#'
#' Classical constrained index families. Each imposes a linear constraint on the
#' expected genetic response and maximises merit subject to it, in closed form.
#'
#' @details
#' # Which family answers which question
#'
#' \describe{
#'   \item{`"kempthorne_nordskog"`}{*Hold these traits still.* Maximises the
#'     correlation with net merit subject to zero expected response in the
#'     restricted traits. Use when a trait must not move: a quality standard, a
#'     maturity window, a plant height ceiling.}
#'   \item{`"tallis"`, `"mallard"`, `"harville"`}{*Move these traits in this
#'     ratio.* These published predetermined-proportional-gain formulations are
#'     algebraically equivalent up to an arbitrary scale of the coefficient
#'     vector. They constrain response proportions, not unattainable absolute
#'     response magnitudes. The separate method names retain citation and
#'     reporting provenance.}
#'   \item{`"restricted_smith_hazel"`}{Smith-Hazel with the restriction applied
#'     by projection rather than by Lagrange multipliers. Algebraically
#'     equivalent to `"kempthorne_nordskog"` and retained because the projection
#'     form is what most textbooks present.}
#' }
#'
#' # The algebra
#'
#' Let \eqn{\mathbf{C}} be the \eqn{k \times p} matrix picking out the
#' constrained traits. Expected response is proportional to
#' \eqn{\mathbf{G}\mathbf{b}}, so a zero-response restriction is
#' \eqn{\mathbf{C}\mathbf{G}\mathbf{b} = \mathbf{0}}. Maximising
#' \eqn{\mathbf{b}^\mathsf{T}\mathbf{G}\mathbf{a}} subject to that and to a
#' scale convention gives
#'
#' \deqn{\mathbf{b} = \left[\mathbf{I} - \mathbf{P}^{-1}\mathbf{G}
#' \mathbf{C}^\mathsf{T}(\mathbf{C}\mathbf{G}\mathbf{P}^{-1}\mathbf{G}
#' \mathbf{C}^\mathsf{T})^{-1}\mathbf{C}\mathbf{G}\right]
#' \mathbf{P}^{-1}\mathbf{G}\mathbf{a}}{
#' b = [I - P^-1 G C' (C G P^-1 G C')^-1 C G] P^-1 G a}
#'
#' which is Kempthorne and Nordskog (1959). The proportional-gain methods
#' replace trait selectors with contrasts that are zero precisely when the
#' responses have the requested proportions.
#'
#' # Why these are worth having beside `run_dgsi()`
#'
#' `run_dgsi()` reaches a constrained solution by searching, which introduces a
#' dependence on the number of replicates and a selection optimism that has to
#' be corrected for. Where the constraint is linear, these families give the
#' same answer exactly. Comparing the two is a useful check: a large divergence
#' means the search has not converged, not that the constraint is unusual.
#'
#' @param values Data frame or matrix of candidate trait values.
#' @param trait_cols Trait column names.
#' @param method Index family; see Details.
#' @param G,P Genetic and phenotypic covariance matrices, named by trait.
#' @param economic_weights Named economic weights defining net merit. Required
#'   by zero-restriction families and used to choose among any free response
#'   directions left by a proportionality constraint.
#' @param restricted_traits Traits held to zero expected response, for
#'   `"kempthorne_nordskog"` and `"restricted_smith_hazel"`.
#' @param target_gains Named responses for the constrained traits. Proportions
#'   for `"tallis"`, `"mallard"` and `"harville"`.
#' @param constraint_matrix Optional numeric contrast matrix with one column per
#'   trait. For a zero-restriction method it enforces
#'   `constraint_matrix %*% expected_response = 0`; this supports, for example,
#'   stability contrasts returned by [expand_environments()].
#' @param penalty Deprecated. Finite soft penalties were never Harville's
#'   published method and are rejected; use an exact published constraint.
#' @param lower_is_better Traits for which smaller original values are
#'   favourable.
#' @param id_col Optional column of `values` holding candidate identifiers.
#' @param center_traits,scale_traits,scale_by Passed through to the same
#'   transformation [selection_index()] applies.
#' @param n_select Number of candidates selected.
#' @param selection_intensity Optional standardised selection intensity.
#' @param main_trait Trait against which relative efficiency is computed.
#'
#' @return An object of class `desiredgainr_index`, as returned by
#'   [selection_index()], with an added `constraint` element recording the
#'   restriction and how nearly it was met.
#'
#' @references
#' Kempthorne, O. and Nordskog, A.W. (1959) Restricted selection indices.
#' *Biometrics* 15, 10-19.
#'
#' Tallis, G.M. (1962) A selection index for optimum genotype.
#' *Biometrics* 18, 120-122.
#'
#' Mallard, J. (1972) La theorie et le calcul des index de selection avec
#' restrictions: synthese critique. *Biometrics* 28, 713-735.
#'
#' Kemp, C.D. and Harville, D.A. (1975) Index selection with proportionality
#' constraints.
#' *Biometrics* 31, 223-225.
#' \doi{10.2307/2529722}
#'
#' @examples
#' set.seed(1)
#' traits <- c("yield", "height", "protein")
#' G <- matrix(
#'   c(
#'     1.0, 0.4, 0.2,
#'     0.4, 0.8, 0.1,
#'     0.2, 0.1, 0.5
#'   ), 3,
#'   dimnames = list(traits, traits)
#' )
#' P <- G + diag(c(1.2, 1.0, 0.9))
#' dimnames(P) <- list(traits, traits)
#' values <- as.data.frame(matrix(
#'   stats::rnorm(90),
#'   ncol = 3, dimnames = list(paste0("g", 1:30), traits)
#' ))
#'
#' # Improve yield and protein while holding height still.
#' fit <- restricted_index(
#'   values, traits,
#'   method = "kempthorne_nordskog",
#'   G = G, P = P,
#'   economic_weights = c(yield = 2, height = 0, protein = 1),
#'   restricted_traits = "height", n_select = 6
#' )
#' fit$constraint$achieved_response
#'
#' @seealso [selection_index()], [run_dgsi()], [gain_feasibility()]
#' @export
restricted_index <- function(
  values,
  trait_cols,
  method = c(
    "kempthorne_nordskog", "tallis", "mallard", "harville",
    "restricted_smith_hazel"
  ),
  G,
  P,
  economic_weights = NULL,
  restricted_traits = NULL,
  target_gains = NULL,
  constraint_matrix = NULL,
  penalty = Inf,
  lower_is_better = NULL,
  id_col = NULL,
  center_traits = TRUE,
  scale_traits = TRUE,
  scale_by = c("sample", "phenotypic"),
  n_select = NULL,
  selection_intensity = NULL,
  main_trait = trait_cols[1L]
) {
  method <- match.arg(method)
  scale_by <- match.arg(scale_by)

  # The transformation, identifier handling, ranking and evaluation are shared
  # with selection_index(). Fitting a base index first and then replacing its
  # coefficients keeps the two families on identical footing, so a comparison
  # between them cannot be confounded by a difference in preprocessing.
  base_weights <- economic_weights
  if (is.null(base_weights)) {
    base_weights <- stats::setNames(rep(1, length(trait_cols)), trait_cols)
  }
  fitted <- selection_index(
    values = values, trait_cols = trait_cols, id_col = id_col,
    method = "smith_hazel", G = G, P = P,
    economic_weights = base_weights,
    lower_is_better = lower_is_better,
    center_traits = center_traits, scale_traits = scale_traits,
    scale_by = scale_by, n_select = n_select,
    selection_intensity = selection_intensity, main_trait = main_trait
  )

  G_a <- fitted$G
  P_a <- fitted$P
  a <- .dgr_named_vector(base_weights, trait_cols, "economic_weights")

  constrained <- .dgr_constrained_coefficients(
    method = method, G = G_a, P = P_a, a = a, trait_cols = trait_cols,
    restricted_traits = restricted_traits, target_gains = target_gains,
    constraint_matrix = constraint_matrix, penalty = penalty,
    selection_intensity = fitted$selection_intensity
  )
  result <- .dgr_refit_index(
    fitted, constrained$coefficients, a, constrained$constraint
  )
  proportional <- method %in% c("tallis", "mallard", "harville")
  if (proportional && setequal(names(target_gains), trait_cols)) {
    P_inverse <- .dgr_inverse(P_a, "P")$inverse
    estimated_bv <- fitted$scaled_values %*% P_inverse %*% G_a
    colnames(estimated_bv) <- trait_cols
    result$restricted_breeding_values <- restricted_breeding_values(
      estimated_bv, G_a, direction = target_gains
    )
    result$constraint$satoh_response <- evaluate_restricted_response(
      result$evaluation$expected_response, target_gains, G_a
    )
    result$constraint$note <- paste(
      result$constraint$note,
      "Satoh's beta and restricted breeding values are reported because the",
      "proportional direction covers every trait."
    )
  }
  result
}

#' Solve for coefficients under a linear response constraint
#' @noRd
.dgr_constrained_coefficients <- function(
  method, G, P, a, trait_cols, restricted_traits, target_gains,
  constraint_matrix, penalty, selection_intensity
) {
  p <- length(trait_cols)
  P_inverse <- .dgr_inverse(P, "P")$inverse
  unrestricted <- as.numeric(P_inverse %*% G %*% a)

  build_selector <- function(names_used) {
    unknown <- setdiff(names_used, trait_cols)
    if (length(unknown)) {
      stop("Unknown constrained traits: ", paste(unknown, collapse = ", "),
        call. = FALSE
      )
    }
    selector <- matrix(0,
      nrow = length(names_used), ncol = p,
      dimnames = list(names_used, trait_cols)
    )
    for (i in seq_along(names_used)) selector[i, names_used[i]] <- 1
    selector
  }

  proportional <- method %in% c("tallis", "mallard", "harville")
  if (identical(method, "harville") && is.finite(penalty)) {
    stop(
      "Finite penalty values are not Harville's published proportionality ",
      "method. Supply target_gains and leave penalty = Inf, or use an exact ",
      "zero-response restriction.",
      call. = FALSE
    )
  }

  if (!proportional) {
    if (!is.null(constraint_matrix)) {
      C <- as.matrix(constraint_matrix)
      storage.mode(C) <- "double"
      if (!nrow(C) || ncol(C) != p || any(!is.finite(C))) {
        stop("constraint_matrix must be a finite matrix with one column per ",
          "trait.",
          call. = FALSE
        )
      }
      if (!is.null(colnames(C))) {
        absent <- setdiff(trait_cols, colnames(C))
        if (length(absent)) {
          stop("constraint_matrix is missing traits: ",
            paste(absent, collapse = ", "),
            call. = FALSE
          )
        }
        C <- C[, trait_cols, drop = FALSE]
      } else {
        colnames(C) <- trait_cols
      }
      if (is.null(rownames(C))) {
        rownames(C) <- paste0("constraint_", seq_len(nrow(C)))
      }
    } else {
      if (!length(restricted_traits)) {
        stop("restricted_traits or constraint_matrix is required for method = '",
          method, "'.",
          call. = FALSE
        )
      }
      if (length(restricted_traits) >= p) {
        stop("At least one trait must be left free; restricting all ", p,
          " traits leaves nothing to select on.",
          call. = FALSE
        )
      }
      C <- build_selector(restricted_traits)
    }
    right_hand <- rep(0, nrow(C))
    names(right_hand) <- rownames(C)
  } else {
    if (is.null(target_gains) || !length(target_gains)) {
      stop("target_gains is required for method = '", method, "'.",
        call. = FALSE
      )
    }
    constrained_names <- names(target_gains)
    if (is.null(constrained_names)) {
      stop("target_gains must be a named vector.", call. = FALSE)
    }
    C <- build_selector(constrained_names)
    right_hand <- as.numeric(target_gains)
    names(right_hand) <- constrained_names
    # Only ratios are constrained. These k - 1 contrasts vanish exactly when
    # the expected responses are proportional to target_gains. Because the
    # standardising factor i/sd(I) is common to every trait, imposing this on
    # G b is equivalent to imposing it on the actual expected response.
    if (length(right_hand) < 2L) {
      stop("A proportional-gain method needs target_gains for at least two ",
        "traits.",
        call. = FALSE
      )
    }
    if (any(right_hand == 0)) {
      stop("A zero target is a zero-response restriction; express it with ",
        "restricted_traits or constraint_matrix instead.",
        call. = FALSE
      )
    }
    difference <- matrix(
      0,
      nrow = length(right_hand) - 1L, ncol = nrow(C),
      dimnames = list(
        paste0("proportion_", seq_len(length(right_hand) - 1L)),
        names(right_hand)
      )
    )
    for (i in seq_len(nrow(difference))) {
      difference[i, i] <- 1 / right_hand[i]
      difference[i, i + 1L] <- -1 / right_hand[i + 1L]
    }
    C <- difference %*% C
    right_hand <- rep(0, nrow(C))
    names(right_hand) <- rownames(C)
  }

  CG <- C %*% G
  middle <- CG %*% P_inverse %*% t(CG)
  middle <- (middle + t(middle)) / 2
  middle_inverse <- .dgr_inverse(middle, "C G P^-1 G C'")$inverse

  correction <- P_inverse %*% t(CG) %*%
    (middle_inverse %*% (CG %*% unrestricted - right_hand))
  coefficients <- as.numeric(unrestricted - correction)
  names(coefficients) <- trait_cols

  index_sd <- sqrt(as.numeric(crossprod(coefficients, P %*% coefficients)))
  expected <- if (is.finite(selection_intensity) && index_sd > 0) {
    as.numeric(selection_intensity * G %*% coefficients / index_sd)
  } else {
    as.numeric(G %*% coefficients / index_sd)
  }
  names(expected) <- trait_cols
  constraint_value <- as.numeric(C %*% expected)
  names(constraint_value) <- rownames(C)
  achieved <- if (proportional) {
    expected[names(target_gains)]
  } else if (!is.null(constraint_matrix)) {
    constraint_value
  } else {
    expected[restricted_traits]
  }
  list(
    coefficients = coefficients,
    constraint = list(
      method = method,
      restricted_traits = restricted_traits,
      target_gains = target_gains,
      penalty = penalty,
      constraint_matrix = C,
      required_response = right_hand,
      achieved_response = achieved,
      constraint_value = constraint_value,
      largest_violation = max(abs(constraint_value)),
      note = switch(method,
        kempthorne_nordskog = paste(
          "Expected response in the restricted traits is held at zero. The",
          "cost is a lower correlation with net merit than the unrestricted",
          "Smith-Hazel index; compare R_HI between the two."
        ),
        restricted_smith_hazel = paste(
          "Algebraically identical to Kempthorne-Nordskog; retained because",
          "the projection form is the one most texts present."
        ),
        tallis = paste(
          "Responses in the constrained traits are held in the supplied",
          "proportions. The magnitude is not constrained, because attainable",
          "magnitude is set by selection intensity rather than by the index."
        ),
        mallard = paste(
          "Mallard's predetermined-gain formulation is used as a response",
          "proportionality constraint. Absolute magnitude is not claimed,",
          "because multiplying every coefficient by a constant leaves the",
          "selected ranking and standardized expected response unchanged."
        ),
        harville = paste(
          "Kemp and Harville's published proportionality constraint is used.",
          "This is an exact contrast constraint, not a locally invented soft",
          "penalty."
        )
      )
    )
  )
}

#' Rebuild a fitted index around replacement coefficients
#' @noRd
.dgr_refit_index <- function(fitted, coefficients, aggregate_weights,
                             constraint) {
  X <- fitted$scaled_values
  trait_cols <- fitted$trait_cols
  scores <- as.numeric(X %*% coefficients)
  candidate_id <- fitted$candidate_id
  n_select <- fitted$n_select

  selection <- rep(FALSE, nrow(X))
  if (!is.null(n_select)) {
    selection[order(-scores, candidate_id)[seq_len(n_select)]] <- TRUE
  }
  n_selected <- sum(selection)
  intensity <- if (!n_selected) NA_real_ else .dgr_intensity(n_selected / nrow(X))

  ranking <- data.table::data.table(
    id = candidate_id,
    score = scores,
    rank = data.table::frank(-scores, ties.method = "min"),
    selected = selection
  )
  data.table::setorderv(ranking, c("score", "id"), c(-1L, 1L))

  observed <- NULL
  if (any(selection)) {
    observed <- data.table::data.table(
      Trait = trait_cols,
      Mean_all = colMeans(X),
      Mean_selected = colMeans(X[selection, , drop = FALSE]),
      Differential = colMeans(X[selection, , drop = FALSE]) - colMeans(X)
    )
  }

  fitted$method <- constraint$method
  fitted$coefficients <- coefficients
  fitted$aggregate_weights <- aggregate_weights
  fitted$ranking <- ranking
  fitted$selected <- ranking[selected == TRUE]
  fitted$n_selected <- n_selected
  fitted$selection_intensity <- intensity
  fitted$observed_differential <- observed
  fitted$objective <- aggregate_weights
  fitted$constraint <- constraint
  fitted$evaluation <- evaluate_index(
    coefficients = coefficients,
    G = fitted$G,
    P = if (is.null(fitted$P)) stats::cov(X) else fitted$P,
    aggregate_weights = aggregate_weights,
    selection_intensity = intensity,
    main_trait = fitted$evaluation$main_trait,
    scores = scores
  )
  fitted$effective_weights <- effective_weights(
    coefficients, fitted$G, fitted$P
  )
  fitted
}
