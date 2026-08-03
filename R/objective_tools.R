# Tools for defining and interrogating a multi-trait breeding objective.
#
# These functions do not invent an objective. They translate between the two
# ways breeders state one, test whether a stated objective is attainable, and
# quantify how sensitive a decision is to it.

#' Solve a symmetric linear system and report its conditioning
#'
#' @param M Symmetric matrix.
#' @param name Object name used in messages.
#' @param warn_condition Reciprocal-condition threshold below which a warning
#'   is issued. The default corresponds to a condition number of one hundred
#'   million, the conventional point at which roughly half the available
#'   double precision has been lost to the inversion.
#'
#' @return A list with the inverse, the reciprocal condition number, and the
#'   eigenvalues.
#' @noRd
.dgr_inverse <- function(M, name, warn_condition = 1e-8) {
  decomposition <- eigen(M, symmetric = TRUE)
  values <- decomposition$values
  largest <- max(abs(values))
  reciprocal_condition <- if (largest > 0) min(abs(values)) / largest else 0
  if (reciprocal_condition < warn_condition) {
    warning(
      sprintf(
        paste(
          "%s is numerically ill-conditioned (reciprocal condition number",
          "%.3g). Coefficients obtained by inverting it are unreliable.",
          "Standardising the traits before analysis usually resolves this."
        ),
        name, reciprocal_condition
      ),
      call. = FALSE
    )
  }
  if (min(values) <= 0) {
    stop(name, " must be positive definite to be inverted.", call. = FALSE)
  }
  inverse <- decomposition$vectors %*%
    (t(decomposition$vectors) / values)
  inverse <- (inverse + t(inverse)) / 2
  dimnames(inverse) <- dimnames(M)
  list(
    inverse = inverse,
    reciprocal_condition = reciprocal_condition,
    eigenvalues = values
  )
}

#' Standardised selection intensity for a selected proportion
#'
#' @param proportion Selected proportion in `(0, 1]`.
#'
#' @return The standardised selection intensity under normality.
#' @noRd
.dgr_intensity <- function(proportion) {
  if (proportion >= 1) {
    return(0)
  }
  stats::dnorm(stats::qnorm(1 - proportion)) / proportion
}

#' Selected proportion delivering a given selection intensity
#'
#' @param intensity Positive standardised selection intensity.
#'
#' @return The selected proportion, or `NA_real_` when the intensity exceeds
#'   what truncation selection can deliver in a finite population.
#' @noRd
.dgr_proportion_for_intensity <- function(intensity) {
  if (!is.finite(intensity) || intensity <= 0) {
    return(NA_real_)
  }
  objective <- function(p) .dgr_intensity(p) - intensity
  lower <- 1e-12
  if (objective(lower) < 0) {
    return(NA_real_)
  }
  stats::uniroot(objective, interval = c(lower, 1 - 1e-12))$root
}

#' Convert a desired-gain vector to original trait units
#'
#' @param d Named numeric vector.
#' @param G,P Covariance matrices in original trait units.
#' @param units One of `"trait"`, `"genetic_sd"`, or `"phenotypic_sd"`.
#'
#' @return The desired gains expressed in original trait units.
#' @noRd
.dgr_gain_to_trait_units <- function(d, G, P, units) {
  switch(units,
    trait = d,
    genetic_sd = d * sqrt(diag(G)),
    phenotypic_sd = d * sqrt(diag(P))
  )
}

#' Convert a vector from original trait units to a requested scale
#'
#' @param d Named numeric vector in original trait units.
#' @param G,P Covariance matrices in original trait units.
#' @param units One of `"trait"`, `"genetic_sd"`, or `"phenotypic_sd"`.
#'
#' @return The vector expressed in the requested units.
#' @noRd
.dgr_gain_from_trait_units <- function(d, G, P, units) {
  switch(units,
    trait = d,
    genetic_sd = d / sqrt(diag(G)),
    phenotypic_sd = d / sqrt(diag(P))
  )
}

#' Orientation matrix placing every trait in a larger-is-better space
#'
#' The rest of the package accepts `lower_is_better`, and requiring the
#' objective-setting functions alone to be handed pre-oriented matrices is an
#' easy way to obtain a confidently wrong answer.
#'
#' @param trait_cols Trait names.
#' @param lower_is_better Traits for which smaller values are favourable.
#'
#' @return A named vector of 1 and -1.
#' @noRd
.dgr_direction <- function(trait_cols, lower_is_better) {
  direction <- rep(1, length(trait_cols))
  names(direction) <- trait_cols
  if (length(lower_is_better)) {
    unknown <- setdiff(lower_is_better, trait_cols)
    if (length(unknown)) {
      stop("Unknown lower_is_better traits: ",
        paste(unknown, collapse = ", "),
        call. = FALSE
      )
    }
    direction[lower_is_better] <- -1
  }
  direction
}

#' Rotate a covariance matrix into the favourable-direction space
#'
#' @param M Covariance matrix.
#' @param direction Named vector of 1 and -1.
#'
#' @return The oriented covariance matrix.
#' @noRd
.dgr_orient_covariance <- function(M, direction) {
  D <- diag(direction, nrow = length(direction))
  oriented <- D %*% M %*% D
  dimnames(oriented) <- dimnames(M)
  oriented
}

#' Define acceptable intervals for desired genetic gains
#'
#' Breeders are often more confident about an acceptable range than about one
#' exact desired-gain vector. This function records those ranges in familiar
#' raw trait directions and converts them to the favourable-direction
#' convention used by DesiredGainR. Thus a requested reduction of 0.25 to 1.0
#' genetic standard deviations in disease severity can be entered as
#' `lower = -1` and `upper = -0.25` while declaring the disease trait in
#' `lower_is_better`.
#'
#' The resulting intervals constrain the *relative desired-gain directions*
#' searched by [optimize_desired_gains()]. They are not promises that every
#' value in the interval is attainable. Use [gain_feasibility()] to test an
#' absolute point and multi-cycle simulation to compare admissible directions.
#'
#' @param lower,upper Named numeric vectors giving the lower and upper bounds
#'   of acceptable genetic change in the raw trait direction. Both vectors
#'   must contain the same traits. For a trait in `lower_is_better`, reductions
#'   are normally negative; for other traits, improvements are normally
#'   positive.
#' @param lower_is_better Traits for which a reduction is favourable.
#' @param gain_units Units of the bounds: genetic standard deviations,
#'   phenotypic standard deviations, or original trait units. Genetic standard
#'   deviations are recommended for eliciting intervals across traits.
#' @param horizon_cycles Optional number of selection cycles over which the
#'   bounds are intended to be achieved. Recording the horizon prevents a
#'   one-cycle target from being silently evaluated as a five-cycle target.
#'
#' @return A `desiredgainr_gain_intervals` data frame containing favourable-
#'   direction lower and upper bounds, with the raw bounds and units retained
#'   as attributes.
#'
#' @examples
#' intervals <- define_desired_gain_intervals(
#'   lower = c(Yield = 0.5, Disease = -1.0, Quality = 0.1),
#'   upper = c(Yield = 1.5, Disease = -0.25, Quality = 0.8),
#'   lower_is_better = "Disease",
#'   gain_units = "genetic_sd"
#' )
#' intervals
#'
#' @export
define_desired_gain_intervals <- function(
  lower,
  upper,
  lower_is_better = NULL,
  gain_units = c("genetic_sd", "phenotypic_sd", "trait"),
  horizon_cycles = NULL
) {
  gain_units <- match.arg(gain_units)
  if (!is.null(horizon_cycles)) {
    horizon_cycles <- .dgr_positive_integer(horizon_cycles, "horizon_cycles")
  }
  validate_bound <- function(x, name) {
    if (!is.numeric(x) || !length(x) || is.null(names(x)) ||
      any(!nzchar(names(x))) || anyDuplicated(names(x)) ||
      any(!is.finite(x))) {
      stop(name, " must be a finite, uniquely named numeric vector.",
        call. = FALSE
      )
    }
    x
  }
  lower <- validate_bound(lower, "lower")
  upper <- validate_bound(upper, "upper")
  if (!setequal(names(lower), names(upper))) {
    stop("lower and upper must contain the same named traits.",
      call. = FALSE
    )
  }
  upper <- upper[names(lower)]
  if (any(lower > upper)) {
    stop("Every raw lower bound must be less than or equal to its upper bound.",
      call. = FALSE
    )
  }
  unknown <- setdiff(lower_is_better, names(lower))
  if (length(unknown)) {
    stop("Unknown lower_is_better traits: ", paste(unknown, collapse = ", "),
      call. = FALSE
    )
  }
  orientation <- stats::setNames(rep(1, length(lower)), names(lower))
  orientation[lower_is_better] <- -1
  endpoint_one <- lower * orientation
  endpoint_two <- upper * orientation
  result <- data.frame(
    trait = names(lower),
    lower = pmin(endpoint_one, endpoint_two),
    upper = pmax(endpoint_one, endpoint_two),
    row.names = names(lower),
    stringsAsFactors = FALSE
  )
  attr(result, "raw_lower") <- lower
  attr(result, "raw_upper") <- upper
  attr(result, "lower_is_better") <- lower_is_better
  attr(result, "gain_units") <- gain_units
  attr(result, "horizon_cycles") <- horizon_cycles
  class(result) <- c("desiredgainr_gain_intervals", "data.frame")
  result
}

#' @export
print.desiredgainr_gain_intervals <- function(x, ...) {
  cat("<desiredgainr_gain_intervals>\n")
  cat("  Units:", attr(x, "gain_units"), "\n")
  if (!is.null(attr(x, "horizon_cycles"))) {
    cat("  Planning horizon:", attr(x, "horizon_cycles"), "cycles\n")
  }
  cat("  Bounds shown as favourable genetic change (larger is better):\n")
  print.data.frame(unclass(x)[c("trait", "lower", "upper")],
    row.names = FALSE
  )
  lower_is_better <- attr(x, "lower_is_better")
  if (length(lower_is_better)) {
    cat(
      "  Lower-is-better traits were entered as raw reductions:",
      paste(lower_is_better, collapse = ", "), "\n"
    )
  }
  invisible(x)
}

#' Translate desired gains into the economic weights they imply
#'
#' The Smith-Hazel economic index and the Pesek-Baker desired-gain index are
#' two parameterisations of one linear index. Setting
#' \eqn{b = P^{-1}Gw} equal to \eqn{b = P^{-1}G(GP^{-1}G)^{-1}d} gives
#' \deqn{w = G^{-1}PG^{-1}d.}
#' Hence every desired-gain vector corresponds to exactly one economic-weight
#' vector, and [implied_desired_gains()] performs the reverse translation.
#'
#' This is the most direct answer to a breeder who cannot state economic
#' weights but can state relative desired gains, and, more usefully, to a
#' breeder who proposes weights and needs to see what response those weights
#' actually request.
#'
#' @section Negative implied weights are expected, not an error:
#' An implied weight can be negative for a trait the breeder wants to improve.
#' This is correct and interpretable. Where a trait receives a large favourable
#' correlated response from the rest of the objective, achieving only the
#' modest gain requested for it requires the index to hold it back, so its
#' implied weight turns negative.
#'
#' When `lower_is_better` is supplied the weights are returned in the
#' favourable-direction space, in which larger is better for every trait.
#' A negative weight therefore never means that a trait should decrease; the
#' direction has already been applied.
#'
#' @section Weight magnitudes are not comparable across traits:
#' An economic weight carries inverse trait units, so a trait measured on a
#' small scale receives a numerically large weight for the same breeding
#' emphasis. Comparing the raw values across traits therefore misleads: the
#' largest number in the vector can correspond to one of the smallest effects.
#'
#' Multiply each weight by that trait's genetic standard deviation before
#' comparing them, which expresses every trait per unit of the genetic
#' variation actually available. [effective_weights()] performs the equivalent
#' calculation for index coefficients.
#'
#' Consequently these coefficients should not be read as statements of
#' biological or economic importance. Covarrubias-Pazaran (2021) puts the point
#' plainly: the weights lack meaning, especially under strong genetic
#' correlations, and the desired response is the only decision of interest. Use
#' the implied weights to drive an index or to compare objectives, and use
#' [effective_weights()] when the question is which traits the index is really
#' acting on.
#'
#' @references
#' Covarrubias-Pazaran G (2021). *Practical implementation of selection
#' indices.* CGIAR Excellence in Breeding.
#'
#' @param desired_gains Named numeric vector of desired gains. When
#'   `lower_is_better` is supplied these are magnitudes in the
#'   favourable-direction space, so a positive value always means improvement.
#' @param G Genetic variance-covariance matrix, named by trait, in original
#'   trait units.
#' @param P Phenotypic variance-covariance matrix, named by trait, in original
#'   trait units.
#' @param lower_is_better Traits for which smaller original values are
#'   favourable. Supplying this orients `G` and `P` internally, so that
#'   `desired_gains` and the returned weights are both expressed as
#'   improvements. Omitting it means every quantity is interpreted in the raw
#'   trait direction, where improving a trait such as disease severity requires
#'   a negative desired gain.
#' @param gain_units Units of `desired_gains`: `"trait"` for original trait
#'   units, `"genetic_sd"` for genetic standard deviations, or
#'   `"phenotypic_sd"` for phenotypic standard deviations.
#'
#' @return A named numeric vector of implied economic weights, carrying a
#'   `provenance` attribute.
#'
#' @examples
#' traits <- c("yield", "disease")
#' G <- matrix(c(0.60, -0.15, -0.15, 0.40), 2, dimnames = list(traits, traits))
#' P <- matrix(c(1.10, -0.20, -0.20, 0.90), 2, dimnames = list(traits, traits))
#'
#' # Disease severity should fall, so it is declared rather than signed by hand.
#' implied_economic_weights(
#'   c(yield = 0.5, disease = 0.3), G, P,
#'   lower_is_better = "disease"
#' )
#'
#' @seealso [implied_desired_gains()], [gain_feasibility()]
#' @export
implied_economic_weights <- function(
  desired_gains,
  G,
  P,
  lower_is_better = NULL,
  gain_units = c("trait", "genetic_sd", "phenotypic_sd")
) {
  gain_units <- match.arg(gain_units)
  trait_cols <- names(desired_gains)
  if (is.null(trait_cols) || anyDuplicated(trait_cols)) {
    stop("desired_gains must be a numeric vector with unique trait names.",
      call. = FALSE
    )
  }
  G <- .dgr_covariance(G, trait_cols, "G")
  P <- .dgr_covariance(P, trait_cols, "P")
  direction <- .dgr_direction(trait_cols, lower_is_better)
  G <- .dgr_orient_covariance(G, direction)
  P <- .dgr_orient_covariance(P, direction)
  d <- .dgr_named_vector(desired_gains, trait_cols, "desired_gains")
  d <- .dgr_gain_to_trait_units(d, G, P, gain_units)

  G_inverse <- .dgr_inverse(G, "G")$inverse
  weights <- as.numeric(G_inverse %*% P %*% G_inverse %*% d)
  names(weights) <- trait_cols
  attr(weights, "provenance") <- paste(
    "Implied by the supplied desired gains through w = G^-1 P G^-1 d;",
    "not an independently estimated economic value.",
    if (any(direction < 0)) {
      paste(
        "Expressed in the favourable-direction space, so a positive weight",
        "always means the trait matters."
      )
    } else {
      "Expressed in the raw trait direction."
    }
  )
  weights
}

#' Translate economic weights into the desired gains they imply
#'
#' Inverse of [implied_economic_weights()]. For economic weights \eqn{w},
#' the Smith-Hazel index requests the response
#' \deqn{d = GP^{-1}Gw,}
#' up to the scalar set by selection intensity. Reading these implied gains is
#' the quickest way for a breeder to discover that a proposed weight vector
#' asks for something other than what was intended.
#'
#' @param economic_weights Named numeric vector of economic weights, in the
#'   same space as `lower_is_better` implies.
#' @param G Genetic variance-covariance matrix, named by trait, in original
#'   trait units.
#' @param P Phenotypic variance-covariance matrix, named by trait, in original
#'   trait units.
#' @param lower_is_better Traits for which smaller original values are
#'   favourable. See [implied_economic_weights()].
#' @param gain_units Units in which the implied gains are returned. Set this to
#'   match the units in which the weights were derived, so that a round trip
#'   through [implied_economic_weights()] returns the vector it started from.
#'
#' @return A named numeric vector of implied desired gains, carrying a
#'   `provenance` attribute recording the units.
#'
#' @examples
#' traits <- c("yield", "disease")
#' G <- matrix(c(0.60, -0.15, -0.15, 0.40), 2, dimnames = list(traits, traits))
#' P <- matrix(c(1.10, -0.20, -0.20, 0.90), 2, dimnames = list(traits, traits))
#'
#' # The translation is exactly invertible when the units match.
#' d <- c(yield = 0.5, disease = 0.3)
#' w <- implied_economic_weights(
#'   d, G, P,
#'   lower_is_better = "disease", gain_units = "genetic_sd"
#' )
#' implied_desired_gains(
#'   w, G, P,
#'   lower_is_better = "disease", gain_units = "genetic_sd"
#' )
#'
#' @seealso [implied_economic_weights()]
#' @export
implied_desired_gains <- function(
  economic_weights,
  G,
  P,
  lower_is_better = NULL,
  gain_units = c("trait", "genetic_sd", "phenotypic_sd")
) {
  gain_units <- match.arg(gain_units)
  trait_cols <- names(economic_weights)
  if (is.null(trait_cols) || anyDuplicated(trait_cols)) {
    stop("economic_weights must be a numeric vector with unique trait names.",
      call. = FALSE
    )
  }
  G <- .dgr_covariance(G, trait_cols, "G")
  P <- .dgr_covariance(P, trait_cols, "P")
  direction <- .dgr_direction(trait_cols, lower_is_better)
  G <- .dgr_orient_covariance(G, direction)
  P <- .dgr_orient_covariance(P, direction)
  w <- .dgr_named_vector(economic_weights, trait_cols, "economic_weights")

  P_inverse <- .dgr_inverse(P, "P")$inverse
  gains <- as.numeric(G %*% P_inverse %*% G %*% w)
  names(gains) <- trait_cols
  gains <- .dgr_gain_from_trait_units(gains, G, P, gain_units)
  attr(gains, "provenance") <- paste(
    "Implied by the supplied economic weights through d = G P^-1 G w,",
    sprintf("expressed in %s units", gain_units),
    "and up to the scalar set by selection intensity."
  )
  gains
}

#' Test whether a desired-gain vector is attainable
#'
#' For any linear index the response vector is
#' \eqn{R = iGb/\sqrt{b^\mathsf{T}Pb}}, so the attainable responses lie on the
#' ellipsoid \eqn{R^\mathsf{T}G^{-1}PG^{-1}R = i^{2}}. Therefore a target
#' response \eqn{d} requires exactly
#' \deqn{i_{\mathrm{required}} = \sqrt{d^\mathsf{T}G^{-1}PG^{-1}d},}
#' and the selected proportion delivering that intensity follows from the
#' normal-truncation relationship.
#'
#' Two consequences deserve emphasis. First, the classical Pesek-Baker index
#' honours only the *direction* of the desired-gain vector; scaling every
#' element by a constant leaves the index unchanged, because the attainable
#' magnitude is fixed by selection intensity. Second, an antagonistic
#' correlation structure can make a modest-looking target unreachable at any
#' practical intensity.
#'
#' This function replaces the external feasibility check that
#' Covarrubias-Pazaran (2021) recommends performing in separate software, and
#' should be run before an optimisation is attempted rather than after it
#' fails.
#'
#' @param desired_gains Named numeric vector of target responses. When
#'   `lower_is_better` is supplied these are magnitudes in the
#'   favourable-direction space, so a positive value always means improvement.
#' @param G Genetic variance-covariance matrix, named by trait, in original
#'   trait units.
#' @param P Phenotypic variance-covariance matrix, named by trait, in original
#'   trait units.
#' @param lower_is_better Traits for which smaller original values are
#'   favourable. See [implied_economic_weights()].
#' @param n_candidates Number of selection candidates available.
#' @param n_select Number of candidates to be selected. Supply either this or
#'   `selection_proportion`.
#' @param selection_proportion Proportion to be selected, in `(0, 1]`.
#' @param gain_units Units of `desired_gains`, as in
#'   [implied_economic_weights()].
#'
#' @return An object of class `desiredgainr_feasibility` reporting: (i) the
#'   required selection intensity and the proportion that delivers it, (ii)
#'   whether the target is attainable in a population of `n_candidates`, (iii)
#'   the attainable response in the requested direction at the planned
#'   intensity, and (iv) the per-trait shortfall.
#'
#'   The attainable response and the shortfall are returned twice:
#'   `attainable_response` and `shortfall` are in original trait units, whereas
#'   `attainable_response_input_units` and `shortfall_input_units` are in
#'   whatever `gain_units` the target was stated in. Comparing a request made
#'   in standard deviations against an answer given in trait units is an easy
#'   way to misread the result.
#'
#' @examples
#' traits <- c("yield", "disease")
#' G <- matrix(c(0.60, -0.15, -0.15, 0.40), 2, dimnames = list(traits, traits))
#' P <- matrix(c(1.10, -0.20, -0.20, 0.90), 2, dimnames = list(traits, traits))
#' gain_feasibility(
#'   desired_gains = c(yield = 1.0, disease = -0.8),
#'   G = G, P = P, n_candidates = 500, n_select = 50,
#'   gain_units = "genetic_sd"
#' )
#'
#' @references
#' Covarrubias-Pazaran G (2021). *Practical implementation of selection
#' indices.* CGIAR Excellence in Breeding.
#'
#' @export
gain_feasibility <- function(
  desired_gains,
  G,
  P,
  n_candidates,
  n_select = NULL,
  selection_proportion = NULL,
  lower_is_better = NULL,
  gain_units = c("trait", "genetic_sd", "phenotypic_sd")
) {
  gain_units <- match.arg(gain_units)
  trait_cols <- names(desired_gains)
  if (is.null(trait_cols) || anyDuplicated(trait_cols)) {
    stop("desired_gains must be a numeric vector with unique trait names.",
      call. = FALSE
    )
  }
  if (!is.numeric(n_candidates) || length(n_candidates) != 1L ||
    !is.finite(n_candidates) || n_candidates < 2) {
    stop("n_candidates must be a single number of at least 2.", call. = FALSE)
  }
  if (is.null(n_select) == is.null(selection_proportion)) {
    stop("Supply exactly one of n_select and selection_proportion.",
      call. = FALSE
    )
  }
  planned_proportion <- if (is.null(selection_proportion)) {
    n_select / n_candidates
  } else {
    selection_proportion
  }
  if (!is.finite(planned_proportion) || planned_proportion <= 0 ||
    planned_proportion > 1) {
    stop("The planned selected proportion must lie in (0, 1].", call. = FALSE)
  }

  G <- .dgr_covariance(G, trait_cols, "G")
  P <- .dgr_covariance(P, trait_cols, "P")
  direction <- .dgr_direction(trait_cols, lower_is_better)
  G <- .dgr_orient_covariance(G, direction)
  P <- .dgr_orient_covariance(P, direction)
  d_input <- .dgr_named_vector(desired_gains, trait_cols, "desired_gains")
  d <- .dgr_gain_to_trait_units(d_input, G, P, gain_units)

  G_inverse <- .dgr_inverse(G, "G")$inverse
  quadratic_form <- as.numeric(crossprod(d, G_inverse %*% P %*% G_inverse %*% d))
  if (quadratic_form <= 0) {
    stop("The desired-gain vector has no positive required intensity.",
      call. = FALSE
    )
  }
  required_intensity <- sqrt(quadratic_form)
  required_proportion <- .dgr_proportion_for_intensity(required_intensity)
  planned_intensity <- .dgr_intensity(planned_proportion)

  # The strongest intensity a finite population can deliver is the intensity
  # at which exactly one candidate is retained.
  maximum_intensity <- .dgr_intensity(1 / n_candidates)
  attainable <- as.numeric(planned_intensity * d / required_intensity)
  names(attainable) <- trait_cols
  # The caller stated the target in `gain_units`, so the attainable response is
  # reported in those units as well as in original trait units. Returning only
  # trait units invites the reader to compare a request in standard deviations
  # against an answer in tonnes per hectare.
  attainable_input <- .dgr_gain_from_trait_units(attainable, G, P, gain_units)
  names(attainable_input) <- trait_cols

  result <- list(
    desired_gains = d,
    desired_gains_input = d_input,
    gain_units = gain_units,
    required_intensity = required_intensity,
    required_proportion = required_proportion,
    # A required proportion below one candidate cannot be realised by any
    # integer selection, so reporting a rounded-up count of one would overstate
    # what the population can deliver.
    required_n_select = if (is.finite(required_proportion) &&
      required_proportion * n_candidates >= 1) {
      ceiling(required_proportion * n_candidates)
    } else {
      NA_real_
    },
    planned_proportion = planned_proportion,
    planned_intensity = planned_intensity,
    maximum_intensity = maximum_intensity,
    n_candidates = n_candidates,
    feasible_at_planned_intensity = planned_intensity >= required_intensity,
    feasible_in_this_population = maximum_intensity >= required_intensity,
    attainable_response = attainable,
    attainable_response_input_units = attainable_input,
    shortfall = d - attainable,
    shortfall_input_units = d_input - attainable_input,
    attainable_fraction = planned_intensity / required_intensity,
    interpretation = paste(
      "Pesek-Baker honours the direction of the desired-gain vector, not its",
      "magnitude. The attainable response is the requested direction scaled",
      "to the intensity actually applied."
    )
  )
  class(result) <- c("desiredgainr_feasibility", "list")
  result
}

#' @export
print.desiredgainr_feasibility <- function(x, ...) {
  cat("<desiredgainr_feasibility>\n")
  cat(sprintf("  Required selection intensity: %.4f\n", x$required_intensity))
  if (is.finite(x$required_proportion) && is.finite(x$required_n_select)) {
    cat(sprintf(
      "  Requires selecting the top %.3f%% (%s of %s candidates)\n",
      100 * x$required_proportion,
      format(x$required_n_select), format(x$n_candidates)
    ))
  } else if (is.finite(x$required_proportion)) {
    cat(sprintf(
      "  Requires the top %.4f%%, which is fewer than one of %s candidates\n",
      100 * x$required_proportion, format(x$n_candidates)
    ))
  } else {
    cat("  Required proportion: unattainable under normal truncation\n")
  }
  cat(sprintf(
    "  Planned intensity: %.4f (top %.1f%%)\n",
    x$planned_intensity, 100 * x$planned_proportion
  ))
  cat(sprintf(
    "  Feasible at planned intensity: %s\n",
    if (isTRUE(x$feasible_at_planned_intensity)) "yes" else "no"
  ))
  cat(sprintf(
    "  Feasible anywhere in this population: %s\n",
    if (isTRUE(x$feasible_in_this_population)) "yes" else "no"
  ))
  cat(sprintf(
    "  Attainable fraction of the requested gain: %.1f%%\n",
    100 * x$attainable_fraction
  ))
  invisible(x)
}

#' Recover the weights implied by past selection decisions
#'
#' Where a programme has been selecting for years without a formal index, the
#' selection differentials it has achieved already encode its objective.
#' Covarrubias-Pazaran recovers the corresponding linear rule as
#' \deqn{b = P^{-1}s,}
#' where \eqn{s} contains the differentials between the selected group and the
#' full population. This is how a selection index was introduced at both the
#' International Maize and Wheat Improvement Center (CIMMYT) maize programmes
#' and the International Rice Research Institute (IRRI) programmes.
#'
#' The recovered weights answer the question "what linear rule best
#' approximates the decisions already made?". They do not answer "what should
#' this programme select for next?". Therefore treat the result as a starting
#' point, inspect it with [weight_sensitivity()], and adjust the emphasis
#' deliberately before adopting it.
#'
#' @param selected_values Numeric matrix or data frame of trait values for the
#'   candidates that were historically selected, with one column per trait.
#' @param population_values Numeric matrix or data frame of trait values for
#'   the full population from which they were selected.
#' @param trait_cols Character vector naming and ordering the trait columns.
#' @param P Optional phenotypic variance-covariance matrix. When `NULL`, it is
#'   estimated from `population_values` and reported as such.
#'
#' @return An object of class `desiredgainr_retrospective` containing the
#'   recovered coefficients, the achieved selection differentials in both
#'   original and standard-deviation units, and the provenance of `P`.
#'
#' @references
#' Covarrubias-Pazaran G. *Bringing a selection index into the CIMMYT-Maize
#' programs* and *Bringing a selection index into the IRRI programs.* CGIAR
#' Excellence in Breeding.
#'
#' @export
retrospective_weights <- function(
  selected_values,
  population_values,
  trait_cols,
  P = NULL
) {
  if (!is.character(trait_cols) || !length(trait_cols) ||
    anyDuplicated(trait_cols)) {
    stop("trait_cols must contain unique trait names.", call. = FALSE)
  }
  extract <- function(x, name) {
    x <- as.matrix(as.data.frame(x)[, trait_cols, drop = FALSE])
    storage.mode(x) <- "double"
    if (!nrow(x) || any(!is.finite(x))) {
      stop(name, " must contain finite values for every trait.", call. = FALSE)
    }
    x
  }
  selected <- extract(selected_values, "selected_values")
  population <- extract(population_values, "population_values")
  if (nrow(population) < 2L) {
    stop("population_values must contain at least two candidates.",
      call. = FALSE
    )
  }

  P_was_estimated <- is.null(P)
  if (P_was_estimated) {
    P <- stats::cov(population)
    dimnames(P) <- list(trait_cols, trait_cols)
    P_source <- "estimated from population_values"
  } else {
    P <- .dgr_covariance(P, trait_cols, "P")
    P_source <- "user supplied"
  }

  differential <- colMeans(selected) - colMeans(population)
  names(differential) <- trait_cols
  population_sd <- apply(population, 2L, stats::sd)
  population_sd[!is.finite(population_sd) | population_sd == 0] <- NA_real_

  P_inverse <- .dgr_inverse(P, "P")$inverse
  coefficients <- as.numeric(P_inverse %*% differential)
  names(coefficients) <- trait_cols
  # A coefficient carries inverse trait units, so a trait measured on a small
  # scale receives a large number for the same emphasis. Multiplying by the
  # trait standard deviation puts every trait on one comparable footing, which
  # is the form to read when deciding which traits the past decisions actually
  # favoured.
  coefficients_per_sd <- coefficients * population_sd
  names(coefficients_per_sd) <- trait_cols

  result <- list(
    coefficients = coefficients,
    coefficients_per_sd = coefficients_per_sd,
    selection_differential = differential,
    selection_differential_sd = differential / population_sd,
    n_selected = nrow(selected),
    n_population = nrow(population),
    selected_proportion = nrow(selected) / nrow(population),
    P = P,
    P_source = P_source,
    caution = paste(
      "These coefficients reproduce past decisions, including any bias in",
      "them. Use them to initialise an index, then adjust the emphasis",
      "deliberately against the target product profile."
    )
  )
  class(result) <- c("desiredgainr_retrospective", "list")
  result
}

#' @export
print.desiredgainr_retrospective <- function(x, ...) {
  cat("<desiredgainr_retrospective>\n")
  cat(sprintf(
    "  %d selected from %d candidates (%.1f%%)\n",
    x$n_selected, x$n_population, 100 * x$selected_proportion
  ))
  cat("  P:", x$P_source, "\n")
  cat("  Recovered coefficients, per trait standard deviation:\n")
  print(round(x$coefficients_per_sd, 4L))
  cat(
    "  (Raw coefficients are in $coefficients; they carry inverse trait",
    "units and are not comparable across traits.)\n"
  )
  invisible(x)
}

#' Report the effective contribution of each trait to an index
#'
#' Index coefficients are not comparable across traits measured on different
#' scales. Crosbie et al. (1980) observed that assigning equal weights to
#' unstandardised traits places most of the selection pressure on whichever
#' trait carries the largest genetic variance. Hence the interpretable
#' quantity is not \eqn{b_j} but \eqn{b_j\sigma_{g_j}}, the contribution of
#' trait \eqn{j} per unit of the genetic variation actually available.
#'
#' Use this diagnostic whenever traits are analysed on their original scales,
#' and treat a single trait dominating the effective weights as a signal to
#' standardise before proceeding.
#'
#' @param coefficients Named numeric vector of index coefficients.
#' @param G Genetic variance-covariance matrix, named by trait.
#' @param P Optional phenotypic variance-covariance matrix. When supplied, the
#'   phenotypic effective weights are reported alongside the genetic ones.
#' @param warn Whether to warn when one trait exceeds `dominance_threshold`.
#'   The default is `FALSE`, because a concentrated effective weight can arise
#'   from a deliberately asymmetric objective as readily as from mismatched
#'   trait scales, and warning on the former would be noise.
#' @param dominance_threshold Share of the total absolute effective weight
#'   above which a warning is issued for a single trait. The default adapts to
#'   the number of traits, since with two traits a share of 0.75 is merely a
#'   three-to-one emphasis whereas with eight traits it is extreme.
#'
#' @return A `data.table` with one row per trait giving the coefficient, the
#'   genetic and phenotypic effective weights, and their shares. Shares are
#'   `NA` when the effective weights sum to zero.
#'
#' @references
#' Crosbie TM, Mock JJ, Smith OS (1980), as discussed in Guimaraes PHR et al.
#' (2021) *Euphytica* 217:95.
#'
#' @export
effective_weights <- function(
  coefficients,
  G,
  P = NULL,
  warn = FALSE,
  dominance_threshold = NULL
) {
  trait_cols <- names(coefficients)
  if (is.null(trait_cols) || anyDuplicated(trait_cols)) {
    stop("coefficients must be a numeric vector with unique trait names.",
      call. = FALSE
    )
  }
  G <- .dgr_covariance(G, trait_cols, "G")
  b <- .dgr_named_vector(coefficients, trait_cols, "coefficients")
  p <- length(trait_cols)
  if (is.null(dominance_threshold)) {
    # Halfway between an equal share and total dominance, so the rule tightens
    # as traits are added rather than firing on any asymmetry.
    dominance_threshold <- 0.5 * (1 / p) + 0.5
  }

  share_of <- function(x) {
    total <- sum(abs(x))
    if (!is.finite(total) || total <= 0) {
      return(rep(NA_real_, length(x)))
    }
    abs(x) / total
  }
  genetic <- b * sqrt(diag(G))
  genetic_share <- share_of(genetic)
  phenotypic <- rep(NA_real_, p)
  phenotypic_share <- phenotypic
  if (!is.null(P)) {
    P <- .dgr_covariance(P, trait_cols, "P")
    phenotypic <- b * sqrt(diag(P))
    phenotypic_share <- share_of(phenotypic)
  }

  if (isTRUE(warn) && any(is.finite(genetic_share))) {
    dominant <- which.max(genetic_share)
    if (length(dominant) == 1L &&
      is.finite(genetic_share[dominant]) &&
      genetic_share[dominant] > dominance_threshold) {
      warning(
        sprintf(
          paste(
            "Trait '%s' carries %.0f%% of the total effective weight, against",
            "an equal share of %.0f%%. Equal desired gains on unstandardised",
            "traits concentrate selection on the trait with the largest",
            "genetic variance; consider standardising."
          ),
          trait_cols[dominant], 100 * genetic_share[dominant], 100 / p
        ),
        call. = FALSE
      )
    }
  }

  data.table::data.table(
    Trait = trait_cols,
    Coefficient = as.numeric(b),
    Genetic_effective_weight = as.numeric(genetic),
    Genetic_share = as.numeric(genetic_share),
    Phenotypic_effective_weight = as.numeric(phenotypic),
    Phenotypic_share = as.numeric(phenotypic_share)
  )
}

#' Assess how sensitive a selection decision is to the stated objective
#'
#' Economic weights and desired gains are estimates, and Guimaraes et al.
#' (2021) demonstrated that a poorly chosen weight vector can eliminate the
#' gain an index would otherwise deliver. However, a decision is often far
#' less sensitive to the objective than breeders fear. This function
#' quantifies which of the two situations applies.
#'
#' Weight vectors are perturbed multiplicatively on the log scale, the index is
#' rebuilt for each draw, and the resulting selected sets are compared with the
#' set obtained from the stated objective.
#'
#' @param economic_weights Named numeric vector of economic weights.
#' @param values Numeric matrix or data frame of candidate trait values.
#' @param G Genetic variance-covariance matrix, named by trait.
#' @param P Phenotypic variance-covariance matrix, named by trait.
#' @param n_select Number of candidates selected.
#' @param trait_cols Character vector naming and ordering the trait columns.
#'   Defaults to the names of `economic_weights`.
#' @param relative_sd Standard deviation of the log-scale perturbation applied
#'   to each weight.
#' @param n_draws Number of perturbed weight vectors evaluated.
#' @param agreement_threshold Selected-set overlap above which a draw is
#'   counted as reproducing the original decision.
#' @param seed Random seed. The caller's random number generator state is
#'   restored on exit.
#'
#' @return An object of class `desiredgainr_sensitivity` giving the
#'   distribution of selected-set agreement, the rank correlation with the
#'   original index, the stability proportion, and the per-trait weight ratios
#'   that most strongly drive disagreement.
#'
#' @section Influence is not the same as contribution:
#' `weight_influence` measures how strongly perturbing each weight disturbs the
#' selected set. It is not the share of the index that the trait contributes,
#' which is what [effective_weights()] reports, and the two routinely disagree.
#'
#' A trait carrying only a small part of the index can still be the largest
#' lever on the decision, because it moves the index in a direction the
#' remaining traits do not already cover, and the candidates near the selection
#' threshold are therefore reordered by it. A trait that dominates the index
#' may conversely be a weak lever, because the ranking already follows it and
#' scaling it further changes little. Read the two diagnostics together: the
#' first says which weights are worth arguing about, the second says which
#' traits the index is acting on.
#'
#' @export
weight_sensitivity <- function(
  economic_weights,
  values,
  G,
  P,
  n_select,
  trait_cols = names(economic_weights),
  relative_sd = 0.25,
  n_draws = 200L,
  agreement_threshold = 0.9,
  seed = 42L
) {
  if (is.null(trait_cols) || anyDuplicated(trait_cols)) {
    stop("trait_cols must contain unique trait names.", call. = FALSE)
  }
  G <- .dgr_covariance(G, trait_cols, "G")
  P <- .dgr_covariance(P, trait_cols, "P")
  w <- .dgr_named_vector(economic_weights, trait_cols, "economic_weights")
  X <- as.matrix(as.data.frame(values)[, trait_cols, drop = FALSE])
  storage.mode(X) <- "double"
  if (any(!is.finite(X))) {
    stop("values must contain only finite trait values.", call. = FALSE)
  }
  n_select <- .dgr_positive_integer(n_select, "n_select")
  if (n_select > nrow(X)) {
    stop("n_select cannot exceed the number of candidates.", call. = FALSE)
  }
  n_draws <- .dgr_positive_integer(n_draws, "n_draws")
  seed <- .dgr_seed(seed)

  if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    stats::runif(1L)
  }
  entry_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit(assign(".Random.seed", entry_seed, envir = globalenv()), add = TRUE)
  set.seed(seed)

  P_inverse <- .dgr_inverse(P, "P")$inverse
  index_scores <- function(weights) {
    as.numeric(X %*% (P_inverse %*% G %*% weights))
  }
  top_set <- function(score) {
    order(-score, seq_along(score))[seq_len(n_select)]
  }

  reference_score <- index_scores(w)
  reference_set <- top_set(reference_score)

  agreement <- numeric(n_draws)
  rank_correlation <- numeric(n_draws)
  draws <- matrix(NA_real_, n_draws, length(trait_cols),
    dimnames = list(NULL, trait_cols)
  )
  for (draw in seq_len(n_draws)) {
    perturbed <- w * exp(stats::rnorm(length(w), sd = relative_sd))
    draws[draw, ] <- perturbed
    score <- index_scores(perturbed)
    candidate_set <- top_set(score)
    agreement[draw] <- length(intersect(reference_set, candidate_set)) /
      length(union(reference_set, candidate_set))
    rank_correlation[draw] <- suppressWarnings(
      stats::cor(reference_score, score, method = "spearman")
    )
  }

  stable <- agreement >= agreement_threshold
  ratio_influence <- vapply(seq_along(trait_cols), function(j) {
    log_ratio <- log(abs(draws[, j])) - log(abs(w[j]))
    suppressWarnings(abs(stats::cor(log_ratio, agreement)))
  }, numeric(1))
  names(ratio_influence) <- trait_cols

  result <- list(
    reference_weights = w,
    n_draws = n_draws,
    relative_sd = relative_sd,
    n_select = n_select,
    agreement = agreement,
    agreement_quantiles = stats::quantile(
      agreement, c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1)
    ),
    rank_correlation_median = stats::median(rank_correlation, na.rm = TRUE),
    stability_proportion = mean(stable),
    agreement_threshold = agreement_threshold,
    weight_influence = sort(ratio_influence, decreasing = TRUE),
    interpretation = paste(
      "A high stability proportion means the selected set is robust to the",
      "stated weights, and further effort refining them is unlikely to change",
      "the decision. A low proportion means the objective, not the index,",
      "is the binding uncertainty."
    )
  )
  class(result) <- c("desiredgainr_sensitivity", "list")
  result
}

#' @export
print.desiredgainr_sensitivity <- function(x, ...) {
  cat("<desiredgainr_sensitivity>\n")
  cat(sprintf(
    "  %d draws, log-scale perturbation SD %.2f, selecting %d\n",
    x$n_draws, x$relative_sd, x$n_select
  ))
  cat(sprintf(
    "  Median selected-set agreement: %.3f\n",
    stats::median(x$agreement)
  ))
  cat(sprintf(
    "  Decisions reproducing the original set (agreement >= %.2f): %.1f%%\n",
    x$agreement_threshold, 100 * x$stability_proportion
  ))
  cat(sprintf(
    "  Median rank correlation with the stated objective: %.3f\n",
    x$rank_correlation_median
  ))
  cat("  Weights the decision is most sensitive to:\n")
  print(round(utils::head(x$weight_influence, 3L), 3L))
  cat(
    "  (Sensitivity of the decision, not share of the index; see",
    "effective_weights().)\n"
  )
  invisible(x)
}
