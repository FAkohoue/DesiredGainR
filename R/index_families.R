# Classical linear selection index families, plus the two non-index selection
# strategies used as comparators throughout the selection-index literature.

#' Rank-based index score, favourable direction assumed
#'
#' @param X Candidate-by-trait matrix in favourable-direction space.
#' @param weights Optional non-negative trait weights applied to the ranks.
#'
#' @return Numeric vector of summed ranks, smaller being better.
#' @noRd
.dgr_rank_sum <- function(X, weights = NULL) {
  ranks <- apply(-X, 2L, function(column) rank(column, ties.method = "average"))
  if (is.null(dim(ranks))) ranks <- matrix(ranks, nrow = nrow(X))
  if (is.null(weights)) return(rowSums(ranks))
  as.numeric(ranks %*% weights)
}

#' Construct a classical multi-trait selection index
#'
#' `selection_index()` builds any of the classical linear selection indices
#' from a common interface, so that alternative objectives and alternative
#' index families can be compared on identical data. All traits are first
#' oriented so that larger values are favourable.
#'
#' @details
#' # Available methods
#'
#' \describe{
#'   \item{`"smith_hazel"`}{The optimum economic index, \eqn{b = P^{-1}Ga},
#'     maximising the correlation between the index and the aggregate genotype
#'     \eqn{H = a^\mathsf{T}g}. Requires `economic_weights`.}
#'   \item{`"base"`}{The base index of Brim, Cockerham and Clark,
#'     \eqn{b = a}. It requires neither `P` nor `G` and is included because
#'     Rahimi and Debnath (2023) found it to match the optimum index almost
#'     exactly in their maize data. Where the two agree, the effort spent
#'     estimating covariance matrices has bought nothing, which is itself worth
#'     reporting.}
#'   \item{`"pesek_baker"`}{The original desired-gain index of Pesek and Baker
#'     (1969), \eqn{b = G^{-1}d}, applied to genotypic values. Requires
#'     `desired_gains`.}
#'   \item{`"yamada"`}{The desired-gain index of Yamada, Yokouchi and Nishida
#'     (1975) for phenotypic selection criteria,
#'     \eqn{b = P^{-1}G(GP^{-1}G)^{-1}d}. This is the formulation used by
#'     Joukhadar et al. (2024) and by [run_dgsi()].}
#'   \item{`"mulamba_mock"`}{The rank-sum index of Mulamba and Mock (1978).
#'     Candidates are ranked for each trait in the favourable direction and the
#'     ranks are summed. It requires neither economic weights nor covariance
#'     matrices, and Guimaraes et al. (2021) found it to give the most balanced
#'     multi-trait response of the methods they compared.}
#'   \item{`"independent_culling"`}{Not an index. Candidates must exceed a
#'     threshold for every trait. Included as a comparator because it is what
#'     most programmes actually do.}
#'   \item{`"tandem"`}{Not an index. Candidates are selected sequentially, one
#'     trait at a time. Included as a comparator.}
#' }
#'
#' # Two formulations share the name Pesek-Baker, and they coincide
#'
#' The literature applies the name to two expressions. Pesek and Baker (1969)
#' give \eqn{b = G^{-1}d}, which is `method = "pesek_baker"` here and what
#' Rahimi and Debnath (2023) implement. Yamada et al. (1975) give
#' \eqn{b = P^{-1}G(GP^{-1}G)^{-1}d}, which is `method = "yamada"` and what
#' Joukhadar et al. (2024) and [run_dgsi()] use.
#'
#' For a square invertible \eqn{G} these are the same index, because
#' \deqn{P^{-1}G(GP^{-1}G)^{-1}d = P^{-1}G(G^{-1}PG^{-1})d = G^{-1}d.}
#' Both are therefore retained as documented routes to one estimand rather than
#' as competing methods.
#'
#' Two situations make the distinction real. First, when \eqn{G} is singular or
#' rank deficient, the direct inverse does not exist and only the Yamada route
#' is available. Second, the two routes are not numerically equivalent: one
#' inverts \eqn{G}, the other inverts \eqn{P} and then \eqn{GP^{-1}G}, so an
#' ill-conditioned matrix can send them apart. That is the practical reason to
#' state which route produced a published result, and the reason
#' [matrix_diagnostics()] reports conditioning.
#'
#' @param values Numeric matrix or data frame of candidate trait values, with
#'   one column per trait. Candidate identifiers are taken from `id_col` when
#'   supplied, and from the row names otherwise.
#' @param trait_cols Character vector naming and ordering the trait columns.
#' @param id_col Optional name of a column in `values` holding the candidate
#'   identifiers, matching the argument of the same name in [run_dgsi()] and
#'   [run_qgsi()]. When `NULL`, row names are used. A data frame carrying
#'   identifiers in a column rather than in its row names will otherwise be
#'   labelled by row number, so a warning is issued when that appears to have
#'   happened.
#' @param method Index family. See Details.
#' @param G Genetic variance-covariance matrix, named by trait. Required by
#'   every method except `"base"`, `"mulamba_mock"`, `"independent_culling"`,
#'   and `"tandem"`.
#' @param P Phenotypic variance-covariance matrix, named by trait. Required by
#'   `"smith_hazel"` and `"yamada"`.
#' @param economic_weights Named non-negative economic weights in the
#'   favourable-direction trait space. Required by `"smith_hazel"` and
#'   `"base"`, and optionally used to weight ranks in `"mulamba_mock"`.
#' @param desired_gains Named non-negative desired gains in the
#'   favourable-direction trait space. Required by `"pesek_baker"` and
#'   `"yamada"`.
#' @param lower_is_better Traits for which smaller original values are
#'   favourable.
#' @param center_traits Whether to subtract the trait means before indexing.
#'   Centring does not change the ranking, because it shifts every score by the
#'   same constant, but it does place the index mean at zero and therefore
#'   makes the index coefficient of variation undefined. Set this to `FALSE`
#'   when reproducing published results that report `CV_I` on an index built
#'   from raw trait values.
#' @param scale_traits Whether to divide traits by their population standard
#'   deviations before indexing. Covarrubias-Pazaran (2021) recommends this,
#'   and a desired gain of 1 then means one standard deviation of progress.
#' @param culling_thresholds Named thresholds in the favourable-direction trait
#'   space, required by `"independent_culling"`.
#' @param tandem_order Character vector giving the order in which traits are
#'   selected, required by `"tandem"`.
#' @param n_select Number of candidates selected.
#' @param selection_intensity Optional standardised selection intensity used
#'   for the expected-response calculations. When `NULL`, it is derived from
#'   `n_select` under normal truncation.
#' @param main_trait Trait against which relative efficiency is computed.
#'   Defaults to the first trait.
#'
#' @return An object of class `desiredgainr_index` containing the coefficients,
#'   the candidate scores and ranking, the evaluation criteria from
#'   [evaluate_index()], the effective weights, and the provenance of every
#'   input.
#'
#' @examples
#' set.seed(1)
#' traits <- c("yield", "disease")
#' values <- matrix(
#'   c(stats::rnorm(40), stats::rnorm(40)), ncol = 2,
#'   dimnames = list(paste0("G", 1:40), traits)
#' )
#' G <- matrix(c(0.60, -0.15, -0.15, 0.40), 2, dimnames = list(traits, traits))
#' P <- matrix(c(1.10, -0.20, -0.20, 0.90), 2, dimnames = list(traits, traits))
#'
#' fit <- selection_index(
#'   values, traits, method = "smith_hazel",
#'   G = G, P = P,
#'   economic_weights = c(yield = 1, disease = 0.5),
#'   lower_is_better = "disease", n_select = 4
#' )
#' fit
#'
#' @references
#' Brim CA, Cockerham HW, Clark C (1959). *Agronomy Journal* 51:42-46.
#'
#' Guimaraes PHR, Melo PGS, Cordeiro ACC, Torga PP, Rangel PHN, de Castro AP
#' (2021). *Euphytica* 217:95. \doi{10.1007/s10681-021-02819-7}
#'
#' Joukhadar R, Li Y, Thistlethwaite R, Forrest KL, Tibbits JF, Trethowan R,
#' Hayden MJ (2024). *Frontiers in Plant Science* 15:1337388.
#' \doi{10.3389/fpls.2024.1337388}
#'
#' Mulamba NN, Mock JJ (1978). *Egyptian Journal of Genetics and Cytology*
#' 7:40-51.
#'
#' Pesek J, Baker RJ (1969). *Canadian Journal of Plant Science* 49:803-804.
#' \doi{10.4141/cjps69-137}
#'
#' Rahimi M, Debnath S (2023). *Scientific Reports* 13:18977.
#' \doi{10.1038/s41598-023-46368-6}
#'
#' Smith HF (1936). *Annals of Eugenics* 7:240-250.
#'
#' Yamada Y, Yokouchi K, Nishida A (1975). *Japanese Journal of Genetics*
#' 50:33-41. \doi{10.1266/jjg.50.33}
#'
#' @seealso [evaluate_index()], [gain_feasibility()], [run_dgsi()]
#' @export
selection_index <- function(
    values,
    trait_cols,
    id_col = NULL,
    method = c(
      "smith_hazel", "base", "pesek_baker", "yamada",
      "mulamba_mock", "independent_culling", "tandem"
    ),
    G = NULL,
    P = NULL,
    economic_weights = NULL,
    desired_gains = NULL,
    lower_is_better = NULL,
    center_traits = TRUE,
    scale_traits = TRUE,
    culling_thresholds = NULL,
    tandem_order = NULL,
    n_select = NULL,
    selection_intensity = NULL,
    main_trait = trait_cols[1L]
) {
  method <- match.arg(method)
  if (!is.character(trait_cols) || !length(trait_cols) ||
      anyDuplicated(trait_cols)) {
    stop("trait_cols must contain unique trait names.", call. = FALSE)
  }
  frame <- as.data.frame(values)
  absent <- setdiff(trait_cols, names(frame))
  if (length(absent)) {
    stop("values is missing trait columns: ",
         paste(absent, collapse = ", "), call. = FALSE)
  }
  X_raw <- as.matrix(frame[, trait_cols, drop = FALSE])
  storage.mode(X_raw) <- "double"
  if (!nrow(X_raw) || any(!is.finite(X_raw))) {
    stop("values must contain finite values for every trait.", call. = FALSE)
  }
  if (!is.null(id_col)) {
    if (!id_col %in% names(frame)) {
      stop("id_col '", id_col, "' was not found in values.", call. = FALSE)
    }
    candidate_id <- as.character(frame[[id_col]])
    if (anyNA(candidate_id) || anyDuplicated(candidate_id)) {
      stop("Candidate identifiers in '", id_col,
           "' must be unique and non-missing.", call. = FALSE)
    }
  } else {
    candidate_id <- rownames(X_raw)
    if (is.null(candidate_id)) {
      candidate_id <- paste0("row_", seq_len(nrow(X_raw)))
    }
    # A data frame that carries identifiers in a column rather than in its row
    # names would otherwise be labelled by row number, and the result would
    # name candidates that cannot be traced back to the breeding programme.
    looks_positional <- identical(
      candidate_id, as.character(seq_len(nrow(X_raw)))
    )
    character_columns <- setdiff(
      names(frame)[vapply(frame, function(x) {
        is.character(x) || is.factor(x)
      }, logical(1L))],
      trait_cols
    )
    if (looks_positional && length(character_columns)) {
      warning(
        "Candidate identifiers were taken from row numbers, because id_col ",
        "was not supplied and values has no row names. Column(s) ",
        paste(sprintf("'%s'", utils::head(character_columns, 3L)),
              collapse = ", "),
        " look like identifiers; pass one through id_col to label the result.",
        call. = FALSE
      )
    }
  }

  direction <- rep(1, length(trait_cols))
  names(direction) <- trait_cols
  if (length(lower_is_better)) {
    unknown <- setdiff(lower_is_better, trait_cols)
    if (length(unknown)) {
      stop("Unknown lower_is_better traits: ",
           paste(unknown, collapse = ", "), call. = FALSE)
    }
    direction[lower_is_better] <- -1
  }
  centre <- if (isTRUE(center_traits)) {
    colMeans(X_raw)
  } else {
    rep(0, length(trait_cols))
  }
  names(centre) <- trait_cols
  scale_factor <- if (isTRUE(scale_traits)) {
    sd_values <- apply(X_raw, 2L, stats::sd)
    if (any(!is.finite(sd_values) | sd_values <= 0)) {
      stop("Every trait must vary when scale_traits = TRUE.", call. = FALSE)
    }
    sd_values
  } else {
    rep(1, length(trait_cols))
  }
  names(scale_factor) <- trait_cols
  X <- sweep(sweep(X_raw, 2L, centre, "-"), 2L, scale_factor, "/")
  X <- sweep(X, 2L, direction, "*")

  transform <- diag(direction / scale_factor, nrow = length(trait_cols))
  dimnames(transform) <- list(trait_cols, trait_cols)
  prepare_matrix <- function(M, name) {
    if (is.null(M)) return(NULL)
    M <- .dgr_covariance(M, trait_cols, name)
    M <- transform %*% M %*% transform
    dimnames(M) <- list(trait_cols, trait_cols)
    .dgr_check_psd(M, name)
    M
  }
  G <- prepare_matrix(G, "G")
  P <- prepare_matrix(P, "P")

  require_matrix <- function(M, name) {
    if (is.null(M)) {
      stop(name, " is required for method = '", method, "'.", call. = FALSE)
    }
    invisible(TRUE)
  }
  objective_vector <- function(x, name) {
    if (is.null(x)) {
      stop(name, " is required for method = '", method, "'.", call. = FALSE)
    }
    x <- .dgr_named_vector(x, trait_cols, name)
    if (any(x < 0)) {
      stop(
        name, " must contain non-negative magnitudes in the ",
        "favourable-direction space; use lower_is_better to declare traits ",
        "that should decrease.",
        call. = FALSE
      )
    }
    x
  }

  coefficients <- NULL
  aggregate_weights <- NULL
  scores <- NULL
  strategy <- "index"

  if (method == "smith_hazel") {
    require_matrix(G, "G"); require_matrix(P, "P")
    a <- objective_vector(economic_weights, "economic_weights")
    coefficients <- as.numeric(.dgr_inverse(P, "P")$inverse %*% G %*% a)
    aggregate_weights <- a
  } else if (method == "base") {
    a <- objective_vector(economic_weights, "economic_weights")
    coefficients <- as.numeric(a)
    aggregate_weights <- a
  } else if (method == "pesek_baker") {
    require_matrix(G, "G")
    d <- objective_vector(desired_gains, "desired_gains")
    coefficients <- as.numeric(.dgr_inverse(G, "G")$inverse %*% d)
    aggregate_weights <- d
  } else if (method == "yamada") {
    require_matrix(G, "G"); require_matrix(P, "P")
    d <- objective_vector(desired_gains, "desired_gains")
    P_inverse <- .dgr_inverse(P, "P")$inverse
    middle <- crossprod(G, P_inverse %*% G)
    middle <- (middle + t(middle)) / 2
    dimnames(middle) <- list(trait_cols, trait_cols)
    coefficients <- as.numeric(
      P_inverse %*% G %*% (.dgr_inverse(middle, "G P^-1 G")$inverse %*% d)
    )
    aggregate_weights <- d
  } else if (method == "mulamba_mock") {
    rank_weights <- if (is.null(economic_weights)) {
      NULL
    } else {
      objective_vector(economic_weights, "economic_weights")
    }
    scores <- -.dgr_rank_sum(X, rank_weights)
    strategy <- "rank_sum"
    aggregate_weights <- rank_weights
  } else if (method == "independent_culling") {
    thresholds <- objective_vector(culling_thresholds, "culling_thresholds")
    passes <- rowSums(sweep(X, 2L, thresholds, FUN = "<")) == 0L
    scores <- as.numeric(passes)
    strategy <- "culling"
  } else if (method == "tandem") {
    if (is.null(tandem_order) || !all(tandem_order %in% trait_cols)) {
      stop("tandem_order must name traits present in trait_cols.",
           call. = FALSE)
    }
    if (is.null(n_select)) {
      stop("n_select is required for method = 'tandem'.", call. = FALSE)
    }
    retained <- seq_len(nrow(X))
    stages <- length(tandem_order)
    per_stage <- ceiling(nrow(X) * (n_select / nrow(X))^(seq_len(stages) / stages))
    for (stage in seq_len(stages)) {
      trait <- tandem_order[stage]
      keep <- min(per_stage[stage], length(retained))
      ordering <- retained[order(-X[retained, trait], candidate_id[retained])]
      retained <- ordering[seq_len(keep)]
    }
    scores <- rep(0, nrow(X))
    scores[retained] <- rev(seq_along(retained))
    strategy <- "tandem"
  }

  if (is.null(scores)) {
    names(coefficients) <- trait_cols
    scores <- as.numeric(X %*% coefficients)
  }

  selection <- rep(FALSE, nrow(X))
  if (!is.null(n_select)) {
    n_select <- .dgr_positive_integer(n_select, "n_select")
    if (n_select > nrow(X)) {
      stop("n_select cannot exceed the number of candidates.", call. = FALSE)
    }
    selection[order(-scores, candidate_id)[seq_len(n_select)]] <- TRUE
  }

  if (is.null(selection_intensity)) {
    selection_intensity <- if (is.null(n_select)) {
      NA_real_
    } else {
      .dgr_intensity(n_select / nrow(X))
    }
  }

  evaluation <- NULL
  effective <- NULL
  if (strategy == "index" && !is.null(G)) {
    evaluation <- evaluate_index(
      coefficients = coefficients,
      G = G,
      P = if (is.null(P)) stats::cov(X) else P,
      aggregate_weights = aggregate_weights,
      selection_intensity = selection_intensity,
      main_trait = main_trait,
      scores = scores
    )
    effective <- effective_weights(coefficients, G, P)
  }

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

  result <- list(
    method = method,
    strategy = strategy,
    call = match.call(),
    trait_cols = trait_cols,
    coefficients = coefficients,
    aggregate_weights = aggregate_weights,
    ranking = ranking,
    # Filter on the column, not on the free-standing logical vector. The
    # vector is in the original candidate order whereas ranking has been
    # sorted by score, so using it here would select positionally wrong rows.
    selected = ranking[selected == TRUE],
    evaluation = evaluation,
    effective_weights = effective,
    observed_differential = observed,
    selection_intensity = selection_intensity,
    n_select = n_select,
    n_candidates = nrow(X),
    transformation = list(
      centre = centre, scale = scale_factor, direction = direction,
      centred = isTRUE(center_traits), scaled = isTRUE(scale_traits)
    ),
    G = G,
    P = P
  )
  class(result) <- c("desiredgainr_index", "list")
  result
}

#' @export
print.desiredgainr_index <- function(x, ...) {
  cat("<desiredgainr_index>\n")
  cat("  Method:", x$method, "\n")
  cat("  Candidates:", x$n_candidates,
      " Traits:", length(x$trait_cols), "\n")
  cat("  Traits standardised:",
      if (isTRUE(x$transformation$scaled)) "yes" else "no", "\n")
  if (!is.null(x$n_select)) {
    cat(sprintf(
      "  Selected: %d (%.1f%%), intensity %.3f\n",
      x$n_select, 100 * x$n_select / x$n_candidates, x$selection_intensity
    ))
  }
  if (!is.null(x$coefficients)) {
    cat("  Coefficients:\n")
    print(round(x$coefficients, 4L))
  }
  if (!is.null(x$evaluation)) {
    cat(sprintf(
      "  R_HI %.4f  dH %.4f  RE %.4f\n",
      x$evaluation$R_HI, x$evaluation$delta_H, x$evaluation$RE
    ))
    if (is.finite(x$evaluation$CV_I)) {
      cat(sprintf("  CV_I %.2f%%\n", x$evaluation$CV_I))
    } else {
      cat("  CV_I undefined for a centred index\n")
    }
  }
  invisible(x)
}
