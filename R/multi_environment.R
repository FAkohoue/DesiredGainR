# Multi-environment and genotype-by-environment structure.
#
# A multi-environment index is a single-environment index on an expanded trait
# set: every trait-in-environment is its own trait, and the covariance between
# them carries the genotype-by-environment structure. Writing it that way means
# the whole existing index layer applies unchanged, and the only new work is
# building the expanded covariance correctly and defining what merit means when
# the same trait appears in several environments.

#' Expand a trait-by-environment structure into an index problem
#'
#' Builds the covariance matrices and economic weights for an index over
#' trait-environment combinations, so that a multi-environment problem can be
#' handed to [selection_index()], [restricted_index()] or [evaluate_index()]
#' without further preparation.
#'
#' @details
#' # The model
#'
#' For trait \eqn{j} in environment \eqn{e}, the genetic value is
#'
#' \deqn{g_{je} = a_j + (ge)_{je}}{g_je = a_j + (ge)_je}
#'
#' with \eqn{a_j} the environment-invariant main effect and \eqn{(ge)_{je}} the
#' interaction. The covariance between the same trait in two environments is
#' then \eqn{\mathbf{G}_{jj}}, and between different traits and environments
#' \eqn{\mathbf{G}_{jk}} scaled by the environmental correlation. This gives
#' the separable structure
#'
#' \deqn{\mathrm{Cov}(g_{je}, g_{kf}) = \mathbf{G}_{jk} \, \mathbf{C}_{ef}}{
#' Cov(g_je, g_kf) = G_jk * C_ef}
#'
#' where \eqn{\mathbf{C}} is the between-environment genetic correlation
#' matrix. A \eqn{\mathbf{C}} with off-diagonals near one means the trait ranks
#' genotypes almost identically everywhere and a single-environment index would
#' have sufficed; values near zero mean the environments are effectively
#' different traits.
#'
#' Separability is an assumption, and it is the one to question first if the
#' result looks wrong. Where a fitted unstructured multi-environment model is
#' available, pass its covariance directly through
#' [import_covariance()] instead of building one here.
#'
#' # Environment-specific economic weights
#'
#' A trait need not be worth the same everywhere: yield in a marginal
#' environment serving a small area contributes less to net merit than yield in
#' the main target environment. `environment_weights` scales each environment's
#' contribution, and the economic weight of trait \eqn{j} in environment
#' \eqn{e} becomes \eqn{w_j \times v_e}.
#'
#' # Stability as an objective
#'
#' Setting `include_stability = TRUE` returns, for each trait, response
#' contrasts equal to the difference between the trait in each environment and
#' its mean across environments. Pass `stability$constraint_matrix` to
#' [restricted_index()] to impose zero response in those contrasts. They are
#' not appended as duplicate index variables, which would make the covariance
#' matrix singular because every contrast is a linear combination of the
#' original trait-environment variables.
#'
#' Stability is not free. Constraining it costs response in the mean, and
#' whether the trade is worth making is exactly the kind of question
#' [gain_feasibility()] and [restricted_index()] exist to answer.
#'
#' @param trait_cols Trait names.
#' @param environments Environment names.
#' @param G Genetic covariance between traits, named by trait.
#' @param P Phenotypic covariance between traits, named by trait.
#' @param environment_correlation Between-environment genetic correlation
#'   matrix, named by environment. A single number applies that correlation to
#'   every pair. Values near 1 mean little interaction.
#' @param residual_correlation Between-environment residual correlation.
#'   Defaults to 0, which is right when environments are separate trials.
#' @param economic_weights Named economic weights per trait.
#' @param environment_weights Named relative value of each environment,
#'   normalised to sum to one.
#' @param include_stability Whether to return per-trait stability response
#'   contrasts and a `constraint_matrix` consumable by [restricted_index()].
#'
#' @return A list of class `desiredgainr_multi_environment` holding the
#'   expanded `G`, `P`, economic weights and labels.
#'
#' @examples
#' traits <- c("yield", "protein")
#' environments <- c("irrigated", "rainfed")
#' G <- matrix(c(1.0, 0.2, 0.2, 0.5), 2, dimnames = list(traits, traits))
#' P <- matrix(c(2.5, 0.4, 0.4, 1.2), 2, dimnames = list(traits, traits))
#'
#' expanded <- expand_environments(
#'   trait_cols = traits, environments = environments,
#'   G = G, P = P,
#'   environment_correlation = 0.6,
#'   economic_weights = c(yield = 2, protein = 1),
#'   environment_weights = c(irrigated = 0.7, rainfed = 0.3)
#' )
#' dim(expanded$G)
#' expanded$economic_weights
#'
#' @seealso [selection_index()], [import_covariance()], [restricted_index()]
#' @export
expand_environments <- function(
  trait_cols,
  environments,
  G,
  P,
  environment_correlation,
  residual_correlation = 0,
  economic_weights = NULL,
  environment_weights = NULL,
  include_stability = FALSE
) {
  if (!is.character(trait_cols) || anyDuplicated(trait_cols)) {
    stop("trait_cols must contain unique trait names.", call. = FALSE)
  }
  if (!is.character(environments) || anyDuplicated(environments) ||
    length(environments) < 2L) {
    stop("environments must contain at least two unique names.", call. = FALSE)
  }
  p <- length(trait_cols)
  e <- length(environments)

  G <- .dgr_covariance(G, trait_cols, "G")
  P <- .dgr_covariance(P, trait_cols, "P")
  .dgr_check_compatible(G, P, "expand_environments()")

  as_correlation <- function(x, name) {
    if (is.matrix(x)) {
      x <- .dgr_covariance(x, environments, name)
      if (max(abs(diag(x) - 1)) > 1e-8) {
        stop(name, " must be a correlation matrix with a unit diagonal.",
          call. = FALSE
        )
      }
      return(x)
    }
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < -1 || x > 1) {
      stop(name, " must be a correlation matrix or a single value in [-1, 1].",
        call. = FALSE
      )
    }
    out <- matrix(x, e, e, dimnames = list(environments, environments))
    diag(out) <- 1
    out
  }
  C_genetic <- as_correlation(
    environment_correlation,
    "environment_correlation"
  )
  C_residual <- as_correlation(residual_correlation, "residual_correlation")
  if (min(eigen(C_genetic, symmetric = TRUE, only.values = TRUE)$values) <=
    -1e-8) {
    stop("environment_correlation is not a valid correlation matrix: it has ",
      "a negative eigenvalue.",
      call. = FALSE
    )
  }

  labels <- as.vector(outer(trait_cols, environments, paste, sep = "_"))
  # kronecker() gives exactly the separable structure Cov(g_je, g_kf) =
  # G_jk * C_ef, with the environment varying slowest so the labels above line
  # up with the rows.
  G_expanded <- kronecker(C_genetic, G)
  dimnames(G_expanded) <- list(labels, labels)

  # The residual is expanded separately and P reassembled from the two, so
  # that P - G stays a covariance matrix by construction. Expanding P directly
  # with the genetic environmental correlation would not guarantee that.
  E <- P - G
  E_expanded <- kronecker(C_residual, E)
  P_expanded <- G_expanded + E_expanded
  dimnames(P_expanded) <- list(labels, labels)
  .dgr_check_compatible(
    G_expanded, P_expanded, "the expanded multi-environment covariance"
  )

  weights <- NULL
  if (!is.null(economic_weights)) {
    economic_weights <- .dgr_named_vector(
      economic_weights, trait_cols, "economic_weights"
    )
    if (is.null(environment_weights)) {
      environment_weights <- stats::setNames(rep(1 / e, e), environments)
    } else {
      environment_weights <- .dgr_named_vector(
        environment_weights, environments, "environment_weights"
      )
      if (any(environment_weights < 0)) {
        stop("environment_weights must be non-negative.", call. = FALSE)
      }
      total <- sum(environment_weights)
      if (total <= 0) {
        stop("environment_weights must not sum to zero.", call. = FALSE)
      }
      environment_weights <- environment_weights / total
    }
    weights <- as.vector(outer(economic_weights, environment_weights, "*"))
    names(weights) <- labels
  }

  stability <- NULL
  if (isTRUE(include_stability)) {
    # One contrast per trait per environment, measuring departure from that
    # trait's mean across environments. The set is rank-deficient by
    # construction, since the deviations sum to zero within a trait, so the
    # last environment is dropped for each trait.
    contrast <- matrix(0, nrow = p * (e - 1L), ncol = p * e)
    contrast_names <- character(p * (e - 1L))
    row <- 0L
    for (trait_index in seq_len(p)) {
      for (environment_index in seq_len(e - 1L)) {
        row <- row + 1L
        columns <- trait_index + p * (seq_len(e) - 1L)
        contrast[row, columns] <- -1 / e
        contrast[row, trait_index + p * (environment_index - 1L)] <-
          1 - 1 / e
        contrast_names[row] <- paste0(
          trait_cols[trait_index], "_stability_",
          environments[environment_index]
        )
      }
    }
    stability <- list(
      contrast = contrast,
      constraint_matrix = contrast,
      labels = contrast_names,
      G = contrast %*% G_expanded %*% t(contrast),
      P = contrast %*% P_expanded %*% t(contrast),
      note = paste(
        "Each row is the deviation of a trait in one environment from that",
        "trait's mean across environments. A genotype with deviations of zero",
        "performs at its own mean everywhere, which is what stability means",
        "here. Selecting for stability costs response in the mean; use",
        "restricted_index() to price that trade rather than assume it."
      )
    )
    dimnames(stability$G) <- list(contrast_names, contrast_names)
    dimnames(stability$P) <- list(contrast_names, contrast_names)
  }

  interaction_share <- 1 - mean(C_genetic[upper.tri(C_genetic)])
  result <- list(
    trait_cols = trait_cols,
    environments = environments,
    labels = labels,
    G = G_expanded,
    P = P_expanded,
    economic_weights = weights,
    environment_weights = environment_weights,
    environment_correlation = C_genetic,
    stability = stability,
    interaction_share = interaction_share,
    interaction_share_is_heuristic = TRUE,
    structure = paste(
      "Separable: Cov(g_je, g_kf) = G_jk * C_ef. The reported interaction",
      "share is the heuristic 1 - mean environmental correlation, not a",
      "general variance-component identity. The mean between-environment",
      "genetic correlation is",
      format(mean(C_genetic[upper.tri(C_genetic)]), digits = 3),
      if (interaction_share < 0.1) {
        paste(
          "- the environments rank genotypes almost identically, so a",
          "single-environment index on the trait means would give nearly the",
          "same answer at a fraction of the complexity."
        )
      } else if (interaction_share > 0.6) {
        paste(
          "- the environments rank genotypes very differently. Treating them",
          "as one target population may not be appropriate; consider separate",
          "programmes."
        )
      } else {
        "- interaction is material and worth modelling."
      }
    )
  )
  class(result) <- c("desiredgainr_multi_environment", "list")
  result
}

#' @export
print.desiredgainr_multi_environment <- function(x, ...) {
  cat("<desiredgainr_multi_environment>\n")
  cat(sprintf(
    "  %d traits x %d environments = %d index variables\n",
    length(x$trait_cols), length(x$environments), length(x$labels)
  ))
  cat("  Environments:", paste(x$environments, collapse = ", "), "\n")
  cat(sprintf("  Interaction share: %.2f\n", x$interaction_share))
  if (!is.null(x$environment_weights)) {
    cat("  Environment weights:\n")
    print(round(x$environment_weights, 3L))
  }
  if (!is.null(x$stability)) {
    cat(sprintf(
      "  Stability contrasts: %d\n", length(x$stability$labels)
    ))
  }
  cat("\n  ", x$structure, "\n")
  invisible(x)
}

#' Reshape long-format multi-environment trial data for indexing
#'
#' Converts genotype-by-environment records in long format into the wide
#' trait-by-environment matrix that [expand_environments()] describes.
#'
#' @param data A data frame with one row per genotype, environment and trait
#'   combination.
#' @param id_col,environment_col Column names identifying the genotype and the
#'   environment.
#' @param trait_cols Trait column names.
#' @param environments Optional subset and ordering of environments.
#' @param missing_policy What to do with genotypes not present in every
#'   environment: `"error"`, `"drop"` to keep only complete genotypes, or
#'   `"mean_impute"`.
#'
#' @return A data frame with one row per genotype and one column per
#'   trait-environment combination, named to match `expand_environments()`.
#'
#' @examples
#' long <- expand.grid(
#'   geno = paste0("g", 1:5),
#'   env = c("irrigated", "rainfed"),
#'   stringsAsFactors = FALSE
#' )
#' set.seed(1)
#' long$yield <- stats::rnorm(nrow(long))
#' long$protein <- stats::rnorm(nrow(long))
#' widen_environments(long, "geno", "env", c("yield", "protein"))
#'
#' @seealso [expand_environments()]
#' @export
widen_environments <- function(
  data, id_col, environment_col, trait_cols,
  environments = NULL,
  missing_policy = c("error", "drop", "mean_impute")
) {
  missing_policy <- match.arg(missing_policy)
  frame <- as.data.frame(data)
  required <- c(id_col, environment_col, trait_cols)
  absent <- setdiff(required, names(frame))
  if (length(absent)) {
    stop("data is missing columns: ", paste(absent, collapse = ", "),
      call. = FALSE
    )
  }
  frame[[environment_col]] <- as.character(frame[[environment_col]])
  if (is.null(environments)) {
    environments <- sort(unique(frame[[environment_col]]))
  }
  frame <- frame[frame[[environment_col]] %in% environments, , drop = FALSE]
  keys <- paste(as.character(frame[[id_col]]), frame[[environment_col]],
    sep = "\r"
  )
  if (anyDuplicated(keys)) {
    duplicated_keys <- unique(keys[duplicated(keys) | duplicated(keys,
      fromLast = TRUE
    )])
    example <- strsplit(duplicated_keys[1L], "\r", fixed = TRUE)[[1L]]
    stop(
      "Each genotype-environment key must be unique; duplicate records were ",
      "found, including genotype '", example[1L], "' in environment '",
      example[2L], "'. Aggregate biological replicates explicitly before ",
      "widening so the chosen summary is auditable.",
      call. = FALSE
    )
  }
  ids <- unique(as.character(frame[[id_col]]))

  wide <- data.frame(id = ids, stringsAsFactors = FALSE)
  names(wide) <- id_col
  for (environment in environments) {
    block <- frame[frame[[environment_col]] == environment, , drop = FALSE]
    matched <- match(ids, as.character(block[[id_col]]))
    for (trait in trait_cols) {
      wide[[paste(trait, environment, sep = "_")]] <- block[[trait]][matched]
    }
  }

  value_columns <- setdiff(names(wide), id_col)
  incomplete <- !stats::complete.cases(wide[, value_columns, drop = FALSE])
  if (any(incomplete)) {
    if (identical(missing_policy, "error")) {
      stop(
        sum(incomplete), " genotype(s) were not observed in every ",
        "environment. A genotype missing from an environment has no value ",
        "for that index variable, so choose an explicit missing_policy.",
        call. = FALSE
      )
    }
    if (identical(missing_policy, "drop")) {
      wide <- wide[!incomplete, , drop = FALSE]
    } else {
      for (column in value_columns) {
        replacement <- mean(wide[[column]], na.rm = TRUE)
        if (!is.finite(replacement)) {
          stop("Column ", column, " has no observed values to impute from.",
            call. = FALSE
          )
        }
        wide[[column]][is.na(wide[[column]])] <- replacement
      }
    }
  }
  rownames(wide) <- NULL
  wide
}
