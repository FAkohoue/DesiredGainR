#' Fit an optimised desired-gain selection index
#'
#' `run_dgsi()` fits a desired-gain selection index (DGSI) from candidate
#' genetic values. All traits are first oriented so that larger values are
#' favourable. Independent stochastic replicates are run internally, and the
#' replicate with the smallest objective is returned automatically.
#'
#' @param init_data Data frame containing candidate identifiers and metadata.
#' @param cand_data Data frame containing candidate identifiers and one column
#'   per trait.
#' @param trait_cols Character vector naming the trait columns.
#' @param dg Named numeric vector of desired gains, expressed in **candidate
#'   standard-deviation units of the favourable-direction trait space**. This
#'   is true regardless of `scale_traits`: the realised response that `dg` is
#'   compared against is always divided by the candidate column standard
#'   deviations, so `dg = c(yield = 0.5)` requests a half-standard-deviation
#'   shift in the selected mean, never half a tonne per hectare. To express a
#'   target in original trait units, divide it by that trait's candidate
#'   standard deviation before passing it here. See Details.
#' @param P Phenotypic or index-variable covariance matrix. If `NULL`, an
#'   empirical working covariance matrix is estimated from `ref_data`; its
#'   provenance is reported and it is not described as a genetic covariance.
#' @param G Genetic covariance matrix or a complete object returned by
#'   [estimate_genetic_covariance()]. This input is required. Passing the
#'   complete estimate object preserves its method, estimand, assumptions, and
#'   diagnostics in the DGSI result.
#' @param ref_data Optional reference data used for scaling and, when needed,
#'   estimation of the working `P` matrix.
#' @param id_col Candidate identifier column.
#' @param scale_traits Whether to centre and scale traits using `ref_data`.
#' @param lower_is_better Traits for which smaller original values are
#'   favourable.
#' @param select_mode Either `"top_n"` or `"eligible_top_n"`. In
#'   `"eligible_top_n"` mode, thresholds define eligibility and the index ranks
#'   eligible candidates; thresholds never replace index ranking.
#' @param n_select Maximum number of candidates retained by the optimisation.
#' @param trait_min Named favourable-direction eligibility thresholds in the
#'   transformed analysis scale.
#' @param empty_eligibility Action when no candidate is eligible.
#' @param missing_policy Missing-value policy. The default rejects missing
#'   values. Imputation is never performed silently.
#' @param n_iter Search iterations per replicate.
#' @param n_rep Number of independent stochastic replicates.
#' @param sd_scale Perturbation scale for sampled desired-gain vectors.
#' @param seed Random seed.
#' @param ridge_P,ridge_M Non-negative ridge constants used in matrix solves.
#' @param objective_weights Optional named non-negative weights for the
#'   trait-specific squared response deviations.
#' @param plateau_window Number of final iterations used to diagnose a plateau.
#' @param plateau_tolerance Maximum relative improvement in the plateau window.
#' @param validation_data Optional independent candidate values. With
#'   `replicate_selection = "holdout"`, these data choose the stochastic
#'   replicate and are never used to fit its coefficients. When absent, an
#'   internal candidate holdout is split off before optimisation.
#' @param return_all_reps Whether to retain full replicate results.
#' @param allow_incompatible_estimated_P When `P` is estimated from `ref_data`
#'   and falls below `G` by more than sampling error explains, the two matrices
#'   cannot describe the same population and the call stops. Set this to `TRUE`
#'   to continue anyway; the override and the offending eigenvalue are recorded
#'   in `covariance_provenance$compatibility`.
#' @param replicate_selection How the winning stochastic search is chosen.
#'   `"holdout"`, the default, fits every replicate without the held-out
#'   candidates and chooses among them on that untouched set. `"training"`
#'   reproduces the earlier
#'   behaviour of taking the lowest training objective, which is a minimum over
#'   `n_rep` draws and is biased downward.
#' @param holdout_fraction Fraction of candidates reserved for choosing among
#'   replicates when `replicate_selection = "holdout"`.
#' @param debug Whether to print progress messages.
#'
#' @details
#' # Units of `dg` and of `realised_response`
#'
#' For each trait the realised response is
#' \deqn{r_j = \frac{\bar{x}_{j,\mathrm{selected}} -
#' \bar{x}_{j,\mathrm{all}}}{s_j},}
#' where \eqn{s_j} is the standard deviation of trait \eqn{j} across all
#' candidates in the favourable-direction analysis space. Both `dg` and
#' `realised_response` are therefore in candidate standard-deviation units.
#'
#' This holds whatever `scale_traits` is set to. The coefficient solve happens
#' in the units of `G` and `P`, which are raw trait units when
#' `scale_traits = FALSE`, so `dg` is converted from standard-deviation units
#' to that space before the solve and the coefficients are invariant to the
#' units the traits happen to be recorded in. Releases up to 0.5.0 passed the
#' standard-deviation vector into the solve directly, which asked for
#' \eqn{d_j} *raw* units per trait; with candidate standard deviations of 10
#' and 1 and equal requested gains, that delivered standardised response in the
#' ratio 1:10 rather than 1:1.
#'
#' # Empirical differential against transmitted response
#'
#' `realised_response` is a property of the candidates that were selected: the
#' standardised differential between the selected group and the whole set. It
#' is not what the next generation inherits.
#'
#' `theoretical_response` reports the model-based expectation
#' \deqn{\Delta G = i \frac{\mathbf{G}\mathbf{b}}
#' {\sqrt{\mathbf{b}^\mathsf{T}\mathbf{P}\mathbf{b}}}}{
#' i * G b / sqrt(b' P b)}
#' in the analysis space, in the original trait units, and standardised. This
#' is the criterion on which a desired-gain index can be compared with the
#' classical families, which report the same quantity through
#' [evaluate_index()]. Keep the two apart when reporting: a large differential
#' among candidates with a small transmitted response means selection is
#' acting mostly on non-heritable variation.
#'
#' # Objective function
#'
#' The search minimises
#' \deqn{\sum_j v_j \left(\frac{r_j - d_j}{\max(|d_j|,\,0.25)}\right)^2,}
#' with `objective_weights` \eqn{v_j}. The floor of `0.25` on the denominator
#' prevents traits with a near-zero desired gain from dominating the objective
#' through division by a vanishing scale factor. Proposal standard deviations
#' are `sd_scale * pmax(abs(centre), 0.05)`, where the `0.05` floor keeps the
#' search from stalling when a component of the desired-gain vector is at or
#' near zero. Both floors are fixed constants in this release.
#'
#' # Optimisation and reproducibility
#'
#' The search is a stochastic hill climb over perturbed desired-gain vectors:
#' each proposal is mapped to coefficients by
#' \eqn{b = P^{-1}G(G^\mathsf{T}P^{-1}G)^{-1}d} and accepted only if it lowers
#' the objective. The objective is a step function of \eqn{b}, because it
#' depends on \eqn{b} only through the identity of the selected set, so a
#' derivative-free search is used. With the default holdout rule, the split is
#' made before optimisation and the replicate is chosen on candidates that did
#' not contribute to fitting. `validation_data`, when supplied, is preferred
#' to an internal split and can therefore select the replicate as well as
#' evaluate it. The final coefficients are applied to all candidates only
#' after that choice; they are not refitted. `run_dgsi()` seeds the RNG from
#' `seed` and restores the caller's RNG state before returning.
#'
#' Ties in the index score are broken by ascending candidate identifier, so
#' the selected set does not depend on input row order.
#'
#' @return An object of class `desired_gain_index`. The `best_replicate`
#'   component identifies the replicate selected automatically. The
#'   `replicate_diagnostics`, `rank_correlation`, `coefficient_stability`, and
#'   `selected_set_agreement` components describe optimisation stability.
#'
#' @examples
#' set.seed(3)
#' traits <- c("yield", "disease")
#' candidates <- data.frame(
#'   GenoID = paste0("G", seq_len(40)),
#'   yield = rnorm(40),
#'   disease = rnorm(40)
#' )
#' # In practice G comes from a fitted multi-trait genetic model, or from
#' # estimate_genetic_covariance(); it is not the covariance of raw phenotypes.
#' G <- matrix(c(1.0, -0.3, -0.3, 0.8), 2, dimnames = list(traits, traits))
#' P <- G + diag(c(0.8, 0.6))
#' dimnames(P) <- list(traits, traits)
#'
#' fit <- run_dgsi(
#'   init_data = candidates["GenoID"],
#'   cand_data = candidates,
#'   trait_cols = traits,
#'   dg = c(yield = 0.6, disease = 0.4),
#'   G = G, P = P,
#'   lower_is_better = "disease",
#'   n_select = 8,
#'   n_iter = 50,
#'   n_rep = 3,
#'   seed = 3
#' )
#' fit$coefficients
#' fit$realised_response
#' fit$replicate_diagnostics
#'
#' @export
run_dgsi <- function(
  init_data,
  cand_data,
  trait_cols,
  dg,
  P = NULL,
  G,
  ref_data = NULL,
  id_col = "GenoID",
  scale_traits = FALSE,
  lower_is_better = NULL,
  select_mode = c("top_n", "eligible_top_n"),
  n_select = 100L,
  trait_min = NULL,
  empty_eligibility = c("error", "fallback_top_n"),
  missing_policy = c("error", "complete_cases", "mean_impute"),
  n_iter = 1000L,
  n_rep = 20L,
  sd_scale = 1,
  seed = 42L,
  ridge_P = 1e-6,
  ridge_M = 1e-6,
  objective_weights = NULL,
  plateau_window = 100L,
  plateau_tolerance = 1e-4,
  validation_data = NULL,
  return_all_reps = TRUE,
  allow_incompatible_estimated_P = FALSE,
  replicate_selection = c("holdout", "training"),
  holdout_fraction = 0.3,
  debug = FALSE
) {
  replicate_selection <- match.arg(replicate_selection)
  select_mode <- match.arg(select_mode)
  empty_eligibility <- match.arg(empty_eligibility)
  missing_policy <- match.arg(missing_policy)

  if (!length(trait_cols) || anyDuplicated(trait_cols)) {
    stop("trait_cols must contain unique trait names.", call. = FALSE)
  }
  n_select <- .dgr_positive_integer(n_select, "n_select")
  seed <- .dgr_seed(seed)
  for (control in c("ridge_P", "ridge_M", "plateau_tolerance")) {
    value <- get(control)
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value < 0) {
      stop(control, " must be a single non-negative finite number.",
        call. = FALSE
      )
    }
  }
  n_iter <- .dgr_positive_integer(n_iter, "n_iter")
  n_rep <- .dgr_positive_integer(n_rep, "n_rep")
  plateau_window <- .dgr_positive_integer(plateau_window, "plateau_window")
  if (!is.numeric(sd_scale) || length(sd_scale) != 1L ||
    !is.finite(sd_scale) || sd_scale <= 0) {
    stop("sd_scale must be a positive finite number.", call. = FALSE)
  }
  if (any(!is.finite(c(ridge_P, ridge_M))) ||
    any(c(ridge_P, ridge_M) < 0)) {
    stop("ridge_P and ridge_M must be non-negative finite numbers.",
      call. = FALSE
    )
  }

  prep <- .dgr_prepare_values(
    init_data = init_data,
    cand_data = cand_data,
    ref_data = ref_data,
    validation_data = validation_data,
    trait_cols = trait_cols,
    id_col = id_col,
    lower_is_better = lower_is_better,
    scale_traits = scale_traits,
    centre_traits = scale_traits,
    missing_policy = missing_policy
  )
  X <- prep$candidate_matrix
  Xref <- prep$reference_matrix
  p <- length(trait_cols)

  # Covarrubias-Pazaran (2021) standardises adjusted means before indexing, and
  # Crosbie et al. (1980) showed that equal desired gains on unstandardised
  # traits concentrate selection on whichever trait carries the largest
  # variance. Warn only when the trait scales actually differ enough for this
  # to matter, so that analyses on comparable scales stay quiet.
  if (!isTRUE(scale_traits) && p > 1L) {
    candidate_sd <- apply(X, 2L, stats::sd)
    candidate_sd <- candidate_sd[is.finite(candidate_sd) & candidate_sd > 0]
    if (length(candidate_sd) > 1L &&
      max(candidate_sd) / min(candidate_sd) > 5) {
      warning(
        sprintf(
          paste(
            "Trait standard deviations differ by a factor of %.0f and",
            "scale_traits = FALSE. Desired gains are interpreted in candidate",
            "standard-deviation units regardless, so the index will be",
            "dominated by the highest-variance trait. Consider",
            "scale_traits = TRUE."
          ),
          max(candidate_sd) / min(candidate_sd)
        ),
        call. = FALSE
      )
    }
  }

  dg <- .dgr_named_vector(dg, trait_cols, "dg")
  if (is.null(objective_weights)) {
    objective_weights <- rep(1, p)
    names(objective_weights) <- trait_cols
  } else {
    objective_weights <- .dgr_named_vector(
      objective_weights, trait_cols, "objective_weights"
    )
    if (any(objective_weights < 0) || !any(objective_weights > 0)) {
      stop("objective_weights must be non-negative with at least one positive value.",
        call. = FALSE
      )
    }
  }

  if (missing(G) || is.null(G)) {
    stop(
      paste(
        "G is required. Supply a genetic covariance matrix from an appropriate",
        "model, or explicitly construct a working approximation with",
        "estimate_genetic_covariance()."
      ),
      call. = FALSE
    )
  }
  if (inherits(G, "desiredgainr_covariance_estimate")) {
    G_provenance <- list(
      source = G$provenance,
      method = G$method,
      estimand = G$estimand,
      assumptions = G$assumptions,
      diagnostics = G$diagnostics
    )
    G <- G$G
  } else {
    G_provenance <- list(
      source = "user supplied genetic covariance matrix",
      method = "external",
      estimand = "defined by the upstream model",
      assumptions = NULL,
      diagnostics = NULL
    )
  }
  G <- .dgr_covariance(G, trait_cols, "G")

  # Argument validation runs before any covariance work. Estimating P can emit
  # warnings, and a call that is going to fail on a missing argument should
  # fail on that argument rather than first warning about a matrix the caller
  # will never use.
  if (select_mode == "eligible_top_n") {
    if (is.null(trait_min)) {
      stop("trait_min is required when select_mode = 'eligible_top_n'.",
        call. = FALSE
      )
    }
    trait_min <- .dgr_named_vector(trait_min, trait_cols, "trait_min")
  }

  P_was_estimated <- is.null(P)
  P_rank <- NA_integer_
  if (P_was_estimated) {
    # A covariance estimated from fewer records than traits is singular, and a
    # singular matrix passes a positive-semidefinite check unchanged. The ridge
    # then makes it invertible, so the failure is silent and the resulting
    # coefficients are noise.
    if (nrow(Xref) <= p) {
      stop(
        sprintf(
          paste(
            "P was not supplied and cannot be estimated: %d reference record(s)",
            "for %d traits gives a singular covariance matrix. Supply P from",
            "the fitted model, or provide at least %d reference records."
          ),
          nrow(Xref), p, p + 1L
        ),
        call. = FALSE
      )
    }
    P <- stats::cov(Xref)
    P_rank <- sum(
      eigen(P, symmetric = TRUE, only.values = TRUE)$values >
        1e-10 * max(1, max(abs(diag(P))))
    )
    if (nrow(Xref) < 5L * p) {
      warning(
        sprintf(
          paste(
            "P was estimated from only %d reference records for %d traits.",
            "A covariance matrix estimated from so few records is poorly",
            "determined, and the index coefficients inherit that noise."
          ),
          nrow(Xref), p
        ),
        call. = FALSE
      )
    }
    P_source <- sprintf(
      "empirical working covariance estimated from %d reference records",
      nrow(Xref)
    )
  } else {
    P <- .dgr_covariance(P, trait_cols, "P")
    P_source <- "user supplied"
  }

  transform <- diag(prep$direction / prep$scale, p)
  if (!P_was_estimated) {
    P <- transform %*% P %*% transform
  }
  G <- transform %*% G %*% transform
  dimnames(P) <- dimnames(G) <- list(trait_cols, trait_cols)
  .dgr_check_psd(P, "P")
  .dgr_check_psd(G, "G")
  # The admissibility of P - G is recorded either way, so that a result can be
  # audited without re-deriving it. A supplied P is held to the strict
  # standard; an estimated P is a sample covariance and is given a sampling
  # allowance, but that exception is now part of the API rather than an
  # undocumented leniency.
  residual <- (P - G + t(P - G)) / 2
  smallest_residual_eigenvalue <- min(
    eigen(residual, symmetric = TRUE, only.values = TRUE)$values
  )
  compatibility <- list(
    smallest_residual_eigenvalue = smallest_residual_eigenvalue,
    smallest_standardised_residual_eigenvalue = NA_real_,
    sampling_threshold = NA_real_,
    status = NA_character_,
    override = FALSE
  )
  if (P_was_estimated) {
    # Work in phenotypic-SD units. A raw-eigenvalue allowance based on the
    # largest variance changes when a trait is expressed in different units;
    # the congruence below is invariant to any diagonal unit conversion.
    phenotypic_scale <- sqrt(diag(P))
    standardise <- diag(1 / phenotypic_scale, p)
    standardised_residual <- standardise %*% residual %*% standardise
    smallest_standardised <- min(eigen(
      standardised_residual,
      symmetric = TRUE, only.values = TRUE
    )$values)
    sampling_scale <- sqrt(p / nrow(Xref))
    compatibility$smallest_standardised_residual_eigenvalue <-
      smallest_standardised
    compatibility$sampling_threshold <- -2 * sampling_scale
    compatibility$status <- if (smallest_standardised >= 0) {
      "admissible"
    } else if (smallest_standardised >= -2 * sampling_scale) {
      "within sampling allowance"
    } else {
      "inadmissible"
    }
    if (identical(compatibility$status, "inadmissible")) {
      compatibility$override <- isTRUE(allow_incompatible_estimated_P)
      message_text <- sprintf(
        paste(
          "The working P estimated from %d reference records is smaller than G",
          "in at least one direction by more than sampling error explains:",
          "the smallest eigenvalue of the variance-standardised P - G is",
          "%.3g against a sampling threshold of %.3g. The two matrices",
          "cannot both describe the same",
          "population, so heritabilities and accuracies derived from them are",
          "not interpretable. Supply P from the fitted model that produced G."
        ),
        nrow(Xref), smallest_standardised, -2 * sampling_scale
      )
      if (isTRUE(allow_incompatible_estimated_P)) {
        warning(message_text, " Continuing because ",
          "allow_incompatible_estimated_P = TRUE.",
          call. = FALSE
        )
      } else {
        stop(message_text, "\n  Set allow_incompatible_estimated_P = TRUE to ",
          "continue anyway; the override is recorded in the result.",
          call. = FALSE
        )
      }
    }
  } else {
    .dgr_check_compatible(G, P, "run_dgsi()")
    compatibility$status <- "admissible"
  }

  eligible <- rep(TRUE, nrow(X))
  if (select_mode == "eligible_top_n") {
    eligible <- rowSums(sweep(X, 2L, trait_min, FUN = "<")) == 0L
    if (!any(eligible)) {
      if (empty_eligibility == "error") {
        stop("No candidate meets all favourable-direction eligibility thresholds.",
          call. = FALSE
        )
      }
      eligible[] <- TRUE
    }
  }

  # Ties are broken by candidate identifier so that the selected set does not
  # depend on the row order of the input data.
  candidate_key <- as.character(prep$init_data[[id_col]])
  select_rows <- function(score, eligibility = eligible,
                          key = candidate_key) {
    pool <- which(eligibility)
    pool[order(-score[pool], key[pool])][
      seq_len(min(n_select, length(pool)))
    ]
  }
  # Column means and standard deviations do not change during the search, so
  # they are computed once per matrix instead of once per iteration.
  .column_stats <- function(matrix) {
    base_sd <- apply(matrix, 2L, stats::sd)
    base_sd[!is.finite(base_sd) | base_sd == 0] <- 1
    list(centre = colMeans(matrix), sd = base_sd)
  }
  X_stats <- .column_stats(X)
  realised_response <- function(score, matrix = X, eligibility = eligible,
                                key = candidate_key, stats_cache = X_stats) {
    keep <- select_rows(score, eligibility, key)
    response <- (colMeans(matrix[keep, , drop = FALSE]) -
      stats_cache$centre) / stats_cache$sd
    names(response) <- trait_cols
    list(response = response, keep = keep)
  }

  P_r <- P + diag(ridge_P, p)
  Xpg <- solve(P_r, G)
  middle <- crossprod(G, Xpg)
  middle <- (middle + t(middle)) / 2 + diag(ridge_M, p)

  # `dg` and `realised_response` are in candidate standard-deviation units, but
  # the Yamada map operates in the units of G and P, which are the analysis
  # space. With scale_traits = TRUE the two coincide, because the analysis
  # space is already standardised. With scale_traits = FALSE the analysis space
  # is raw trait units, and passing an SD-unit vector into the solve silently
  # asks for a response of d_j RAW units per trait instead of d_j standard
  # deviations.
  #
  # The consequence is a real misdirection, not a labelling nicety. For two
  # traits with candidate standard deviations 10 and 1 and equal requested SD
  # gains, the solve targets raw response proportional to (1, 1), which is
  # standardised response proportional to (0.1, 1): ten times too little of the
  # first trait. Converting here makes the coefficients invariant to the units
  # the traits happen to be measured in.
  gain_scale <- X_stats$sd
  coefficient <- function(d) {
    value <- as.numeric(Xpg %*% solve(middle, d * gain_scale))
    names(value) <- trait_cols
    value
  }
  objective <- function(response) {
    denominator <- pmax(abs(dg), 0.25)
    sum(objective_weights * ((response - dg) / denominator)^2)
  }

  # Restore the caller's random number generator state on exit. run_dgsi()
  # seeds the RNG for reproducibility, but a package should not leave the
  # user's stream permanently advanced or reseeded.
  if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    stats::runif(1L)
  }
  .dgr_entry_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit(
    assign(".Random.seed", .dgr_entry_seed, envir = globalenv()),
    add = TRUE
  )

  set.seed(seed)

  # Decide the fitting/comparison split before any coefficient search. An
  # internal holdout that is sampled after fitting is not a holdout: changing
  # those observations can already have changed every fitted replicate.
  fit_rows <- seq_len(nrow(X))
  comparison <- NULL
  if (identical(replicate_selection, "holdout") && n_rep > 1L) {
    if (!is.numeric(holdout_fraction) || length(holdout_fraction) != 1L ||
      !is.finite(holdout_fraction) ||
      holdout_fraction <= 0 || holdout_fraction >= 1) {
      stop("holdout_fraction must lie strictly between 0 and 1.",
        call. = FALSE
      )
    }
    if (!is.null(prep$validation_matrix)) {
      comparison <- list(
        matrix = prep$validation_matrix,
        eligibility = rep(TRUE, nrow(prep$validation_matrix)),
        key = sprintf("validation_%08d", seq_len(nrow(prep$validation_matrix))),
        source = "external validation data",
        rows = NULL
      )
    } else {
      n_holdout <- max(p + 1L, floor(nrow(X) * holdout_fraction))
      minimum_training <- max(p + 1L, n_select)
      if (nrow(X) - n_holdout < minimum_training) {
        warning(
          "Too few candidates for a pre-fit holdout while retaining at least ",
          minimum_training, " training candidates. The winning replicate was ",
          "chosen on the training objective, which is optimistically biased.",
          call. = FALSE
        )
      } else {
        # Sample identifiers in a canonical order so an input-row permutation
        # cannot change the holdout under the same seed.
        holdout_keys <- sample(sort(candidate_key), n_holdout)
        holdout_rows <- sort(match(holdout_keys, candidate_key))
        fit_rows <- setdiff(seq_len(nrow(X)), holdout_rows)
        comparison <- list(
          matrix = X[holdout_rows, , drop = FALSE],
          eligibility = eligible[holdout_rows],
          key = candidate_key[holdout_rows],
          source = "internal pre-fit holdout",
          rows = holdout_rows
        )
      }
    }
  }
  fit_matrix <- X[fit_rows, , drop = FALSE]
  fit_eligible <- eligible[fit_rows]
  fit_key <- candidate_key[fit_rows]
  fit_stats <- .column_stats(fit_matrix)

  replicates <- vector("list", n_rep)
  for (rr in seq_len(n_rep)) {
    centre <- dg
    initial_b <- coefficient(centre)
    initial_eval <- realised_response(
      as.numeric(fit_matrix %*% initial_b),
      matrix = fit_matrix,
      eligibility = fit_eligible, key = fit_key, stats_cache = fit_stats
    )
    best <- list(
      objective = objective(initial_eval$response),
      d = centre,
      b = initial_b,
      response = initial_eval$response,
      keep = initial_eval$keep,
      iteration = 0L
    )
    objective_trace <- numeric(n_iter)
    best_trace <- numeric(n_iter)

    for (iteration in seq_len(n_iter)) {
      proposal_sd <- sd_scale * pmax(abs(centre), 0.05)
      proposal <- stats::rnorm(p, mean = centre, sd = proposal_sd)
      names(proposal) <- trait_cols
      b <- coefficient(proposal)
      evaluation <- realised_response(
        as.numeric(fit_matrix %*% b),
        matrix = fit_matrix,
        eligibility = fit_eligible, key = fit_key, stats_cache = fit_stats
      )
      value <- objective(evaluation$response)
      objective_trace[iteration] <- value
      if (is.finite(value) && value < best$objective) {
        best <- list(
          objective = value,
          d = proposal,
          b = b,
          response = evaluation$response,
          keep = evaluation$keep,
          iteration = iteration
        )
        centre <- proposal
      }
      best_trace[iteration] <- best$objective
    }

    score <- as.numeric(X %*% best$b)
    full_evaluation <- realised_response(score)
    window_start <- max(1L, n_iter - plateau_window + 1L)
    relative_change <- (
      best_trace[window_start] - best_trace[n_iter]
    ) / max(abs(best_trace[window_start]), 1e-12)
    replicates[[rr]] <- c(
      best[c("d", "b", "iteration")],
      list(
        objective = best$objective,
        training_objective = best$objective,
        training_response = best$response,
        response = full_evaluation$response,
        keep = full_evaluation$keep,
        full_candidate_objective = objective(full_evaluation$response),
        replicate = rr,
        score = score,
        objective_trace = objective_trace,
        best_trace = best_trace,
        plateau = relative_change <= plateau_tolerance,
        final_window_relative_improvement = relative_change
      )
    )
  }

  objectives <- vapply(replicates, function(x) x$objective, numeric(1))

  # Replicates fitted without the comparison observations can now be compared
  # honestly on that untouched matrix. External validation data are preferred
  # because they leave every candidate available for fitting.
  holdout_selection <- NULL
  if (!is.null(comparison)) {
    comparison_stats <- .column_stats(comparison$matrix)
    comparison_objectives <- vapply(replicates, function(replicate) {
      evaluation <- realised_response(
        as.numeric(comparison$matrix %*% replicate$b),
        matrix = comparison$matrix,
        eligibility = comparison$eligibility,
        key = comparison$key,
        stats_cache = comparison_stats
      )
      objective(evaluation$response)
    }, numeric(1L))
    holdout_selection <- list(
      n_holdout = nrow(comparison$matrix),
      holdout_rows = comparison$rows,
      source = comparison$source,
      objectives = comparison_objectives,
      chosen = which.min(comparison_objectives),
      training_rows = fit_rows
    )
  }

  best_replicate <- if (is.null(holdout_selection)) {
    which.min(objectives)
  } else {
    holdout_selection$chosen
  }
  best <- replicates[[best_replicate]]
  score_matrix <- do.call(cbind, lapply(replicates, function(x) x$score))
  coefficient_matrix <- do.call(
    cbind, lapply(replicates, function(x) x$b)
  )
  selected_matrix <- vapply(
    replicates,
    function(x) seq_len(nrow(X)) %in% x$keep,
    logical(nrow(X))
  )

  rank_correlation <- if (n_rep == 1L) {
    matrix(1, 1L, 1L, dimnames = list("replicate_1", "replicate_1"))
  } else {
    suppressWarnings(stats::cor(score_matrix, method = "spearman"))
  }
  selected_agreement <- .dgr_jaccard_matrix(selected_matrix)
  coefficient_stability <- data.table::data.table(
    Trait = trait_cols,
    Mean = rowMeans(coefficient_matrix),
    SD = apply(coefficient_matrix, 1L, stats::sd),
    CV = apply(coefficient_matrix, 1L, function(x) {
      stats::sd(x) / max(abs(mean(x)), 1e-12)
    })
  )
  replicate_diagnostics <- data.table::rbindlist(lapply(replicates, function(x) {
    data.table::data.table(
      Replicate = x$replicate,
      Objective = x$objective,
      Iteration_of_best = x$iteration,
      Selected = length(x$keep),
      Plateau = x$plateau,
      Final_window_relative_improvement =
        x$final_window_relative_improvement
    )
  }))
  replicate_diagnostics[, Chosen := Replicate == best_replicate]

  ranked <- data.table::copy(prep$init_data)
  ranked[, SelectionIndex := best$score]
  ranked[, Eligible := eligible]
  ranked[, Selected := seq_len(.N) %in% best$keep]
  ranked[, Rank := data.table::frank(-SelectionIndex, ties.method = "average")]
  data.table::setorder(ranked, -SelectionIndex)

  validation <- NULL
  if (!is.null(prep$validation_matrix)) {
    validation_score <- as.numeric(prep$validation_matrix %*% best$b)
    validation_eval <- realised_response(
      validation_score,
      matrix = prep$validation_matrix,
      eligibility = rep(TRUE, nrow(prep$validation_matrix)),
      key = seq_len(nrow(prep$validation_matrix)),
      stats_cache = .column_stats(prep$validation_matrix)
    )
    validation <- list(
      score = validation_score,
      realised_response = validation_eval$response,
      selected_rows = validation_eval$keep,
      objective = objective(validation_eval$response)
    )
  }

  result <- list(
    method = "DGSI",
    call = match.call(),
    desired_gain = dg,
    P = P,
    G = G,
    covariance_provenance = list(
      P = P_source,
      P_numerical_rank = P_rank,
      P_was_estimated = P_was_estimated,
      compatibility = compatibility,
      G = G_provenance
    ),
    # Keep training, comparison and full-candidate objectives separate. Only
    # the comparison objective is out-of-sample under the holdout rule.
    optimism = local({
      values <- objectives
      held_out <- !is.null(holdout_selection)
      selection_value <- if (held_out) {
        holdout_selection$objectives[best_replicate]
      } else {
        values[best_replicate]
      }
      list(
        selection_rule = if (held_out) {
          holdout_selection$source
        } else {
          "training"
        },
        chosen_objective = selection_value,
        chosen_training_objective = values[best_replicate],
        full_candidate_objective = best$full_candidate_objective,
        minimum_objective = min(values),
        mean_objective = mean(values),
        median_objective = stats::median(values),
        gap = mean(values) - min(values),
        relative_gap = if (mean(values) > 0) {
          (mean(values) - min(values)) / mean(values)
        } else {
          NA_real_
        },
        n_replicates = length(values),
        holdout = holdout_selection,
        note = if (held_out) {
          paste(
            "The replicate was chosen on", holdout_selection$n_holdout,
            "observations from", holdout_selection$source, "that were not",
            "used to fit any replicate. Its comparison objective is",
            "out-of-sample. The coefficients were then applied without",
            "refitting to the full candidate set."
          )
        } else {
          paste(
            "The replicate was chosen on the same candidates it then selects,",
            "so its objective is a minimum over replicates and is biased",
            "downward. Use the default holdout rule with enough candidates or",
            "supply independent validation_data to remove this selection bias."
          )
        }
      )
    }),
    transformation = prep$transformation,
    missing_data = prep$missing_data,
    eligibility = list(
      mode = select_mode,
      thresholds = trait_min,
      eligible_n = sum(eligible),
      selected_n = length(best$keep)
    ),
    optimised_d = best$d,
    coefficients = best$b,
    realised_response = best$response,
    # The empirical differential above is a property of the candidates that
    # happened to be selected. The transmitted genetic response is a different
    # quantity, is what the next generation actually inherits, and is the
    # criterion on which DGSI can be compared with the classical index
    # families. Reporting only the first invites the two to be confused.
    theoretical_response = local({
      b <- best$b
      index_sd <- sqrt(as.numeric(crossprod(b, P %*% b)))
      if (!is.finite(index_sd) || index_sd <= 0) {
        NULL
      } else {
        intensity <- .dgr_intensity(length(best$keep) / nrow(X))
        analysis <- as.numeric(intensity * (G %*% b) / index_sd)
        names(analysis) <- trait_cols
        # Back out of the direction-and-scale transform to give the same
        # response in the units the traits were supplied in.
        original <- analysis * prep$scale * prep$direction
        names(original) <- trait_cols
        list(
          analysis_units = analysis,
          original_units = original,
          standardised = analysis / gain_scale,
          selection_intensity = intensity,
          index_sd = index_sd,
          equation = "i * G b / sqrt(b' P b)",
          note = paste(
            "Expected transmitted genetic response under truncation selection",
            "at the achieved selected fraction. This is a model-based",
            "prediction from G and P, not the differential observed among the",
            "selected candidates, which is reported as realised_response."
          )
        )
      }
    }),
    objective = if (is.null(holdout_selection)) {
      best$objective
    } else {
      holdout_selection$objectives[best_replicate]
    },
    # Joukhadar et al. (2024) benchmark the iterative solution against the
    # classical index obtained by substituting the desired gains directly.
    # That comparator costs nothing, because it is the search starting point.
    non_iterated = local({
      baseline_b <- coefficient(dg)
      baseline <- realised_response(as.numeric(X %*% baseline_b))
      list(
        coefficients = baseline_b,
        realised_response = baseline$response,
        objective = objective(baseline$response),
        note = paste(
          "Classical Yamada solution using the supplied desired gains without",
          "iteration. Because d enters as a scalar direction, this is",
          "invariant to multiplying every desired gain by a constant."
        )
      )
    }),
    # A diagnostic must never be able to break the fit, and the scale warning
    # is issued separately and conditionally above.
    effective_weights = tryCatch(
      effective_weights(best$b, G, P, warn = FALSE),
      error = function(e) NULL
    ),
    best_replicate = best_replicate,
    replicate_diagnostics = replicate_diagnostics,
    rank_correlation = rank_correlation,
    coefficient_stability = coefficient_stability,
    selected_set_agreement = selected_agreement,
    ranked_geno = ranked,
    selected_geno = ranked[Selected == TRUE],
    non_selected_geno = ranked[Selected == FALSE],
    validation = validation
  )
  if (isTRUE(return_all_reps)) result$all_reps <- replicates
  # "desired_gain_index" is retained as the first class because dependent
  # packages test for it, and renaming it would break them. The consistent
  # name is added alongside so that new code has one, and both dispatch to the
  # same methods.
  class(result) <- c("desired_gain_index", "desiredgainr_dgsi", "list")
  .desiredgainr_dbg(
    debug,
    "Selected replicate %d automatically (objective %.6g).",
    best_replicate,
    best$objective
  )
  result
}

#' @export
print.desired_gain_index <- function(x, ...) {
  cat("<desired_gain_index>\n")
  cat(sprintf(
    "  Traits: %d   Candidates: %d   Selected: %d\n",
    length(x$trait_cols), nrow(x$ranked_geno), x$eligibility$selected_n
  ))
  cat("  Selection rule:", x$eligibility$mode, "\n")
  cat("  Coefficients:\n")
  print(round(x$coefficients, 4L))
  cat("  Realised response against the desired gains:\n")
  print(round(x$realised_response, 4L))
  if (!is.null(x$optimism)) {
    cat(sprintf(
      "  Objective %.6g, chosen from %d replicates.\n",
      x$optimism$chosen_objective, x$optimism$n_replicates
    ))
    if (identical(x$optimism$selection_rule, "training")) {
      cat("  This is the minimum training objective and is biased downward.\n")
    } else {
      cat("  Replicate selection used:", x$optimism$selection_rule, "\n")
      cat("  The comparison observations were excluded before fitting.\n")
    }
  }
  cat("  P provenance:", x$covariance_provenance$P, "\n")
  invisible(x)
}

#' Validate a count
#'
#' Integrality is confirmed before any coercion. Calling `as.integer()` first
#' would silently truncate 1.9 to 1, and would return NA with only a warning
#' for values above the integer limit, so both checks are made on the supplied
#' value while it is still a double.
#'
#' @noRd
.dgr_positive_integer <- function(x, name) {
  if (length(x) != 1L) {
    stop(name, " must be a single number, and was of length ", length(x),
      ".",
      call. = FALSE
    )
  }
  # NA is checked before the type, because a bare NA is logical and reporting
  # "must be numeric, and was logical" would describe the wrong problem.
  if (is.atomic(x) && is.na(x)) {
    stop(name, " must be a positive integer and was NA.", call. = FALSE)
  }
  if (!is.numeric(x)) {
    stop(name, " must be a single number, and was ", class(x)[1L], ".",
      call. = FALSE
    )
  }
  if (!is.finite(x)) {
    stop(name, " must be a positive integer and was ", x, ".", call. = FALSE)
  }
  if (x < 1) {
    stop(name, " must be a positive integer and was ", x,
      ". Zero and negative counts are not meaningful here.",
      call. = FALSE
    )
  }
  if (x > .Machine$integer.max) {
    stop(name, " is larger than the largest representable integer.",
      call. = FALSE
    )
  }
  if (x != trunc(x)) {
    stop(name, " must be a whole number and was ", x,
      ". Round it explicitly if that is what was intended.",
      call. = FALSE
    )
  }
  as.integer(x)
}

#' Require a genetic covariance to be compatible with its phenotypic covariance
#'
#' A genetic covariance matrix is a component of the phenotypic one, so the
#' residual \eqn{\mathbf{P} - \mathbf{G}} is itself a covariance matrix and
#' must be positive semidefinite. Nothing in the package enforced this, and
#' every consequence of violating it is a number that looks valid and is not: a
#' heritability above one, an accuracy above one, a squared correlation above
#' one, or a negative mean squared prediction error.
#'
#' Checking the diagonal alone is not enough. `diag(G) <= diag(P)` for every
#' trait is necessary but not sufficient, because a linear combination of
#' traits can have a genetic variance exceeding its phenotypic variance while
#' no single trait does. The eigenvalues of the difference are the right test.
#'
#' @param G,P Covariance matrices in the same trait space and order.
#' @param context Where the failure arose, used in the message.
#' @param g_name,p_name Names to use when describing the matrices.
#' @param tolerance Relative tolerance for the smallest eigenvalue, allowing
#'   for floating-point noise in an otherwise admissible pair.
#'
#' @return `TRUE` invisibly, or an error.
#' @noRd
.dgr_check_compatible <- function(
  G, P, context, g_name = "G", p_name = "P", tolerance = 1e-8
) {
  if (is.null(G) || is.null(P)) {
    return(invisible(TRUE))
  }
  if (!identical(dim(G), dim(P))) {
    stop(g_name, " and ", p_name, " must have the same dimensions in ",
      context, ".",
      call. = FALSE
    )
  }
  residual <- P - G
  residual <- (residual + t(residual)) / 2
  values <- eigen(residual, symmetric = TRUE, only.values = TRUE)$values
  scale <- max(1, max(abs(diag(P))))
  if (min(values) >= -tolerance * scale) {
    return(invisible(TRUE))
  }

  # Name the offending combination, because "the matrices are incompatible" is
  # not actionable and the single-trait check often passes.
  offending <- eigen(residual, symmetric = TRUE)$vectors[, which.min(values)]
  names(offending) <- rownames(P)
  worst_trait <- names(which.max(diag(G) / diag(P)))
  ratio <- max(diag(G) / diag(P))
  stop(
    p_name, " - ", g_name, " is not positive semidefinite in ", context,
    ", so the two matrices cannot both be true.\n",
    "  Its smallest eigenvalue is ", format(min(values), digits = 3),
    ", meaning some combination of traits has more genetic variance than ",
    "phenotypic variance.\n",
    "  Largest single-trait ratio diag(", g_name, ")/diag(", p_name, "): ",
    format(ratio, digits = 3), " for ", worst_trait,
    if (ratio <= 1) {
      " (no single trait exceeds 1, so the problem is in the correlations)"
    } else {
      ""
    },
    ".\n",
    "  Re-estimate the pair jointly, or repair the residual with ",
    "bend_covariance() and rebuild ", p_name, " as ", g_name, " + residual.",
    call. = FALSE
  )
}

#' Clamp a correlation that has drifted marginally outside its range
#'
#' Applied only after `.dgr_check_compatible()` has confirmed the inputs are
#' admissible, so anything beyond floating-point noise is a real error and is
#' allowed to propagate rather than being quietly hidden. The link form
#' `[...]` is not used here because an internal function has no manual page for
#' it to resolve to.
#'
#' @noRd
.dgr_clamp_correlation <- function(x, name, tolerance = 1e-6) {
  if (!is.finite(x)) {
    return(x)
  }
  if (x > 1 + tolerance || x < -1 - tolerance) {
    stop(name, " evaluated to ", format(x, digits = 6),
      ", which is outside the range of a correlation by more than ",
      "rounding error. This indicates inconsistent covariance inputs ",
      "rather than a numerical artefact.",
      call. = FALSE
    )
  }
  min(1, max(-1, x))
}

#' Validate a count that may legitimately be zero
#' @noRd
.dgr_non_negative_integer <- function(x, name) {
  if (length(x) != 1L) {
    stop(name, " must be a single number, and was of length ", length(x), ".",
      call. = FALSE
    )
  }
  if (is.atomic(x) && is.na(x)) {
    stop(name, " must be a non-negative integer and was NA.", call. = FALSE)
  }
  if (!is.numeric(x) || !is.finite(x)) {
    stop(name, " must be a single finite non-negative integer.", call. = FALSE)
  }
  if (x < 0) {
    stop(name, " must be non-negative and was ", x, ".", call. = FALSE)
  }
  if (x > .Machine$integer.max) {
    stop(name, " is larger than the largest representable integer.",
      call. = FALSE
    )
  }
  if (x != trunc(x)) {
    stop(name, " must be a whole number and was ", x, ".", call. = FALSE)
  }
  as.integer(x)
}

#' Validate a random seed
#' @noRd
.dgr_seed <- function(x, name = "seed") {
  if (is.null(x)) {
    return(NULL)
  }
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
    x != trunc(x) || abs(x) > .Machine$integer.max) {
    stop(name, " must be a single whole number within the integer range.",
      call. = FALSE
    )
  }
  as.integer(x)
}

.dgr_named_vector <- function(x, trait_cols, name) {
  if (!is.numeric(x) || is.null(names(x))) {
    stop(name, " must be a named numeric vector.", call. = FALSE)
  }
  missing_names <- setdiff(trait_cols, names(x))
  if (length(missing_names)) {
    stop(name, " is missing: ", paste(missing_names, collapse = ", "),
      call. = FALSE
    )
  }
  x <- as.numeric(x[trait_cols])
  names(x) <- trait_cols
  if (any(!is.finite(x))) {
    stop(name, " must contain finite values.", call. = FALSE)
  }
  x
}

.dgr_covariance <- function(x, trait_cols, name) {
  if (!is.matrix(x) || any(dim(x) != length(trait_cols))) {
    stop(name, " must be a square matrix matching trait_cols.", call. = FALSE)
  }
  if (!is.null(rownames(x)) && !is.null(colnames(x))) {
    if (!all(trait_cols %in% rownames(x)) ||
      !all(trait_cols %in% colnames(x))) {
      stop(name, " dimnames must contain every trait.", call. = FALSE)
    }
    x <- x[trait_cols, trait_cols, drop = FALSE]
  }
  storage.mode(x) <- "double"
  if (any(!is.finite(x))) {
    stop(name, " must contain finite values.", call. = FALSE)
  }
  x <- (x + t(x)) / 2
  dimnames(x) <- list(trait_cols, trait_cols)
  x
}

.dgr_check_psd <- function(x, name, tolerance = 1e-8) {
  minimum <- min(eigen(x, symmetric = TRUE, only.values = TRUE)$values)
  scale <- max(1, max(abs(diag(x))))
  if (minimum < -tolerance * scale) {
    stop(
      name, " must be positive semidefinite; its smallest eigenvalue is ",
      format(minimum, digits = 3), ".\n",
      "  Multi-trait REML routinely returns such a matrix when the traits ",
      "outnumber the information available to estimate their covariances.\n",
      "  Repair it with bend_covariance() and report the adjustment, or ",
      "supply a matrix estimated on fewer traits.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.dgr_prepare_values <- function(
  init_data,
  cand_data,
  ref_data,
  validation_data,
  trait_cols,
  id_col,
  lower_is_better,
  scale_traits,
  centre_traits = scale_traits,
  missing_policy
) {
  init <- data.table::as.data.table(data.table::copy(init_data))
  cand <- data.table::as.data.table(data.table::copy(cand_data))
  ref <- data.table::as.data.table(data.table::copy(
    if (is.null(ref_data)) cand_data else ref_data
  ))
  if (!id_col %in% names(init) || !id_col %in% names(cand)) {
    stop("id_col must occur in init_data and cand_data.", call. = FALSE)
  }
  if (anyDuplicated(init[[id_col]]) || anyDuplicated(cand[[id_col]])) {
    stop("Candidate identifiers must be unique.", call. = FALSE)
  }
  for (object_name in c("cand", "ref")) {
    object <- get(object_name)
    absent <- setdiff(trait_cols, names(object))
    if (length(absent)) {
      stop(object_name, " is missing trait columns: ",
        paste(absent, collapse = ", "),
        call. = FALSE
      )
    }
    for (trait in trait_cols) {
      converted <- suppressWarnings(as.numeric(object[[trait]]))
      introduced <- is.na(converted) & !is.na(object[[trait]])
      if (any(introduced)) {
        stop("Non-numeric values occur in ", object_name, "$", trait, ".",
          call. = FALSE
        )
      }
      data.table::set(object, j = trait, value = converted)
    }
    assign(object_name, object)
  }

  common <- intersect(init[[id_col]], cand[[id_col]])
  if (!length(common)) stop("No candidate identifiers overlap.", call. = FALSE)
  cand <- cand[match(common, cand[[id_col]])]
  init <- init[match(common, init[[id_col]])]

  missing_counts <- list(
    candidate = colSums(is.na(cand[, ..trait_cols])),
    reference = colSums(is.na(ref[, ..trait_cols]))
  )
  if (missing_policy == "error" &&
    (any(missing_counts$candidate > 0) ||
      any(missing_counts$reference > 0))) {
    stop(
      "Missing trait values detected. Choose an explicit missing_policy to continue.",
      call. = FALSE
    )
  }
  if (missing_policy == "complete_cases") {
    keep_cand <- stats::complete.cases(cand[, ..trait_cols])
    keep_ref <- stats::complete.cases(ref[, ..trait_cols])
    cand <- cand[keep_cand]
    init <- init[keep_cand]
    ref <- ref[keep_ref]
  }
  if (missing_policy == "mean_impute") {
    for (trait in trait_cols) {
      ref_mean <- mean(ref[[trait]], na.rm = TRUE)
      if (!is.finite(ref_mean)) {
        stop("Trait ", trait, " has no finite reference values.", call. = FALSE)
      }
      data.table::set(
        cand, which(is.na(cand[[trait]])), trait, ref_mean
      )
      data.table::set(
        ref, which(is.na(ref[[trait]])), trait, ref_mean
      )
    }
  }
  if (!nrow(cand) || nrow(ref) < 2L) {
    stop("Insufficient complete candidate or reference records.", call. = FALSE)
  }

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
  ref_raw <- as.matrix(ref[, ..trait_cols])
  cand_raw <- as.matrix(cand[, ..trait_cols])
  centre <- if (isTRUE(centre_traits)) colMeans(ref_raw) else rep(0, length(trait_cols))
  scale <- if (isTRUE(scale_traits)) apply(ref_raw, 2L, stats::sd) else rep(1, length(trait_cols))
  if (any(!is.finite(scale)) || any(scale <= 0)) {
    stop("Every scaled trait must have positive finite reference variation.",
      call. = FALSE
    )
  }
  transform_matrix <- function(x) {
    x <- sweep(x, 2L, centre, "-")
    x <- sweep(x, 2L, scale, "/")
    sweep(x, 2L, direction, "*")
  }

  validation_matrix <- NULL
  if (!is.null(validation_data)) {
    validation <- data.table::as.data.table(data.table::copy(validation_data))
    absent <- setdiff(trait_cols, names(validation))
    if (length(absent)) {
      stop("validation_data is missing: ", paste(absent, collapse = ", "),
        call. = FALSE
      )
    }
    validation_matrix <- as.matrix(validation[, ..trait_cols])
    storage.mode(validation_matrix) <- "double"
    if (anyNA(validation_matrix)) {
      stop("validation_data must be complete; it is never imputed.",
        call. = FALSE
      )
    }
    validation_matrix <- transform_matrix(validation_matrix)
  }

  list(
    init_data = init,
    candidate_matrix = transform_matrix(cand_raw),
    reference_matrix = transform_matrix(ref_raw),
    reference_ids = if (id_col %in% names(ref)) ref[[id_col]] else NULL,
    validation_matrix = validation_matrix,
    direction = direction,
    scale = scale,
    transformation = list(
      centre = stats::setNames(centre, trait_cols),
      scale = stats::setNames(scale, trait_cols),
      direction = direction
    ),
    missing_data = c(list(policy = missing_policy), missing_counts)
  )
}

.dgr_jaccard_matrix <- function(selected) {
  n <- ncol(selected)
  result <- matrix(1, n, n)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      union <- sum(selected[, i] | selected[, j])
      result[i, j] <- if (union == 0L) {
        1
      } else {
        sum(selected[, i] & selected[, j]) / union
      }
    }
  }
  labels <- paste0("replicate_", seq_len(n))
  dimnames(result) <- list(labels, labels)
  result
}
