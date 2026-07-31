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
#' @param validation_data Optional independent candidate values used only to
#'   evaluate the winning coefficients after optimisation.
#' @param return_all_reps Whether to retain full replicate results.
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
#' derivative-free search is used. The replicate with the lowest objective is
#' returned; because that choice is made on the same candidates that are then
#' selected, the reported objective is optimistically biased, and the result
#' depends on `n_rep`. Use `validation_data` for an unbiased evaluation of the
#' winning coefficients. `run_dgsi()` seeds the RNG from `seed` and restores
#' the caller's RNG state before returning.
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
#'
#' fit <- run_dgsi(
#'   init_data = candidates["GenoID"],
#'   cand_data = candidates,
#'   trait_cols = traits,
#'   dg = c(yield = 0.6, disease = 0.4),
#'   G = G,
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
    debug = FALSE
) {
  select_mode <- match.arg(select_mode)
  empty_eligibility <- match.arg(empty_eligibility)
  missing_policy <- match.arg(missing_policy)

  if (!length(trait_cols) || anyDuplicated(trait_cols)) {
    stop("trait_cols must contain unique trait names.", call. = FALSE)
  }
  if (!is.numeric(n_select) || length(n_select) != 1L ||
      !is.finite(n_select) || n_select < 1) {
    stop("n_select must be a positive integer.", call. = FALSE)
  }
  n_select <- as.integer(n_select)
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
         call. = FALSE)
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
           call. = FALSE)
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
  P_was_estimated <- is.null(P)
  if (P_was_estimated) {
    P <- stats::cov(Xref)
    P_source <- "empirical working covariance estimated from ref_data"
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

  if (select_mode == "eligible_top_n") {
    if (is.null(trait_min)) {
      stop("trait_min is required when select_mode = 'eligible_top_n'.",
           call. = FALSE)
    }
    trait_min <- .dgr_named_vector(trait_min, trait_cols, "trait_min")
  }

  eligible <- rep(TRUE, nrow(X))
  if (select_mode == "eligible_top_n") {
    eligible <- rowSums(sweep(X, 2L, trait_min, FUN = "<")) == 0L
    if (!any(eligible)) {
      if (empty_eligibility == "error") {
        stop("No candidate meets all favourable-direction eligibility thresholds.",
             call. = FALSE)
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
  coefficient <- function(d) {
    value <- as.numeric(Xpg %*% solve(middle, d))
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

  set.seed(as.integer(seed))
  replicates <- vector("list", n_rep)
  for (rr in seq_len(n_rep)) {
    centre <- dg
    initial_b <- coefficient(centre)
    initial_eval <- realised_response(as.numeric(X %*% initial_b))
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
      b <- coefficient(proposal)
      evaluation <- realised_response(as.numeric(X %*% b))
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
    window_start <- max(1L, n_iter - plateau_window + 1L)
    relative_change <- (
      best_trace[window_start] - best_trace[n_iter]
    ) / max(abs(best_trace[window_start]), 1e-12)
    replicates[[rr]] <- c(
      best,
      list(
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
  best_replicate <- which.min(objectives)
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
      G = G_provenance
    ),
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
    objective = best$objective,
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
  class(result) <- c("desired_gain_index", "list")
  .desiredgainr_dbg(
    debug,
    "Selected replicate %d automatically (objective %.6g).",
    best_replicate,
    best$objective
  )
  result
}

.dgr_positive_integer <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < 1 || x != as.integer(x)) {
    stop(name, " must be a positive integer.", call. = FALSE)
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
         call. = FALSE)
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
    stop(name, " must be positive semidefinite.", call. = FALSE)
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
           paste(absent, collapse = ", "), call. = FALSE)
    }
    for (trait in trait_cols) {
      converted <- suppressWarnings(as.numeric(object[[trait]]))
      introduced <- is.na(converted) & !is.na(object[[trait]])
      if (any(introduced)) {
        stop("Non-numeric values occur in ", object_name, "$", trait, ".",
             call. = FALSE)
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
           paste(unknown, collapse = ", "), call. = FALSE)
    }
    direction[lower_is_better] <- -1
  }
  ref_raw <- as.matrix(ref[, ..trait_cols])
  cand_raw <- as.matrix(cand[, ..trait_cols])
  centre <- if (isTRUE(centre_traits)) colMeans(ref_raw) else rep(0, length(trait_cols))
  scale <- if (isTRUE(scale_traits)) apply(ref_raw, 2L, stats::sd) else rep(1, length(trait_cols))
  if (any(!is.finite(scale)) || any(scale <= 0)) {
    stop("Every scaled trait must have positive finite reference variation.",
         call. = FALSE)
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
           call. = FALSE)
    }
    validation_matrix <- as.matrix(validation[, ..trait_cols])
    storage.mode(validation_matrix) <- "double"
    if (anyNA(validation_matrix)) {
      stop("validation_data must be complete; it is never imputed.",
           call. = FALSE)
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
      result[i, j] <- if (union == 0L) 1 else
        sum(selected[, i] & selected[, j]) / union
    }
  }
  labels <- paste0("replicate_", seq_len(n))
  dimnames(result) <- list(labels, labels)
  result
}
