#' Run DGSI, QGSI, or both workflows
#'
#' This convenience wrapper keeps the two breeding objectives separate.
#' `dg` contains desired responses for [run_dgsi()], whereas
#' `qgsi_linear_weights` and `W` contain the linear and quadratic economic
#' weights for [run_qgsi()]. A desired-gain vector is never reused as a QGSI
#' economic-weight vector.
#'
#' @param mode One of `"both"`, `"dg"`, or `"qgsi"`.
#' @param init_data Candidate identifiers and metadata.
#' @param trait_cols Trait columns shared by the requested workflows.
#' @param id_col Candidate identifier column.
#' @param dg Named desired-gain vector required for DGSI modes.
#' @param cand_data Candidate genetic values required for DGSI modes.
#' @param ref_data Optional DGSI reference values.
#' @param P Optional DGSI phenotypic or index-variable covariance matrix.
#' @param G Genetic covariance matrix required for DGSI modes.
#' @param select_mode DGSI selection mode.
#' @param n_select DGSI selection count.
#' @param trait_min_sd Optional favourable-direction DGSI eligibility
#'   thresholds.
#' @param fallback_to_top_n Whether DGSI falls back when no candidate is
#'   eligible.
#' @param n_iter,n_rep,sd_scale,seed,ridge_P,ridge_M DGSI optimisation
#'   controls passed to [run_dgsi()].
#' @param dg_scale_traits Whether DGSI values are centred and scaled.
#' @param gebv_data Candidate GEBVs required for QGSI modes.
#' @param qgsi_linear_weights Named QGSI linear economic weights.
#' @param W Symmetric QGSI squared/cross-product economic-weight matrix.
#' @param qgsi_reference_gebv_data Optional reference GEBVs for estimating
#'   `Gamma`.
#' @param Gamma Optional genomic covariance matrix.
#' @param relationship_matrix Optional named genomic relationship matrix.
#' @param true_G Optional true genetic covariance matrix for simulation-based
#'   QGSI accuracy and MSPE calculations.
#' @param qgsi_center_traits,qgsi_scale_traits QGSI trait transformations.
#' @param qgsi_missing_policy QGSI missing-value policy.
#' @param qgsi_n_select Optional QGSI selection count.
#' @param qgsi_selection_proportion Optional QGSI selected proportion.
#' @param lower_is_better Traits for which smaller original values are
#'   favourable.
#' @param merge_outputs Whether to merge rankings when `mode = "both"`.
#' @param compare_sort_by Sorting rule passed to [compare_dg_and_qgsi()].
#' @param debug Whether to print progress messages.
#'
#' @return A list containing the requested `dg_result`, `qgsi_result`, and,
#'   when requested, `comparison_result`.
#'
#' @examples
#' set.seed(7)
#' traits <- c("yield", "height")
#' values <- data.frame(
#'   GenoID = paste0("G", 1:20),
#'   yield = rnorm(20),
#'   height = rnorm(20)
#' )
#' G <- stats::cov(values[traits])
#' dimnames(G) <- list(traits, traits)
#' W <- matrix(
#'   c(0.05, 0, 0, -0.03), 2,
#'   dimnames = list(traits, traits)
#' )
#' result <- run_dgsi_qgsi_pipeline(
#'   mode = "both",
#'   init_data = values["GenoID"],
#'   trait_cols = traits,
#'   dg = c(yield = 0.5, height = 0.2),
#'   cand_data = values,
#'   G = G,
#'   n_select = 4,
#'   n_iter = 20,
#'   n_rep = 2,
#'   gebv_data = values,
#'   qgsi_linear_weights = c(yield = 1, height = 0.2),
#'   W = W,
#'   qgsi_n_select = 4,
#'   debug = FALSE
#' )
#'
#' @export
run_dgsi_qgsi_pipeline <- function(
    mode = c("both", "dg", "qgsi"),
    init_data,
    trait_cols,
    id_col = "GenoID",
    dg = NULL,
    cand_data = NULL,
    ref_data = NULL,
    P = NULL,
    G = NULL,
    select_mode = c("top_n", "trait_thresholds"),
    n_select = 100L,
    trait_min_sd = NULL,
    fallback_to_top_n = TRUE,
    n_iter = 1000L,
    n_rep = 20L,
    sd_scale = 1,
    seed = 42L,
    ridge_P = 1e-6,
    ridge_M = 1e-6,
    dg_scale_traits = FALSE,
    gebv_data = NULL,
    qgsi_linear_weights = NULL,
    W = NULL,
    qgsi_reference_gebv_data = NULL,
    Gamma = NULL,
    relationship_matrix = NULL,
    true_G = NULL,
    qgsi_center_traits = TRUE,
    qgsi_scale_traits = FALSE,
    qgsi_missing_policy = c("error", "complete_cases", "mean_impute"),
    qgsi_n_select = NULL,
    qgsi_selection_proportion = NULL,
    lower_is_better = NULL,
    merge_outputs = TRUE,
    compare_sort_by = c("DG_rank", "QGSI_rank", "DG", "QGSI", "none"),
    debug = TRUE
) {
  mode <- match.arg(mode)
  select_mode <- match.arg(select_mode)
  qgsi_missing_policy <- match.arg(qgsi_missing_policy)
  compare_sort_by <- match.arg(compare_sort_by)
  out <- list()

  if (mode %in% c("dg", "both")) {
    if (is.null(cand_data)) {
      stop("cand_data is required when mode includes 'dg'.", call. = FALSE)
    }
    if (is.null(dg)) {
      stop("dg is required when mode includes 'dg'.", call. = FALSE)
    }
    if (is.null(G)) {
      stop("G is required when mode includes 'dg'.", call. = FALSE)
    }
    out$dg_result <- run_dgsi(
      init_data = init_data,
      cand_data = cand_data,
      trait_cols = trait_cols,
      dg = dg,
      P = P,
      G = G,
      ref_data = ref_data,
      id_col = id_col,
      scale_traits = dg_scale_traits,
      lower_is_better = lower_is_better,
      select_mode = if (select_mode == "trait_thresholds") {
        "eligible_top_n"
      } else {
        "top_n"
      },
      n_select = n_select,
      trait_min = trait_min_sd,
      empty_eligibility = if (isTRUE(fallback_to_top_n)) {
        "fallback_top_n"
      } else {
        "error"
      },
      n_iter = n_iter,
      n_rep = n_rep,
      sd_scale = sd_scale,
      seed = seed,
      ridge_P = ridge_P,
      ridge_M = ridge_M,
      debug = debug
    )
  }

  if (mode %in% c("qgsi", "both")) {
    if (is.null(gebv_data)) {
      stop("gebv_data is required when mode includes 'qgsi'.",
           call. = FALSE)
    }
    if (is.null(qgsi_linear_weights)) {
      stop(
        paste(
          "qgsi_linear_weights is required when mode includes 'qgsi';",
          "dg is not reused as an economic-weight vector."
        ),
        call. = FALSE
      )
    }
    if (is.null(W)) {
      stop("W is required when mode includes 'qgsi'.", call. = FALSE)
    }
    out$qgsi_result <- run_qgsi(
      init_data = init_data,
      gebv_data = gebv_data,
      trait_cols = trait_cols,
      linear_weights = qgsi_linear_weights,
      W = W,
      id_col = id_col,
      reference_gebv_data = qgsi_reference_gebv_data,
      Gamma = Gamma,
      relationship_matrix = relationship_matrix,
      true_G = true_G,
      lower_is_better = lower_is_better,
      center_traits = qgsi_center_traits,
      scale_traits = qgsi_scale_traits,
      missing_policy = qgsi_missing_policy,
      n_select = qgsi_n_select,
      selection_proportion = qgsi_selection_proportion,
      debug = debug
    )
  }

  if (mode == "both" && isTRUE(merge_outputs)) {
    out$comparison_result <- compare_dg_and_qgsi(
      dg_result = out$dg_result,
      qgsi_result = out$qgsi_result,
      id_col = id_col,
      sort_by = compare_sort_by,
      debug = debug
    )
  }
  out
}
