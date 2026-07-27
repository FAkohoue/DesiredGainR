#' @rdname run_qgsi
#' @param dg Deprecated compatibility name for `linear_weights`. It is
#'   interpreted as linear economic weights, not desired genetic gains.
#' @param W_d Deprecated compatibility name for `W`.
#' @param quadratic_diag_weights Optional explicit diagonal economic weights
#'   used to create `W_d` only when a full matrix is absent.
#' @param impute_missing Deprecated logical missing-value switch.
#' @param return_components Whether to return candidate-specific contributions.
#' @export
run_qgsi_desired_gain <- function(
    init_data,
    gebv_data,
    trait_cols,
    id_col = "GenoID",
    dg,
    W_d = NULL,
    quadratic_diag_weights = NULL,
    lower_is_better = NULL,
    center_traits = FALSE,
    scale_traits = FALSE,
    impute_missing = FALSE,
    return_components = TRUE,
    debug = TRUE
) {
  warning(
    paste(
      "run_qgsi_desired_gain() is deprecated because desired gains are not",
      "QGSI economic weights. Use run_qgsi(linear_weights = ..., W = ...).",
      "For compatibility, dg is interpreted as linear_weights in this call."
    ),
    call. = FALSE
  )
  if (is.null(W_d)) {
    if (is.null(quadratic_diag_weights)) {
      stop(
        paste(
          "Supply W_d. A quadratic economic-weight matrix is not",
          "constructed from desired gains."
        ),
        call. = FALSE
      )
    }
    diagonal <- .dgr_named_vector(
      quadratic_diag_weights, trait_cols, "quadratic_diag_weights"
    )
    W_d <- diag(diagonal, nrow = length(trait_cols))
    dimnames(W_d) <- list(trait_cols, trait_cols)
  }
  result <- run_qgsi(
    init_data = init_data,
    gebv_data = gebv_data,
    trait_cols = trait_cols,
    linear_weights = dg,
    W = W_d,
    id_col = id_col,
    lower_is_better = lower_is_better,
    center_traits = center_traits,
    scale_traits = scale_traits,
    missing_policy = if (isTRUE(impute_missing)) "mean_impute" else "error",
    return_contributions = return_components,
    debug = debug
  )

  # Compatibility aliases for scripts written against DesiredGainR 0.2.0.
  result$ranked_geno[, LinearDGPart := LinearPart]
  result$ranked_geno[, QuadraticDGPart := QuadraticPart]
  result$ranked_geno[, QGSI_DG := QGSI]
  result$ranked_geno[, Rank_QGSI_DG := Rank]
  result$dg <- result$linear_weights
  result$W_d <- result$W
  result
}
