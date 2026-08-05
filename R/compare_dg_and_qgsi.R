#' Compare an iteratively optimised DGSI with a QGSI
#'
#' Merges results from [run_dgsi()] and [run_qgsi()] by candidate identifier.
#' It reports scores, ranks, and the selection flag from each fit. DGSI means
#' desired-gain selection index. The `run_dgsi()` object contains the established
#' desired-gain index after the iterative search applied by Joukhadar et al.
#' (2024). QGSI means quadratic genomic selection index.
#'
#' The comparison is descriptive. DGSI targets a desired response direction.
#' QGSI represents linear, squared, and cross-product economic value. Agreement
#' does not validate either objective and disagreement does not identify a
#' winner. The breeder must first decide which objective represents the
#' programme.
#'
#' @param dg_result A result returned by [run_dgsi()].
#' @param qgsi_result A result returned by [run_qgsi()].
#' @param id_col Candidate identifier column.
#' @param include_metadata Whether to retain DGSI metadata columns.
#' @param compute_rank_differences Whether to add signed and absolute rank
#'   differences.
#' @param add_correlation_summary Whether to calculate Pearson score and
#'   Spearman rank correlations.
#' @param sort_by Output sorting rule.
#' @param debug Whether to print progress messages.
#'
#' @return A list with `comparison_table`, `correlation_summary`, and
#'   `decision_summary`. The last component reports the common candidate count,
#'   selected counts, selected overlap, and Jaccard similarity.
#' @export
compare_dg_and_qgsi <- function(
  dg_result,
  qgsi_result,
  id_col = "GenoID",
  include_metadata = TRUE,
  compute_rank_differences = TRUE,
  add_correlation_summary = TRUE,
  sort_by = c("DG_rank", "QGSI_rank", "DG", "QGSI", "none"),
  debug = FALSE
) {
  sort_by <- match.arg(sort_by)
  if (is.null(dg_result$ranked_geno)) {
    stop("dg_result does not contain ranked_geno.", call. = FALSE)
  }
  if (is.null(qgsi_result$ranked_geno)) {
    stop("qgsi_result does not contain ranked_geno.", call. = FALSE)
  }

  dg <- data.table::as.data.table(data.table::copy(dg_result$ranked_geno))
  qg <- data.table::as.data.table(data.table::copy(qgsi_result$ranked_geno))
  if (!id_col %in% names(dg) || !id_col %in% names(qg)) {
    stop("id_col must occur in both ranked_geno tables.", call. = FALSE)
  }
  if (anyDuplicated(dg[[id_col]]) || anyDuplicated(qg[[id_col]])) {
    stop("Candidate identifiers must be unique in both results.",
      call. = FALSE
    )
  }

  if (!"Rank" %in% names(dg) && "SelectionIndex" %in% names(dg)) {
    dg[, Rank := data.table::frank(-SelectionIndex, ties.method = "average")]
  }
  legacy_qg <- c(
    LinearDGPart = "LinearPart",
    QuadraticDGPart = "QuadraticPart",
    QGSI_DG = "QGSI",
    Rank_QGSI_DG = "Rank"
  )
  for (old in names(legacy_qg)) {
    new <- unname(legacy_qg[[old]])
    if (!new %in% names(qg) && old %in% names(qg)) {
      data.table::setnames(qg, old, new)
    }
  }

  required_dg <- c("SelectionIndex", "Rank")
  required_qg <- c("LinearPart", "QuadraticPart", "QGSI", "Rank")
  if (length(setdiff(required_dg, names(dg)))) {
    stop("dg_result lacks canonical score or rank columns.", call. = FALSE)
  }
  if (length(setdiff(required_qg, names(qg)))) {
    stop("qgsi_result lacks canonical component, score, or rank columns.",
      call. = FALSE
    )
  }

  dg_keep <- if (isTRUE(include_metadata)) {
    names(dg)
  } else {
    intersect(c(id_col, "SelectionIndex", "Rank", "Selected"), names(dg))
  }
  dg_sub <- dg[, dg_keep, with = FALSE]
  data.table::setnames(
    dg_sub,
    intersect(c("SelectionIndex", "Rank", "Selected"), names(dg_sub)),
    paste0(
      "DG_",
      intersect(c("SelectionIndex", "Rank", "Selected"), names(dg_sub))
    )
  )

  qg_keep <- intersect(
    c(id_col, "LinearPart", "QuadraticPart", "QGSI", "Rank", "Selected"),
    names(qg)
  )
  qg_sub <- qg[, qg_keep, with = FALSE]
  data.table::setnames(
    qg_sub,
    intersect(
      c("LinearPart", "QuadraticPart", "Rank", "Selected"),
      names(qg_sub)
    ),
    paste0(
      "QGSI_",
      intersect(
        c("LinearPart", "QuadraticPart", "Rank", "Selected"),
        names(qg_sub)
      )
    )
  )

  merged <- merge(dg_sub, qg_sub, by = id_col, all = TRUE, sort = FALSE)
  if (isTRUE(compute_rank_differences)) {
    merged[, RankDiff_DG_minus_QGSI := DG_Rank - QGSI_Rank]
    merged[, AbsRankDiff_DG_vs_QGSI := abs(RankDiff_DG_minus_QGSI)]
  }

  correlation_summary <- NULL
  if (isTRUE(add_correlation_summary)) {
    score_ok <- is.finite(merged$DG_SelectionIndex) &
      is.finite(merged$QGSI)
    rank_ok <- is.finite(merged$DG_Rank) &
      is.finite(merged$QGSI_Rank)
    rows <- list()
    if (sum(score_ok) >= 3L) {
      rows[[length(rows) + 1L]] <- data.table::data.table(
        Comparison = "DGSI score vs QGSI score",
        Measure = "Pearson",
        Correlation = stats::cor(
          merged$DG_SelectionIndex[score_ok],
          merged$QGSI[score_ok]
        ),
        N = sum(score_ok)
      )
    }
    if (sum(rank_ok) >= 3L) {
      rows[[length(rows) + 1L]] <- data.table::data.table(
        Comparison = "DGSI rank vs QGSI rank",
        Measure = "Spearman",
        Correlation = stats::cor(
          merged$DG_Rank[rank_ok],
          merged$QGSI_Rank[rank_ok],
          method = "spearman"
        ),
        N = sum(rank_ok)
      )
    }
    if (length(rows)) {
      correlation_summary <- data.table::rbindlist(rows)
    }
  }

  if (sort_by == "DG") data.table::setorder(merged, -DG_SelectionIndex)
  if (sort_by == "QGSI") data.table::setorder(merged, -QGSI)
  if (sort_by == "DG_rank") data.table::setorder(merged, DG_Rank)
  if (sort_by == "QGSI_rank") data.table::setorder(merged, QGSI_Rank)

  has_selection <- all(c("DG_Selected", "QGSI_Selected") %in% names(merged))
  decision_summary <- if (has_selection) {
    common <- !is.na(merged$DG_Selected) & !is.na(merged$QGSI_Selected)
    dg_selected <- common & merged$DG_Selected
    qgsi_selected <- common & merged$QGSI_Selected
    overlap <- sum(dg_selected & qgsi_selected)
    union_size <- sum(dg_selected | qgsi_selected)
    data.table::data.table(
      Common_candidates = sum(common),
      DGSI_selected = sum(dg_selected),
      QGSI_selected = sum(qgsi_selected),
      Selected_overlap = overlap,
      Selected_Jaccard = if (union_size > 0L) overlap / union_size else NA_real_
    )
  } else {
    data.table::data.table(
      Common_candidates = sum(
        is.finite(merged$DG_Rank) & is.finite(merged$QGSI_Rank)
      ),
      DGSI_selected = NA_integer_,
      QGSI_selected = NA_integer_,
      Selected_overlap = NA_integer_,
      Selected_Jaccard = NA_real_
    )
  }
  .desiredgainr_dbg(
    debug, "Compared %d candidates across DGSI and QGSI.", nrow(merged)
  )
  list(
    comparison_table = merged,
    correlation_summary = correlation_summary,
    decision_summary = decision_summary,
    interpretation = paste(
      "Agreement is descriptive because DGSI targets desired responses",
      "whereas QGSI predicts a quadratic economic merit."
    )
  )
}
