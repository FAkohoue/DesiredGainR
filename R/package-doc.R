#' DesiredGainR: Auditable Multi-trait Selection Indices
#'
#' Package providing tools for:
#' \itemize{
#'   \item desired-gain selection indices, including iterative optimisation of
#'     the input desired-gain direction,
#'   \item quadratic genomic selection index computation from trait genomic
#'     estimated breeding values, explicit economic weights, genomic covariance
#'     estimation, and published response diagnostics,
#'   \item cross-family comparison of classical, restricted, general,
#'     desired-gain, and quadratic genomic candidate decisions,
#'   \item a high-level wrapper for full pipeline execution.
#' }
#'
#' @import data.table
#' @importFrom utils globalVariables
#' @importFrom stats predict
#' @keywords internal
"_PACKAGE"
