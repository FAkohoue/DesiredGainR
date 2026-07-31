# Internal utilities for DesiredGainR
#
# Non-standard-evaluation globals are declared once in R/globals.R.

#' Internal debug message helper
#'
#' @param debug Logical.
#' @param ... Arguments passed to `sprintf()`.
#'
#' @noRd
.desiredgainr_dbg <- function(debug, ...) {
  if (isTRUE(debug)) message(sprintf(...))
}

#' Internal z-score helper
#'
#' Not called by any exported function. Retained only so that external code
#' calling `DesiredGainR:::.desiredgainr_z()` does not break. Scheduled for
#' removal once no reverse dependency relies on it.
#'
#' @param x Numeric vector.
#'
#' @return Numeric vector.
#' @noRd
.desiredgainr_z <- function(x) {
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

#' Validate a square matrix
#'
#' Superseded internally by `.dgr_covariance()` and
#' `.dgr_qg_symmetric_matrix()`, which also check names, finiteness, and
#' symmetry. Retained only so that external code calling
#' `DesiredGainR:::.validate_square_matrix()` does not break. Scheduled for
#' removal once no reverse dependency relies on it.
#'
#' @param M Matrix.
#' @param p Expected dimension.
#' @param nm Object name for messages.
#'
#' @keywords internal
.validate_square_matrix <- function(M, p, nm = "matrix") {
  if (!is.matrix(M)) stop(sprintf("%s must be a matrix.", nm), call. = FALSE)
  if (nrow(M) != p || ncol(M) != p) {
    stop(sprintf("%s must be a %d x %d matrix.", nm, p, p), call. = FALSE)
  }
  invisible(TRUE)
}
